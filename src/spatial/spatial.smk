from dataclasses import dataclass

import openpyxl

# The spatial analysis gets its own subdirectory in ./reports and within the data (sub)directories
# Raw data is in data/vectra

###################
# Sample metadata #
###################

# This is unique to the spatial analyses for now, but could be generalized
# to all data types and moved to the main Snakefile.

@dataclass(frozen=True)
class Sample():
    sample_id: str
    patient_id: str

    def from_dict(x):
        if x["sample_type"] == "slide":
            return VectraSample(
                sample_id = x["sample_ID"],
                patient_id = x["patient_ID"],
                batch_HALO = x["batch_HALO"],
                panel = x["experimental_platform"],
            )
        else:
            return Sample(
                sample_id = x["sample_ID"],
                patient_id = x["patient_ID"],
            )
        

@dataclass(frozen=True)
class VectraSample(Sample):
    batch_HALO: str
    panel: str

def read_sample_xlsx(fn):
    # Avoid pandas magic by using openpyxl and Apache Arrow
    # First row should contain a header
    worksheet = openpyxl.load_workbook(fn)["sample_data"]
    rows = worksheet.values
    columns = next(rows)
    samples = list()
    for row in rows:
        row = dict(zip(columns, row))
        # Only include samples that should be processed such that we don't have to exclude them
        # in every rule
        if row["Included"] == 1:
            samples.append(Sample.from_dict(row))
    if len(samples) == 0:
        raise Exception('No samples included in the analysis')
    if len(samples) != len(set(samples)):
        import collections
        c = {s: n for s, n in collections.Counter(samples).items() if n > 1}
        print(c)
        raise Exception('Duplicated samples')
    return samples

samples = read_sample_xlsx("data/metadata/PPBC_metadata.xlsx")
vectra_samples = {s.sample_id: s for s in samples if isinstance(s, VectraSample)}

#################
# Preprocessing #
#################

# Create a manageable data structure
rule organize_vectra:
  input:
    # Exact date and time of Vectra staining
    batch_info = "data/metadata/spatial/MPIF26en27 batches gekleurd.xlsx", 
    example_halo_archive = 
        "data/vectra/halo/Batch 1/MPIF26/T21-60303/Halo archive 2021-06-24 11-33 - v3.2.1851",
    metadata = "data/metadata/PPBC_metadata.xlsx",
    rmd="src/spatial/organize_vectra_samples.Rmd",
    script="src/rmarkdown.R"
  output: "data/metadata/spatial/00_file_location_dictionary.csv"
  log: "src/spatial/organize_vectra_samples.html"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{log}"
    " --vectra_dir data/vectra  --metadata '{input.metadata}'  --batch_info '{input.batch_info}'"
    " --example_halo_archive '{input.example_halo_archive}'"
    " --out '{output}'"

#QC checks based on aggregated summary files
rule summary_qc:
  input:
    filedict="data/metadata/spatial/00_file_location_dictionary.csv",
    rmd="reports/spatial/01_summary_QC.Rmd",
    script="src/rmarkdown.R"
  output:
    mpif26_summary="data/vectra/interim/summaries/01_MPIF26_batch1_summaries.csv",
    mpif27_summary="data/vectra/interim/summaries/01_MPIF27_batch1_summaries.csv",
    html="reports/spatial/01_summary_QC.html"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"

#Convert object results and summary to .nc format
#Do not invoke directly: use all_objects or object_QC instead
#On harris, command must read python3 instead of python
rule load_objects:
  input:
    script = "src/spatial/import_halo_objects.py",
    objects="data/vectra/raw/objects/{t_number}_{panel}_{batch}_object_results.csv",
    summary="data/vectra/raw/summary/{t_number}_{panel}_{batch}_summary_results.csv",
  output:
    objects="data/vectra/interim/objects/{t_number}_{panel}_{batch}.nc",
  shell:
    "export OMP_NUM_THREADS=1\n"
    "python3 {input.script} {input.summary} {input.objects} {output.objects}"
    

rule all_objects:
  input:
    ["data/vectra/interim/objects/{s.sample_id}_{s.panel}_{s.HALO_batch}.nc"
     for s in vectra_samples.values()]

#Generate an overview of marker combination abundance by panel
rule count_cells:
  input:
    "src/spatial/read_cells.R",
    objects = "data/vectra/interim/objects/{t_number}_{panel}_{batch}.nc",
    script="src/spatial/counts_cells.R",
  output:
    nc="results/spatial/marker_cell_counts/{t_number}_{panel}_{batch}.csv",
  shell:
    "mkdir -p results/spatial/marker_cell_counts\n"
    "Rscript {input.script} {input.objects} {output.csv}"


#QC reports based on object files rather than summary files
rule object_QC:
  input:
    rmd="reports/spatial/02_object_QC.Rmd",
    script="src/rmarkdown.R",
    mpif26_objects="data/vectra/interim/MPIF26_df.Rds",
    mpif27_objects="data/vectra/interim/MPIF27_df.Rds",
  output:
    html="reports/spatial/02_object_QC.html",
  shell:
    "Rscript {input.script} {input.rmd} $(realpath -s {output.html})"      

#Identify and visualize problematic marker combinations
#TODO rename this notebook: actual correction now takes place elsewhere
rule report_marker_correction:
  input:
    mpif26="data/vectra/interim/MPIF26_df.Rds",
    mpif27="data/vectra/interim/MPIF27_df.Rds",
    cell_count_by_marker="results/spatial/cell_counts_by_marker.csv",
    rmd="reports/spatial/03a_marker_correction.Rmd",
    script="src/rmarkdown.R"
  output:
    html="reports/spatial/03a_marker_correction.html",
  shell:
    "Rscript {input.script} {input.rmd} $(realpath -s {output.html})"

#Correct double positivity, deal with disallowed marker pairs,
#and assign a "cell type" based on marker positivity
#Do not call this rule alone: it exists to provide input for downstream rules
rule call_cell_types:
  input:
    script="src/spatial/call_cell_types.R",
    objects="data/vectra/interim/{panel}_df.Rds",
  output:
    "data/vectra/processed/objects_{panel}.Rds",
  shell:
    "Rscript {input.script} {input.objects} {wildcards.panel} {output}"

#Ensure disallowed marker pairs no longer exist
rule test_marker_correction_results:
  input:
    mpif26="data/vectra/processed/objects_MPIF26.Rds",
    mpif27="data/vectra/processed/objects_MPIF27.Rds",
    rmd="reports/spatial/03b_test_marker_correction.Rmd",
    script="src/rmarkdown.R"
  output:
    html="reports/spatial/03b_test_marker_correction.html",
  shell:
    "Rscript {input.script} {input.rmd} $(realpath -s {output.html})"

#Quantify cell types and the marker combinations they represent
rule define_cell_types_report:
  input:
    mpif26="data/vectra/processed/objects_MPIF26.Rds",
    mpif27="data/vectra/processed/objects_MPIF27.Rds",
    rmd="reports/spatial/04_define_cell_types.Rmd",
    script="src/rmarkdown.R"
  output:
    html="reports/spatial/04_define_cell_types.html"
  shell:
    "Rscript {input.script} {input.rmd} $(realpath -s {output.html})"
    
#Model overall marker combination density by panel
#Use all_densities instead of invoking this rule directly
rule cell_type_density:
  input:
    script="src/spatial/model_density.R",
    objects="data/vectra/processed/objects_{panel}.Rds",
  output:
    tsv="results/spatial/density/{panel}.tsv",
  shell:
    "mkdir -p results/spatial/density\n"
    "Rscript {input.script} {input.objects} {output}"

panel_densities = expand(
  "results/spatial/density/{panel}.tsv",
  panel = ['MPIF26', 'MPIF27'])

rule all_densities:
  input: panel_densities

rule density_report:
  input:
    panel_densities,
    rmd="reports/spatial/05_density.Rmd",
    script="src/rmarkdown.R",
  output:
    html="reports/spatial/05_density.html"
  shell:
    "Rscript {input.script} {input.rmd} $(realpath -s {output.html})"

rule ppbc_density:
  input:
    mpif26="results/spatial/density/MPIF26.tsv",
    mpif27="results/spatial/density/MPIF27.tsv",
    meta="data/metadata/spatial/00_vectra_metadata.csv",
    rmd="reports/spatial/06_density_outcome.Rmd",
    script="src/rmarkdown.R"
  output:
    html="reports/spatial/06_density_outcome.html"
  shell:
    "Rscript {input.script} {input.rmd} $(realpath -s {output.html})"
