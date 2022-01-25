from dataclasses import dataclass

import openpyxl

# The spatial analysis gets its own subdirectory in ./reports and within the data (sub)directories
# Raw data is in data/vectra


wildcard_constraints:
    # T numbers, with a possible block ID appended
    t_number="T[0-9]+-[0-9]+(_[I0-9]+)?",
    panel="MPIF[0-9]+",
    batch="batch[0-9]+"


###################
# Sample metadata #
###################

# This is unique to the spatial analyses for now, but could be generalized
# to all data types and moved to the main Snakefile.

@dataclass(frozen=True)
class Sample():
    sample_id: str
    patient_id: str
    included: bool

    def from_dict(x):
        if x["sample_type"] == "slide":
            return VectraSample(
                sample_id = x["sample_ID"],
                patient_id = x["patient_ID"],
                included = x["Included"] == 1,
                batch_HALO = x["batch_HALO"].lower().replace(" ", "") if x["batch_HALO"] is not None else None,
                panel = x["experimental_platform"],
            )
        else:
            assert x["sample_type"] in ["RNA"]
            return Sample(
                sample_id = x["sample_ID"],
                patient_id = x["patient_ID"],
                included = x["Included"] == 1,
            )
        

@dataclass(frozen=True)
class VectraSample(Sample):
    batch_HALO: str
    panel: str

def read_sample_xlsx(fn):
    # Avoid pandas magic by using openpyxl
    # First row should contain a header
    worksheet = openpyxl.load_workbook(fn)["sample_data"]
    rows = worksheet.values
    columns = next(rows)
    samples = list()
    for row in rows:
        row = dict(zip(columns, row))
        samples.append(Sample.from_dict(row))
    if len(samples) == 0:
        raise Exception('No samples')
    if len(samples) != len(set(samples)):
        import collections
        c = {s: n for s, n in collections.Counter(samples).items() if n > 1}
        print(c)
        raise Exception('Duplicated samples')
    return samples

samples = read_sample_xlsx("data/metadata/PPBC_metadata.xlsx")

# Only include samples that should be processed such that we don't have to exclude them
# in every rule
vectra_samples = [s for s in samples if isinstance(s, VectraSample) and s.included]
if len(vectra_samples) == 0:
    raise Exception('No vectra samples included in the analysis')

#FIXME: Check this sample, something goes wrong with calling
vectra_samples = [s for s in vectra_samples if s.sample_id != "T20-62429"]

#################
# Preprocessing #
#################

# Create a manageable data structure
rule organize_vectra:
  input:
    # Exact date and time of Vectra staining
    batch_info = "data/metadata/spatial/MPIF26en27 batches gekleurd.xlsx", 
    example_halo_archive = 
        "data/vectra/halo/Batch 1/MPIF26/T21-60304/Halo archive 2021-06-21 12-17 - v3.2.1851/",
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
rule summary_QC:
  input:
    filedict="data/metadata/spatial/00_file_location_dictionary.csv",
    rmd="reports/spatial/01_summary_QC.Rmd",
    script="src/rmarkdown.R"
  output:
    html="reports/spatial/01_summary_QC.html"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"

#Convert object results and summary to .nc format
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
    
#Create a single data frame that includes object info for all samples
#May be extremely slow to load due to size
rule aggregrate_objects:
  input:
    object_files = [f"data/vectra/interim/objects/{s.sample_id}_{s.panel}_{s.batch_HALO}.nc"
                    for s in vectra_samples],
    script = "src/spatial/read_cells.R"
  output:
    "data/vectra/interim/objects.Rds"
  shell:
    "Rscript {input.script} {input.object_files} {output}"

#QC reports based on object files rather than summary files
rule object_QC:
  input:
    #For Tycho: something like this?
    #nc = lambda wildcards: glob(`data/vectra/interim/objects/.*{wildcards.batch}.nc`)
    rmd = "reports/spatial/02_object_QC.Rmd"
  output:
    html="reports/spatial/object_qc_by_batch/02_object_QC_{batch}.html"
  shell:
    "mkdir -p reports/spatial/object_qc_by_batch\n"
    "Rscript -e \"rmarkdown::render('{input.rmd}'," 
    "params=list(batch='{wildcards.batch}'),"
    "output_file = here::here('{output.html}'))\""
    
#Generate an overview of marker combination abundance
rule count_marker_combos:
  input:
    "src/spatial/read_cells.R",
    objects = "data/vectra/interim/objects/{t_number}_{panel}_{batch}.nc",
    script="src/spatial/counts_cells.R",
  output: "results/spatial/marker_cell_counts/{t_number}_{panel}_{batch}.csv",
  shell:
    "mkdir -p results/spatial/marker_cell_counts\n"
    "Rscript {input.script} {input.objects} {output}"
    
#Aggregate marker combinations by panel, and report on them
rule aggregate_markers:
  input:
    counts = [f"results/spatial/marker_cell_counts/{s.sample_id}_{s.panel}_{s.batch_HALO}.csv"
                    for s in vectra_samples],
    rmd = "reports/spatial/03_aggregate_marker_combos.Rmd",
    lib="src/spatial/aggregate_counts.R",
    script="src/rmarkdown.R"
  output: 
    marker_region="results/spatial/marker_combos_by_region.csv",
    html = "reports/spatial/03_aggregate_marker_combos.html"
  shell: 
    "Rscript {input.script} {input.rmd} $(realpath -s {output.html})"

#Visualize problematic marker combinations by batch
rule batch_marker_coexpression:
  input:
    #FIXME
    #nc = lambda wildcards: glob(`data/vectra/interim/objects/.*{wildcards.batch}.csv`),
    rmd="reports/spatial/03_batch_marker_coexpression.Rmd"
  output:
    html="reports/spatial/batch_marker_viz/03_{batch}_marker_coexpression.html",
  shell:
    "mkdir -p reports/spatial/batch_marker_viz\n"
    "Rscript -e \"rmarkdown::render('{input.rmd}'," 
    "params=list(batch='{wildcards.batch}'),"
    "output_file = here::here('{output.html}'))\""

#Correct double positivity, deal with disallowed marker pairs,
#and assign a "cell type" based on marker positivity
rule call_cell_types:
  input:
    script="src/spatial/call_cell_types.R",
    objects = "data/vectra/interim/objects/{t_number}_{panel}_{batch}.nc",
    lib="src/spatial/read_cells.R",
  output:
    "data/vectra/processed/objects/{t_number}_{panel}_{batch}.Rds",
  shell:
    "Rscript {input.script} {input.objects} {wildcards.panel} {output}"

#Quantify cell types and the marker combinations they represent
rule report_cell_types:
  input:
    cells = [f"data/vectra/processed/objects/{s.sample_id}_{s.panel}_{s.batch_HALO}.Rds"
                    for s in vectra_samples],
    rmd="reports/spatial/04_report_cell_types.Rmd",
    script="src/rmarkdown.R"
  output:
    html="reports/spatial/04_report_cell_types.html"
  shell:
    "Rscript {input.script} {input.rmd} $(realpath -s {output.html})"
    
#Model overall marker combination density by panel
#Use aggregrate_densities instead of invoking this rule directly
rule cell_type_density:
  input:
    script="src/spatial/model_density.R",
    objects="data/vectra/processed/objects/{t_number}_{panel}_{batch}.Rds",
  output:
    tsv="results/spatial/density/{t_number}_{panel}_{batch}.tsv",
  shell:
    "mkdir -p results/spatial/density\n"
    "Rscript {input.script} {input.objects} {output}"

#combine densities for reporting
rule aggregrate_densities:
  input:
    tsv=[f"results/spatial/density/{s.sample_id}_{s.panel}_{s.batch_HALO}.tsv"
     for s in vectra_samples],
    script="src/spatial/cat_tsv.R",
  output: 'results/spatial/density.tsv'
  shell:
    "Rscript {input.script} {input.tsv} {output}"

#Histograms, scatterplots and heatmaps of cell types by tumor/stroma and panel
rule density_report:
  input:
    "results/spatial/density.tsv",
    rmd="reports/spatial/05_density.Rmd",
    script="src/rmarkdown.R",
  output:
    html="reports/spatial/05_density.html"
  shell:
    "Rscript {input.script} {input.rmd} $(realpath -s {output.html})"

rule ppbc_density:
  input:
    densities = "results/spatial/density.tsv",
    meta="data/metadata/PPBC_metadata.xlsx",
    rmd="reports/spatial/06_ppbc_density.Rmd",
    script="src/rmarkdown.R"
  output:
    html="reports/spatial/06_ppbc_density.html"
  shell:
    "Rscript {input.script} {input.rmd} $(realpath -s {output.html})"