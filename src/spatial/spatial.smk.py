from dataclasses import dataclass

import openpyxl

# The spatial analysis gets its own subdirectory in ./reports and within the data (sub)directories
# Raw data is in data/vectra

wildcard_constraints:
    # T numbers, with a possible block ID appended
    t_number="T[0-9]+-[0-9]+(_[I0-9]+)?",
    panel="MPIF[0-9]+",
    batch="batch[0-9]+",
    tissue="stroma|tumor"


###################
# Sample metadata #
###################

# This is unique to the spatial analyses for now, but could be generalized
# to all data types and moved to the main Snakefile.

@dataclass(frozen=True)
class Sample():
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
        elif x["sample_type"] == "RNA": 
            return RNASample(
                fastq = x["sample_ID"],
                patient_id = x["patient_ID"],
                included = x["Included"] == 1,
            )
        else:
            # Should probably not refer to CD38 staining as "HE", request change to metadata
            assert x["sample_type"] in ["HE"]
            return HESample(
                patient_id = x["patient_ID"],
                included = x["Included"] == 1,
                scoring = x["experimental_platform"],
            )

@dataclass(frozen=True)
class VectraSample(Sample):
    sample_id: str
    batch_HALO: str
    panel: str
    
@dataclass(frozen=True)
class RNASample(Sample):
    fastq: str

@dataclass(frozen=True)
class HESample(Sample):
    scoring: str    

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
    if len(samples) != len(set(samples)): # Set allows no duplicates
        print("Total samples read:" + str(len(samples)))
        print("Unique samples read" + len(set(samples)))
        import collections
        c = {s: n for s, n in collections.Counter(samples).items() if n > 1}
        print(c)
        #raise Exception('Duplicated samples')
    return samples

samples = read_sample_xlsx("data/external/PPBC_metadata_20220811.xlsx")

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

# Note: sometimes Snakemake tries to write too many files at the same time
# Try rerunning rules that fail the first time, they often succeed the second time

# Create a manageable data structure
rule organize_vectra:
  input:
    # Exact date and time of Vectra staining
    batch_info = "data/external/MPIF26en27 batches gekleurd.xlsx", 
    example_halo_archive = 
        "data/vectra/halo/Batch 1/MPIF26/T21-60304/Halo archive 2021-06-21 12-17 - v3.2.1851/",
    metadata = "data/external/PPBC_metadata_20220811.xlsx",
    rmd="src/spatial/organize_vectra_samples.Rmd",
    script="src/utils/rmarkdown.R"
  output: 
    # Snakemake complains if not all output files and logs have the same wildcards
    # For simplicity, track sample-specific files, aggregate files as shell text
    xmls = "data/vectra/raw/annotations/{t_number}_{panel}_{batch}_annotations.xml",
    objects = "data/vectra/raw/objects/{t_number}_{panel}_{batch}_object_results.csv",
    summary = "data/vectra/raw/summary/{t_number}_{panel}_{batch}_summary_results.csv"
  shell:
    "Rscript {input.script} {input.rmd} $PWD'/src/spatial/organize_vectra_samples.html'"
    " --vectra_dir data/vectra  --metadata '{input.metadata}' "
    " --batch_info '{input.batch_info}'"
    " --example_halo_archive '{input.example_halo_archive}'"
    " --out '{output}'"
    " --filedict 'data/vectra/metadata/file_location_dictionary.csv'"

#QC checks based on aggregated summary files
rule summary_QC:
  input:
    summaries=[f"data/vectra/raw/summary/{s.sample_id}_{s.panel}_{s.batch_HALO}_summary_results.csv"
                    for s in vectra_samples],
    # produced by organize_vectra but can't be tracked due to wildcards                
    # filedict="data/vectra/metadata/file_location_dictionary.csv",
    rmd="reports/spatial/01_summary_QC.Rmd",
    script="src/utils/rmarkdown.R"
  params:
    dir="data/vectra/raw/summary"
  output:
    html="reports/spatial/01_summary_QC.html"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"
    " --dir '{params.dir}' "

rule load_annotations:
    input:
        script = "src/spatial/import_halo_annotation_wkb.py",
        xml = "data/vectra/raw/annotations/{t_number}_{panel}_{batch}_annotations.xml",
    output:
        "data/vectra/interim/annotations/{t_number}_{panel}_{batch}_tumor.wkb"
    shell:
        "python3 {input.script} {input.xml} '.*[tT]umor.*' {output}"

#Convert object results and summary to .nc format
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
    

#QC reports based on object files rather than summary files
#Generate batchwise output by adding this to rule all:
#expand("reports/spatial/object_qc_by_batch/02_object_QC_{batch}.html", 
# batch = ['batch' + str(i) for i in range (1, 8)])
rule batch_object_QC:
  input:
    object_files = [f"data/vectra/interim/objects/{s.sample_id}_{s.panel}_{s.batch_HALO}.nc"
                    for s in vectra_samples],
    rmd = "reports/spatial/02_object_QC.Rmd",
    script="src/utils/rmarkdown.R"
  output:
    html="reports/spatial/object_qc_by_batch/02_object_QC_{batch}.html"
  shell:
    "mkdir -p reports/spatial/object_qc_by_batch\n"
    "Rscript {input.script} {input.rmd} $(realpath -s {output.html})"
    " --batch {wildcards.batch}"


#####################
# Cell Type Calling #
#####################
    
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
    script="src/utils/rmarkdown.R"
  output: 
    marker_region="results/spatial/marker_combos_by_region.csv",
    html = "reports/spatial/03_aggregate_marker_combos.html"
  shell: 
    "Rscript {input.script} {input.rmd} $(realpath -s {output.html})"

#Visualize problematic marker combinations by batch
rule batch_marker_coexpression:
  input:
    counts = [f"results/spatial/marker_cell_counts/{s.sample_id}_{s.panel}_{s.batch_HALO}.csv"
                    for s in vectra_samples],
    rmd="reports/spatial/04_batch_marker_coexpression.Rmd"
  output:
    html="reports/spatial/batch_marker_viz/04_{batch}_marker_coexpression.html",
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
    rmd="reports/spatial/05_report_cell_types.Rmd",
    script="src/utils/rmarkdown.R"
  output:
    html="reports/spatial/05_report_cell_types.html"
  shell:
    "Rscript {input.script} {input.rmd} $(realpath -s {output.html})"

#######################
# Tissue Segmentation #
#######################
 
# Decided against tumor vs stroma tissue segmentation based on
# the results of the this analysis. Remove?

rule segment_tissue_density:
  input:
    cells = "data/vectra/processed/objects/{sample}_{panel}_{batch}.Rds",
    config = "src/spatial/segment_tissue_density_config.yaml",
    script = "src/spatial/segment_tissue_density.R",
  output: "data/vectra/processed/segmentation/{sample}_{panel}_{batch}.Rds"
  shell:
    "mkdir -p data/vectra/processed/segmentation/\n"
    "Rscript {input.script} {input.cells} {input.config} {output}"

rule plot_segmentations:
  input:
    segmentations = [
        f"data/vectra/processed/segmentation/{s.sample_id}_{s.panel}_{s.batch_HALO}.Rds"
        for s in vectra_samples],
    script = "src/spatial/plot_segmentation.R",
  output: "figures/spatial/segmentation.pdf"
  shell:
    "mkdir -p figures/spatial/\n"
    "Rscript {input.script} {input.segmentations} {output}"

rule report_segmentation:
  input:
    "data/vectra/processed/objects/T21-60303_MPIF26_batch2.Rds",
    "data/vectra/processed/objects/T21-60303_MPIF27_batch2.Rds",
    "src/spatial/segment_tissue_density_config.yaml",
    rmd="reports/spatial/11_tissue_segmentation.Rmd",
    script="src/utils/rmarkdown.R"
  output:
    html="reports/spatial/11_tissue_segmentation.html"
  shell:
    "Rscript {input.script} {input.rmd} $(realpath -s {output.html})"

############################################
# Spatial density of cell types in tissues #
############################################
    
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
    rmd="reports/spatial/06_density.Rmd",
    script="src/utils/rmarkdown.R",
  output:
    html="reports/spatial/06_density.html"
  shell:
    "Rscript {input.script} {input.rmd} $(realpath -s {output.html})"
    
#Link PPBC study group and clinical outcomes to cell densities
rule process_density:
  input:
    densities = "results/spatial/density.tsv",
    meta="data/external/PPBC_metadata_20220811.xlsx",
    script="src/spatial/process_density_outcomes.R",
    lib="src/utils/parse_args.R"
  output:
    density_outcome = "data/vectra/processed/density_ppbc.Rds"
  shell:
    "Rscript {input.script} --densities {input.densities}"
    " --meta {input.meta}"

#Kruskal wallis tests and beehive plots for PPBC association with cell density
rule kruskal_density:
  input:
    density_outcome = "data/vectra/processed/density_ppbc.Rds",
    rmd = "reports/spatial/07_kruskal_density.Rmd",
    script = "src/utils/rmarkdown.R",
    lib = "src/spatial/kruskal_density.R" 
  params:
    min_cell_count = 20000
  output:
    html = "reports/spatial/07_kruskal_density.html"
  shell:
    "Rscript {input.script} {input.rmd} $(realpath -s {output.html})"
    " --min_cell_count '{params.min_cell_count}'"
    " --density_outcome '{input.density_outcome}'"

#Relationship between involution and breastfeeding duration and immune density
rule invbf_time_density:
  input:
    density_outcome = "data/vectra/processed/density_ppbc.Rds",
    rmd="reports/spatial/08_inv_time_density.Rmd",
    script="src/utils/rmarkdown.R"
  params:
    min_cell_count = 20000
  output:
    html="reports/spatial/08_inv_time_density.html"
  shell:
    "Rscript {input.script} {input.rmd} $(realpath -s {output.html})"
    " --min_cell_count '{params.min_cell_count}'"
    " --density_outcome '{input.density_outcome}'"
    
#Cox regressions for cell densities, OS and DRS
rule cox_density:
  input:
    density_outcome = "data/vectra/processed/density_ppbc.Rds",
    rmd="reports/spatial/09_cox_density.Rmd",
    script="src/utils/rmarkdown.R",
    lib="src/spatial/cox_spatial.R"
  params:
    min_cell_count = 20000
  output:
    html="reports/spatial/09_cox_total_density_{outcome}.html"
  shell:
    "Rscript {input.script} {input.rmd} $(realpath -s {output.html})"
    " --min_cell_count '{params.min_cell_count}'"
    " --density_outcome '{input.density_outcome}'"
    
#Check whether CD20 intensity in Vectra is correlated with CD20 RNAseq expression    
rule cd20_clusters:
  input:
    density_outcome = "data/vectra/processed/density_ppbc.Rds",
    inv_clusters = "results/rnaseq/clustering/11_inv_clusters.xlsx",
    rmd="reports/spatial/10_ig_clusters_cd20.Rmd",
    script="src/utils/rmarkdown.R"
  output:
    html="reports/spatial/10_ig_clusters_cd20.html"
  shell:
    "Rscript {input.script} {input.rmd} $(realpath -s {output.html})"
    " --density_outcome '{input.density_outcome}'"
    " --inv_clusters '{input.inv_clusters}'"

#################################
# Second order spatial measures #
#################################

rule model_l:
  "L measure of the immune cell types"
  input:
    script="src/spatial/model_l.R", 
    objects="data/vectra/processed/objects/{t_number}_{panel}_{batch}.Rds",
    annotation="data/vectra/interim/annotations/{t_number}_{panel}_{batch}_tumor.wkb",
  output:
    tsv="results/spatial/l/{t_number}_{panel}_{batch}.tsv",
  shell:
    "mkdir -p results/spatial/l\n"
    "Rscript {input.script} {input.objects} {input.annotation} {output} F"


rule model_lcross_panck:
  "L cross measure between panCK and immune cell types"
  input:
    script="src/spatial/model_lcross_panck.R", 
    objects="data/vectra/processed/objects/{t_number}_{panel}_{batch}.Rds",
    annotation="data/vectra/interim/annotations/{t_number}_{panel}_{batch}_tumor.wkb",
  output:
    tsv="results/spatial/lcross_panck/{t_number}_{panel}_{batch}.tsv",
  shell:
    "mkdir -p results/spatial/lcross_panck\n"
    "Rscript {input.script} {input.objects} {input.annotation} {output} F"

rule model_lcross_immune:
  "L cross measure between immune cell types"
  input:
    script="src/spatial/model_lcross_immune_cells.R",
    objects="data/vectra/processed/objects/{t_number}_{panel}_{batch}.Rds",
    annotation="data/vectra/interim/annotations/{t_number}_{panel}_{batch}_tumor.wkb",
  output:
    tsv="results/spatial/lcross_immune/{t_number}_{panel}_{batch}.tsv",
  shell:
    "mkdir -p results/spatial/lcross_immune\n"
    "Rscript {input.script} {input.objects} {input.annotation} {output} F"

rule report_spatstat:
  input:
    [f"results/spatial/l/{s.sample_id}_{s.panel}_{s.batch_HALO}.tsv"
         for s in vectra_samples],
    [f"results/spatial/lcross_panck/{s.sample_id}_{s.panel}_{s.batch_HALO}.tsv"
         for s in vectra_samples],
    [f"results/spatial/lcross_immune/{s.sample_id}_{s.panel}_{s.batch_HALO}.tsv"
         for s in vectra_samples],
    lib="src/spatial/outcome.R",
    rmd="reports/spatial/12_spatstat_overview.Rmd",
    metadata="data/external/PPBC_metadata_20220811.xlsx", 
    density="data/vectra/processed/density_ppbc.Rds",
    script="src/utils/rmarkdown.R"
  output:
    html="reports/spatial/12_spatstat_overview.html",
    spatstat="data/vectra/processed/11_spatstat.Rds",
    spatstat_l="data/vectra/processed/11_spatstat_l.Rds",
    spatstat_lcross_immune="data/vectra/processed/11_spatstat_lcross_immune.Rds",
    spatstat_lcross_panck="data/vectra/processed/11_spatstat_lcross_panck.Rds"
  shell:
    "Rscript {input.script} {input.rmd} $(realpath -s {output.html})"
    " --metadata {input.metadata}"
    " --density {input.density}"
    " --lib {input.lib}"
