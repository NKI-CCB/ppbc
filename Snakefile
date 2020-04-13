#Snakemake integration post alignment with Salmon

rule process_metadata:
    input:
        raw_metadata="data/external/Hercoderingslijst_v09032020_KM.xlsx"
    output:
        report="reports/01_process_metadata_tximport.html",
        multiple_patient_fastqs="data/metadata/01_patients_with_multiple_fastqs.csv",
        excluded_samples="data/metadata/01_pre_excluded_samples.csv",
        sample_annotation="data/metadata/01_sample_annot.tsv",
        gene_annotation="data/metadata/01_tx_annot.tsv",
        tx="data/Rds/01_tx.Rds",
        rdata="reports/01_process_metadata_tximport.RData"
    script:
        "reports/01_process_metadata_tximport.Rmd"