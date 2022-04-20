library(here)
library(tidyverse)

#Prepare cell densities for survival analysis by linking them to PPBC subtypes and clinical outcomes.

## Load data

density <- read_tsv(here("results/spatial/density.tsv"), show_col_types = F)
#head(density)

density %>%
  group_by(panel) %>%
  count()

#The metadata is separated into three sheets.
read_sheets_to_list <- function(xlsfile){
  readxl::excel_sheets(xlsfile) %>% 
    set_names() %>% 
    map(readxl::read_excel, path = xlsfile)
}

meta <- read_sheets_to_list(here("data/metadata/PPBC_metadata.xlsx"))

#Here the t_number can be mapped to the patient_ID (via sample_ID).
#head(meta$sample_data)

#Survival data can be found here.
#head(meta$patient_data)

#Explanations are in the codebook.
#head(meta$codebook)


## Survival outcomes in PPBC

# Join t_numbers with survival, link to density.

dens <- meta$sample_data %>%
  filter(experimental_platform != "RNAseq" & !is.na(sample_ID)) %>%
  select(t_number = sample_ID, patient_ID) %>% distinct() %>%
  left_join(., meta$patient_data, by = "patient_ID") %>%
  left_join(density, ., by = "t_number")

stopifnot(nrow(filter(dens, is.na(death)))==0)

#head(dens)

#Four possible possible PPBC categories exist : 
  
#* nulliparous breast cancer(non_prbc)
#* pregnant breast cancer (prbc)
#* post-partum breast cancer: lactating (ppbc_lac)
#* post-partum breast cancer: involuting (ppbc_inv)

unique(dens$study_group)

#Format PPBC and survival categories.
#Grade is sometimes missing
#dens %>% group_by(study_group, database, grade) %>% count()

dens <- dens %>%
  mutate(classifier_label = factor(classifier_label,
                                   levels = c("Stroma", "Tumor", "Total")),
         panel = factor(panel, levels = c("MPIF26", "MPIF27")),
         #Shorter group names, similar to study_group but clearer 
         #and with no potential delimiters
         group = case_when(
           study_group == "ppbc_inv" ~ "inv",
           study_group == "non_prbc" ~ "nonprbc",
           study_group == "ppbc_lac" ~ "lac",
           TRUE ~ study_group) 
  ) %>%
  mutate(group = factor(group, levels = c("nonprbc", "prbc", "lac", "inv"))) %>%
  relocate(group, .after = study_group) %>%
  mutate(stage = factor(stage, levels = c("stage I", "stage II",
                                          "stage III", "stage IV"))) %>%
  mutate(grade = factor(grade, levels = c("grade I", "grade II", "grade III")))

dens %>%
  select(study_group, group) %>%
  distinct()

## Write data
saveRDS(dens, here("data/vectra/processed/density_ppbc.Rds"))

sessionInfo()
