library(here)
library(tidyverse)

# Prepare cell densities for survival analysis by linking them to PPBC subtypes and clinical outcomes.

source(here("src/utils/parse_args.R"))
argdf <- retrieve_args()

print(argdf)

## Load data

density <- read_tsv(here(filter(argdf, argname == "densities")$argval),
                    show_col_types = F)
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

meta <- read_sheets_to_list(here(filter(argdf, argname == "meta")$argval))

#Here the t_number can be mapped to the patient_ID (via sample_ID).
#head(meta$sample_data)

#Survival data can be found here.
#head(meta$patient_data)

#Explanations are in the codebook.
#head(meta$codebook)

## Survival outcomes in PPBC

# Join t_numbers with survival, link to density.

dens <- meta$sample_data %>%
  filter(sample_type == "slide" & Included == 1) %>%
  select(t_number = sample_ID, patient_ID) %>% distinct() %>%
  left_join(., meta$patient_data, by = "patient_ID") %>%
  left_join(density, ., by = "t_number")

stopifnot(nrow(filter(dens, is.na(death)))==0)

unique(dens$study_group)

if(nrow(filter(dens, is.na(grade))) > 0){
  print("Grade is sometimes missing")
  dens %>% group_by(study_group, database, grade) %>%
    count() %>% print(n=Inf)
}

if(nrow(filter(dens, is.na(stage))) > 0){
  print("Stage is sometimes missing")
  dens %>% group_by(study_group, database, stage) %>%
    count() %>% print(n=Inf)
}

#Format PPBC and survival categories.
dens <- dens %>%
  mutate(classifier_label = factor(classifier_label,
                                   levels = c("Stroma", "Tumor", "Total")),
         panel = factor(panel, levels = c("MPIF26", "MPIF27"))
  ) %>%
  mutate(study_group = factor(study_group, levels = c("npbc", "prbc",
                                                      "ppbcdl", "ppbcpw"))) %>%
  mutate(PPBC = factor(PPBC, levels = c("nulliparous", "pregnant",
                                               "lactating", "involuting"))) %>%
  mutate(stage = factor(stage, levels = c("stage I", "stage II",
                                          "stage III", "stage IV"))) %>%
  mutate(grade = factor(grade, levels = c("grade I", "grade II", "grade III")))

dens %>%
  select(study_group, PPBC, stage, grade) %>%
  distinct() %>% print(n=Inf)

## Write data
saveRDS(dens, here("data/vectra/processed/density_ppbc.Rds"))

sessionInfo()
