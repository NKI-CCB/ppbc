#' Read patient data
#' 
#' @description Read the patient data with enforcement of column types. 
#' Patient data can subsequently be merged with the relevant sample data. Some 
#' columns reflect manual calculations in Excel. These are excluded during the
#' read-in step, and will be recalculated later.
#'
#' @param path Character, path to patient data tsv file.
#' @param verbose Logical, whether to print information about column processing
#'
#' @return A tibble containing processed patient data
#' @export
#'
#' @examples read_patient_data(here("data/external/patient_data.tsv"))
read_patient_data <- function(path, verbose = T){
  
  pd <- read_tsv(path, show_col_types = F, na = c("", "NA"))
  
  stopifnot(pd$Ref == pd$patient_ID)
  colnames_before <- colnames(pd)
  
  # parse_factor is similar to factor(), but will generate warnings if 
  # elements of x are not found in levels.
  if(verbose){
    print("Recoding no_ppbc in months_involution and months_breastfeeding to NA")
    print("Renaming TILprecent to TILpercent")
    print("Omitting manual Excel calculations")
  }
  #return(pd)
  pd <- pd %>%
    dplyr::transmute( # Transmute discards columns not mentioned
      patient_ID, sample_name,
      study_group = readr::parse_factor(
        study_group, levels = c("npbc", "prbc", "ppbcdl", "ppbcpw"), include_na=F
      ),
      PPBC = readr::parse_factor(
        PPBC, levels = c("nulliparous", "pregnant", "lactating", "involuting"), include_na=F
      ),
      clin_subtype = readr::parse_factor(
        clin_subtype, levels = c("TripleNegative", "LumA", "Luminal B (HER2neg)",
                                 "Luminal B (HER2pos)", "HER2pos (non-luminal)"), include_na=F
      ),
      ER = as.integer(ER),
      PR = as.integer(PR),
      HER2 = as.integer(HER2),
      time_OS_months = as.double(time_OS_months),
      death = as.integer(death),
      time_DRS_months = as.double(time_DRS_months),
      distant_recurrence = as.integer(distant_recurrence),
      stage = readr::parse_factor(
        stage, levels = c("stage I", "stage II", "stage III", "stage IV"), include_na=F
      ),
      grade = readr::parse_factor(
        grade, levels = c("grade I", "grade II", "grade III"), include_na=F
      ),
      country = readr::parse_factor(
        country, levels = c("BE", "IT", "GER", "CA", "NL"), include_na=F
      ),
      # Remove no_ppbc to parse as dbl
      months_involution = ifelse(months_involution == "no_ppbc",
                                 NA, months_involution),
      months_involution = as.double(months_involution),
      months_breastfeeding = ifelse(months_breastfeeding == "no_ppbc",
                                    NA, months_breastfeeding),
      months_breastfeeding = as.double(months_breastfeeding),
      surgery = as.integer(surgery),
      BE = as.integer(BE),
      ME = as.integer(ME),
      SN = as.integer(SN),
      OE = as.integer(OE),
      RT = as.integer(RT),
      HT = as.integer(HT),
      HT_type = readr::parse_factor(
        HT_type, levels = c(NA, "no_HT", "adj_HT", "neoadj_HT", "combo_HT"), include_na=F
      ),
      CT = as.integer(CT),
      CT_type = readr::parse_factor(
        CT_type, levels = c(NA, "adj_CT", "no_CT", "neoadj_CT", "combo_CT"), include_na=F
      ),
      herceptin = as.integer(herceptin),
      herceptin_type = readr::parse_factor(
        herceptin_type, levels = c(NA, "no_herceptin", "adj_herceptin",
                                   "neoadj_herceptin",  "combo_herceptin"), include_na=F
      ),
      database = readr::parse_factor(database, levels = c("MBC", "INCIP"), include_na=F),
      year_birth = as.double(year_birth),
      year_diagnosis = as.double(year_diagnosis),
      FU_year = as.double(FU_year),
      FU_time_months = as.double(FU_time_months),
      # Stromal CD38 scores according to pathologist
      # Is reportedly a percentage (dbl check requested)
      SI_CD38_total = as.double(SI_CD38_total),
      # SI_CD38_total is high if above or equal to the median, otherwise low
      # was calculated without a formula in Excel, recalculate downstream
      # SI_CD38,
      # Intratumoral CD38 scores according to pathologist
      TI_CD38_total = as.double(TI_CD38_total),
      # TI_CD38_total is high if above or equal to the median, otherwise low
      # also redo this manual Excel calculation downstream
      TILpercent = as.double(TILprecent),
      # Manual calculation of TIL tertiles, redo downstream
      # TIL,
      # Manual recoding of medium or high TIL teritles as "high"
      # Redo downstream
      # TIL_short,
      # Pathologist's TAPC score of 0-3
      TAPC_score,
      # Manual recoding of TAPC 0-1 as "low" and 2-3 as "high"
      # TAPC,
      reason_death = readr::parse_factor(reason_death, levels = c("alive", "BC"), include_na=F)
    )
  
  if(verbose){
    print(paste("Columns kept:", paste(colnames(pd), collapse = ", ")))
    print(paste("Columns discarded:",
                paste(setdiff(colnames_before, colnames(pd)),
                      collapse = ", ")))
  }
  pd
}