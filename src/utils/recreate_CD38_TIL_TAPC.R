# Recreate the TAPC, TIL, and CD38 median/tertile categories that were 
# originally calculated manually in Excel

recreate_HE_categories <- function(df, verbose = T){
  print(paste("Median SI_CD38:", median(pd$SI_CD38_total, na.rm = T)))
  print(paste("Median TI_CD38:", median(pd$TI_CD38_total, na.rm = T)))
  
  df <- df %>% 
    mutate(SI_CD38 = ifelse(SI_CD38_total >= median(SI_CD38_total, na.rm=T),
                     "high", "low"), .after = SI_CD38_total) %>%
    mutate(TI_CD38 = ifelse(TI_CD38_total >= median(TI_CD38_total, na.rm=T),
                            "high", "low"), .after = TI_CD38_total) %>%
    mutate(TIL = dplyr::ntile(TILprecent, 3)) %>%
    mutate(TIL = case_when(
      TIL == 1 ~ "low",
      TIL == 2 ~ "medium",
      TIL == 3 ~ "high"
    ))
  
  
}