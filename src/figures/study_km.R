#' Kaplan meier for study group and covariate
#'
#' @description Will return a ggplot KM facetted around the covariate in question.
#' Filtering is performed for time > 20 months and death due to breast cancer.
#'
#' @param df A data frame
#' @param covariate A column in df that contains a factor covariate for facetting
#' @param surv_formula The survival formula to fit
#' @param event_type OS or DRS (overall survival vs distant recurrence)
#' @param study_groups A PPBC study group
#' @param returndata Logical, whether to return the df used for plotting
#' @param returnfit 
#'
#' @return A ggplot
#' @export
#'
#' @exampleskm_TAPC_OS <- km(figuredata, covariate = "TAPC", event_type = "OS",
#' surv_formula=as.formula("Surv(time=time_OS_months, event=death) ~ study_group + TAPC"))
km <- function(df, covariate, surv_formula, event_type = c("OS", "DRS"),
               study_groups = c("PP-BCpw", "Pr-BC", "NP-BC", "PP-BCdl"),
               returndata = F, returnfit = F){
  
  df <- df[!is.na(df[,covariate]),]
  
  
  df <- df %>% dplyr::filter(study_group %in% {{study_groups}}) %>%
    droplevels()
  
  # For debugging or manual exploration of relevant columns
  if(returndata){
    return(df[,c(covariate, "study_group", "sample_name", "time_OS_months", "death",
                 "reason_death", "time_DRS_months", "distant_recurrence")])
  }
  
  # Perform filtering depending on time variable
  if(event_type == "OS"){
    df <- df %>% filter(!is.na(death), !is.na(time_OS_months)) %>%
      filter(time_OS_months > 20 | reason_death == "BC") %>%
      as.data.frame()
    y <- "Overall survival"
  } else {
    df <- df %>% filter(!is.na(distant_recurrence), !is.na(time_DRS_months)) %>%
      filter(time_DRS_months > 20 | reason_death == "BC") %>%
      as.data.frame()
    y <- "Distant recurrence"
  }
  
  # This is a weird workaround to avoid "object of type 'symbol' is not subsettable" errors
  # Which would occur if you used fit <- survfit(surv_formula, data = df)
  # See https://github.com/kassambara/survminer/issues/252
  fit <- do.call(survfit, args = list(formula = surv_formula, data = df))
  if(returnfit){return(fit)}
  
  ggsurvplot(
    fit, data = df,
    xlab = "Months", 
    ylab = paste(y,"probability"),
    title = paste(y, covariate),
    palette = study_colors[names(study_colors) %in% unique(df$study_group)],
    ggtheme = theme_bw(),
    pval = T,
    facet.by = covariate,
    legend.title="Study group:"
  )
}