# Cox regression-related functions for linking cell spatial density with clinical outcome

library(survival)
library(dplyr)
library(broom)

#' Tidy Cox results
#'
#' @description Performs both an LRT between the full and reduced Cox models and
#' extracts estimates and HR from the coefficient of interest, then tidies results
#' 
#' @param reduced_fit A coxph fit of the reduced formula for anova,
#'  i.e. Surv(time, event) ~ 1
#' @param full_fit A coxph fit of the full formula for anova, i.e.
#' Surv(time, event) ~ 1 + Total_density
#' @param cell_type Character, the cell density to be modelled
#' @param panel Character, the panel in which that marker is found.
#'
#' @return A tibble containing the summary information of the anova between
#' the full and reduced fit coxph objects
#' @export
#'
#' @examples 
#' tidy_cox(reduced_fit = coxph("Surv(time, event) ~ 1", data = df),
#' full_fit = coxph("Surv(time, event) ~ 1 + total_density", data = df),
#' cell_type = "CD20+", panel = "MPIF26")
tidy_cox <- function(reduced_fit, full_fit, cell_type, panel){
  lrt <- anova(reduced_fit, full_fit) %>% broom::tidy()
  full_fit %>%
    broom::tidy(conf.int = T) %>%
    mutate(HR = exp(estimate),.after=estimate) %>%
    mutate(cell_type = {{cell_type}},
           panel = {{panel}},
           .after=term) %>%
    mutate(lrt.pval = lrt$p.value[[2]],
           lrt.logLik = lrt$logLik[[2]],
           lrt.statistic =  lrt$statistic[[2]],
           .after = panel) %>%
    rename(coef.estimate = estimate,
           coef.statistic = statistic,
           coef.pval = p.value)
}

#' Perform Cox regressions on on cell type density
#' 
#' @description Performs an anova between coxph fit produced by formula `f` and
#' `f + Total_density` of a given cell type and panel 
#'
#' @param df A data frame containing outcome columns matching formula f
#' @param panel Character, the IF panel
#' @param cell_type Character, the cell type density to model
#' @param f Character, the reduced Surv formula against which to test a full formula
#' which includes a "Total_density" term
#' @param tidied Logical, whether to return a tibble via tidy_cox(), or print 
#' the coxph summary directly
#'
#' @return A data frame or printed coxph summary
#' @export
#'
#' @examples cox_density(cox_dens, "MPIF27", "CD3+CD8+", f = univ_form, tidied = T)
cox_density <- function(df, panel, cell_type, f, tidied = T){
  
  #Subset by panel and celltype
  df <- df %>% filter(panel == {{panel}}, cell_type == {{cell_type}}) %>%
    rename(Total_density = density)
  
  #If t_numbers are duplicated, something is wrong
  stopifnot(nrow(filter(df, duplicated(t_number)))==0)
  
  # manual check that the tidy output is what it should be
  if(!tidied){
    reduced_fit = coxph(as.formula(f), data = as.data.frame(df))
    full_fit = coxph(as.formula(paste(f, "Total_density", sep = "+")),
                     data = as.data.frame(df))
    return(
      list(full = coxph(as.formula(paste(f, "Total_density", sep = "+")),
                        data = as.data.frame(df)),
           anova = anova(reduced_fit, full_fit)
      ))
  }
  
  #Perform Cox regression and LRT for tumor and stroma
  res <- tidy_cox(
    reduced_fit = coxph(as.formula(f), data = as.data.frame(df)),
    full_fit = coxph(as.formula(paste(f, "Total_density", sep = "+")),
                     data = as.data.frame(df)),
    cell_type, panel
  ) %>%
    mutate(term = str_replace(term, "stage", "tumor "))
  
  
  res
}
