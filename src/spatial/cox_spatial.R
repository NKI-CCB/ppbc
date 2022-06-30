#### Cox regression-related functions for linking cell spatial density with clinical outcome ----

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
  
  # allow manual check that the tidy output is what it should be
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

#' Cox regressions for all cell types while auto-detecting the panel
#'
#' @description This is a wrapper function for cox_density that auto-detects in which 
#' panel the provided cell type can be found, then runs cox_density for that panel.
#' If the cell type can be found in multiple panels, the panels are run separately.
#'
#' @param df A data frame containing outcome columns matching formula f
#' @param cell_type Character, the cell type density to model
#' @param f Character, the reduced Surv formula against which to test a full formula
#' which includes a "Total_density" term. See `cox_density()` for details.
#'
#' @return A data frame containing the results of the Cox regression, separated by panel
#' @export
#'
#' @examples lapply(relevant_cells, function(x){
#' all_cox_density(cox_dens, x,
#'                 f = "Surv(time = time_OS_months, event = death) ~ 1")
#'                 })
all_cox_density <- function(df, cell_type, f){
  
  df <- df %>% filter(cell_type == {{cell_type}})
  
  mpif26_cells <- df %>% filter(panel == "MPIF26") %>% pull(cell_type) %>%
    unique()
  
  mpif27_cells <- df %>% filter(panel == "MPIF27") %>% pull(cell_type) %>%
    unique()
  
  #Skip regression if cell type is absent in that panel
  if(cell_type %in% mpif26_cells & cell_type %in% mpif27_cells){
    res <- bind_rows(
      cox_density(df, "MPIF26", cell_type, f),
      cox_density(df, "MPIF27", cell_type, f)
    )
  } else if (cell_type %in% mpif26_cells) {
    res <- cox_density(df, "MPIF26", cell_type, f)
  } else if (cell_type %in% mpif27_cells){
    res <- cox_density(df, "MPIF27", cell_type, f)
  } else {
    stop("Provide a valid cell type")
  }
  
  #Remove non-density entries from output
  res <- res %>% filter(str_detect(term, "density"))
  res
}

#' Kaplan-meier plots for density
#'
#' @param df A data frame containing columns which match formula `f`
#' @param panel  panel Character, the IF panel
#' @param cell_type Character, the cell type density to model
#' @param f Character, the Surv formula containing at least a quantile term. I.e.,
#' "Surv(time = time_OS_months, event = death) ~ quantile"
#' @param outcome Character, either "OS" or "DRS"
#' @param ngroups Integer, the number of ntiles to plot. Must be 2, 3 or 4.
#' @param pal Character, a color palette to be passed to ggsurvplot()
#'
#' @return A ggsurvplot
#' @export
#'
#' @examples km_plot(cox_dens, panel = "MPIF27", cell_type = "CD20+", outcome,
#' f = paste0(survres, " ~ quantile"))
km_plot <- function(df, panel, cell_type, f, outcome, ngroups = 2, pal = "npg"){
  
  stopifnot(ngroups %in% c(2,3,4))
  
  labels = switch(as.character(ngroups),
                  "2" = c("low", "high"),
                  "3" = c("low", "medium", "high"),
                  "4" = c("Q1", "Q2", "Q3", "Q4")
  )
  
  df <- df %>% filter(panel == {{panel}}, cell_type == {{cell_type}}) %>%
    rename(Total_density = density)
  
  #Add density ntiles
  df <- df %>%
    mutate(ntile = dplyr::ntile(Total_density, ngroups), .after=Total_density) %>%
    mutate(quantile = factor(ntile, labels = labels), .after=ntile)
  
  survminer::ggsurvplot(
    fit = survminer::surv_fit(as.formula(f), data = as.data.frame(df)), 
    facet.by = "study_group",
    xlab = "Months", 
    ylab = paste(outcome, "probability"),
    title = paste(ptitle, ": ", cell_type, " total density in ", panel),
    pval = T,
    pval.method = T,
    palette = pal
  )
}

#testing


#### Deprecated regional analyses ----

#' Reshape regions for Cox regressions
#' 
#' @description Separate tumor and stroma densities to prevent unintentional sample duplication
#' @param df 
#'
#' @return
#' @export
#'
#' @examples
cox_reshape <- function(df){
  
  tum_df <- df %>%
    select(-n, -area) %>%
    filter(classifier_label == "Tumor") %>%
    pivot_wider(names_from = classifier_label, values_from = density,
                values_fill = 0, names_glue = "{classifier_label}_density") %>%
    relocate(Tumor_density, .before = everything())
  
  str_df <- df %>%
    select(-n, -area) %>%
    filter(classifier_label == "Stroma") %>%
    pivot_wider(names_from = classifier_label, values_from = density,
                values_fill = 0, names_glue = "{classifier_label}_density") %>%
    relocate(Stroma_density, .before = everything())
  
  return(list(tumor = tum_df, stroma = str_df))
}

cox_regional_density <- function(df, panel, cell_type, f, testing = F){
  
  #Subset by panel and celltype
  df <- df %>% filter(panel == {{panel}}, cell_type == {{cell_type}})
  
  #Reshape the data frame for Cox regression
  tum_df <- cox_reshape(df)$tumor
  str_df <- cox_reshape(df)$stroma
  
  #If t_numbers are duplicated, something is wrong
  stopifnot(nrow(filter(tum_df, duplicated(t_number)))==0 &
              nrow(filter(str_df, duplicated(t_number)))==0)
  
  # for comparing broom output to original call
  # also allows checking effect of performing an anova instead of 
  # taking the relevant covariate p value in a full rank model
  if(testing){
    reduced_fit = coxph(as.formula(f), data = as.data.frame(tum_df))
    full_fit = coxph(as.formula(paste(f, "Tumor_density", sep = "+")),
                     data = as.data.frame(tum_df))
    return(
      list(fulltumor = coxph(as.formula(paste(f, "Tumor_density", sep = "+")),
                             data = as.data.frame(tum_df)),
           anovatumor = anova(reduced_fit, full_fit)
      ))
  }
  
  #Perform Cox regression and LRT for tumor and stroma
  tum <- tidy_cox(
    reduced_fit = coxph(as.formula(f), data = as.data.frame(tum_df)),
    full_fit = coxph(as.formula(paste(f, "Tumor_density", sep = "+")),
                     data = as.data.frame(tum_df)),
    cell_type, panel
  ) %>%
    mutate(term = str_replace(term, "stage", "tumor "))
  
  strom <- tidy_cox(
    reduced_fit = coxph(as.formula(f), data = as.data.frame(str_df)),
    full_fit = coxph(as.formula(paste(f, "Stroma_density", sep = "+")),
                     data = as.data.frame(str_df)),
    cell_type, panel
  ) %>%
    mutate(term = str_replace(term, "stage", "tumor "))
  
  bind_rows(tum, strom)
}

#' Title
#'
#' @param df 
#' @param panel 
#' @param cell_type 
#' @param f 
#' @param outcome 
#' @param ngroups 
#' @param pal 
#'
#' @return
#' @export
#'
#' @examples  km_plot(cox_dens, panel = "MPIF27", cell_type = "CD20+", outcome,f = paste0(survres, " ~ quantile"))
km_regional_plot <- function(df, panel, cell_type, f, outcome, ngroups = 2, pal = "npg"){
  
  stopifnot(ngroups %in% c(2,3,4))
  
  labels = switch(as.character(ngroups),
                  "2" = c("low", "high"),
                  "3" = c("low", "medium", "high"),
                  "4" = c("Q1", "Q2", "Q3", "Q4")
  )
  
  df <- df %>% filter(panel == {{panel}}, cell_type == {{cell_type}})
  
  #Separate tumor and stroma densities to prevent unintentional sample duplication
  tum_df <- cox_reshape(df)$tumor
  str_df <- cox_reshape(df)$stroma
  
  #Add density ntiles
  tum_df <- tum_df %>%
    mutate(ntile = dplyr::ntile(Tumor_density, ngroups), .after=Tumor_density) %>%
    mutate(quantile = factor(ntile, labels = labels), .after=ntile)
  
  str_df <- str_df %>%
    mutate(ntile = dplyr::ntile(Stroma_density, ngroups), .after=Stroma_density) %>%
    mutate(quantile = factor(ntile, labels = labels), .after=ntile)
  
  list(
    tumor = survminer::ggsurvplot(
      fit = survminer::surv_fit(as.formula(f), data = as.data.frame(tum_df)), 
      facet.by = "study_group",
      xlab = "Months", 
      ylab = paste(outcome, "probability"),
      title = paste(ptitle, ": ", cell_type, " tumor density in ", panel),
      pval = T,
      pval.method = T,
      palette = pal
    ),  
    stroma = survminer::ggsurvplot(
      fit = survminer::surv_fit(as.formula(f), data = as.data.frame(str_df)), 
      facet.by = "study_group",
      xlab = "Months", 
      ylab = paste(outcome, "probability"),
      title = paste0(ptitle, ": ", cell_type, " stroma density in ", panel),
      pval = T,
      pval.method = T,
      palette = pal
    )
  )
}

cox_regional_interaction <- function(df, panel, cell_type, baseform, testing=F){
  
  df <- df %>% filter(panel == {{panel}}, cell_type == {{cell_type}})
  
  tum_df <- cox_reshape(df)$tumor
  str_df <- cox_reshape(df)$stroma
  
  stopifnot(nrow(filter(tum_df, duplicated(t_number)))==0 &
              nrow(filter(str_df, duplicated(t_number)))==0)
  f_tumor <- as.formula(paste(baseform, paste("group", "Tumor_density",sep=":"), sep = "+"))
  f_stroma <- as.formula(paste(baseform, paste("group", "Stroma_density",sep=":"), sep = "+"))
  
  if(testing){
    reduced_fit = coxph(as.formula(baseform), data = as.data.frame(tum_df))
    full_fit = coxph(as.formula(f_tumor),
                     data = as.data.frame(tum_df))
    return(
      list(fulltumor = coxph(as.formula(f_tumor),
                             data = as.data.frame(tum_df)),
           anovatumor = anova(reduced_fit, full_fit)
      ))
  }
  
  #Perform Cox regression and LRT for tumor and stroma
  tum <- tidy_cox(
    reduced_fit = coxph(as.formula(baseform), data = as.data.frame(tum_df)),
    full_fit = coxph(as.formula(f_tumor),
                     data = as.data.frame(tum_df)),
    cell_type, panel
  ) %>%
    mutate(term = str_replace(term, "stage", "tumor "))
  
  strom <- tidy_cox(
    reduced_fit = coxph(as.formula(baseform), data = as.data.frame(str_df)),
    full_fit = coxph(as.formula(f_stroma),
                     data = as.data.frame(str_df)),
    cell_type, panel
  ) %>%
    mutate(term = str_replace(term, "stage", "tumor "))
  
  bind_rows(tum, strom) %>%
    mutate(term = str_remove(term, "group"))
}
#testing
cox_interaction(cox_dens, "MPIF27", "CD20+", baseform = multiv_subform, testing = T)
cox_interaction(cox_dens, "MPIF27", "CD20+", baseform = multiv_subform)