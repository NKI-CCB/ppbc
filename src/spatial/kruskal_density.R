
#' Tidy pairwise wilcox tests
#'
#' @description Reshape results of pairwise wilcox tests between groups, so that
# every cell type can be reported in a single row
#'
#' @param pairwise_res Object of class "pairwise.htest"
#'
#' @return A tibble in tidy format
#' @export
#'
#' @examples pairwise.wilcox.test(df$density, df$study_group, 
#' p.adjust.method = "BH", exact = F) %>% tidy_pairwise()
tidy_pairwise <- function(pairwise_res){
  pairwise_res$p.value %>%
    as.data.frame() %>% rownames_to_column("study_group") %>%
    pivot_longer(cols = c(-study_group),names_to = "study_group2", values_to = "wilcox.bh") %>%
    filter(!is.na(wilcox.bh)) %>%
    mutate(pairwise.wilcox = paste(study_group, study_group2, sep = "_vs_")) %>%
    select(-study_group, -study_group2) %>%
    pivot_wider(names_from = pairwise.wilcox, values_from = wilcox.bh,
                names_glue = "{pairwise.wilcox}_bh")
}

#Perform both kruskal wallis and pairwise wilcox tests on each cell type and panel
#for tumor and stroma separately
#' Density statistics between groups
#'
#' @description Perform Kruskal-Wallis and pairwise Wilcoxon tests between PPBC study groups
#'
#' @param df A data frame containing columns `study_group` (character, PPBC classes),
#' `density` (numeric), `cell_type`, and `panel`
#' @param panel Character, the immune panel. Results are calculated separately for
#' cell types present in more than one immune panel.
#' @param cell_type Character, the cell type for which to perform the test
#' @param tidied Logical, whether to return tidy output or pairwise.htest object
#'
#' @return A data frame containing kruskal wallis tests between study_groups
#' @export
#'
#' @examples stat_density(df = dens, panel = "MPIF26", cell_type = "CD3+FoxP3-", tidied = T)
stat_density <- function(df, panel = c("MPIF26", "MPIF27"), cell_type, tidied = T){
  
  stopifnot(cell_type %in% unique(df$cell_type))
  
  df <- filter(df, panel == {{panel}})
  df <- df %>% filter(cell_type == {{cell_type}})
  #return(df)
  
  #Test between all study_group groups together for tumor
  krusk <- tibble(
    cell_type = {{cell_type}},
    panel = {{panel}},
    kruskal_p = kruskal.test(density ~ study_group, data = df)$p.value
  )
  
  pairwise <- pairwise.wilcox.test(df$density, df$study_group,
                                   p.adjust.method = "BH", exact = F)
  
  if(!tidied){return(list(krusk, pairwise))}
  
  pairwise <- tidy_pairwise(pairwise)
  
  res <- cbind(krusk, pairwise)
  res
}

#' Auto-detect appropriate cell types per panel for density statistics
#'
#' @description This is a wrapper functions that allows looping over all cell types
#' without encountering errors for cell types that only exist in one panel
#' 
#' @param df A data frame suitable for `stat_density`, with column "panel"
#' that contains exclusively values "MPIF26" and "MPIF27"
#' @param cell_type Character, the cell type to perform the KW tests and
#' pairwise wilcoxon tests
#'
#' @return A data frame with density statistics
#' @export
#'
#' @examples relevant_cells <- unique(dens$cell_type)
#' relevant_cells <- relevant_cells[!relevant_cells %in% c("PanCK+", "Other")]
#' density_by_ppbc <- lapply(relevant_cells,
#'                           function(x){multi_dens(dens, x)})
multi_dens <- function(df, cell_type){
  
  # separate cells by panel
  mpif26_cells <- df %>% filter(panel == "MPIF26") %>% pull(cell_type) %>%
    unique()
  
  mpif27_cells <- df %>% filter(panel == "MPIF27") %>% pull(cell_type) %>%
    unique()
  
  # Calculate density statistics, skipping panels if the cell type is absent
  if(cell_type %in% mpif26_cells & cell_type %in% mpif27_cells){
    #print("Both panels")
    
    res <- bind_rows(
      stat_density(df, panel = "MPIF26", cell_type = cell_type),
      stat_density(df, panel = "MPIF27", cell_type = cell_type)
    )
    
  } else if (cell_type %in% mpif26_cells) {
    
    res <- stat_density(df, panel = "MPIF26", cell_type = cell_type)
    
  } else if (cell_type %in% mpif27_cells){
    
    res <- stat_density(df, panel = "MPIF27", cell_type = cell_type)
    
  } else {
    stop("Provide a valid cell type")
  }
  res %>%
    mutate(classifier_label = "Total_density", .after = panel)
  
}

#' Beehive plots for density between PPBC study groups
#'
#' @param df A data frame suitable for `stat_density()`, plus additional columns
#' for `colorby`
#' @param panel Character, the immune panel
#' @param cell_type Character, the cell type to be visualized
#' @param colorby Character, the name of the column used for color aesthetics
#' @param debug Logical, returns the data frame being plotted instead of the plot
#'
#' @return A ggplot object
#' @export
#'
#' @examples  plot_density(panel = "MPIF26", cell_type = "CD20+",
#'  colorby = "distant_recurrence")
plot_density <- function(df = dens,
                         panel = c("MPIF26", "MPIF27"),
                         cell_type, colorby, debug=F){
  
  stopifnot(cell_type %in% unique(df$cell_type))
  
  df <- filter(df, panel == {{panel}})
  df <- df %>% filter(cell_type == {{cell_type}})
  df <- df %>% mutate(
    death = factor(death, levels = c(0,1)),
    distant_recurrence = factor(distant_recurrence, levels = c(0,1))
  )
  
  #DRS missing for 3 samples
  df <- df %>%
    filter(!is.na(.data[[colorby]]))
  
  if(debug){return(df)}
  
  comps <- list(
    c("ppbcpw", "npbc")
  )
  
  df %>%
    mutate(classifier_label = ifelse(classifier_label == "Total",
                                     "Total density",
                                     classifier_label)) %>%
    ggplot(aes(x = study_group, y = density)) +
    geom_boxplot(alpha = 0) +
    geom_jitter(aes(color = get(colorby)), height = 0 ) +
    facet_wrap(~classifier_label) +
    labs(color = str_replace_all(colorby, "_", " "),
         y = paste({{cell_type}}, "density"),
         x = "study group") +
    stat_compare_means() +#KW test all samples
    stat_compare_means(comparisons = comps) +#pairwise
    ggtitle(paste({{panel}}, {{cell_type}}, "density in PPBC"))
}

