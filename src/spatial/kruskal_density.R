#Reshape results of pairwise wilcox tests between groups, so that
#every cell type can be reported in a single row
tidy_pairwise <- function(pairwise_res){
  pairwise_res$p.value %>%
    as.data.frame() %>% rownames_to_column("group") %>%
    pivot_longer(cols = c(-group),names_to = "group2", values_to = "wilcox.bh") %>%
    filter(!is.na(wilcox.bh)) %>%
    mutate(pairwise.wilcox = paste(group, group2, sep = "_vs_")) %>%
    select(-group, -group2) %>%
    pivot_wider(names_from = pairwise.wilcox, values_from = wilcox.bh,
                names_glue = "{pairwise.wilcox}_bh")
}

#Perform both kruskal wallis and pairwise wilcox tests on each cell type and panel
#for tumor and stroma separately
stat_density <- function(df, panel = c("MPIF26", "MPIF27"), cell_type, tidied = T){
  
  stopifnot(cell_type %in% unique(df$cell_type))
  
  df <- filter(df, panel == {{panel}})
  df <- df %>% filter(cell_type == {{cell_type}})
  #return(df)
  
  #Test between all study_group groups together for tumor
  krusk <- tibble(
    cell_type = {{cell_type}},
    panel = {{panel}},
    kruskal_p = kruskal.test(density ~ group, data = df)$p.value
  )
  
  pairwise <- pairwise.wilcox.test(df$density, df$group,
                                   p.adjust.method = "BH", exact = F)
  
  if(!tidied){return(list(krusk, pairwise))}
  
  pairwise <- tidy_pairwise(pairwise)
  
  res <- cbind(krusk, pairwise)
  res
}

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
    c("inv", "nonprbc")#,
    #Not really enough samples for these
    #c("inv", "prbc"),
    #c("inv", "lac")
  )
  
  df %>%
    mutate(classifier_label = ifelse(classifier_label == "Total",
                                     "Total density",
                                     classifier_label)) %>%
    ggplot(aes(x = group, y = density)) +
    geom_boxplot(alpha = 0) +
    geom_jitter(aes(color = get(colorby)), height = 0 ) +
    facet_wrap(~classifier_label) +
    labs(color = colorby, y = paste({{cell_type}}, "density")) +
    stat_compare_means() +#KW test all samples
    stat_compare_means(comparisons = comps) +#pairwise
    ggtitle(paste({{panel}}, {{cell_type}}, "density in PPBC"))
}

