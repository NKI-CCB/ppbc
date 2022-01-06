library(tidyverse)
library(spatstat)
library(dbscan)

# Snippet to compute areas on your slide

foci_data_with_assigned_locations <- read_delim(
            '~/nabucco/results/spatial_analysis/kde_estimations/foci_data_with_assigned_locations___kde_sigma___bw_ppl.tsv',
            delim='\t') # dataframe with x/y coordinates
sigmas_ppl <- read_delim(
            '~/nabucco/results/spatial_analysis/kde_estimations/kde__per_foci/optimized_kde_bandwith__bw_ppl.tsv',
            delim='\t') # dataframe with optimized bandwiths KDE (computed in the 'classify' script/snippet)

# Create DF to store areas 
# This df will contain Area of Tumor and Area of Stroma, but only area of Tumor will be used
# (Afterwards area of Tumor+Stroma will be calculated, and Stroma area will be infered as Area(Tumor+Stroma)-Area(Tumor) )
# Threshold to compute areas will be 0.1 (different thresholds are computed for exploratory purposes)
areas_df <- data.frame(c('a'), c('b'), c(1), c(1))
colnames(areas_df) <- c('foci_sample','tumor_stroma', 'threshold','area')
areas_df <- areas_df %>% slice(-1)

for(focisample in unique(foci_data_with_assigned_locations %>% pull(sample))){
  for(tumor_stroma in c('tumor','stroma')){
    pl <- list()
    # Iterate through different thresholds
    for(threshold in c(0.05, 0.1, 0.2)){
      t_s_cell <- ifelse(tumor_stroma == 'tumor', 'PanCK+','negative')
      t_s_compartment <- ifelse(tumor_stroma=='tumor','T','S')
      
      # Create point pattern of either tumor or stromal cells
      x <- foci_data_with_assigned_locations %>% 
        filter(sample == focisample) %>% 
        filter(phenotype == t_s_cell)
      spdat_tumor <- ppp(
        x = x$Xcenter,
        y = x$Ycenter,
        window = owin(
          xrange = c(min(x$Xcenter), max(x$Xcenter)),
          yrange = c(min(x$Ycenter), max(x$Ycenter))
        ),
        marks = x$phenotype
      )
      Window(spdat_tumor) <- ripras(spdat_tumor)
      
      # Compute KDE
      k1 <- density(spdat_tumor, edge=TRUE, 
                    sigma=sigmas_ppl %>% 
                      filter(foci_sample == focisample) %>% 
                      filter(`tumor_stroma` == t_s_compartment) %>% pull(sigma))

      k1 <- k1 / max(k1) # normalize KDE by max

      pl[threshold] <- ggplot() + 
        geom_point(data=data.frame(k1) %>%  mutate(value = ifelse(value < threshold, NA, value)), 
                   aes(x=x, y=y, color=value)) + 
        geom_point(data=data.frame(spdat_tumor), aes(x=x, y=y),size=0.1, color='white') +
        ggtitle(paste('White points = cells\nColoured area = pixels considered for area calculation',
                      '\n', focisample, tumor_stroma, 'threshold=',threshold))
      pl[threshold]
      ggsave(filename = paste(focisample,'__', tumor_stroma, '_AAA', threshold, '.png'),
             width=7, height=5,
             path = '~/nabucco/results/spatial_analysis/kde_estimations/foci_area_calculations/considered_areas_maps/',
             # useDingbats=FALSE, 
             dpi=150)
      
      data.frame(k1) %>% ggplot(aes(x=value)) + geom_histogram() + theme_bw() +
        geom_vline(xintercept=0.05) + geom_vline(xintercept=0.1) + geom_vline(xintercept=0.2) +
        ggtitle(paste('KDE estimates histogram',
                      '\n', focisample, tumor_stroma))
      ggsave(filename = paste(focisample,'__', tumor_stroma, '.pdf'),
             width=4, height=3,
             path = '~/nabucco/results/spatial_analysis/kde_estimations/foci_area_calculations/histograms/',
             # useDingbats=FALSE, 
             dpi=150)
      
      # Compute area of a pixel, and multiply by number of pixels that 'survive' after filtering by threshold to estimate area
      x_vector <- data.frame(k1) %>% pull(x)
      x_dists <- lapply(2:(length(x_vector)), function(x) x_vector[x] - x_vector[x-1]) %>% as.numeric
      
      x_distance_mean <- x_dists[x_dists > 0] %>% mean
      
      y_vector <- data.frame(k1) %>% pull(y)
      y_dists <- lapply(2:(length(y_vector)), function(x) y_vector[x] - y_vector[x-1]) %>% as.numeric
      y_distance_mean <- y_dists[y_dists > 0] %>% mean
      
      area <- data.frame(k1) %>%  mutate(value = ifelse(value < threshold, NA, value)) %>%
        filter(!is.na(value)) %>%
        nrow() * x_distance_mean * y_distance_mean
      area <- area / 1000000 # Convert area to squared microns     
      
      # Append calculated area
      areas_df <- bind_rows(areas_df, 
                            data.frame(foci_sample = focisample,
                                       tumor_stroma = tumor_stroma, 
                                       threshold = threshold, 
                                       area = area)
      )
    }
  }
}

# This df will contain Area of Tumor and Area of Stroma, but only area of Tumor will be used
# (Afterwards area of Tumor+Stroma will be calculated, and Stroma area will be infered as Area(Tumor+Stroma)-Area(Tumor) )
# Threshold to compute areas will be 0.1 (different thresholds are computed for exploratory purposes)
write_delim(areas_df,
            '~/nabucco/results/spatial_analysis/kde_estimations/foci_area_calculations/areas_df__separated_tumor_stroma.tsv', delim='\t')


# Calculate total area of tissue
kdecombined_df <- data.frame(x=c(1),y=c(1),value=c(1),optional=c(TRUE),foci_sample=c('a'))
kdecombined_df <- kdecombined_df %>% slice(-1)
# This df will contain Area of Tumor+Stroma (So, KDE estimated altogether using Tumor AND Stroma cells)
# (Afterwards area of Tumor+Stroma will be calculated, and Stroma area will be infered as Area(Tumor+Stroma)-Area(Tumor) )
# Threshold to compute areas will be 0.1 (different thresholds are computed for exploratory purposes)
areas_total_df <- data.frame(c('a'), c('b'), c(1), c(1))
colnames(areas_total_df) <- c('foci_sample','tumor_stroma', 'threshold','area')
areas_total_df <- areas_total_df %>% slice(-1)
for(focisample in unique(foci_data_with_assigned_locations %>% pull(sample))){
  for(threshold in c(0.05, 0.1, 0.2)){
    tumor_stroma <- 'tumor'
    
    t_s_cell <- ifelse(tumor_stroma == 'tumor', 'PanCK+','negative')
    t_s_compartment <- ifelse(tumor_stroma=='tumor','T','S')
    
    x_fulldata <- foci_data_with_assigned_locations %>% 
      filter(sample == focisample) %>% 
      filter(phenotype %in% c('PanCK+','negative'))
    
    full_window <-   ppp(
      x = x_fulldata$Xcenter,
      y = x_fulldata$Ycenter,
      window = owin(
        xrange = c(min(x_fulldata$Xcenter), max(x_fulldata$Xcenter)),
        yrange = c(min(x_fulldata$Ycenter), max(x_fulldata$Ycenter))
      ),
      marks = x_fulldata$phenotype
    )
    Window(full_window) <- ripras(full_window)
    full_window_owin <- Window(full_window)
    
    x <- foci_data_with_assigned_locations %>% 
      filter(sample == focisample) %>% 
      filter(phenotype == t_s_cell)
    
    spdat_tumor <- ppp(
      x = x$Xcenter,
      y = x$Ycenter,
      window = owin(
        xrange = c(min(foci_data_with_assigned_locations %>% filter(sample == focisample) %>% pull(Xcenter)), 
                   max(foci_data_with_assigned_locations %>% filter(sample == focisample) %>% pull(Xcenter))),
        yrange = c(min(foci_data_with_assigned_locations %>% filter(sample == focisample) %>% pull(Ycenter)),
                   max(foci_data_with_assigned_locations %>% filter(sample == focisample) %>% pull(Ycenter)))
      ),
      marks = x$phenotype
    )
    Window(spdat_tumor) <- full_window_owin # implement full window 
    
    k1 <- density(spdat_tumor, edge=TRUE, 
                  sigma=sigmas_ppl %>% 
                    filter(foci_sample == focisample) %>% 
                    filter(`tumor_stroma` == t_s_compartment) %>% pull(sigma))
    k1 <- k1 / max(k1) # normalize KDE

    sptumor <- spdat_tumor
    k1tumor <- k1
    
    tumor_stroma <- 'stroma'
    
    t_s_cell <- ifelse(tumor_stroma == 'tumor', 'PanCK+','negative')
    t_s_compartment <- ifelse(tumor_stroma=='tumor','T','S')
    
    x <- foci_data_with_assigned_locations %>% 
      filter(sample == focisample) %>% 
      filter(phenotype == t_s_cell)
    spdat_tumor <- ppp(
      x = x$Xcenter,
      y = x$Ycenter,
      window = owin(
        xrange = c(min(foci_data_with_assigned_locations %>% filter(sample == focisample) %>% pull(Xcenter)), 
                   max(foci_data_with_assigned_locations %>% filter(sample == focisample) %>% pull(Xcenter))),
        yrange = c(min(foci_data_with_assigned_locations %>% filter(sample == focisample) %>% pull(Ycenter)),
                   max(foci_data_with_assigned_locations %>% filter(sample == focisample) %>% pull(Ycenter)))
      ),
      marks = x$phenotype
    )
    Window(spdat_tumor) <- full_window_owin
    
    k1 <- density(spdat_tumor, edge=TRUE, 
                  sigma=sigmas_ppl %>% 
                    filter(foci_sample == focisample) %>% 
                    filter(`tumor_stroma` == t_s_compartment) %>% pull(sigma))
    # sigma=sigmas_diggle %>% filter(foci_sample == 'T18-23748 6') %>% filter(`T/S` == 'T') %>% pull(sigma))
    
    k1 <- k1 / max(k1)
    # plot(k1, main=NULL, las=1)
    # contour(k1, add=TRUE)
    
    spstroma <- spdat_tumor
    k1stroma<- k1
    
    kdecombined_df <- bind_rows(kdecombined_df, 
                                data.frame(k1tumor + k1stroma) %>% 
                                  mutate(foci_sample = rep(focisample, nrow(data.frame(k1tumor + k1stroma))))
    )
    
    data.frame(k1tumor + k1stroma) %>% 
      ggplot(aes(x=value)) + geom_histogram() + theme_bw() +
      geom_vline(xintercept=0.05) + geom_vline(xintercept=0.1) + geom_vline(xintercept=0.2) +
      ggtitle(paste('KDE Tumor + KDE stroma, estimates histogram',
                    '\n', focisample))
    
    ggsave(filename = paste(focisample,'__', tumor_stroma, '.pdf'),
           width=4, height=3,
           path = '~/nabucco/results/spatial_analysis/kde_estimations/foci_area_calculations/histogram_kdetumorpluskdestroma/',
           # useDingbats=FALSE, 
           dpi=150)
    
    
    plot <- ggplot() + 
      geom_point(data=data.frame(k1tumor+k1stroma) %>%  mutate(value = ifelse(value < threshold, NA, value)), 
                 aes(x=x, y=y, color=value)) + 
      geom_point(data=bind_rows(data.frame(sptumor), data.frame(spstroma)), aes(x=x, y=y),size=0.01, color='white') +
      ggtitle(paste('White points = cells\nColoured area = pixels considered for area calculation',
                    '\n', focisample, tumor_stroma, 'threshold=',threshold))
    plot
    ggsave(filename = paste(focisample,'__', '___threshold_', threshold, '.png'),
           width=7, height=5,
           path = '~/nabucco/results/spatial_analysis/kde_estimations/foci_area_calculations/areas_maps_tumor_plus_stroma/',
           # useDingbats=FALSE, 
           dpi=150)
    
    # Compute area of a pixel, and multiply by number of pixels that 'survive' after filtering by threshold to estimate area
    x_vector <- data.frame(k1tumor + k1stroma) %>% pull(x)
    x_dists <- lapply(2:(length(x_vector)), function(x) x_vector[x] - x_vector[x-1]) %>% as.numeric
    
    x_distance_mean <- x_dists[x_dists > 0] %>% mean
    
    y_vector <- data.frame(k1tumor + k1stroma) %>% pull(y)
    y_dists <- lapply(2:(length(y_vector)), function(x) y_vector[x] - y_vector[x-1]) %>% as.numeric
    y_distance_mean <- y_dists[y_dists > 0] %>% mean
    
    area <- data.frame(k1tumor + k1stroma) %>%  
      mutate(value = ifelse(value < threshold, NA, value)) %>%
      filter(!is.na(value)) %>%
      nrow() * x_distance_mean * y_distance_mean
    area <- area / 1000000 # Compute area in squared microns     
    
    # Append calculated area
    areas_total_df <- bind_rows(areas_total_df, 
                                data.frame(foci_sample = focisample,
                                           tumor_stroma = 'Tumor+Stroma', 
                                           threshold = threshold, 
                                           area = area)
    )
    
  }
}

# This df  contains Area of Tumor +Stroma
# (Afterwards area of Tumor+Stroma will be calculated, and Stroma area will be infered as Area(Tumor+Stroma)-Area(Tumor) )
# Threshold to compute areas will be 0.1 (different thresholds are computed for exploratory purposes)
write_delim(areas_total_df,
            '~/nabucco/results/spatial_analysis/kde_estimations/foci_area_calculations/areas_total_df__separated_tumor_stroma.tsv', delim='\t')


# Compute the areas 
# Estimated by summing pixels with normalized KDE >= 0.1
# Area Tumor --> From KDE Tumor
# Area Stroma --> Area (Tumor+Stroma) - Area (Tumor)
estimated_areas_foci <- left_join(areas_df %>% filter(tumor_stroma == 'tumor') %>%
                                    filter(threshold == 0.1) %>% 
                                    dplyr::select(-threshold) %>%
                                    spread(key='tumor_stroma', value='area'),
                                  areas_total_df %>%
                                    filter(threshold == 0.1) %>% 
                                    dplyr::select(-threshold) %>%
                                    spread(key='tumor_stroma', value='area')
            ) %>% 
      mutate(Stroma = `Tumor+Stroma` - `tumor`) %>% 
      mutate(tumor_perc = (`tumor` / `Tumor+Stroma`)*100) %>% 
      arrange(desc(tumor_perc)) 

write_delim(estimated_areas_foci,
            '~/nabucco/results/spatial_analysis/kde_estimations/foci_area_calculations/estimated_areas_foci.tsv', delim='\t')
