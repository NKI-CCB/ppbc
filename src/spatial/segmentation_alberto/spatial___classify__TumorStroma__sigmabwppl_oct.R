library(tidyverse)
library(spatstat)
library(dbscan)

# function to assign Tumor/Stroma using KDEtumor, KDEstroma, and input coordinates
tumor_or_stroma <- function(kde_tumor, kde_stroma, location){
  
  tumor_value <- kde_tumor %>% data.frame %>%
    mutate(dist=(x-location[1])^2+(y-location[2])^2) %>%
    arrange(dist) %>% slice(1) %>% pull(value)
  
  stroma_value <- kde_stroma %>% data.frame %>%
    mutate(dist=(x-location[1])^2+(y-location[2])^2) %>%
    arrange(dist) %>% slice(1) %>% pull(value)
  
  assigned_loc <- ifelse(tumor_value > stroma_value, 'T', 'S')
  
  return(assigned_loc)
}

# Load data
# This if your df with the x/y coordinates
## it should contain at least: tnumber, Xcenter, Ycenter, phenotype
foci_data <- read_rds("/DATA/share/VECTRA/RMKDN/cells_areas_clinical.rds")$cellobjects
odir <- '~/nabucco/results/spatial/'# define output directory

#### 1] Compute "foci's" (each tumor island will be assigned an ID)
#### (only for biopsy material, when you have 'tumor islands'. )
for(tnumber_sel in unique(foci_data$tnumber)){
  # compute clusters
  dbclusters <- dbscan(as.matrix(foci_data %>% filter(tnumber == tnumber_sel) %>% dplyr::select(Xcenter,Ycenter)), 
                       eps = 300, minPts = 50)
  # Assign clusters / foci's
  foci_data[foci_data$tnumber == tnumber_sel,'ClusterID'] <- dbclusters$cluster %>% as.character
}
foci_data <- foci_data %>% mutate(sample = paste(tnumber, ClusterID))


#### 2] Compute the smoothing bandwith KDE by likelihood cross-validation for the kernel estimation of point process intensity. 
#### (will be done for each foci, using either tumor cells or stromal/negative cells)
# DF to store optimal bandwiths:
sigmas_ppl <- data.frame(matrix(nrow=0,ncol=3))
colnames(sigmas_ppl) <- c('foci_sample','sigma', 'tumor_stroma')
sigmas_ppl$foci_sample <- as.character(sigmas_ppl$foci_sample)
sigmas_ppl$sigma <- as.numeric(sigmas_ppl$sigma)
sigmas_ppl$tumor_stroma <- as.character(sigmas_ppl$tumor_stroma)

# Do the computations for the tumor (PanCK+) or Stroma (negative) data 
for(foci_sample in samples){
  # Compute point pattern tumor cells
  x <- foci_data %>% filter(sample == foci_sample) %>% filter(phenotype == 'PanCK+')
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
  
  # Compute point pattern stromal cells
  x <- foci_data %>% filter(sample == foci_sample) %>% filter(phenotype == 'negative')
  spdat_stroma <- ppp(
    x = x$Xcenter,
    y = x$Ycenter,
    window = owin(
      xrange = c(min(x$Xcenter), max(x$Xcenter)),
      yrange = c(min(x$Ycenter), max(x$Ycenter))
    ),
    marks = x$phenotype
  )
  Window(spdat_stroma) <- ripras(spdat_stroma)
  
  # Compute optimal sigma tumor
  optimal_sigma <- bw.ppl(spdat_tumor)
  sigmas_ppl <- bind_rows(sigmas_ppl, 
                          data.frame(foci_sample=foci_sample, sigma=as.numeric(optimal_sigma), tumor_stroma='T'))
  
  # Compute optimal sigma stroma
  optimal_sigma <- bw.ppl(spdat_stroma)
  sigmas_ppl <- bind_rows(sigmas_ppl, 
                          data.frame(foci_sample=foci_sample, sigma=as.numeric(optimal_sigma), tumor_stroma='S'))
}

colnames(sigmas_ppl) <- c('foci_sample','sigma','tumor_stroma')
write_delim(sigmas_ppl, 
            paste(odir, 'optimized_kde_bandwith__bw_ppl.tsv',sep=''),
            delim='\t')


#### 3] Classify Tumor / Stroma region for each cell position
####
message('Assigning locations')
overall_assigned_locs <- list()
for(foci in samples){
  message(foci)
  test <- foci_data %>% filter(sample == foci) %>% dplyr::select(Xcenter, Ycenter)  %>% as.matrix
  ## Compute point pattern for tumor cells
  x <- foci_data %>% filter(sample == foci) %>% filter(phenotype == 'PanCK+')
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
  
  ## Compute point pattern for stromal cells
  x <- foci_data %>% filter(sample == foci) %>% filter(phenotype == 'negative')
  spdat_stroma <- ppp(
    x = x$Xcenter,
    y = x$Ycenter,
    window = owin(
      xrange = c(min(x$Xcenter), max(x$Xcenter)),
      yrange = c(min(x$Ycenter), max(x$Ycenter))
    ),
    marks = x$phenotype
  )
  
  # Compute KDE tumor
  sigma_tumor <- sigmas_ppl %>% filter(foci_sample == foci) %>% 
    filter(get('tumor_stroma') == 'T') %>% pull(sigma) %>% as.numeric
  k1_tumor_raw <-     density(spdat_tumor, edge=TRUE, 
                              sigma=sigma_tumor)
  # Compute KDE stroma
  sigma_stroma <- sigmas_ppl %>% filter(foci_sample == foci) %>% 
    filter(get('tumor_stroma') == 'S') %>% pull(sigma) %>% as.numeric
  k1_stroma_raw <- density(spdat_stroma, edge=TRUE, 
                           sigma=sigma_stroma)
  
  # Classify Tumor / Stroma region for all cell positions from the sample
  ## KDE's are normalized by maximum value of KDE
  assigned_locs <- sapply(1:nrow(foci_data %>% filter(sample == foci) ), 
                          function(x) tumor_or_stroma(k1_tumor_raw / max(k1_tumor_raw), 
                                                      k1_stroma_raw / max(k1_stroma_raw), 
                                                      c(test[x,1], test[x,2])))
  overall_assigned_locs[[foci]] <- assigned_locs
}


foci_data_with_assigned_locations <- foci_data %>% slice(0) %>%
  mutate(ClusterID=as.character(ClusterID))

# Add assigned locations to each cell from your sample
for(foci in samples){
  foci_data_with_assigned_locations <- bind_rows(foci_data_with_assigned_locations,
                                                 foci_data %>% filter(sample == foci) %>%
                                                              mutate(assigned_loc = overall_assigned_locs[[foci]])
                                                 )
}

write_delim(foci_data_with_assigned_locations, 
            paste(odir, 'classified_cells.tsv',sep=''),
            delim='\t')
