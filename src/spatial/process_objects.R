#Standardize object file columns

library(here)
library(tidyverse)

outDir <- here("data/vectra/interim/objects")
dir.create(outDir, showWarnings = F)

objectfiles <- list.files(here("data/vectra/symlinks/objects"))
objectdir <- here("data/vectra/symlinks/objects")

# `rename_all(make.names)` will result in nonsensical output
# manually create syntactically valid columns

obj_name_repair <- function(cols){
  cols %>%
    tolower() %>%
    str_remove("\\(opal.*\\) ") %>%
    str_remove("\\:") %>%
    str_replace("\\(µm²\\)", "um2") %>%
    str_replace("\\(µm\\)", "um") %>%
    str_replace_all(" ", "_")
} 

# extract identifier into column and remove unnecessary columns
clean_obj <- function(df){
  df %>%
    separate(filename, into = c("t_number","panel","batch"), sep = "_", remove=F, extra = "drop") %>%
    relocate(t_number, panel, batch, .before = everything()) %>%
    select(-image_file_name)
} 

# read object and write processed object
process_object <- function(obj_file){
  
  obj <- read_csv(file.path(objectdir, obj_file)) %>%
    mutate(filename = obj_file)
  
  colnames(obj) <- obj_name_repair(colnames(obj)) 
  obj <- clean_obj(obj)
  
  outFile <- file.path(outDir, str_replace(obj$filename[1], "object_results.csv", "processed_object.csv"))
  
  obj <- obj %>% select(-filename)
  
  write_csv(obj, outFile)
}

#lapply(objectfiles[1:2], process_object)

for (f in objectfiles){process_object(f)}

