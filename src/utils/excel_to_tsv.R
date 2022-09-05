library(here)
library(dplyr)
library(readxl)
library(readr)

write_sheets_to_tsv <- function(inFile, outDir){
  
  sheets <- readxl::excel_sheets(inFile)
  print(paste("Converting", inFile, "to", "tsv"))
  
  for(sheet in sheets){
    df <- readxl::read_excel(inFile, sheet = sheet)
    outFile <- paste0(file.path(outDir, sheet), ".tsv")
    print(paste("Writing sheet", sheet, "to", outFile))
    readr::write_tsv(df, outFile)
  }
}

# sys.nframe == 0
# Allows sourcing this file to load the functions without running the code below
# If the script is called directly, the code below is run
if(sys.nframe()==0){
  
  # Parse the command line arguments
  source(here("src/utils/parse_args.R"))
  argdf <- retrieve_args()
  print(argdf)
  inFile <- dplyr::filter(argdf, argname == "metadata")$argval
  outDir <- dplyr::filter(argdf, argname == "outDir")$argval
  write_sheets_to_tsv(here(inFile), here(outDir))
}
