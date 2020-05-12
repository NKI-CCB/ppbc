#' Retrieve basic entrez information for a given gene
#'
#' @param id An entrez gene ID
#'
#' @return A list containing the following elements:
#' Name: The gene name
#' Description: Gene description (usually name written out)
#' Other aliases: Other known names for that gene as a single vector
#' Entrez summary: A text summary describing the most common functions of that gene
#' @export
#'
#' @examples get_entrez_summary(id = 931) #CD20/MS4A1
get_entrez_summary <- function(id){
  
  gene_summary <- rentrez::entrez_summary(db = "gene", id = id)
  
  res <- list()
  
  res$name <- rentrez::extract_from_esummary(gene_summary, "name")
  res$description <- rentrez::extract_from_esummary(gene_summary, "description")
  res$otheraliases <- rentrez::extract_from_esummary(gene_summary, "otheraliases")
  #res$otherdesignations <- rentrez::extract_from_esummary(gene_summary, "otherdesignations")
  res$entrez_summary <- rentrez::extract_from_esummary(gene_summary, "summary")
  return(res)
}


#' Retrieve the uniprot summary for a given gene by webscraping
#'
#' @param id A uniprot ID
#'
#' @return A single text string containing the uniprot summary for that gene
#' @export
#'
#' @examples get_uniprot_summary("O75534")  #Csde1
get_uniprot_summary <- function(id){
  
  base_url = "https://www.uniprot.org/uniprot/"
  url = paste0(base_url, id)
  
  uniprot_summary <- tryCatch(
  xml2::read_html(url) %>%
    rvest::html_node(xpath='/html/body/main/div/div[3]/div[3]/div[1]') %>%
    rvest::html_text(),
    error = function(e){NA} 
  )
  
  uniprot_summary <- uniprot_summary %>% str_split("[:digit:] Publications")
  uniprot_summary <- uniprot_summary[[1]][[1]]
  return(uniprot_summary)
}






