webscrape_genecards <- function(gene_id){
  
  base_url = "http://www.genecards.org/cgi-bin/carddisp.pl?gene="
  url = paste0(base_url,gene_id)
  page = xml2::read_html(url)
  
  res <- list()
  
  res$entrez_summary <- page %>%
    rvest::html_node(xpath='/html/body/div[1]/div[3]/div/div/main/div[2]/div/div/section[2]/div[1]/ul/li/p') %>%
    rvest::html_text()
  
  res$civic_summary <- page %>%
    rvest::html_node(xpath='/html/body/div[1]/div[3]/div/div/main/div[2]/div/div/section[2]/div[2]/ul/li/p') %>%
    rvest::html_text()
  
  res$genecard_summary <- page %>%
    rvest::html_node(xpath='/html/body/div[1]/div[3]/div/div/main/div[2]/div/div/section[2]/div[3]/p') %>%
    rvest::html_text()
  
  res$tocris_summary <- page %>%
    rvest::html_node(xpath='/html/body/div[1]/div[3]/div/div/main/div[2]/div/div/section[2]/div[5]/ul/li/p') %>%
    rvest::html_text()
  
  #Doesn't work?
  res$uniprot_summary <- page %>%
    rvest::html_node(xpath="/html/body/div[1]/div[3]/div/div/main/div[2]/div/div/section[2]/div[4]/ul/li/div") %>%
    rvest::html_text()
  
  return(res)
}

#Works
webscrape_genecards("TP53")

#Doesn't entirely work
webscrape_genecards("MS4A1")
