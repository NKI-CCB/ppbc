library(httr)
library(XML)

#### Entrez gene summary  ----
#CD20 = id 931

#Option 1
#Retrieve the page as an xml document using httr
#Then parse into data frame with XML
#Not all info available via this route (no additional aliases)
page <- httr::GET("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id=931&retmode=xml")
xml <- httr::content(page, "parsed")
df <- XML::xmlToDataFrame(xmlParse(xml))

colnames(df)

#Summary info
as.character(df$Entrezgene_summary)


#Option 2
#Using rentrez
#Allows both summary and other aliases
library(rentrez)
gene_summary <- rentrez::entrez_summary(db = "gene", id = 931)
rentrez::extract_from_esummary(gene_summary, "otheraliases")
rentrez::extract_from_esummary(gene_summary, "otherdesignations")
rentrez::extract_from_esummary(gene_summary, "summary")

#### Uniprot Summary ----
#Non functional without the magick framework, which we lack permissions to install
#install.packages("UniprotR")
#GetProteinFunction("O14520")

#Genes accessible as XMLs using this format:
#https://www.uniprot.org/uniprot/P11836.xml #CD20
#df <- XML::xmlToDataFrame(xmlParse(xml)) #Duplicate subscripts
#l <- XML::xmlToList(xmlParse(xml)) #Yields a completely unintelligible list

#Or via text using this format
#url <- "https://www.uniprot.org/uniprot/P11836.txt"
#page <- httr::GET(url)
#txt <- httr::content(page, "parsed")

con = url("https://www.uniprot.org/uniprot/P11836.txt")
lines = readLines(con)

#Lines starting with CC contain the necessary info but are annoying to parse together
cc <- lines[str_detect(lines, "^CC")]

#Or just web scrape it

library(rvest)
url = "https://www.uniprot.org/uniprot/P11836"
url = "https://www.uniprot.org/uniprot/O75534"
page = xml2::read_html(url)


uniprot_summary <- page %>%
  rvest::html_node(xpath='/html/body/main/div/div[3]/div[3]/div[1]') %>%
  rvest::html_text()

uniprot_summary <- uniprot_summary %>% str_split("[:digit:] Publications")
uniprot_summary[[1]][[1]]  
