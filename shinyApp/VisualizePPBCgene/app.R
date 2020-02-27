#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyalert)
library(shinycssloaders)
library(here)
library(survival)
library(survminer)
library(tidyverse)

appDir <- "/DATA/share/postpartumbc/shinyApp/VisualizePPBCgene"
dataDir <- file.path(appDir, "data")

#Gene annotation data
gx_annot <- readRDS(file.path(dataDir,"12e_gx_annot.Rds"))

#Genewise overall survival
os <- readxl::read_excel(path = file.path(dataDir, "12_cox_allgenes.xlsx"),
                         sheet = "multi_cox_surv")
#Genewise distant recurrence
drs <- readxl::read_excel(path = file.path(dataDir, "12_cox_allgenes.xlsx"),
                          sheet = "multi_cox_dr")

#Interaction model overall survival
inv_int_os <- readxl::read_excel(file.path(dataDir, "12d_os_gene:inv_newformula.xlsx"))

#Interaction model distant recurrence
inv_int_drs <- readxl::read_excel(file.path(dataDir, "12d_drs_gene:inv_newformula.xlsx"))

#Differential expression results
res_list <- readRDS(file.path(dataDir, "12e_diffex_res_list.Rds"))

#Sample metadata
sample_data <- readRDS(file.path(dataDir,"12e_survival_sample_data.Rds"))

#TMM/log normalized gene expression matrices, with ensembl ids and gene symbols
sym_mat <- readRDS(file.path(dataDir, "12e_symbol_tmmnorm_genesxsample.Rds"))
ens_mat <- readRDS(file.path(dataDir, "12e_ensembl_tmmnorm_genesxsample.Rds"))

#Functions for plotting
source(here("src", "survival_tools.R"))
source(here("src", "retrieve_gene_summaries.R"))

#### UI ----
ui <- fluidPage(
    
    # Application title
    titlePanel("PPBC gene explorer"),
    
    #### Sidebar elements ----
    sidebarPanel(
        # Input gene ID
        textInput("id", "Gene ID", "MS4A1"),
        
        # Select symbol type
        selectInput("id_type", label = h3("Select id type"), 
                    choices = list("Gene name/symbol" = "symbol", "Ensembl ID" = "ensembl")),
        fluidRow(column(6,"Gene name:"), column(6, textOutput("gn"))),
        fluidRow(column(6,"Ensembl gene ID:"), column(6, textOutput("ens"))),
        fluidRow(column(6,"Entrez ID:"), column(6, textOutput("entrez"))),
        fluidRow(column(6,"Uniprot ID:"), column(6, textOutput("uniprot"))),
        fluidRow(column(6,"Gene biotype:"), column(6, textOutput("gene_type"))),
        fluidRow(column(6,"Description:"), column(6, textOutput("description")))
    ),
    

    #### Output data ----
    mainPanel(
        tabsetPanel(
            tabPanel("Gene summary",
                     useShinyalert(),
                     hr(),
                     h4("Entrez summary:"),
                     textOutput("entrez_summary") %>% shinycssloaders::withSpinner(),
                     br(),
                     hr(),
                     h4("Uniprot summary:"),
                     textOutput("uniprot_summary") %>% shinycssloaders::withSpinner(),
                     hr()
            ),
            tabPanel("Survival",
                     dataTableOutput("survival_summary"),
                     plotOutput("ntile_os") %>% shinycssloaders::withSpinner(),
                     br(),
                     br(),
                     plotOutput("ntile_drs") %>% shinycssloaders::withSpinner(),
                     br(),
                     br()
            ),
            tabPanel("Differential expression",
                     #textOutput("entrez_summary"),
                     dataTableOutput("diffex_report") %>% shinycssloaders::withSpinner(),
                     plotOutput("beehive") %>% shinycssloaders::withSpinner(),
                     br(),
                     br()
            )
        )
    )
)


#### Server ----
server <- function(input, output) {
    
    #### Show a warning message if the gene name is mapped to multiple IDs----

  observeEvent(input$id, {
    if (nrow(gene_lookup(input$id, id_type = input$id_type, dictionary = gx_annot)) > 1) {
    shinyalert(title = "Warning",
               text = paste0(
                 paste0("More than one ensembl ID detected for ", input$id, ": "),
                 "\n",
                 paste(gene_lookup(input$id, id_type = input$id_type, dictionary = gx_annot)$ensembl_gene_id, collapse = ", "),
                 "\n",
                 paste("Showing results for first in series:",
                       gene_lookup(input$id, id_type = input$id_type, dictionary = gx_annot)$ensembl_gene_id[1])),
               type = "warning")
  }})
    #Will only show if multiple ensembl IDs retrieved
    #output$warning <- renderText({
     #   paste("Multiple ensembl ids found for", input$id, ":",
      #        paste(gene_lookup(input$id, id_type = input$id_type, dictionary = gx_annot)$ensembl_gene_id, collapse = ", ")
       # )
    #})
    
    #output$thisens <- renderText({
     #   paste("Showing results for first in list: ", gene_lookup(input$id, id_type = input$id_type, dictionary = gx_annot)$ensembl_gene_id)[1]
    #})
    
    #An empty reactive element
    #e <- reactive({
     #   str_remove_all(input$id, ".*")
    #})
    
    #### Get the various IDs as reactive objects----
    gn <- reactive({
        head(gene_lookup(input$id, id_type = input$id_type, dictionary = gx_annot)$gene_name, 1)
    })
    output$gn <- renderText({gn()})
    
    ens <- reactive({
        head(gene_lookup(input$id, id_type = input$id_type, dictionary = gx_annot)$ensembl_gene_id, 1)
    })
    output$ens <- renderText({ens()})
    
    entrez <- reactive({
        head(gene_lookup(input$id, id_type = input$id_type, dictionary = gx_annot)$entrez_id, 1)
    })
    output$entrez <- renderText({entrez()})
    
    uniprot <- reactive({
        head(gene_lookup(input$id, id_type = input$id_type, dictionary = gx_annot)$uniprot_id, 1)
    })
    output$uniprot <- renderText({uniprot()})
    
    
    output$gene_type <- renderText({ 
        str_replace_all(head(gene_lookup(input$id, id_type = input$id_type, dictionary = gx_annot)$gene_type, 1), "_", " ")
    })
    
    output$description <- renderText({ 
        str_replace_all(head(gene_lookup(input$id, id_type = input$id_type, dictionary = gx_annot)$description, 1), "_", " ")
    })
    
    #### Retrieve gene summaries from entrez and uniprot ----
    
    output$entrez_summary <- renderText({
        get_entrez_summary(id = entrez())$entrez_summary
    })
    
    output$uniprot_summary <- renderText({
        get_uniprot_summary(id = uniprot())
    })
    
    #### Render the overall survival and distant recurrence table results for that gene---- 
    output$survival_summary <- renderDataTable({
        gene_survival(id = ens(),
                      id_type = "ensembl", s = os, d = drs, ios = inv_int_os, idrs = inv_int_drs,
                      gene_dict = gx_annot) %>%
            dplyr::select(-gene_name, -ensembl_gene_id, -gene_type, -description)
    })
    
    #### Render the differential expression table results for that gene----
    output$diffex_report <- renderDataTable({
        diffex_report(ensembl_id = ens(), 
                      list_reports = res_list, pthresh = 0.05, abslogfcthresh = 0.5) %>%
            dplyr::select(comparison, padj, log2FoldChange, sig)
    })
    
    #### Render the ntile kaplan meier os and drs plots for that gene----
    output$ntile_os <- renderPlot({
        km_ntiles_ovr(gene_id =  ens(),
                      id_type = "ensembl", 
                      survival_type = "os", ovr_column = "involution",
                      n = 3, labels = c("low", "medium", "high"), 
                      p_method = "anova", return_list = F, 
                      sampledata = sample_data, gene_dict = gx_annot, geneEx = ens_mat)
    })
    
    output$ntile_drs <- renderPlot({
        km_ntiles_ovr(gene_id = ens(),
                      id_type = "ensembl", 
                      survival_type = "drs", ovr_column = "involution",
                      n = 3, labels = c("low", "medium", "high"), 
                      p_method = "anova", return_list = F, 
                      sampledata = sample_data, gene_dict = gx_annot, geneEx = ens_mat)
    })
    
    #### Render the gene expression boxplot ----
    output$beehive <- renderPlot({
        tmm_plots(id = ens(),ensembl_mat = ens_mat, sampledata = sample_data,
                  id_type = "ensembl")$beehive +
            ggtitle(paste("TMM/log normalized expression of", gn()))
    })
}

#### Run the application ----
shinyApp(ui = ui, server = server)

