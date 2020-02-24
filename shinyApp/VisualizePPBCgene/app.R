#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
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
        textOutput("gn"),
        textOutput("ens"),
        textOutput("gene_type"),
        textOutput("description")
    ),
    
    #### Output data ----
    mainPanel(
        tabsetPanel(
            tabPanel("Results summary",
                     conditionalPanel(
                         condition = "output.len_ens",
                         span(h3("Warning"), style="color:red"),
                         span(textOutput("warning"), style = "color:red"),
                         span(textOutput("thisens"), style = "color:red")),
                     tableOutput("survival_summary"),
                     tableOutput("diffex_report")
            ),
            tabPanel("Kaplan-meier gene ntiles",
                     plotOutput("ntile_os"),
                     plotOutput("ntile_drs")    
            ),
            tabPanel("Boxplot gene expression",
                     plotOutput("beehive")
            )
        )
    )
)


#### Server ----
server <- function(input, output) {
    
    #### Show a warning message if the gene name is mapped to multiple IDs----
    output$len_ens <- reactive({
        nrow(gene_lookup(input$id, id_type = input$id_type, dictionary = gx_annot)) > 1
    })
    #Output values can only be used when they are rendered elsewhere in the ui
    #Workaround is to set outputOptions
    outputOptions(output, "len_ens", suspendWhenHidden = FALSE)
    
    #Will only show if multiple ensembl IDs retrieved
    output$warning <- renderText({
        paste("Multiple ensembl ids found for", input$id, ":",
              paste(gene_lookup(input$id, id_type = input$id_type, dictionary = gx_annot)$ensembl_gene_id, collapse = ", ")
        )
    })
    
    output$thisens <- renderText({
        paste("Showing results for first in list: ", gene_lookup(input$id, id_type = input$id_type, dictionary = gx_annot)$ensembl_gene_id)[1]
    })
    
    #### Get the gene symbol and ensembl ID as reactive objects----
    gn <- reactive({
        head(gene_lookup(input$id, id_type = input$id_type, dictionary = gx_annot)$gene_name, 1)
    })
    
    ens <- reactive({
        head(gene_lookup(input$id, id_type = input$id_type, dictionary = gx_annot)$ensembl_gene_id, 1)
    })
    
    #### Format the above fields alongside gene biotype and description for plotting----
    output$gn <- renderText({ 
        paste("Gene name:", gn())
    })
    
    output$ens <- renderText({
        paste("Ensembl gene id:", ens()) 
    })
    
    output$gene_type <- renderText({ 
        paste("Gene biotype:", str_replace_all(head(gene_lookup(input$id, id_type = input$id_type, dictionary = gx_annot)$gene_type, 1), "_", " "))
    })
    
    output$description <- renderText({ 
        paste("Description:", str_replace_all(head(gene_lookup(input$id, id_type = input$id_type, dictionary = gx_annot)$description, 1), "_", " "))
    })
    
    #### Render the overall survival and distant recurrence table results for that gene---- 
    output$survival_summary <- renderTable({
        gene_survival(id = ens(),
                      id_type = "ensembl", s = os, d = drs, ios = inv_int_os, idrs = inv_int_drs,
                      gene_dict = gx_annot) %>%
            dplyr::select(-gene_name, -ensembl_gene_id, -gene_type, -description)
    })
    
    #### Render the differential expression table results for that gene----
    output$diffex_report <- renderTable({
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
        tmm_plots(id = ens(),
                  id_type = "ensembl")$beehive +
            ggtitle(paste("TMM/log normalized expression of", gn()))
    })
}

#### Run the application ----
shinyApp(ui = ui, server = server)

