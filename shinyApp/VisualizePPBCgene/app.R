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

ui <- fluidPage(
    
    # Application title
    titlePanel("PPBC gene explorer"),
    
    # Input gene ID
    textInput("id", "Gene ID", "MS4A1"),
    
    # Select symbol type
    selectInput("id_type", label = h3("Select id type"), 
                choices = list("Gene name/symbol" = "symbol", "Ensembl ID" = "ensembl")),
    
    # Output data
    mainPanel(
        textOutput("gn"),
        textOutput("ens"),
        textOutput("gene_type"),
        textOutput("description"),
        tableOutput("survival_summary"),
        tableOutput("diffex_report"),
        plotOutput("ntile_os"),
        plotOutput("ntile_drs"),
        plotOutput("beehive")
    )
)


# Define server logic required to draw a histogram
server <- function(input, output) {
    
    
    output$gn <- renderText({ 
        paste("Gene name:", gene_lookup(input$id, id_type = input$id_type, dictionary = gx_annot)$gene_name)
        })
    
    output$ens <- renderText({
        paste("Ensembl gene id:", gene_lookup(input$id, id_type = input$id_type, dictionary = gx_annot)$ensembl_gene_id) 
        })
    
    output$gene_type <- renderText({ 
        paste("Gene biotype:", str_replace_all(gene_lookup(input$id, id_type = input$id_type, dictionary = gx_annot)$gene_type, "_", " "))
        })
    
    output$description <- renderText({ 
        paste("Description:", str_replace_all(gene_lookup(input$id, id_type = input$id_type, dictionary = gx_annot)$description, "_", " "))
    })
    
    output$survival_summary <- renderTable({
        gene_survival(id = gene_lookup(input$id, id_type = input$id_type, dictionary = gx_annot)$ensembl_gene_id,
                      id_type = "ensembl", s = os, d = drs, ios = inv_int_os, idrs = inv_int_drs,
                      gene_dict = gx_annot) %>%
            dplyr::select(-gene_name, -ensembl_gene_id, -gene_type, -description)
    })
    
    output$diffex_report <- renderTable({
        diffex_report(ensembl_id = gene_lookup(input$id, id_type = input$id_type, dictionary = gx_annot)$ensembl_gene_id, 
                      list_reports = res_list, pthresh = 0.05, abslogfcthresh = 0.5) %>%
            dplyr::select(comparison, padj, log2FoldChange, sig)
    })

    output$ntile_os <- renderPlot({
        km_ntiles_ovr(gene_id =  gene_lookup(input$id, id_type = input$id_type, dictionary = gx_annot)$ensembl_gene_id,
                      id_type = "ensembl", 
                      survival_type = "os", ovr_column = "involution",
                      n = 3, labels = c("low", "medium", "high"), 
                      p_method = "anova", return_list = F, 
                      sampledata = sample_data, gene_dict = gx_annot, geneEx = ens_mat)
    })
    
    output$ntile_drs <- renderPlot({
        km_ntiles_ovr(gene_id = gene_lookup(input$id, id_type = input$id_type, dictionary = gx_annot)$ensembl_gene_id,
                      id_type = "ensembl", 
                      survival_type = "drs", ovr_column = "involution",
                      n = 3, labels = c("low", "medium", "high"), 
                      p_method = "anova", return_list = F, 
                      sampledata = sample_data, gene_dict = gx_annot, geneEx = ens_mat)
    })
    
    output$beehive <- renderPlot({
        tmm_plots(id = gene_lookup(input$id, id_type = input$id_type, dictionary = gx_annot)$ensembl_gene_id,
                  id_type = "ensembl")$beehive +
            ggtitle(paste("TMM/log normalized expression of", gene_lookup(input$id, id_type = input$id_type, dictionary = gx_annot)$gene_name))
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

