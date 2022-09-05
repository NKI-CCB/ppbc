library(shiny)
library(shinyalert)
library(shinycssloaders)
library(here)
library(survival)
library(survminer)
library(viridis)
library(rvest)
library(rentrez)
library(tidyverse)

#### Load data ----
{
appDir <- here("shinyApp/VisualizePPBCgene")
dataDir <- file.path(appDir, "data")

#Gene annotation data
gx_annot <- readRDS(file.path(dataDir,"app_gx_annot.Rds"))

#readxl::excel_sheets(file.path(dataDir, "12_cox_allgenes.xlsx"))

#Genewise overall survival
os <- readxl::read_excel(path = file.path(dataDir, "12_cox_allgenes.xlsx"),
                         sheet = "multi_genewise_os")

#Genewise distant recurrence
drs <- readxl::read_excel(path = file.path(dataDir, "12_cox_allgenes.xlsx"),
                          sheet = "multi_genewise_drs")

#Interaction model overall survival
inv_int_os <- readxl::read_excel(path = file.path(dataDir, "12_cox_allgenes.xlsx"),
                                 sheet = "multi_interaction_os")

#Interaction model distant recurrence
inv_int_drs <- readxl::read_excel(path = file.path(dataDir, "12_cox_allgenes.xlsx"),
                                  sheet = "multi_interaction_drs")

#Differential expression results
res_list <- readRDS(file.path(dataDir, "app_diffex_res_list.Rds"))
#names(res_list)

#Sample metadata
sample_data <- readRDS(file.path(dataDir,"app_survival_sample_data.Rds"))

#TMM/log normalized gene expression matrices, with ensembl ids and gene symbols
sym_mat <- readRDS(file.path(dataDir, "app_symbol_tmmnorm_genesxsample.Rds"))
ens_mat <- readRDS(file.path(dataDir, "app_ensembl_tmmnorm_genesxsample.Rds"))
}

#Functions for plotting
source(here("src", "rnaseq", "survival_tools.R"))
source(here("src", "rnaseq", "retrieve_gene_summaries.R"))

#### UI ----
ui <- fluidPage(
    
    # Application title
    titlePanel("PPBC gene explorer"),
    
    #### Sidebar elements ----
    sidebarPanel(
        # Input gene ID
        textInput("id", "Gene ID", "MS4A1"),
        
        # Select symbol type
        selectInput("id_type", label = "Select id type", 
                    choices = list("Gene name" = "symbol", "Ensembl ID" = "ensembl")),
        fluidRow(column(6,"Gene name:"), column(6, textOutput("gn"))),
        fluidRow(column(6,"Ensembl gene ID:"), column(6, textOutput("ens"))),
        fluidRow(column(6,"Entrez ID:"), column(6, textOutput("entrez"))),
        fluidRow(column(6,"Uniprot ID:"), column(6, textOutput("uniprot"))),
        fluidRow(column(6,"Gene biotype:"), column(6, textOutput("gene_type"))),
        downloadButton("report", "Generate report")#,
        #textOutput("checkrender")
    ),
    


    mainPanel(
        tabsetPanel(
          
          #### Survival panel----
            tabPanel("Survival",
                     hr(),
                     h3("Genewise survival results"),
                     tags$p("Model types 'overall survival' and 'distant recurrence' will identify genes associated with OS/DRS in all study groups",
                            "using the following multivariate Cox model:",
                            br(),
                            tags$code("survival ~ age + year of diagnosis + stage + grade + treatment + PAM50 + gene"),
                            br(),
                            "Displayed is the p value for the gene covariate, corrected for multiple testing across the entire dataset.",
                            br(),
                            br(),
                            "By contrast, interaction-type models seek to identify genes that behave differently in the involution patients",
                            "than in patients who are nulliparous, lactating or pregnant. The interaction model is as follows:",
                            br(),
                            tags$code("survival ~ age + year of diagnosis + stage + grade + treatment + PAM50 + gene + involution + involution*gene"),
                            br(),
                            "'Involution' is a binary variable indicating whether that sample belongs to an involuting patient (or not).",
                            "Displayed is the p value for the involution * gene covariate.",
                            br(),
                            "For all models, gene expression input was TMM and log2 normalized.",
                            "Only genes which pass the minimum count threshold were considered."),
                     
                     hr(),
                     tableOutput("survival_summary") %>% shinycssloaders::withSpinner(),
                     hr(),
                     h4("Kaplan-Meier curves for gene ntiles"),
                     hr(),
                     tags$p(
                       "Shown are user-defined quantiles ranging from 2-4 for the selected gene.",
                       "Samples are allocated to a ntile category using", tags$code("dplyr::ntile()"), ".",
                       br(),
                       br(),
                       "Unadjusted (univariate) curves are created using",
                       tags$code("survminer::ggsurvplot()"), ".",
                       "Displayed is the logrank p value.",
                       br(),
                       br(),
                       "Adjusted curves are created via", tags$code("survminer::ggadjustedcurves()"), "using the marginal method.",
                       "The marginal method uses a regression strategy to balance subpopulations.",
                       "A detailed description of the methodology can be found",
                       tags$a(href="https://cran.r-project.org/web/packages/survival/vignettes/adjcurve.pdf", "here"), ".",
                       br(),
                       br(),
                       "For adjusted curves, the p value is calculated via an anova that compares the fit from formula",
                       tags$code("survival ~ clinical covariate + gene ntile"),
                       "to the reduced formula fit", tags$code("survival ~ clinical covariates"), ".",
                       br(),
                       "Hormonetherapy is stratified to account for ER status, but PAM50 is not, as doing so introduces too many NA points."
                       
                     ),
                     selectInput("ntile", label = "Select Ntiles:", 
                                 choices = list("2" = 2,
                                                "3" = 3,
                                                "4" = 4),
                                 selected = "3"),
                     hr(),
                     plotOutput("ntile_os") %>% shinycssloaders::withSpinner(),
                     br(),
                     br(),
                     hr(),
                     plotOutput("ntile_drs") %>% shinycssloaders::withSpinner(),
                     br(),
                     br(),
                     tags$p("Treatment is a binary value with the following possible categories:",
                            br(),
                            "surgery, radiotherapy, hormonetherapy, chemotherapy, and herceptin")
            ),
          
          #### Diffex panel ----
            tabPanel("Differential expression",
                     hr(),
                     h3("Differential expression overview"),
                     tags$p(
                       "The table below shows the results for the selected gene",
                       "from a series of DESeq2 analyses:",
                       tags$ol(
                         tags$li("A likelihood ratio test (LRT), which identifies genes that differ in",
                                 tags$em("at least one"), "of the groups (notebook 6)."),
                         tags$li("A pairwise comparison, which compairs one group vs another (notebook 7). Possible values",
                                 tags$ol(
                                   tags$li("nonprbc = nulliparous"),
                                   tags$li("prbc = pregnant"),
                                   tags$li("lac = lactation"),
                                   tags$li("inv = involution")
                                 )),
                         tags$li("A one-vs-rest comparison, in which one group is",
                                 "paired vs all the rest pooled together ('vs rest', notebook 8).")
                       ),
                       br(),
                       "For all methods, the apeglm method was used for fold change shrinkage.",
                       "Only genes which pass the minimum count threshold were considered."
                     ),
                     hr(),
                     tableOutput("diffex_report") %>% shinycssloaders::withSpinner(),
                     hr(),
                     h4("Boxplot of gene expression"),
                     tags$p(
                       "Because the Cox regression was run on TMM-log normalized counts,",
                       "the same method is used to display expression below.",
                       "Of course, differential expression with DESeq2 was performed",
                       "on raw counts as is the prescribed methodology."
                     ),
                     selectInput("beehive_color", label = h4("Select color variable"), 
                                 choices = list("Overall survival" = "survival",
                                                "Distant recurrence" = "drs",
                                                "PAM50" = "PAM50")),
                     plotOutput("beehive") %>% shinycssloaders::withSpinner(),
                     br(),
                     br()
            ),
          
          #### Gene function panel ----
            tabPanel("Gene function",
                     hr(),
                     h4("Description:"),
                     textOutput("description") %>% shinycssloaders::withSpinner(),
                     hr(),
                     h4("Ensembl ID(s):"),
                     textOutput("ens_all") %>% shinycssloaders::withSpinner(),
                     hr(),
                     h4("Other aliases:"),
                     textOutput("aliases") %>% shinycssloaders::withSpinner(),
                     hr(),
                     h4("Entrez summary:"),
                     textOutput("entrez_summary") %>% shinycssloaders::withSpinner(),
                     br(),
                     hr(),
                     h4("Uniprot summary:"),
                     textOutput("uniprot_summary") %>% shinycssloaders::withSpinner(),
                     hr()
            )
        )
    )
)


#### Server ----
server <- function(input, output) {
    
  #### Multiple ensembl IDs ----
  
  #Show a warning message if the gene name is mapped to multiple IDs

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
   
    #### ID retrieval ----
  
    #Get the various IDs as reactive objects
    gn <- reactive({
        head(gene_lookup(input$id, id_type = input$id_type, dictionary = gx_annot)$gene_name, 1)
    })
    output$gn <- renderText({gn()})
    
    ens <- reactive({
        head(gene_lookup(input$id, id_type = input$id_type, dictionary = gx_annot)$ensembl_gene_id, 1)
    })
    output$ens <- renderText({ens()})
    
    output$ens_all <- renderText({
      paste(gene_lookup(input$id, id_type = input$id_type, dictionary = gx_annot)$ensembl_gene_id, collapse = ", ")
    })
    
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
    
    
    #### Survival----
    
    #Survival results table
    output$survival_summary <- renderTable({
        gene_survival(id = ens(),
                      id_type = "ensembl", s = os, d = drs, ios = inv_int_os, idrs = inv_int_drs,
                      gene_dict = gx_annot) %>%
            dplyr::select(-gene_name, -ensembl_gene_id, -gene_type, -description)
    })
    
    #Render the ntile kaplan meier os and drs plots for that gene
    output$ntile_os <- renderPlot({
      km_ntiles_ovr(gene_id =  ens(),
                    id_type = "ensembl", 
                    survival_type = "os", ovr_column = "involution",
                    n = as.integer(input$ntile),
                    line_colors = viridis::scale_color_viridis(discrete = T),
                    legend_positions = c("none","none", "bottom","bottom"),
                    p_method = "anova", return_list = F, 
                    sampledata = sample_data, gene_dict = gx_annot, geneEx = ens_mat)
    })
    
    output$ntile_drs <- renderPlot({
      km_ntiles_ovr(gene_id = ens(),
                    id_type = "ensembl", 
                    survival_type = "drs", ovr_column = "involution",
                    n = as.integer(input$ntile),
                    line_colors = viridis::scale_color_viridis(discrete = T),
                    legend_positions = c("none","none", "bottom","bottom"),
                    p_method = "anova", return_list = F, 
                    sampledata = sample_data, gene_dict = gx_annot, geneEx = ens_mat)
    })
    
    
    #### Differential expression----
    
    #Render the differential expression table results for that gene
    output$diffex_report <- renderTable({
        diffex_report(ensembl_id = ens(), 
                      list_reports = res_list, pthresh = 0.05, abslogfcthresh = 0.5) %>%
            dplyr::select(comparison, padj, log2FoldChange, sig)
    })
    
    
    # Gene expression boxplot
    output$beehive <- renderPlot({
        tmm_plots(id = ens(),ensembl_mat = ens_mat, sampledata = sample_data,
                  colorby = input$beehive_color, id_type = "ensembl")$beehive +
            ggtitle(paste("TMM/log normalized expression of", gn()))
    })
    
    #### Gene function  ----
    
    # Retrieve other aliases
    output$aliases <- renderText({
      get_entrez_summary(id = entrez())$otheraliases
    })
    
    #Retrieve gene summaries from entrez and uniprot
    output$entrez_summary <- renderText({
      get_entrez_summary(id = entrez())$entrez_summary
    })
    
    output$uniprot_summary <- renderText({
      get_uniprot_summary(id = uniprot())
    })
    
    #### Download handler ----
    
    #output$checkrender <- renderText({
     # if (identical(rmarkdown::metadata$runtime, "shiny")) {
      #  TRUE
      #} else {
      #  FALSE
      #}
    #})
    
    output$report <- downloadHandler(
      filename = "PPBC_gene_report.pdf",
      content = function(file) {
        withProgress(message = 'Rendering, please wait!', {
          # Copy the report file to a temporary directory before processing it, in
          # case we don't have write permissions to the current working dir (which
          # can happen when deployed).
          tempReport <- file.path(tempdir(), "gene_report_template.Rmd")
          file.copy("gene_report_template.Rmd", tempReport, overwrite = TRUE)
          
          # Set up parameters to pass to Rmd document
          params <- list(id = input$id,
                         id_type = input$id_type,
                         ntile = input$ntile)
          # Knit the document, passing in the `params` list, and eval it in a
          # child of the global environment (this isolates the code in the document
          # from the code in this app)
          rmarkdown::render(
            tempReport,
            output_file = file,
            params = params,
            envir = new.env(parent = globalenv())
          )
        })
      }
    )
} 


#### Run the application ----
shinyApp(ui = ui, server = server)

