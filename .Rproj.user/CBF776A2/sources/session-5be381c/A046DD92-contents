#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#


library(shiny)
library(bslib)
library(tidyr)
library(ggplot2)
library(readxl)

options(shiny.maxRequestSize = 1000 * 1024^2)

source("modules/helpers.R")
source("modules/load_data.R")
source("modules/gene_violinplot.R")
source("modules/gene_selector_module.R")
source("modules/filtering_annot_module.R")
source("modules/filtering_data_module.R")
source("modules/housekeeper_module.R")


# Define UI for application that draws a histogram
ui <- fluidPage(
  tags$head(
    # Preconnect & preload Poppins font
    tags$link(rel="preconnect", href="https://fonts.googleapis.com"),
    tags$link(rel="preconnect", href="https://fonts.gstatic.com", crossorigin=NA),
    tags$link(
      href="https://fonts.googleapis.com/css2?family=Poppins:wght@400;600&display=swap",
      rel="stylesheet"
    )#,
    # Custom styles
    #tags$link(rel="stylesheet", href="styles.css")
  ),
  theme = bs_theme(
    version = 4,
    base_font = "Poppins",
    heading_font = "Poppins",
    bootswatch = "flatly",
    bg = "#F4F7FA",
    fg = "#2A4D69",
    primary = "#4B86B4",
    secondary = "#ADCBE3",
    font_scale = 1.0
  ),
  titlePanel("RNA-seq Expression Database Viewer"),
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(
        id = "sidebar_tabs",
        tabPanel(
          "Instructions",
          h4("Instructions"),
          tags$ol(
            tags$li("Upload your RNAseq RDS file and annotation file in the 'Load data' tab."),
            tags$li("Annotation file must include: InternalID, IndividualID, Family, Material, and optionally 'Affected'."),
            tags$li("Switch to 'Gene Expression' to filter samples and select a gene."),
            tags$li("Use dropdowns to filter by material, controls, family, or individual. Download plots as needed."),
            tags$li("See 'Housekeeper Expression' for housekeeper gene plots. Download plots as needed.")
          )
        ),
        tabPanel(
          "Adding Data",
          h4("Adding Data"),
          tags$ol(
            tags$li("Add new samples to your annotation file with required columns: Salmon, InternalID, IndividualID, Family, Material, and optionally 'affected'."),
            tags$li("For new samples, fill in the Salmon column with the path to quantification files."),
            tags$li("Upload your updated RDS and annotation files in 'Load data'."),
            tags$li("Click 'Load missing samples' to add new samples."),
            tags$li("Click 'Save updated data as RDS' to download the updated dataset."),
            tags$li("Reuse the updated annotation and dataset files for future sessions.")
          ),
          helpText("See documentation or contact the app maintainer for more details.")
        ),
        tabPanel(
          "Annotation File Structure",
          h4("Annotation File Structure"),
          tags$p("Your annotation file should be a tab-separated (.tsv) or comma-separated (.csv) text file with the following columns:"),
          tags$ul(
            tags$li(strong("InternalID:"), " Unique identifier for each sample matching column names in the RDS file (e.g. ACC1111111 )."),
            tags$li(strong("IndividualID:"), " Identifier for the individual/sample (e.g. 20201-I-1A)."),
            tags$li(strong("Family:"), " Family or group identifier (e.g. 20201)."),
            tags$li(strong("Material:"), " Sample material type (e.g., blood, fibroblast, muscle)."),
            tags$li(strong("Salmon:"), " (Optional) Path to quantification files for new samples (e.g. data/salmon/blood/ACC1111111_quant.sf)."),
            tags$li(strong("Affected:"), " (Optional) Status of the sample (e.g., affected/unaffected)."),
            tags$li("Additional columns are allowed but not required")
          ),
          tags$p("Make sure column names are spelled exactly as shown.")
        )
      ),
      width = 3
    ),
    mainPanel(
      tabsetPanel(
        id = "tabs",
        tabPanel(
          "Load data", 
          buttons_tab1_UI("load_btns")
        ),
        tabPanel(
          "Gene Expression", 
          # Filter embedded directly in this tab
          annotFilterUI("gene_filters"),
          geneSearchUI("gene_search"),
          gene_plot_UI("gene_plot")
        ),
        tabPanel(
          "Housekeeper Expression",
          housekeeperUI("housekeeper")
        )
      ),
      width = 8
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  # Load data module
  loaded_data <- buttons_tab1_Server("load_btns")
  
  # Gene sync module (pass the reactive data)
  gene_selected <- geneSearchServer("gene_search", loaded_data$loaded_data)
  
  # Annotation filter module for gene expression tab
  gene_filters <- annotFilterServer("gene_filters", loaded_data$annotation_data)
  
  # Data filtering module: filter by material and selected individuals
  filtered_data <- annotationDataFilteringServer(
    id = "gene_data_filter",
    df_filtered_annot = gene_filters$filtered_annot,
    df_data = loaded_data$loaded_data,
    selected_acc_id = gene_filters$selected_acc_id,
    selected_control_type = gene_filters$selected_control_type
  )

  # Filter controls and selected by gene
  ctrl_gene_data <- reactive({
    filter_by_gene(filtered_data$df_filtered_data_ctrl(), gene_selected())
  })
  selected_gene_data <- reactive({
    filter_by_gene(filtered_data$df_filtered_data_selected(), gene_selected())
  })

  #observe({
  #  print("DEBUG: gene_selected()")
  #  print(head(gene_selected()))
  #  print("DEBUG: ctrl_gene_data()")
  #  print(head(ctrl_gene_data()))
  #})
  
  # Pass the gene-filtered control data to the plot module
  gene_plot_Server(
    "gene_plot", 
    ctrl_data = ctrl_gene_data,
    selected_data = selected_gene_data,
    filtered_annot = gene_filters$filtered_annot
  )

  # Table outputs for controls and selected
  output$ctrl_gene_table <- renderTable({
    ctrl_gene_data()
  })
  output$selected_gene_table <- renderTable({
    selected_gene_data()
  })

  # Debug observe (optional)
  #observe({
  #  print("DEBUG: gene_filters$filtered_annot()")
  #  print(gene_filters$filtered_annot()[1:5, ])
  #  print("DEBUG: gene_filters$selected_individual_id()")
  #  print(gene_filters$selected_individual_id())    
  #  print("DEBUG: gene_selected()")
  #  print(gene_selected())
  #  print("Forcing evaluation of df_filtered_data_ctrl")
  #  print(filtered_data$df_filtered_data_ctrl()[2,])
  #  print("Forcing evaluation of df_filtered_data_selected")
  #  print(filtered_data$df_filtered_data_selected()[2,])
  #})

  # Housekeeper module (NEW, replaces old block)
  housekeeperServer(
    "housekeeper",
    filtered_data = filtered_data,
    gene_filters = gene_filters
  )
}

# Run the application 
shinyApp(ui = ui, server = server)
