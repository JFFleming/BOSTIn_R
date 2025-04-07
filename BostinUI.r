library(shiny)
library(phangorn)
library(Biostrings)
library(reticulate)

# Define UI for data upload app ----
ui <- fluidPage(

  # App title ----
  titlePanel("BOSTIN: Broad Overview of Sequence and Topology INcongruence"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(
      img(src="BostinLogo.png",height=150,width=150,style="display: block; margin-left: auto; margin-right: auto; margin-top: 10px;"),
      tags$hr(),
      # Input: Select a file ----
      fileInput("file1", "Choose Fasta File",
                multiple = FALSE,
                accept = c("fasta",".fasta","fst", "faa","fna",".fa",".fas")),

      # Horizontal line ----
      tags$hr(),

      # Input: Checkbox if file has header ----
      checkboxInput(inputId = "comphet", "Compositional Heterogeneity: nRCFV", TRUE),
      checkboxInput(inputId = "sitesat", "Site Saturation: DE-Score for Amino Acids, C-Score for Nucleotides", TRUE),
      checkboxInput(inputId = "branchhet", "Branch Length Heterogeneity: LB-Score", TRUE),
      
      # Input: Select separator ----
      radioButtons(inputId = "d_type", "DataType",
                   choices = c(AminoAcid = "AA",
                               DNA = "DNA"),
                   selected = "AA"),
      actionButton("do","Submit"),
      # Horizontal line ----
      tags$hr(),

    ),

    # Main panel for displaying outputs ----
   mainPanel(tabsetPanel(
      id = "tabsetPanelID",
      type = "tabs",
      tabPanel("Compositional Heterogeneity", tabsetPanel( 
        tabPanel("nRCFV Summary", tableOutput("rcfv")),
        tabPanel("ntsRCFV", tableOutput("tsrcfv")), 
        tabPanel("ntsRCFV Histogram", plotOutput("tshist")),
        tabPanel("ncsRCFV", tableOutput("csrcfv")),
        tabPanel("ncsRCFV Histogram", plotOutput("cshist")),
        tabPanel("Bostin", uiOutput("ch_bostin"))
      )),
      tabPanel("Site Saturation", uiOutput("dynamic_tabs")),
      
      tabPanel("Branch Length Heterogeneity", tabsetPanel(
        tabPanel("NJ Tree", plotOutput("treeplot")),
        tabPanel("LB Score Histogram", plotOutput("histplot")),
        tabPanel("LB Score Summary", tableOutput("summary")),
        tabPanel("LBi Scores", tableOutput("list")),
        tabPanel("Bostin", uiOutput("lb_bostin"))
      ))
    ))


  )
)
