#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
source('convertSymbols.R')
library(shiny)

#source("census-app/helpers.R")
#counties <- readRDS("census-app/data/counties.rds")

# Define UI for data upload app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Uploading Files"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      
      # Input: Checkbox if file has header ----
      #checkboxInput("header", "Header", TRUE),
      
      # Input: Select separator ----
      radioButtons("sep", "Separator",
                   choices = c(Comma = ",",
                               Semicolon = ";",
                               Tab = "\t",
                               Whitespace = " "),
                   selected = " "),
      
      # Input: Select quotes ----
      radioButtons("quote", "Quote",
                   choices = c(None = "",
                               "Double Quote" = '"',
                               "Single Quote" = "'"),
                   selected = ""),
      
      # Horizontal line ----
      tags$hr(),
      
      # Input: Select number of rows to display ----
      radioButtons("disp", "Display",
                   choices = c(Head = "head",
                               All = "all"),
                   selected = "head"),
      
      selectInput("genome", "Select Genome:",
                  c("hg19","mm10")),
      
      selectInput("idtype", "Select Gene ID type:",
                  c("geneid","unigene","refseq","ensembl","name")),
      # Horizontal line ----
      tags$hr(),
      # Input: Select a file ----
      fileInput("file1", "Upload a list of genes",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv"))
      
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Data file ----
      textOutput("geneout"),
      tableOutput("contents"),
      tableOutput("changed"),
      
      tabsetPanel(type = "tabs",
                  tabPanel("PCA", plotOutput("plot")),
                  tabPanel("Heatmap", verbatimTextOutput("summary")),
                  tabPanel("Network Graph", tableOutput("table"))
      )
    )
    
  )
)

# Define server logic to read selected file ----
server <- function(input, output) {
  
  output$contents <- renderTable({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    
    req(input$file1)
    
 
    
    df <- read.delim(input$file1$datapath,
                   #header = input$header,
                   sep = input$sep,
                   quote = input$quote,
                   stringsAsFactors = F)
    
    if(input$disp == "head") {
      return(head(df))
    }
    else {
      return(df)
    }
    

  

    
  })
  output$changed <- renderTable({
    if(is.null(input$file1))  
      return(NULL) 
    
    newgenes <- convertGenes(input$file1$datapath, input$genome, input$idtype)
    return(newgenes)
  })
  output$geneout <- renderText({
    paste("You chose", input$genome)
  })
  
}
# Run the app ----
shinyApp(ui, server)

