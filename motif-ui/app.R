#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
source('convertSymbols.R')
source('subsetGenes.R')
source('genHeatmap.R')
library(shiny)
library(ggplot2)
library(gplots)

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
      #tableOutput("contents"),
     
      
      tabsetPanel(type = "tabs",
                  tabPanel("Input Genes", tableOutput("contents"), tableOutput("changed")),
                  tabPanel("PCA", plotOutput("plot"),
                            downloadButton("save", "Download")),
                  #tabPanel("Heatmap", tableOutput("heatmap")),
                  tabPanel("Heatmap", plotOutput("heatmap"),
                           downloadButton("downloadData", "Download")),
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
    
    newgenes <<- reactive(convertGenes(input$file1$datapath, input$genome, input$idtype))
    return(head(newgenes()))
  })
  
  output$geneout <- renderText({
    paste("You chose", input$genome)
  })
  
  
  
  individualGraph <- reactive({
    if(is.null(input$file1))  
      return(NULL) 
    #newgenes <- convertGenes(input$file1$datapath, input$genome, input$idtype)
    # pcaDF <- generatePCA(input$file1$datapath)
    pcaDF <- generatePCA(newgenes(),input$genome)
    
    #ggplot(pcaDF)
    groupnum <- rownames(pcaDF)
    ## new plot with 150
    xmax <- round(max(pcaDF$PC1))+1
    xmin <- round(min(pcaDF$PC1))-1
    
    ymax <- round(max(pcaDF$PC2))+1
    ymin <- round(min(pcaDF$PC2))-1
    ggplot(pcaDF,aes(PC1 ,PC2)) +
      ggtitle("Genes Grouped by Common Motifs") +
      xlab(paste("PC1",round(pc1var * 100, 0), "%"))+
      ylab(paste("PC2",round(pc2var * 100, 0), "%")) +
      xlim(xmin, xmax) + 
      ylim(ymin,ymax) +
      #geom_text(label= groupnum, size = 3, nudge_x = 1 ,nudge_y = 2) +
      #scale_fill_discrete(name = "Grouping") + 
      geom_point(size = 3) +
      #geom_step() +
      #scale_color_manual(values = c("#7FC97F","#F0027F","#386CB0"))
      theme_classic() 
    # dev.off()
  })
  
  
  output$plot <- renderPlot({
    req(individualGraph())
    individualGraph()
   #  if(is.null(input$file1))  
   #    return(NULL) 
   #  #newgenes <- convertGenes(input$file1$datapath, input$genome, input$idtype)
   # # pcaDF <- generatePCA(input$file1$datapath)
   #  pcaDF <- generatePCA(newgenes(),input$genome)
   #  
   #  #ggplot(pcaDF)
   #  groupnum <- rownames(pcaDF)
   #  ## new plot with 150
   #  xmax <- round(max(pcaDF$PC1))+1
   #  xmin <- round(min(pcaDF$PC1))-1
   #  
   #  ymax <- round(max(pcaDF$PC2))+1
   #  ymin <- round(min(pcaDF$PC2))-1
   #  ggplot(pcaDF,aes(PC1 ,PC2)) +
   #     ggtitle("Genes Grouped by Common Motifs") +
   #     xlab(paste("PC1",round(pc1var * 100, 0), "%"))+
   #     ylab(paste("PC2",round(pc2var * 100, 0), "%")) +
   #     xlim(xmin, xmax) + 
   #     ylim(ymin,ymax) +
   #     #geom_text(label= groupnum, size = 3, nudge_x = 1 ,nudge_y = 2) +
   #     #scale_fill_discrete(name = "Grouping") + 
   #     geom_point(size = 3) +
   #     #geom_step() +
   #     #scale_color_manual(values = c("#7FC97F","#F0027F","#386CB0"))
   #     theme_classic() 
   #  # dev.off()
  })
  
  output$save <- downloadHandler(
    filename = "save.png" ,
    content = function(file) {
      #ggsave(p(), filename = file)
      req(individualGraph())
      ggsave(file, plot = individualGraph(), device = 'png')
      dev.off()
    })
  # output$heatmap <- renderTable({
  #   if(is.null(input$file1))  
  #     return(NULL) 
  #   map <- genHeatmap(newgenes(),input$genome) #should return table
  #   return(map)
  # })
  
  output$heatmap <- renderPlot({
    if(is.null(input$file1))  
      return(NULL) 
    mapmat <<- genHeatmap(newgenes(),input$genome) #should return table
    
    heatmap.2(as.matrix(mapmat), Rowv = T, Colv = T, col = heat.colors, 
              trace = "none", labRow = rownames(targetmatnmum),
              #lhei = c(0.5,5),     
              #lhei = c(0.5,1),     
              #lwid = c(0.5,0.5),
              cexRow=0.15,
              cexCol=0.15,
              hclustfun=function(x) hclust(x, method="ward.D"))
             # invisible(dev.off())
              cat("Heatmap finished!\n")
              
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      #"out.txt"
      gsub(".*\\.","outputgenes.",input$file1)
      #paste(input$file1, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(newgenes(), file, row.names = FALSE, col.names = F, quote = F)
    }
  )

  
  
  
}
# Run the app ----
shinyApp(ui, server)

