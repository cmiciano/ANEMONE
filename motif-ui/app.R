#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
source('convertSymbols.R')
source('subsetGenes.R')
source('genHeatmap.R')
library(shiny)
library(ggplot2)
library(gplots)
library(heatmaply)
library(shinyHeatmaply)


#source("census-app/helpers.R")
#counties <- readRDS("census-app/data/counties.rds")

# Define UI for data upload app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Motif/Gene Clustering Tool"),
  
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
                  tabPanel("Get Started", 
                           h2("Welcome to the tool, put your genes of interest in a text file"),
                           h3("1. Select the genome of interest, hg19 or mm10"),
                           h3("2. Select the current Gene ID type of your data. We will convert it."),
                           h3("3. Upload your file of gene symbols"),
                           h3("4. Navigate to the other tabs to see visualizations of your data")),
                  tabPanel("Input Genes", tableOutput("contents"), tableOutput("changed")),
                  tabPanel("PCA", plotOutput("plot"),
                            downloadButton("save", "Download")),
                  #tabPanel("Heatmap", tableOutput("heatmap")),
                  tabPanel("Heatmap", plotlyOutput("heatmap", width = "800px", height = "800px" ),
                           downloadButton("downloadData", "Download")),
                  tabPanel("Static Heatmap", plotOutput("statmap", width = "800px", height = "800px"),
                           downloadButton("downloadHeat", "Download")),
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
    progressPCA <- shiny::Progress$new()
    on.exit(progressPCA$close())
    progressPCA$set(message = "Making PCA plot", value = 0)
    
    
    # for (i in 1:15) {
    #   progressPCA$set(value = i)
    #   Sys.sleep(0.5)
    # }
    progressPCA$set(message = "Generating expression matrix", value = 0.50)
    pcaList <- generatePCA(newgenes(),input$genome)
    
    pcaDF <- pcaList[[1]]
    pc1var <- pcaList[[2]]
    pc2var <- pcaList[[3]]    
    
    progressPCA$set(message = "Returning plot output", value = 0.9)
    
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

  output$heatmap <- renderPlotly({
     if(is.null(input$file1))
       return(NULL)
    
     progressHeat <- shiny::Progress$new()
     on.exit(progressHeat$close())
     progressHeat$set(message = "Making heatmap matrix", value = 0)
    
     mapmat <<- genHeatmap(newgenes(),input$genome) #should return table
     data("mtcars")
     progressHeat$set(message = "Plotting heatmap", value = 0.50)
     cat("nrows mapmat:", nrow(mapmat))
     #mapmat
     #submat <- mapmat[1:5,1:5]
     #heatmaply(submat)
     #heatmaply(submat, labRow = NA, labCol = NA)
     heatmaply(mapmat, cexRow = 0.5, cexCol = 0.3, hclust_method = "ward.D2")
     
     #heatmaply(mapmat, cexRow = 0.5, cexCol = 0.3, file = "mapmat.pdf")
     #heatmaply(mtcars)
     
     #heatmaply(mtcars, xlab = "Features", ylab = "Cars") 
               ##WORKS
     # heatmap.2(as.matrix(mapmat), Rowv = T, Colv = T, col = heat.colors,
     #           trace = "none", labRow = rownames(targetmatnmum),
     #           #lhei = c(0.5,5),
     #           #lhei = c(0.5,1),
     #           #lwid = c(0.5,0.5),
     #           cexRow=0.15,
     #           cexCol=0.15,
     #           hclustfun=function(x) hclust(x, method="ward.D"))
     #          # invisible(dev.off())
     #           cat("Heatmap finished!\n")
  
   })
  
  output$statmap <- renderPlot({
    progressHeat <- shiny::Progress$new()
    on.exit(progressHeat$close())
    progressHeat$set(message = "Making heatmap matrix", value = 0)
    
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

  
  indHeat <- reactive({
    if(is.null(input$file1))  
      return(NULL) 
    print("in ind heat")
    
    mapmat <<- genHeatmap(newgenes(),input$genome) #should return table
   
    heatmap.2(as.matrix(mapmat), Rowv = T, Colv = T, col = heat.colors, 
              trace = "none", labRow = rownames(targetmatnmum),
              #lhei = c(0.5,5),     
              #lhei = c(0.5,1),     
              #lwid = c(0.5,0.5),
              cexRow=0.15,
              cexCol=0.15,
              hclustfun=function(x) hclust(x, method="ward.D"))
     invisible(dev.off())
    cat("Heatmap finished!\n")
  })
  
  output$downloadHeat <- ({
    downloadHandler(
      filename = function() { 'heatmap.pdf' },
      content = function(file) {
      #png(file, width=1200, height=1200,res=300, pointsize=8)
      pdf(file)
      par(mar=c(10,2,2,3), cex=1.0)
      cat("setup pdf")
      print(indHeat())
      cat("Finished in download")
      dev.off()
        
      })
  })
  # output$heatmap <- renderPlot({
  #   req(indHeat())
  #   indHeat()
  # })
  # 
  output$downloadData <- downloadHandler(
    filename = function() {
      #"out.txt"
      gsub(".*\\.","outputgenes.",input$file1)
      #paste(input$file1, ".csv", sep = "")
    },
    content = function(file) {
      pdf("genemotifctint.pdf", 8, 15)
      par(mar=c(10,2,2,3), cex=1.0)
      req(indHeat())
      indHeat()
      
      dev.off()
    }
  )

  
  
  
}
# Run the app ----
shinyApp(ui, server)

