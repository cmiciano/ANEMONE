#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
source('convertSymbols.R')
source('subsetGenes.R')
source('genHeatmap.R')
source('sig_tf_matrix_motif.R')
library(shiny)
library(ggplot2)
library(gplots)
library(heatmaply)
library(shinyHeatmaply)
library(bsplus)
library(htmltools)

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
      radioButtons("sep", "Separator:",
      #radioButtons("sep", "Specify what delimiter your genes are separated by:",
                                
                   choices = c(Comma = ",",
                               Semicolon = ";",
                               Tab = "\t",
                               Whitespace = " "),
                   selected = " "),
      
      # Input: Select quotes ----
      radioButtons("quote", "Quote:",
      #radioButtons("quote", "Specify if there are quotes surrounding your gene IDs",
                   choices = c(None = "",
                               "Double Quote" = '"',
                               "Single Quote" = "'"),
                   selected = ""),
      
      # Horizontal line ----
      tags$hr(),
      
      # Input: Select number of rows to display ----
      #radioButtons("disp", "Specify whether to display head of converted genes or all converted genes",
      #radioButtons("disp", "Display Head/All",
      #             choices = c(Head = "head",
      #                         All = "all"),
      #             selected = "head"),
      
      
      selectInput("genome", "Select Genome:",
                  c("hg19","mm10"))
      %>%
        shinyInput_label_embed(
          shiny_iconlink() %>%
            bs_embed_popover(
              title = "Select what genome your genes belong to", content = "Choose a favorite", placement = "left"
            )
        ),
      
      selectInput("idtype", "Select Gene ID type:",
                  c("geneid","unigene","refseq","ensembl","name"))
      %>%
        shinyInput_label_embed(
          shiny_iconlink() %>%
            bs_embed_popover(
              title = "Select what Gene ID type your genes are written as", content = "Choose a favorite", placement = "left"
            )
        ),
      # Horizontal line ----
      tags$hr(),
      # Input: Select a file ----
      fileInput("file1", "Upload a list of genes as a text file",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv"))
      
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Data file ----
      #textOutput("geneout"), #You chose
      #tableOutput("contents"),
     
      
      tabsetPanel(type = "tabs",
                  tabPanel("Get Started", 
                           h2("Welcome to the tool, put your genes of interest in a text file"),
                           h4("The goal of this tool is take your differentially expressed genes of interest
                              and cluster these genes based on the motifs that are associated with each of them.
                              Our assumption is that genes that occur near similar motifs will be regulated by
                              similar transcription factors."),
                           h4("In the first tab, we will convert your gene IDs to gene symbols."),
                           h4("Afterwards, we will cluster your genes based on their motifs and output them to
                               a PCA plot and two heatmaps, one that is interactive, and another that is static."),
                           h4("Finally, we can visualize what transcription factors are regulating these genes
                              with a network graph."),
                           
                           
                           h2("Steps"),
                           h4("1. Select the genome of interest, hg19 or mm10"),
                           h4("2. Select the current Gene ID type of your data. We will convert it."),
                           h4("3. Upload your file of gene symbols"),
                           h4("4. Navigate to the other tabs to see visualizations of your data")),
                  tabPanel("Input Genes", 
                           fluidRow(
                             column(1,tableOutput("contents")),
                             column(2,offset = 1, tableOutput("changed"))
                           )
                  ),
                  tabPanel("PCA", plotOutput("plot"),
                            downloadButton("save", "Download")),
                  #tabPanel("Heatmap", tableOutput("heatmap")),
                  #tabPanel("Heatmap", plotlyOutput("heatmap", width = "800px", height = "800px" )),
                  #tabPanel("Static Heatmap", plotOutput("statmap", width = "800px", height = "800px"),
                  #         downloadButton("downloadHeat", "Download")),
                  tabPanel("Heatmaps", uiOutput("heattab")),
                  tabPanel("Network Graph", uiOutput("nettab"))
                           #          visNetworkOutput("net",  width = "700px", height = "700px"),
                           #          downloadButton("downloadNet", "Download"))
                  # tabPanel("Network Graph", 
                  #          visNetworkOutput("net",  width = "700px", height = "700px"),
                  #          downloadButton("downloadNet", "Download")),
                  # tabPanel("Network Graph", 
                  #          visNetworkOutput("circ",  width = "700px", height = "700px"))

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
    return(head(df))
    
    #if(input$disp == "head") {
    #  return(head(df))
    #}
    #else {
    #  return(df)
    #}
    
  })
  

  output$changed <- renderTable({
    if(is.null(input$file1))  
      return(NULL) 
    
    newgenes <<- reactive(convertGenes(input$file1$datapath, input$genome, input$idtype))
    genetab <- head(newgenes())
    
    validate(
      need(nrow(genetab) > 0, "No genes found, make sure you inputted the correct gene ID type")
    )
    genetab
      
    #return(head(newgenes()))
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
      #xlim(-10, 10) + 
      #ylim(-5,5) +
      #geom_text(label= groupnum, size = 3, nudge_x = 1 ,nudge_y = 2) +
      #scale_fill_discrete(name = "Grouping") + 
      geom_point(size = 3) +
      #geom_step() +
      #scale_color_manual(values = c("#7FC97F","#F0027F","#386CB0"))
      theme_classic() 
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
      while (!is.null(dev.list()))  dev.off()
      
      #dev.off()
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
     #heatmaply(mtcars, cexRow = 0.5, cexCol = 0.3,
               
     heatmaply(mapmat, cexRow = 0.5, cexCol = 0.3,
               hclustfun = function(x) hclust(x, method="ward.D"),
               #distfun = dist,
               seriate = "mean",
               row_dend_left = TRUE, 
               plot_method = "plotly")
               #hclust_method = "ward.D2")
  
   })
  
   ##currently displayed heatmap displays output from this heatmap
   output$statmap <- renderPlot({
     progressHeat <- shiny::Progress$new()
     on.exit(progressHeat$close())
     progressHeat$set(message = "Making heatmap matrix", value = 0)
     
     mapmat <<- genHeatmap(newgenes(),input$genome) #should return table
     data("mtcars")
     #indHeat()
     #heatmap.2(as.matrix(mtcars), Rowv = T, Colv = T,  col = viridis(n = 256, alpha = 1, begin
     heatmap.2(as.matrix(mapmat), Rowv = T, Colv = T,  col = viridis(n = 256, alpha = 1, begin
                                                                    = 0, end = 1, option = "viridis"),
                trace = "none", labRow = rownames(mapmat),
                #trace = "none", labRow = rownames(mtcars),
    
                #lhei = c(0.5,5),
                #lhei = c(0.5,1),
                #lwid = c(0.5,0.5),
                cexRow=0.5,
                cexCol=0.3,
                hclustfun=function(x) hclust(x, method="ward.D")
                #distfun = dist
               )
               # invisible(dev.off())
        #        cat("Heatmap finished!\n")
   })

   output$heattab <- renderUI({
     tabsetPanel(id = "heattab", 
                 tabPanel("Interactive Heatmap", 
                          tabPanel("Heatmap", plotlyOutput("heatmap", width = "800px", height = "800px" ))
                 ),
                 tabPanel("Static Heatmap",
                          tabPanel("Static Heatmap", plotOutput("statmap", width = "800px", height = "800px")),
                          downloadButton("downloadHeat", "Download")
                 )
     )
   })
     
  
  indHeat <- reactive({
    if(is.null(input$file1))  
      return(NULL) 
    print("in ind heat")
    
    mapmat <<- genHeatmap(newgenes(),input$genome) #should return table
    data("mtcars")
    #heatmap.2(as.matrix(mtcars), Rowv = T, Colv = T, col = viridis(n = 256, alpha = 1, begin
    heatmap.2(as.matrix(mapmat), Rowv = T, Colv = T, col = viridis(n = 256, alpha = 1, begin
                                                                            = 0, end = 1, option = "viridis"), 
                trace = "none", labRow = rownames(mapmat),
              # trace = "none", labRow = rownames(mtcars),
    
              #lhei = c(0.5,5),     
              #lhei = c(0.5,1),     
              #lwid = c(0.5,0.5),
              cexRow=0.5,#0.15
              cexCol=0.3,#0.15
              hclustfun=function(x) hclust(x, method="ward.D")
              )
     #invisible(dev.off())
    cat("Heatmap finished!\n")
  })
  
  output$downloadHeat <- ({
    downloadHandler(
      filename = function() { 'heatmap.pdf' },
      content = function(file) {
      #req(indHeat())
      #png(file, width=1200, height=1200,res=300, pointsize=8)
      pdf(file)
      par(mar=c(10,2,2,3), cex=1.0)
      cat("setup pdf")
      print(indHeat())
      cat("Finished in download")
      while (!is.null(dev.list()))  dev.off()
      
      })
  })
  # output$statmap <- renderPlot({
  #   progressHeat <- shiny::Progress$new()
  #   on.exit(progressHeat$close())
  #   progressHeat$set(message = "Making heatmap matrix", value = 0)
  #   #   
  #   cat("startind")
  #   req(indHeat())
  #   indHeat()
  #   cat("finind")
  #   
  # })
   
  # output$downloadData <- downloadHandler(
  #   filename = function() {
  #     #"out.txt"
  #     gsub(".*\\.","outputgenes.",input$file1)
  #     #paste(input$file1, ".csv", sep = "")
  #   },
  #   content = function(file) {
  #     pdf("genemotifctint.pdf", 8, 15)
  #     par(mar=c(10,2,2,3), cex=1.0)
  #     req(indHeat())
  #     indHeat()
  #     
  #     dev.off()
  #   }
  # )
  

  
  output$net  <- renderVisNetwork({
  #  ... visOptions(nodesIdSelection = TRUE)
    #visout <- visNetwork(nodes, edges, width = "100%")
    #visout
    networks <<- makenetgraph(newgenes(), input$genome, 0.05)
    networks[[1]]
    #networks[[2]]
  }) # created input$mynetwork_selected
  
  output$circ  <- renderVisNetwork({
    #  ... visOptions(nodesIdSelection = TRUE)
    #visout <- visNetwork(nodes, edges, width = "100%")
    #visout
    cat("In circ")
    networks <<- makenetgraph(newgenes(), input$genome, 0.05)
    networks[[2]]
    #networks[[2]]
  }) 
  
  output$nettab <- renderUI({
    tabsetPanel(id = "subTabPanel1", 
                tabPanel("Network-Style", 
                         visNetworkOutput("net",  width = "700px", height = "500px")#,
                         
                ),
                tabPanel("Circle-Style",
                         visNetworkOutput("circ",  width = "700px", height = "500px")#,
                         #downloadButton("downloadNet", "Download")
                         
                )
    )
    
    
    
  })
 # created input$mynetwork_selected
  
  #visSave(network, file = "network.html")


  output$downloadNet <- ({
    downloadHandler(
      filename = function() { 'network.html' },
      content = function(file) {
        cat("setup net")
        visSave(network, file)
        cat("Finished in download net")
        while (!is.null(dev.list()))  dev.off()
        
      })
  })
  
  
}
# Run the app ----
shinyApp(ui, server)

