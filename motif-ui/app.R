#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
Sys.setenv(RSTUDIO_PANDOC="/var/shiny-server/www/demo/pandoc")

source('convertSymbols.R')
source('subsetgenes.R')
source('genHeatmap.R')
source('sig_tf_matrix_motif.R')
library(shiny) #Loading required package: shiny
library(ggplot2) #Loading required package: plotly
library(gplots)
library(heatmaply) #Attaching package: ‘heatmaply’
library(shinyHeatmaply)
library(bsplus)
library(htmltools)
library(shinythemes)
library(RColorBrewer)

#source("census-app/helpers.R")
#counties <- readRDS("census-app/data/counties.rds")

# Define UI for data upload app ----
ui <- fluidPage(
  theme = shinytheme("cosmo"),
  
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
                           ".csv")),
      
      actionButton("runAnalysesButton","Button 1")
      
      
      
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Data file ----
      #textOutput("geneout"), #You chose
      #tableOutput("contents"),
     
      
      tabsetPanel(type = "tabs",
                  tabPanel("Get Started", 
                           h1(strong("Welcome to the tool")),
                           h4("The goal of this tool is take your differentially expressed genes of interest
                              and cluster these genes based on the motifs that are associated with each of them.
                              Our assumption is that genes that occur near similar motifs will be regulated by
                              similar transcription factors."),
                           h4("In the first tab, we will convert your gene IDs to gene symbols."),
                           h4("Afterwards, we will cluster your genes based on their motifs and output them to
                               a PCA plot and two heatmaps, one that is interactive, and another that is static."),
                
                           
                           h2(strong("Steps")),
                           h4("1. Select the genome of interest, hg19 or mm10"),
                           h4("2. Select the current Gene ID type of your data. We will convert them to 
                              gene symbols."),
                           h4("3. Upload your file of gene symbols"),
                           h4("4. Navigate to the other tabs to see visualizations of your data"),
                           tags$hr(),
                           
                           h4("Below is a text file of gene IDs you can use as input"),
                           downloadButton("downloadEx", "Download Example File")
                  
                           ),
                  tabPanel("Input Genes", 
                           fluidRow(
                             column(1,
                                    h4("Inputted"),
                                    tableOutput("contents"),
                                    downloadButton("downloadGenes", "Download Converted Genes")),
                             column(2,offset = 1, 
                                    h4("Converted"),
                                    tableOutput("changed"))
                           )
                  ),
                  tabPanel("Motif/Gene Table", uiOutput("mattab")),
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
  values <- reactiveValues()
  
  output$runAnalyses <- renderUI({
    if(is.null(input$file1))  
      return(NULL)
    actionButton("runAnalysesButton", "Run analyses")
  })
  
  data <- observeEvent(input$runAnalysesButton, {
    # The observeEvent takes no dependency on button 2, even though we refer to the input in the following line.
    input$button2  
    showModal(modalDialog(
      title = "Button pressed",
      "You pressed one of the buttons!"
    ))
    progressSig <- shiny::Progress$new()
    on.exit(progressSig$close())
    progressSig$set(message = "Reading in genes...", value = 0.10)
    ## Print original list
    #output$contents <- renderTable({
     
       # input$file1 will be NULL initially. After the user selects
       # and uploads a file, head of that data file by default,
       # or all rows if selected, will be shown.
       req(input$file1)
       df <- read.delim(input$file1$datapath,
                      #header = input$header,
                      sep = input$sep,
                      quote = input$quote,
                      stringsAsFactors = F)
       #return(head(df))
       values$origsym <- df
     
   #})
    progressSig$set(message = "Converting genes...", value = 0.20)
       
    ## Converted genes   
    newgenes <- convertGenes(input$file1$datapath, input$genome, input$idtype)
    values$convertedgenes <- newgenes
       
    # output$changed <- renderTable({
    #   if(is.null(input$file1))
    #     return(NULL)
    #   
    #   newgenes <- reactive(convertGenes(input$file1$datapath, input$genome, input$idtype))
    #   values$convertedgenes <- newgenes
    #   genetab <- head(newgenes())
    #   
    #   validate(
    #     need(nrow(genetab) > 0, "No genes found, make sure you inputted the correct gene ID type")
    #   )
    #   genetab
    #   
    #   #return(head(newgenes()))
    # })
    
    
    progressSig$set(message = "Converting gene IDs to gene matrix", value = 0.30)
    
    ## Convert to gene matrix (targetmatnum)
    genesInt <- values$convertedgenes
    mapmat <- genHeatmap(genesInt,input$genome) #should return table
    values$matobj <- mapmat[[1]]
    values$geneobj <- mapmat[[2]]
    
    
    progressSig$set(message = "Generating static heatmap", value = 0.40)
    
    
    ##Generate static heatmap
    
    
    ##Generate PCA
    progressSig$set(message = "Generating PCA", value = 0.50)
    genesPCA <- values$convertedgenes
    pcaList <- generatePCA(genesPCA,input$genome)
    values$pca <- pcaList
    

    #Gen graph
    genGraph()

    
  }) ##observe event bracket
  

   output$contents <- renderTable({
  
     # input$file1 will be NULL initially. After the user selects
     # and uploads a file, head of that data file by default,
     # or all rows if selected, will be shown.
     #req(values$origsym)
     origcont <- values$origsym
     head(origcont)
     
  #   req(input$file1)
  #   df <- read.delim(input$file1$datapath,
  #                  #header = input$header,
  #                  sep = input$sep,
  #                  quote = input$quote,
  #                  stringsAsFactors = F)
  #   return(head(df))
  #
   })


   output$changed <- renderTable({
     if(is.null(input$file1))
       return(NULL)
     newout <- values$convertedgenes
     head(newout)
     
  # 
  #   newgenes <- reactive(convertGenes(input$file1$datapath, input$genome, input$idtype))
  #   values$convertedgenes <- newgenes
  #   genetab <- head(newgenes())
  # 
  #   validate(
  #     need(nrow(genetab) > 0, "No genes found, make sure you inputted the correct gene ID type")
  #   )
  #   genetab
  # 
  #   #return(head(newgenes()))
   })
  
  output$downloadGenes <- downloadHandler(
    filename = function() {
      "convertedGenes.txt"
    },
    content = function(file) {
      finalgenes <- values$convertedgenes
      write.table(finalgenes(), file, row.names = FALSE, quote = F, sep="\t")
    }
  )
  
  output$geneout <- renderText({
    paste("You chose", input$genome)
  })
  
  output$mattab <- renderTable({
    if(is.null(input$file1))
      return(NULL)

    progressSig <- shiny::Progress$new()
    on.exit(progressSig$close())
    progressSig$set(message = "Generating Table of DE genes by Motif", value = 0.50)
    #genesInt <- values$convertedgenes
    #mapmat <- genHeatmap(genesInt(),input$genome) #should return table
    #values$matobj <- mapmat[[1]]
    #values$geneobj <- mapmat[[2]]
    mattable <- values$geneobj
    progressSig$set(message = "Rendering Table", value = 0.80)
    
    mattable


  })
  

  
  individualGraph <- reactive({
    if(is.null(input$file1))  
      return(NULL) 
    #newgenes <- convertGenes(input$file1$datapath, input$genome, input$idtype)
    # pcaDF <- generatePCA(input$file1$datapath)
    progressPCA <- shiny::Progress$new()
    on.exit(progressPCA$close())
    progressPCA$set(message = "Making PCA plot", value = 0)
    progressPCA$set(message = "Generating expression matrix", value = 0.50)
    genesPCA <- values$convertedgenes
    #validate(
    #  need(nrow(genesPCA) > 0, 'Check that genes is not empty')
    #)
    #validate(
    #  need(nrow(genetab) > 0, "No genes found, make sure you inputted the correct gene ID type")
    #)
    #pcaList <- generatePCA(genesPCA,input$genome)
    
    pcaOut <- values$pca
    pcaDF <- pcaOut[[1]]
    pc1var <- pcaOut[[2]]
    pc2var <- pcaOut[[3]]    
    
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
      geom_point(size = 3) +
      theme_classic() 
  })
  
  
  output$plot <- renderPlot({
    req(individualGraph())
    individualGraph()
  })
  
  output$save <- downloadHandler(
    filename = "save.png" ,
    content = function(file) {
      req(individualGraph())
      ggsave(file, plot = individualGraph(), device = 'png')
      while (!is.null(dev.list()))  dev.off()
      
      #dev.off()
    })
  

  output$heatmap <- renderPlotly({
     #if(is.null(input$file1))
     # return(NULL)
     progressHeat <- shiny::Progress$new()
     on.exit(progressHeat$close())
     progressHeat$set(message = "Making interactive heatmap", value = 0)
    
     genesInt <- values$convertedgenes
     #mapmat <- genHeatmap(genesInt(),input$genome) #should return table
     #values$matobj <- mapmat[[1]]
     targetnum <- values$matobj 
     #data("mtcars")
     progressHeat$set(message = "Plotting heatmap", value = 0.50)
     #cat("nrows mapmat:", nrow(mapmat))

     dendint <- values$dendstat
     #yb<-colorRampPalette(c("blue","white","red"))(100)
     #hmcol<-brewer.pal(11,"RdBu")
     
     #heatmaply(mtcars, cexRow = 0.5, cexCol = 0.3,
     heatmaply(targetnum, cexRow = 0.5, cexCol = 0.3,
               #colors = viridis(n = 256, alpha = 1, begin = 0, end = 1, option = "viridis"), 
               #colors = cm.colors(256),
               #colors = hmcol,
               colors=c("#009999", "#0000FF"), #blue and teal
               hclustfun = function(x) hclust(x, method="ward.D"),
               Colv = rev(dendint),
               seriate = "mean",
               row_dend_left = TRUE, 
               plot_method = "plotly")

   })
  
   ##currently displayed heatmap displays output from this heatmap
   output$statmap <- renderPlot({
     progressHeat <- shiny::Progress$new()
     on.exit(progressHeat$close())
     progressHeat$set(message = "Making heatmap matrix", value = 0)
     genesStat <- values$convertedgenes
     #mapmat <- genHeatmap(genesStat(),input$genome) #should return table
     targetnum <- values$matobj
     #statout <- values$stat
     statobj <- heatmap.2(as.matrix(targetnum), Rowv = T, Colv = T,
                          col = viridis(n = 256, alpha = 1, begin = 0, end = 1, option = "viridis"),
                trace = "none", #labRow = rownames(mapmat),
                labRow = rownames(targetnum),
                cexRow=0.5,
                cexCol=0.3,
                hclustfun=function(x) hclust(x, method="ward.D")
                #distfun = dist
               )
     statobj
     statdend <- statobj[["colDendrogram"]]
     values$dendstat <- statdend
     
   })



   
   
   output$heattab <- renderUI({
     tabsetPanel(id = "heattab", 
                 tabPanel("Static Heatmap",
                          tabPanel("Static Heatmap", plotOutput("statmap", width = "800px", height = "600px")),
                          downloadButton("downloadHeat", "Download")
                 ),
                 tabPanel("Interactive Heatmap", 
                          tabPanel("Heatmap", plotlyOutput("heatmap", width = "800px", height = "600px" )),
                          downloadButton("downloadHeatmaply", "Download")
                          
                 )
             
     )
   })
     
  
  indHeat <- reactive({
    if(is.null(input$file1))  
      return(NULL) 
    print("in ind heat")
    genesInd <- values$convertedgenes
    mapmat <- genHeatmap(genesInd(),input$genome) #should return table
    targetnum <- mapmat[[1]]
    data("mtcars")
    #heatmap.2(as.matrix(mtcars), Rowv = T, Colv = T, col = viridis(n = 256, alpha = 1, begin
    heatmap.2(as.matrix(targetnum), Rowv = T, Colv = T, col = viridis(n = 256, alpha = 1, begin
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
  
  
  indHeatmaply <- reactive({
    if(is.null(input$file1))  
      return(NULL) 
    print("in ind heat")
    indmap <- values$matobj
    
    data("mtcars")
    #heatmap.2(as.matrix(mtcars), Rowv = T, Colv = T, col = viridis(n = 256, alpha = 1, begin
    heatmaply(indmap, cexRow = 0.5, cexCol = 0.3,
              #colors = viridis(n = 256, alpha = 1, begin = 0, end = 1, option = "viridis"), 
              #colors = cm.colors(256),
              #colors = hmcol,
              colors=c("#009999", "#0000FF"), #blue and teal
              hclustfun = function(x) hclust(x, method="ward.D"),
              Colv = rev(dendint),
              seriate = "mean",
              row_dend_left = TRUE, 
              plot_method = "plotly",
              file = "heatmaply_plot.html"
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
  
  output$downloadHeatmaply <- ({
    downloadHandler(
      filename = function() { 'heatmaply.html' },
      content = function(outfile) {
        #req(indHeat())
        #png(file, width=1200, height=1200,res=300, pointsize=8)
        #pdf(file)
        #par(mar=c(10,2,2,3), cex=1.0)
        cat("setup heatmaply")
        maplyout <- values$matobj
        #indHeatmaply()
        dendint <- values$dendstat
        #dir.create("folder")
        #data("mtcars")
        tmp <- heatmaply(maplyout, cexRow = 0.5, cexCol = 0.3,
                   #colors = viridis(n = 256, alpha = 1, begin = 0, end = 1, option = "viridis"), 
                   #colors = cm.colors(256),
                   #colors = hmcol,
                   colors=c("#009999", "#0000FF"), #blue and teal
                   hclustfun = function(x) hclust(x, method="ward.D"),
                   Colv = rev(dendint),
                   seriate = "mean",
                   row_dend_left = TRUE, 
                   plot_method = "plotly",
                   file = outfile
        )
        #tmp

        #browseURL("folder/heatmaply_plot.html")
        cat("Finished in download")
        while (!is.null(dev.list()))  dev.off()
        
      })
  })

  #values$fdroutput <- 
  
  
  genGraph <- reactive({
    cat("in gen graph")
    genesNet <- values$convertedgenes
    
    if (is.null(input$decimal)) {
      cat("null cond")
      nets <- makenetgraph(genesNet, input$genome, 0.05)
      firstcall <- TRUE
      values$call <- firstcall
    }
    else {
      cat("else cond")
      cat("value of input$decimal", input$decimal ,"\n")
      nets <- makenetgraph(genesNet, input$genome, input$decimal)
      
    }
    values$netsobj <- nets[[1]]
    values$circsobj <- nets[[2]]
    values$sigsTF <- nets[[3]]
  })
  
  output$net  <- renderVisNetwork({
    if (is.numeric(input$decimal)){
    #if (input$decimal != 0.05){
      req(genGraph())
    }
    progressNet <- shiny::Progress$new()
    on.exit(progressNet$close())
    progressNet$set(message = "Making Network-Style Graph", value = 0.50)
    genesNet <- values$convertedgenes
    #networks <- makenetgraph(genesNet, input$genome, input$decimal)
    #values$graphobj <- networks
    progressNet$set(message = "Rendering Graph", value = 0.80)
    
    #networks[[1]]
    #networks[[2]]
    
    #values$netobj <- networks[[1]]
    
    netout <- values$netsobj
    netout
    
    #validate(
    #  need(nrow(networks[[3]]) > 0, "No DE genes based on parameters for fdrout try changing parameters")
    #)
  }) # created input$mynetwork_selected
  
  output$circ  <- renderVisNetwork({
    if (input$decimal != 0.05){
      req(genGraph())
    }
    progressCirc <- shiny::Progress$new()
    on.exit(progressCirc$close())
    progressCirc$set(message = "Making Circle-Style Graph", value = 0.50)
    genesCirc <- values$convertedgenes
    #networks <- makenetgraph(genesCirc, input$genome, input$decimal)
    progressCirc$set(message = "Rendering Graph", value = 0.80)
    #networks[[2]]
    #networks[[2]]
    #values$circobj <- networks[[2]]
    circout <- values$circsobj
    circout
    
    
  }) 
  
  output$nettab <- renderUI({
    tabsetPanel(id = "subTabPanel1", 
                tabPanel("Get Started",
                        h1(strong("TF-Gene Relationship Visualization Tool")),
                        h4("The goal of this tool is take your differentially expressed genes of interest
                           and display a network graph visualizing the transcription factors that are known to 
                          regulate these genes."),
                        h4("This visualization is based on TRRUST, a manually curated database of human and 
                           mouse transcriptional regulatory networks."),
                        h4("TRRUST can be found here (https://www.grnpedia.org/trrust/)."),
                        h4("You can subset the nodes you want to display based on the significance cutoff with the
                           slider below."),
                        h4("This score determines what particular transcription factors in the database 
                           interact with your set of genes based on Fisher's exact test"),
                        h4("Adjusting the slider will regenerate the network graph based on the significance
                           cutoff you specify."),
                        h3(strong("Reading the Network Graph")),
                        h4("The network graph can be visualized in either network-style or circle-style."),
                        h4("Green Nodes represent the differentially expressed genes that you inputted as
                           gene IDs in a text file."),
                        h4("Purple Nodes represent the transcription factors that regulate these genes."),
                        h4("The arrows represent the mode of regulation a transcription factor has on a particular
                           gene. Red means that the transcription factor activates the gene, Blue means that
                           the transcription factor represses the gene. Yellow means that the mode of regulation is 
                           unknown. Grey means that there are multiple modes of regulation."),
                        h4("Hovering over transcription factor displays the genes that the transcription factor
                           regulates and their mode of regulation, the PMIDs describing this relationship in 
                           literature, the number of genes for which this TF regulates in your particular gene set,
                           and the significance value of the overlap between the TF and your particular gene set."),
                        h4("This information is also presented in a table which you can download under the
                           TF table tab"),
                        h3(strong("References")),
                        h4("1. TRRUST v2: an expanded reference database of human and mouse 
                            transcriptional regulatory interactions. Nucleic Acids Research 26 Oct, 2017
                          ")
                        
                        
                        
                        
                     
                ),
                tabPanel("Network-Style", 
                         visNetworkOutput("net",  width = "900px", height = "500px"),
                         downloadButton("downloadNet", "Download")
                         
                         
                ),
                tabPanel("Circle-Style",
                         visNetworkOutput("circ",  width = "900px", height = "500px"),
                         downloadButton("downloadCirc", "Download")
                         
                ),
                tabPanel("TF Table",
                         tableOutput("sigTFs"),
                         downloadButton("downloadsigTFs", "Download")
                         
                ),
                sliderInput("decimal", "Significance Cutoff:",
                            min = 0, max = 1,
                            value = 0.05, step = 0.05)
    )
    
    
    
    
  })
  
  output$sigTFs <- renderTable({
    if (input$decimal != 0.05){
      req(genGraph())
    }
    
    if(is.null(input$file1))  
      return(NULL) 
    
    progressSig <- shiny::Progress$new()
    on.exit(progressSig$close())
    progressSig$set(message = "Generating Table of Significant Transcription Factors", value = 0.50)
    genesSig <- values$convertedgenes
    #networks <- makenetgraph(genesSig, input$genome, input$decimal)
    progressSig$set(message = "Rendering Table", value = 0.80)
    
    #networks[[3]]
    #values$sigTFobj <- networks[[3]]
    
    #values$sigTFobj <- networks[[3]]
    sigout <- values$sigsTF
    sigout
   
  })
 # created input$mynetwork_selected
  
  #visSave(network, file = "network.html")


  output$downloadNet <- ({
    downloadHandler(
      filename = function() { 'network.html' },
      content = function(file) {
        cat("setup net")
        finalnet <- values$netobj
        visSave(finalnet, file)
        cat("Finished in download net")
        
      })
  })
  
  output$downloadCirc <- ({
    downloadHandler(
      filename = function() { 'circ.html' },
      content = function(file) {
        cat("setup net")
        finalcirc <- values$circobj
        visSave(finalcirc, file)
        cat("Finished in download net")
        
      })
  })
  
  output$downloadsigTFs <- downloadHandler(
    filename = function() {
      "sigTFs.txt"
    },
    content = function(file) {
      finalsigTFs <- values$sigTFobj
      write.table(finalsigTFs, file, row.names = FALSE, quote = F, sep="\t")
    }
  )
  
  output$downloadEx <- downloadHandler(
    filename = function() {
      "ex_DE_genes.txt"
    },
    content = function(file) {
      exfile <- read.delim("data/33geneuni.txt", stringsAsFactors=FALSE)
      write.table(exfile, file, row.names = FALSE, quote = F, sep="\t")
      
    }
  )
}
# Run the app ----
shinyApp(ui, server)

