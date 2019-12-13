##!/usr/bin/env Rscript

###############
#Makes a graph of nodes of TFs (in red) and Target genes (in blue) from DE genes 
#######
suppressMessages(library(GeneOverlap))
suppressMessages(library(visNetwork))
suppressMessages(library(igraph))
suppressMessages(library(RColorBrewer))

#Node colors 
green <- "#4DAF4A"
purple <-  "#984EA3"
orange <-  "#FF7F00"
blue <- "#377EB8"
red <- "#E41A1C"
black <-  "#999999"


  
  #Parameters for gene list
#if (infile_type == "list") {
#  if (length(args)!=5) {
#    print("1.  input file <file.txt>")
#    print("2.  specify if input file is homer output or gene list <homer | list>")
#    print("3.  specify if species human or mouse <human | mouse >")
#    print("4.  fdr value for output file <0.05>")
#    print("5.  specify if ID type <geneid | unigene | refseq | ensembl | name>")
#    
#    stop()
#  }
#  fdroutput=args[4]
#  fdroutput=as.numeric(fdroutput)
#  idtype = args[5] #troubleshooting
#  cat("FDRoutput: ", fdroutput, "\n")
#}

cat("\n")
cat("#######################################################################\n")
cat("\n")
######### Reading gene list #######

#cat("Reading in gene list:", infile_name, "\n")
#file <- read.delim("motif-ui/outputgenes.txt", header = F, stringsAsFactors = F)
  
#inputgenesdata <- read.delim("motif-ui/outputgenes.txt", header = F, stringsAsFactors = F)
#inputgenesdata <- read.delim(infile_name, header = F, stringsAsFactors = F)
#names(inputgenesdata) <- "gene_symbol"
  
#inputgenesdata <- data.frame(inputgenesdata)
#names(inputgenesdata) <- "gene_symbol"
  

##################################
### Function to make network graph
### Parameters: list of genes, TF database, name of the file to be outputted, and input file type
makenetgraph<-function(genes, dbinput, fdroutput) {
  # #Get number of unique target genes for a unique TF i column
  green <- "#4DAF4A"
  purple <-  "#984EA3"
  orange <-  "#FF7F00"
  blue <- "#377EB8"
  red <- "#E41A1C"
  black <-  "#999999"
  
  
  #inputgenesdata <- genes
  inputgenesdata <- data.frame(genes)
  names(inputgenesdata) <- "gene_symbol"
  
  
  cat("\n")
  cat("#######################################################################\n")
  cat("\n")
  #cat("Making network graph for:", fname, infile_name, "\n")
  cat("Making network graph","\n")
  
  
  if (dbinput == "mm10") {
    cat("Using mm10 genome\n")
    #generef <- read.delim("motif-genes-app/motif-ui/data/mm10_GeneIDsymbol.tsv", 
    #                      header=TRUE, stringsAsFactors = F) #for troubleshooting
    db <- read.delim("data/trrust_rawdata.mouse.txt", 
                            header=FALSE, stringsAsFactors = F)
    generef <- read.delim("data/mm10_GeneIDsymbol.tsv", 
                          header=TRUE, stringsAsFactors = F)
    
  } else if(dbinput == "hg19") {
    cat("Using hg19 genome\n")
    #generef <- read.delim("motif-genes-app/motif-ui/data/hg19_GeneIDsymbol.tsv", 
    #                      header=TRUE, stringsAsFactors = F)
    db <- read.delim("data/trrust_rawdata.human.txt", 
                            header=FALSE, stringsAsFactors = F)
    generef <- read.delim("data/hg19_GeneIDsymbol.tsv", 
                          header=TRUE, stringsAsFactors = F)
    
  } else {
    print("Not a valid genome")
  }
  
  colnames(db) <- c("TF","Target", "Mode of Regulation", "References (PMID)")
  db$uniq_target_ct <- sapply(db$TF, function(x) length(unique(db[db$TF == x ,]$Target)) )
  
  
  # Subsets TFs that have at least 5 target genes in the raw database
  
  great5tg <- subset(db, db$uniq_target_ct >= 5)
  cat("Number of TFs with at least 5 target genes:", length(unique(great5tg$TF)), "\n") #241 
  
  #Subsets input genes if they're in the db and have at least 5 target genes
  check_in_db <- subset(inputgenesdata, inputgenesdata$gene_symbol %in% great5tg$Target == TRUE)
  query_in_db <- subset(great5tg, great5tg$Target %in% check_in_db$gene_symbol )
  cat("Number of query target genes in database:", length(unique(query_in_db$Target)), "\n") #241 
  
  if (nrow(query_in_db) == 0) {
    print("DE genes not found in database")
    cat("\n")
    return()
  }
  cat("Finding number of overlapped genes...\n")
  cat("\n")
  query_in_db$overlap_ct  <- sapply(query_in_db$TF, function(x) length(unique(query_in_db[query_in_db$TF == x,]$Target)))
  great1tg <- subset(query_in_db, query_in_db$overlap_ct > 1)
  cat("Number of TFs that overlap with more than one query target gene:", length(unique(great1tg$TF)), "\n") #241 
  
  if (nrow(great1tg) == 0) {
    print("No TFs found with more than one overlap with given DE genes")
    cat("\n")
    return()
  }
  great1tg_noct <- great1tg[, -which(names(great1tg) %in% c("uniq_target_ct", "overlap_ct"))]
  # Create list of TF factors that overlap with more than one query target gene
  final_list <- lapply(unique(great1tg_noct$TF),
                       function(x) subset(great1tg_noct, great1tg_noct$TF == x))
  # #2367 unique targets for TF that have at least 5 target genes number of unique targets in db
  numuniqDBtarget <- length(unique(great5tg$Target))
  # num of TFs in list
  uniqTF <- unique(great1tg_noct$TF)
  
  # subset based on TF that are in the final list
  great5tg_all_target_query <- subset(great5tg, great5tg$TF %in% uniqTF)
  # List of all target genes by TF that are based on the unique TFs
  all_target_query <- lapply(unique(great5tg_all_target_query$TF),
                             function(x) subset(great5tg_all_target_query,
                                                great5tg_all_target_query$TF == x))
  cat("Calculating p-value between TF and list of DE genes...\n")
  queryoverlap <- lapply(all_target_query, function(x)
    testGeneOverlap(
      newGeneOverlap(
        x$Target,check_in_db$gene_symbol,genome.size = numuniqDBtarget)))
  TFpval <- sapply(queryoverlap, function(x) x@pval)
  cat("Collapsing target genes and PMIDs to matching TFs...\n")
  symmvec <- lapply(final_list, function(x) paste(x[["Target"]], x[["Mode of Regulation"]], sep = "_"))
  symmodevec <- sapply(symmvec, function(x) paste(unlist(x), collapse  = "|"))
  PMIDvec<- sapply(final_list, function(x) paste(x[["References (PMID)"]], collapse = "|"))
  TFpadj <- p.adjust(TFpval, method = "BH")
  TFoverlap <- subset(great1tg, !duplicated(great1tg$TF))$overlap_ct
  cat("Creating TF matrix...\n")
  
  TFsdbnum <- great5tg[great5tg$TF %in% uniqTF, c("TF", "uniq_target_ct") ]
  TFsdbnumuniq <- TFsdbnum[!duplicated(TFsdbnum$TF), c("uniq_target_ct")]
  TFdf <- data.frame(uniqTF, TFpval, TFpadj, TFsdbnumuniq, TFoverlap, symmodevec, PMIDvec, stringsAsFactors = F)
  ######### Aggregate TFs by mode of regulation and PMID ######
  
  names(great1tg_noct)[3:4] <- c("mode","ref")
  combTFtar <- paste(great1tg_noct$TF, great1tg_noct$Target, sep = "_")
  newfr <- data.frame(combTFtar, great1tg_noct$mode, great1tg_noct$ref, stringsAsFactors = F)
  colnames(newfr)[2:3] <- c("mode","ref")
  result <- aggregate(.~combTFtar, data = newfr, paste, collapse = "|")
  
  sep <- strsplit(result$combTFtar, split = "_")
  tfvec <- sapply(sep, "[", 1)
  targvec <- sapply(sep, "[", 2)
  #paste into new frame that is aggregated
  finalframe <- data.frame(tfvec, targvec, result$mode, result$ref)
  cat("Subsetting TF matrix based on user's fdroutput cutoff...\n")
  sig <- subset(TFdf, TFdf$TFpadj <= fdroutput)
  #sig <- subset(TFdf, TFdf$TFpadj <= 0.4) #troubleshooting
  #sig <- subset(TFdf, TFdf$TFpadj <= 0.6) #mouse
  
  ##### If there are no significant genes based on the user's input return #####
  if (nrow(sig) == 0) {
    print("No DE genes based on parameters for fdrout try changing parameters")
    cat("\n")
    return()
  }
  
  finsig <- subset(finalframe, finalframe$tfvec %in% sig$uniqTF )
  
  # finsig$pval <- sapply(finsig$tfvec, 
  #              function(x) sig$TFpval[which(as.character(x) == sig$uniqTF)])
  
  
  finsig$pvaladj <- sapply(finsig$tfvec, 
                           function(x) sig$TFpadj[which(as.character(x) == sig$uniqTF)])
  cat("\n")
  cat("Creating network graph...\n")
  #creates edges
  sigrelation2 <- data.frame(from = finsig$tfvec, to= finsig$targvec, arrows = c("to"),
                             title =paste(finsig$result.mode, finsig$result.ref, sep = "<br>"))
  
  cat("Number of TFs:", length(unique(sigrelation2$from)), "\n")
  cat("Number of unique target genes:", length(unique(sigrelation2$to)), "\n")
  #creates vertices
  overlapinfo <- paste("# of overlapped genes: ", sig$TFoverlap)
  #pvalinfo <- paste("P-Value: ", sig$TFpval) #Right now pvalue is the TFdf pvalue for gene overlap
  
  pvalinfo <- paste("P-Value: ", sig$TFpadj) #Right now pvalue is the TFdf pvalue for gene overlap
  sigcombmoderef <- paste(sig$symmodevec, sig$PMIDvec, overlapinfo, pvalinfo, sep = "<br>")
  sigallvert <- union(finsig$tfvec, finsig$targvec)
  sigallvertdata <-data.frame(sigallvert)
  cat("Number of total nodes:" , length(sigallvert), "\n")
  cat("\n")
  sigallvertdata$color <- ifelse(sigallvertdata$sigallvert %in% sig$uniqTF,green, purple)
  #adds annotations of mode of regulation, PMID, num overlapped genes, and gene overlap pval
  sigallvertdata$title <- ifelse(sigallvertdata$sigallvert %in% sig$uniqTF, sigcombmoderef, "")
  
  #Scales size of nodes based on how many overlapping genes that TF has
  scale1 <- log(sig$TFoverlap) * 25
  sigallvertdata$value<- ifelse(sigallvertdata$sigallvert %in% sig$uniqTF, scale1, 25)
  
  
  #Count the number of arrows pointing to target node
  lengthuniraw <- c()
  for (i in 1:nrow(sigrelation2)){
    tfsave <- as.character(sigrelation2$to[i])
    #numTraw <- length(unique(relation2[relation2$to == tfsave ,]$to))
    pl <-  nrow(sigrelation2[sigrelation2$to == tfsave ,])
    i = i + 1
    lengthuniraw <- append(lengthuniraw, pl)
  }
  sigrelation2$num <- lengthuniraw
  
  
  dfnodenum <- data.frame(sigrelation2$to, sigrelation2$num)
  
  #Scales the size of target genes based on how many TFs regulated them
  for (i in 1:nrow(sigallvertdata)){
    if ( sigallvertdata$value[i] == 25) {
      tofind <- as.character(sigallvertdata$sigallvert[i])
      val <- unique(dfnodenum[tofind == dfnodenum$sigrelation2.to, 2])
      scalenode <- (log(val)  + 1) * 25
      sigallvertdata$value[i] = scalenode                  
      
    }
    i = i+ 1
    
  }
  
  
  
  
  #Colors the edges of the arrows based on the mode of regulation
  vec <- c()
  for (ele in finsig$result.mode) {
    lecol <- if(ele == "Activation" ){
      vec <- append(vec, red)
      
    }else if (ele == "Repression") {
      
      vec <- append(vec, blue) 
    }else if (ele == "Unknown") {
      vec <- append(vec, orange)
    } else {
      vec <- append(vec, black)
    }
    
  }
  sigrelation2$color <- vec
  
  
  
  
  h <- graph.data.frame(sigrelation2, directed = TRUE, vertices = sigallvertdata)
  #Activates
  #Represses
  #Unknown
  #Multiple
  legendnodes <- data.frame(label = c("TF", "Gene"),
                       shape = c("ellipse"),
                       color = c(green, purple))
  legendedges <- data.frame(label = c("Activates", "Represses", "Unknown", "Multiple"),
                            shape = c("ellipse"),
                            color = c(red, blue, orange, black))
  datavo <- toVisNetworkData(h)
  ### Creates network graph ####
  graph <- visNetwork(nodes = datavo$nodes, edges = datavo$edges, height = "700px", width = "100%") %>%
    visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T)) %>%
    visInteraction(keyboard = TRUE, tooltipDelay = 0, tooltipStay = 10000) %>%
    visLayout(randomSeed = 18) %>%
    visLegend(addNodes = legendnodes, addEdges = legendedges, useGroups = F)
  
  
  graph
  cat("Saving network graph...\n")
  #netname <- paste(fname,  infile_name,fdroutput, "net.html", sep = "-")
  #cat("Network graph written out to", netname, "\n")
  cat("Network graph written out","\n")
  
  #graph %>% visSave(file = netname)
  graph %>% visSave(file = "net.html")
  
  #Creates circle graph
  circgraph <- visNetwork(nodes = datavo$nodes, edges = datavo$edges, main = "A really simple example", height = "700px", width = "100%") %>%
    visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T)) %>%
    visInteraction(keyboard = TRUE, tooltipDelay = 0, tooltipStay = 5000) %>%
    visIgraphLayout(layout = "layout_in_circle") %>%
    visLayout(randomSeed = 321) %>%
    visLegend(addNodes = legendnodes, addEdges = legendedges, useGroups = F)
  circgraph
  cat("Saving circle graph...\n")
  
  #circname <- paste(fname,  infile_name,fdroutput, "circ.html", sep = "-")
  
  #cat("Circle graph written out to", circname, "\n")
  cat("Circle graph written out to", "\n")
  
  #graph %>% visSave(file = circname)
  circgraph %>% visSave(file = "circ.html")
  
  cat("Saving TF matrix...\n")
  names(sig) <- c("TF", "Gene Overlap p-value", "Gene Overlap adjusted p-value", "# of target genes in DB",
                  "# of overlapped genes",
                  "Target Gene and Mode of Regulation", "References (PMID)")
  #txtname <- paste(fname,  infile_name,fdroutput, "sigTFs.txt", sep = "-")
  #### Writes out datafram of sig TFs ####
  #cat("TF matrix written out to", txtname, "\n")
  cat("TF matrix written out", "\n")
  
  #sigTable <- write.table(sig,"sig_TFs_network.txt",sep="\t",row.names=FALSE)
  cat("Finished!\n")
  #return(sig)
  ## fdfdfdsfsf
  #fdsfsf
  graphs <- list(graph, circgraph, sig)
  return(graphs)
}
#Functions to make graphs
#dflist <- makenetgraph(file, speciesdb, "GL", infile_type) 

#fdroutput <- 0.05
#infile_type <- "list"
#infile_name <- "motif-ui/outputgenes.txt"

#fname <- "GL"

#ex_tfgenes <- read.table("~/Documents/Salk/TFMatrix/ex_tfgenes.txt", 
#                        quote="\"", comment.char="", stringsAsFactors=FALSE)
#dflist <- makenetgraph(ex_tfgenes, "hg19", 0.05) 
#makenetgraph<-function(file, db, fname, ftype) {
  