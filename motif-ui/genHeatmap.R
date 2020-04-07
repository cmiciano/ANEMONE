genHeatmap <- function(convgenes,genome) {
  print("in genheat")
  #infile_name=args[1]
  cat("Input file name: ",nrow(convgenes), "\n")
  
  #inputgenes <-read.delim("inputgenes.txt", stringsAsFactors = F, header = F) #works
  #inputgenes <-read.delim("~/Downloads/33geneuni.txt", stringsAsFactors = F, header = F) #works
  # inputgenes <- read.delim("/Users/charlenemiciano/motif-genes-app/dr_fc15.txt")
  # inputgenes <- read.delim("/Users/charlenemiciano/motif-genes-app/ur_fc75.txt")
  
  inputgenes <- convgenes #DONT FORGET TO UNCOMMENT
  cat("inputgenes\n")
  print(class(inputgenes))
  print(dim(inputgenes))
  
  print(head(inputgenes[1:10,1]))
  #inputvec <- readLines(infile_name)
  
  #inputvec <- as.vector(inputgenes)
  colnames(inputgenes) <- "name"
  
  #inputvec <- inputvec[1:10]
  if (genome == "hg19") {
    allTFcounts <- read.delim("data/allTFcounts.txt", stringsAsFactors=FALSE)
    #allTFcounts <- read.delim("motif-ui/data/allTFcounts.txt", stringsAsFactors=FALSE)

  } else if (genome == "mm10") {
    allTFcounts <- read.delim("data/allTFcountsmm10.txt", stringsAsFactors=FALSE)
    
  }
  
  
  ### Subset expression matrix by the representative motifs
  repmotifs <- read.delim("data/163repmotifssummary.txt", stringsAsFactors=FALSE)
  #repmotifs <- read.delim("motif-ui/data/163repmotifssummary.txt", stringsAsFactors=FALSE)
  
  
  
  #View(colnames(allTFcounts))
  colallTF <- colnames(allTFcounts)
  #View(repmotifs$short_rep_motif)
  shortTF <- repmotifs$short_rep_motif
  
  allTFsub <- subset(allTFcounts, colallTF %in% shortTF)
  allTFsub <- allTFcounts[,colallTF %in% shortTF]
  ID <- allTFcounts[,1]
  allTFsub <- cbind(ID, allTFsub) #allTFsub represents the 20k ish genes x 163 rep motifs
  
  
  ### Normalize based on all motif occurrences in genome
  allTFnum <- allTFsub[,2:ncol(allTFsub)]  
  allMotifCt <- colSums(allTFnum)
  allCtNorm <- sweep(allTFnum, 2, allMotifCt, `/`)
  is.nan.data.frame <- function(x) {do.call(cbind, lapply(x, is.nan))}
  
  allCtNorm[is.nan(allCtNorm)] <- 0.00000000
  
  allTFsub[,2:ncol(allTFsub)] <- allCtNorm
  

  #Sample
  set.seed(1)
  #inputgenes <- sample(allTFcounts$ID, 1000, replace = F)
  #write.table(inputgenes, "inputgenes.txt", quote = F, row.names = F, col.names = F)
  cat("Subsetting table based on genes of interest...\n")
  targetgenes <- subset(allTFsub, allTFsub$ID %in% inputgenes$name)
  cat("target genes\n")
  cat(nrow(targetgenes))
  #targetgenes <- subset(allTFcounts, allTFcounts$ID %in% inputvec)
  
  motifnm <- colnames(targetgenes)
  genenm <- targetgenes$ID
  
  targetmatnum <- targetgenes[,2:ncol(targetgenes)] #absolute counts of motif occurrences
  
  ## divided by number of motif occurrences
  # targetsums <- colSums(targetmatnum)
  # targetdiv <- sweep(targetmatnum, 2, targetsums, `/`)
  # is.nan.data.frame <- function(x) {do.call(cbind, lapply(x, is.nan))}
  # 
  # targetdiv[is.nan(targetdiv)] <- 0.00000000
  
  #write.table(targetdiv, "ur_fc75_normcounts.txt", row.names = FALSE, col.names = T, quote = F, sep="\t")
  

  
  cat(head(targetmatnum[1:10,1]))
  rownames(targetmatnum) <- genenm
  #targetmatnum <- targetmatnum+5 #genes by motif #originally natural log what if log10?
  targetmatnum <- log(targetmatnum+5) #genes by motif #originally natural log what if log10?
  #targetmatnum <- targetmatnum[7:10,1:10] #genes by motif
 #temporary for figure use
  
  cat(head(targetmatnum[1:10,1]))
  print("ret gen")
  #heatobj <- heatmaply(targetmatnum)
  targetobj <- list(targetmatnum, targetgenes)
  return(targetobj) ##Return dataframe
  #return(heatobj) ##Return heatmap object
   #fsfsf
  options(digits = 15)
  
  #targetsci <- lapply(targetmatnum, formatC(targetmatnum, format = "e", digits = 6))
  
  # write.table(targetmatnum, "ur_fc75_normcounts_genome.txt", row.names = FALSE, col.names = T, quote = F, sep="\t")
  # 
  # 
  # 
  # statobj <- heatmap.2(as.matrix(targetmatnum), Rowv = T, Colv = T,
  #                      col = viridis(n = 256, alpha = 1, begin = 0, end = 1, option = "viridis"),
  #                      trace = "none", 
  #                      labRow = rownames(targetmatnum),
  #                      cexRow=0.5, #1.5
  #                      cexCol=0.3, # 1.5
  #                      #sepwidth = c(0.01, 0.01),
  #                      #lhei = c(1,5),
  #                      #lwid = c(1,6),
  #                      #margins = c(13,20),
  #                      hclustfun=function(x) hclust(x, method="ward.D")
  # )
  # statobj
  # statdend <- statobj[["colDendrogram"]]
  # values$dendstat <- statdend
  # 
  
  
  
  
  # heatmaply(targetmatnum, cexRow = 0.5, cexCol = 1.0,
  #           colors = viridis(n = 256, alpha = 1, begin = 0, end = 1, option = "viridis"), 
  #           #colors = cm.colors(256),
  #           #colors = hmcol,
  #           #colors=c("#009999", "#0000FF"), #blue and teal
  #           hclustfun = function(x) hclust(x, method="ward.D"),
  #           Colv = rev(statdend),
  #           seriate = "mean",
  #           #fontsize_row = 20,
  #           #fontsize_col = 20,
  #           #cellnote_size = 24,#seems to be size of plot
  #           row_dend_left = TRUE, 
  #           plot_method = "plotly")
  
}

