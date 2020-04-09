genHeatmap <- function(convgenes,genome) {
  cat("Input file name: ",nrow(convgenes), "\n")
  
  #inputgenes <-read.delim("inputgenes.txt", stringsAsFactors = F, header = F) #works
  #inputgenes <-read.delim("~/Downloads/33geneuni.txt", stringsAsFactors = F, header = F) #works
  #inputgenes <- read.delim("/Users/charlenemiciano/motif-genes-app/dr_fc15.txt")
  #inputgenes <- read.delim("/Users/charlenemiciano/motif-genes-app/ur_fc75.txt")
  
  inputgenes <- convgenes #DONT FORGET TO UNCOMMENT

  colnames(inputgenes) <- "name"
  
  if (genome == "hg19") {
    allTFcounts <- read.delim("data/allTFcounts.txt", stringsAsFactors=FALSE)

  } else if (genome == "mm10") {
    allTFcounts <- read.delim("data/allTFcountsmm10.txt", stringsAsFactors=FALSE)
    
  }
  
  
  ### Subset expression matrix by the representative motifs
  repmotifs <- read.delim("data/163repmotifssummary.txt", stringsAsFactors=FALSE)
  
  colallTF <- colnames(allTFcounts)
  shortTF <- repmotifs$short_rep_motif
  
  allTFsub <- subset(allTFcounts, colallTF %in% shortTF)
  allTFsub <- allTFcounts[,colallTF %in% shortTF]
  ID <- allTFcounts[,1]
  allTFsub <- cbind(ID, allTFsub) #allTFsub represents the 20k ish genes x 163 rep motifs
  
  
  #1) Starting with raw counts for each gene and motif (should be between 0 and like 100) 
  #2) Divide by total number of counts in all genes for each motif (normalizes for motif background frequency)
  #3) Multiply by 1000000 (make values larger than 1)
  #4) Add 1 and take the log (normalize for non-normal shape of values)
  
  ### Normalize based on all motif occurrences in genome
  allTFnum <- allTFsub[,2:ncol(allTFsub)]  
  allMotifCt <- colSums(allTFnum)
  allCtNorm <- sweep(allTFnum, 2, allMotifCt, `/`)
  is.nan.data.frame <- function(x) {do.call(cbind, lapply(x, is.nan))}
  allCtNorm[is.nan(allCtNorm)] <- 0.00000000
  
  # multiply by 1mil
  allCtMil <- allCtNorm * 1000000
  
  #add 1 and take log
  allCtLog <- log(allCtMil + 1)


  allTFsub[,2:ncol(allTFsub)] <- allCtLog
  
  cat("Subsetting table based on genes of interest...\n")
  targetgenes <- subset(allTFsub, allTFsub$ID %in% inputgenes$name)
  motifnm <- colnames(targetgenes)
  genenm <- targetgenes$ID
  
  targetmatnum <- targetgenes[,2:ncol(targetgenes)] #only counts of motif occurrences for input genes
  
  #targetmatnum <- log(targetmatnum + 5)

  rownames(targetmatnum) <- genenm
  options(digits = 15)
  print("ret gen")
  
  targetobj <- list(targetmatnum, targetgenes)
  return(targetobj) ##Return dataframe

  #write.table(targetgenes, "ur_fc75_normcounts_id_renorm.txt", row.names = FALSE, col.names = T, quote = F, sep="\t")
  set.seed(1)
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
  #  statobj
  #  statdend <- statobj[["colDendrogram"]]
  #  values$dendstat <- statdend

  

  
  
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

