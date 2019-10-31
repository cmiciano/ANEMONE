genHeatmap <- function(infile_name,genome) {
  print("in genheat")
  #infile_name=args[1]
  cat("Input file name: ",nrow(infile_name), "\n")
  #inputgenes <-read.delim(infile_name, stringsAsFactors = F, header = F)
  
  #inputgenes <-read.delim("inputgenes.txt", stringsAsFactors = F, header = F)
  #inputvec <- readLines("../inputgenes.txt")
  #inputgenes <- read.delim("~/Documents/Salk/TFMatrix/33geneuni.txt", header = F)
  #inputgenes <- read.delim("~/Documents/Salk/TFMatrix/ex_tfgenes.txt", header = F)
  
  #inputgenes <-read.delim("../inputgenes.txt", stringsAsFactors = F, header = F) #works
  inputgenes <- infile_name
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
    
  } else if (genome == "mm10") {
    allTFcounts <- read.delim("data/allTFcountsmm10.txt", stringsAsFactors=FALSE)
    
  }
  
  
  ### Subset expression matrix by the representative motifs
  repmotifs <- read.delim("data/163repmotifssummary.txt", stringsAsFactors=FALSE)
  
  
  
  #View(colnames(allTFcounts))
  colallTF <- colnames(allTFcounts)
  #View(repmotifs$short_rep_motif)
  shortTF <- repmotifs$short_rep_motif
  
  allTFsub <- subset(allTFcounts, colallTF %in% shortTF)
  allTFsub <- allTFcounts[,colallTF %in% shortTF]
  ID <- allTFcounts[,1]
  allTFsub <- cbind(ID, allTFsub) #allTFsub represents the 20k ish genes x 163 rep motifs
  
  
  ### Subset allTFcounts.txt into the representative TFs
  
  #genelist <- c("A1CF", "AACS", "AADACL4", "AADAT")
  
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
  
  targetmatnmum <- targetgenes[,2:ncol(targetgenes)]
  cat(head(targetmatnmum[1:10,1]))
  #targetmatnmum <- expint[,2:ncol(expint)]
  
  rownames(targetmatnmum) <- genenm
  targetmatnmum <- log(targetmatnmum+5) #genes by motif
  
  cat(head(targetmatnmum[1:10,1]))
  
  # ## Generate a heatmap of the counts of each motif occurrence in 
  # ## the promoter region of each gene
  # cat("Generating heatmap...\n")
  # #pdf("genemotifctint.pdf", 8, 15)
  # pdf("genemotifctint.pdf", 4, 4)
  # #par(mar=c(10,2,2,3), cex=1.0)
  # par(mar=c(1,1,1,1), cex=1.0)
  # 
  # heatmap.2(as.matrix(targetmatnmum), Rowv = T, Colv = T, col = heat.colors, 
  #           trace = "none", labRow = rownames(targetmatnmum),
  #           #lhei = c(0.5,5),     
  #           lhei = c(0.5,1),     
  #           #lwid = c(0.5,0.5),
  #           cexRow=0.15,
  #           cexCol=0.15)
  #           #hclustfun=function(x) hclust(x, method="ward.D"))
  # invisible(dev.off())
  # cat("Heatmap finished!\n")
  
  
  print("ret gen")
  #heatobj <- heatmaply(targetmatnmum)
  
  return(targetmatnmum) ##Return dataframe
  #return(heatobj) ##Return heatmap object
  
  
}

