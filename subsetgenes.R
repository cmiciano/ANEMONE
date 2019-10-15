#!/usr/bin/env Rscript

#### Takes a list of input genes and subsets allTFcounts.txt by genes of interest
suppressMessages(library(gplots))

args<-commandArgs(trailingOnly=T)
infile_name=args[1]
cat("Input file name: ",infile_name, "\n")
inputgenes <-read.delim(infile_name, stringsAsFactors = F, header = F)

#inputgenes <-read.delim("inputgenes.txt", stringsAsFactors = F, header = F)
inputvec <- readLines("inputgenes.txt")
#inputvec <- as.vector(inputgenes)
allTFcounts <- read.delim("~/allTFcounts.txt", stringsAsFactors=FALSE)


### Subset expression matrix for 
repmotifs <- read.delim("~/163repmotifssummary.txt", stringsAsFactors=FALSE)



View(colnames(allTFcounts))
colallTF <- colnames(allTFcounts)
View(repmotifs$short_rep_motif)
shortTF <- repmotifs$short_rep_motif

allTFsub <- subset(allTFcounts, colallTF %in% shortTF)
allTFsub <- allTFcounts[,colallTF %in% shortTF]
ID <- allTFcounts$ID
allTFsub <- cbind(ID, allTFsub)


### Subset allTFcounts.txt into the representative TFs

#genelist <- c("A1CF", "AACS", "AADACL4", "AADAT")

#Sample
set.seed(1)
#inputgenes <- sample(allTFcounts$ID, 1000, replace = F)
#write.table(inputgenes, "inputgenes.txt", quote = F, row.names = F, col.names = F)
cat("Subsetting table based on genes of interest...\n")
targetgenes <- subset(allTFcounts, allTFcounts$ID %in% inputvec)

motifnm <- colnames(targetgenes)
genenm <- targetgenes$ID

targetmatnmum <- targetgenes[,2:ncol(targetgenes)]
rownames(targetmatnmum) <- genenm
targetmatnmum <- log(targetmatnmum+5) #genes by sample

## Generate a heatmap of the counts of each motif occurrence in 
## the promoter region of each gene
cat("Generating heatmap...\n")
pdf("genemotifct.pdf", 8, 15)
par(mar=c(10,2,2,3), cex=1.0)
heatmap.2(as.matrix(targetmatnmum), Rowv = T, Colv = T, col = heat.colors, 
          trace = "none", labRow = rownames(targetmatnmum),
          lhei = c(0.5,5),          
          cexRow=0.15,
          cexCol=0.15,
          hclustfun=function(x) hclust(x, method="ward.D2"))
invisible(dev.off())
cat("Heatmap finished!\n")
