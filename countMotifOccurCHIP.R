#!/usr/bin/env Rscript

args<-commandArgs(trailingOnly=T)
infile_name=args[1]
cat("Input file name: ",infile_name, "\n")
filetocount <-read.delim(infile_name, stringsAsFactors = F, header = T)
#filetocount <- read.delim("~/annoOct2(POU,Homeobox).txt", stringsAsFactors=FALSE)
#filetocount <- read.delim("~/p53(p53).m",  header = F, stringsAsFactors=FALSE)
#filetocount <- read.delim("~/annoAP-1(bZIP).txt",  header = T, stringsAsFactors=FALSE)
#filetocount <- read.delim("~/annoc-Myc(bHLH)|LNCAP-cMyc-ChIP-Seq(Unpublished).txt", stringsAsFactors=FALSE)

#filetocount <- read.delim("~/annoRAR:RXR(NR),DR5.txt", stringsAsFactors=FALSE)

#filetcount <- `annoOct2(POU,Homeobox)`
#print(nrow(filetocount))
cat("Rows in", infile_name, ":", nrow(filetocount), "\n")

#filetocount <- read.delim("~/Desktop/annoZSCAN22(Zf).txt",
#                             stringsAsFactors=FALSE)
#filetocount <- read.delim("~/Desktop/annop63(p53).txt",
#                        stringsAsFactors=FALSE)


##### checking mouse files 

#filetocount <- read.delim("~/Desktop/annoAtf2(bZIP).txt",
#                                                 stringsAsFactors=FALSE)
#filetocount <- read.delim("~/Desktop/annoARE(NR).txt",
#                          stringsAsFactors=FALSE)
sub <- subset(filetocount, filetocount$Distance.to.TSS > -1000  
              & filetocount$Distance.to.TSS < 100)


##############################################################
##### subset of filetocount and remove duplicate motif occurrences from different chip-seq exp
#sub$comb <- paste(sub$Chr, sub$Start, sub$End, sub$Focus.Ratio.Region.Size)

#sub <- sub[!duplicated(sub$comb),]

## Extracts name of TF and CHIP-seq to write out to file
fullnm <- as.character(infile_name)
#Remove anno and .txt
#fullnm <- "annoRAR:RXR(NR),DR5|ES-RAR-ChIP-Seq(GSE56893)-12.txt"
fullnm <- gsub(".txt", "",fullnm)
fullnm <- gsub("anno", "",fullnm)

#RAR:RXR(NR),DR5|ES-RAR-ChIP-Seq(GSE56893)-12
#cat("nm to write", fullnm)
#names(sub)[1] <- "Motif"
#namemotif <- sub$Motif
#tfmotif <- sapply(strsplit(namemotif,"[/]"), `[`, 1)
#chipmotif <- sapply(strsplit(namemotif,"[/]"), `[`, 2)
#fullnm <- paste(tfmotif, chipmotif, sep = "|")
#fullnm <- fullnm[1]

#Sort by gene name
sub <- sub[order(sub$Gene.Name, decreasing = F),]
#print(nrow(sub))
cat("Rows in promoter region:", nrow(sub), "\n")

sub$numTF <- sapply(sub$Gene.Name, 
               function(x) nrow(sub[sub$Gene.Name == x,]))


subgenes <- sub$Gene.Name
TFcountpergene <- sub$numTF
TFcountdf <- data.frame(subgenes, TFcountpergene)
TFcountdf <- TFcountdf[!duplicated(subgenes),]

firstgene <- subgenes[1]
if(nchar(firstgene) == 0) {
  TFcountdf <- TFcountdf[-c(1),]
} 
#print(nrow(TFcountdf))
cat("Rows in TFcountdf:", nrow(TFcountdf), "\n")

write.table(TFcountdf, paste(fullnm,"TFcountdf.txt", sep = "-"), row.names = F , col.names = F)
print("Written to file")

