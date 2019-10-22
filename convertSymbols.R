#!/usr/bin/env Rscript

#### Takes in a gene list as a text file, a genome (either mm10 or hg19),
## and ID type (geneID, Unigene, RefSeq, Ensembl, name) and converts
## them to gene symbol
#args<-commandArgs(trailingOnly=T)

#infile_name = args[1]
#cat("Reading in gene list:", infile_name, "\n")

convertGenes <- function(infile_name, genome_type, idtype) {
inputgenesdata <- read.delim(infile_name, header = F, stringsAsFactors = F)
#inputgenesdata <-read.delim("inputgenes.txt", stringsAsFactors = F, header = F)

## Example file 
#GeneID
#inputgenesdata <- read.delim("~/Documents/Salk/TFMatrix/31geneidmm10.txt",
#                             stringsAsFactors = F, header = F)

#Unigene
#inputgenesdata <- read.delim("~/Documents/Salk/TFMatrix/31geneunimm10.txt",
#                             stringsAsFactors = F, header = F)

#RefSeq

#inputgenesdata <- read.delim("~/Documents/Salk/TFMatrix/31generefmm10.txt",
#                            stringsAsFactors = F, header = F)
#Ensembl
#inputgenesdata <- read.delim("~/Documents/Salk/TFMatrix/31geneensmm10.txt",
#                            stringsAsFactors = F, header = F)


names(inputgenesdata) <- "gene_symbol"

if (genome_type == "mm10") {
  cat("Using mm10 genome\n")
  generef <- read.delim("mm10_GeneIDsymbol.tsv", 
                        header=TRUE, stringsAsFactors = F)
  
} else if(genome_type == "hg19") {
  cat("Using hg19 genome\n")
  generef <- read.delim("hg19_GeneIDsymbol.tsv", 
                        header=TRUE, stringsAsFactors = F)

} else {
  print("Not a valid genome")
}




if(idtype == "geneid") {
  cat("Converting input gene ID to gene symbol...\n")
  sub <- subset(generef,  generef$GeneID %in% inputgenesdata$gene_symbol)
  
}else if(idtype == "unigene") {
  cat("Converting input Unigene IDs to gene symbol...\n")
  sub <- subset(generef,  generef$Unigene %in% inputgenesdata$gene_symbol)
  
}else if(idtype == "refseq") {
  cat("Converting input RefSeq IDs to gene symbol...\n")
  sub <- subset(generef,  generef$RefSeq %in% inputgenesdata$gene_symbol)

  } else if(idtype == "ensembl") {
  cat("Converting input Ensembl IDs to gene symbol...\n")
  sub <- subset(generef,  generef$Ensembl %in% inputgenesdata$gene_symbol)

  } else if (idtype == "name") {
  cat("Taking in name IDs...\n")
  sub <- subset(generef,  generef$name %in% inputgenesdata$gene_symbol)
  }

#Converts all forms to ENS
outgenes <- as.character(sub$Ensembl)
outgenes <- data.frame(outgenes)
names(outgenes) <- "ensembl"
#write.table(inputgenesdata, "outputgenes.txt", row.names = F, col.names = F, quote = F)
#cat(nrows(inputgenesdata))
return(outgenes)

}

#outfile <- convertGenes("inputgenes.txt", "hg19", "name")

