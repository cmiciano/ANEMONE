#!/usr/bin/env Rscript
### Script to separate RAR:RXR motif file into separate files
files <- list.files(pattern = "RAR:RXR")
filesm <- grep(".m", files, value = T)
filesm
origmotif <- read.delim(filesm, header = F, stringsAsFactors = F)
colnames(origmotif)[7]<- "Focus.Ratio.Region.Size"
RARsplit <- split(origmotif, nchar(origmotif$Focus.Ratio.Region.Size))


RAR12 <- RARsplit[["12"]]
write.table(RAR12, 
            paste("RAR:RXR(NR),DR5|ES-RAR-ChIP-Seq(GSE56893)","12.m", sep = "-"),
            row.names = F , col.names = F, , quote = F, sep = "\t")

RAR18 <- RARsplit[["18"]]
write.table(RAR18, 
            paste("RAR:RXR(NR),DR5|ES-RAR-ChIP-Seq(GSE56893)","18.m", sep = "-"),
            row.names = F , col.names = F, , quote = F, sep = "\t")
