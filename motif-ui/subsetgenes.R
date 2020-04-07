#!/usr/bin/env Rscript

#### Takes a list of input genes and subsets allTFcounts.txt by genes of interest
suppressMessages(library(gplots))

args<-commandArgs(trailingOnly=T)

generatePCA <- function(infile_name,genome) {

targetmatnum <- infile_name
transformed <-targetmatnum #motifs by genes 
print(head(targetmatnum))
print("transformed dim\n")
print(dim(transformed))
no_zeroes <- transformed[ , apply(transformed, 2, var, na.rm = TRUE) != 0]

pcagenes <- prcomp(no_zeroes, center = T , scale. = T)


#compute standard deviation of each principal component
std_dev <- pcagenes$sdev

#compute variance
pr_var <- std_dev^2

#proportion of variance explained
prop_varex <- pr_var/sum(pr_var)

pc1var <- prop_varex[1]
pc2var <- prop_varex[2]
library(ggplot2)

pcagenesDF <- data.frame(pcagenes[["x"]])
pcagenesDF <- pcagenesDF[,c(1,2,3)]
cat("PCA output \n")
cat(head(pcagenesDF[1:10,1]),"\n")
cat(dim(pcagenesDF))

#View(rownames(pcagenesDF))

groupnum <- rownames(pcagenesDF)
## new plot with 150
xmax <- round(max(pcagenesDF$PC1))+1
xmin <- round(min(pcagenesDF$PC1))-1

ymax <- round(max(pcagenesDF$PC2))+1
ymin <- round(min(pcagenesDF$PC2))-1
# ggplot(pcagenesDF,aes(PC1 ,PC2)) +
#   ggtitle("Genes Grouped by Common Motifs") +
#   xlab(paste("PC1",round(pc1var * 100, 0), "%"))+
#   ylab(paste("PC2",round(pc2var * 100, 0), "%")) +
#   xlim(xmin, xmax) + 
#   ylim(ymin,ymax) +
#   #geom_text(label= groupnum, size = 3, nudge_x = 1 ,nudge_y = 2) +
#   #scale_fill_discrete(name = "Grouping") + 
#   geom_point(size = 3) +
#   #geom_step() +
#   #scale_color_manual(values = c("#7FC97F","#F0027F","#386CB0"))
#   theme_classic() 
# dev.off()
listobj <- list(pcagenesDF, pc1var, pc2var)
return(listobj) ##Return PCA plot to plot in shiny app

#return(pcagenesDF) ##Return PCA plot to plot in shiny app

}


