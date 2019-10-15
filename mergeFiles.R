#!/usr/bin/env Rscript

#Merge input files

# ZNF692 <- read.table("~/Desktop/ZNF692(Zf)-TFcountdf.txt",
#                                      quote="\"", comment.char="",
#                                      stringsAsFactors=FALSE)
# 
# ZSCAN22 <- read.table("~/Desktop/ZSCAN22(Zf)-TFcountdf.txt",
#                                       quote="\"", comment.char="", 
#                                       stringsAsFactors=FALSE)
# 
# ZNF675 <- read.table("~/Desktop/ZNF675(Zf)-TFcountdf.txt",
#                      quote="\"", comment.char="",
#                      stringsAsFactors=FALSE)
# ZNF669 <- read.table("~/Desktop/ZNF669(Zf)-TFcountdf.txt",
#                     quote="\"", comment.char="",
#                     stringsAsFactors=FALSE)
# p63 <- read.table("~/Desktop/p63(p53)-TFcountdf.txt",
#                      quote="\"", comment.char="",
#                      stringsAsFactors=FALSE)


# names <- list.files()
# filenames <- grep("TFcountdf", names, value = T)
# splitname <- gsub("-TFcountdf.txt", "", filenames)

#Gets only the names of the TFdf files
#splitname <- strsplit(filenames, "-")
# onlyTF <- sapply(splitname, function(x) x[1])


# filestomerge <- filenames
# readfirst = FALSE
#while there are files to merge, keep merging
# while (length(filestomerge) != 0) {
#   #print( filestomerge)
#   #print("begin loop")
#   if (readfirst == FALSE) {
#     cat("Merging first two files", filestomerge[1], filestomerge[2], "\n")
#     df1 <- read.table(filestomerge[1],quote="\"", comment.char="", 
#                     stringsAsFactors=FALSE)
#     names(df1)[1] <- "ID"
#     names(df1)[length(names(df1))] <- onlyTF[1]
#     #print(names(df1)[length(names(df1))])
#     #cat("df1 names", names(df1))
#     #print(nrow(df1))
#   
#     df2 <- read.table(filestomerge[2],quote="\"", comment.char="", 
#                     stringsAsFactors=FALSE)
#     names(df2)[1] <- "ID"
#     names(df2)[length(names(df2))] <- onlyTF[2]
#     #cat("df2 names", names(df2))
#     #print(nrow(df2))
#   
#     comb <- merge(df1, df2, all = T)
#     #cat("comb init nrow", nrow(comb), "\n")
#     readfirst = TRUE
#     filestomerge <- filestomerge[c(-1,-2)]
#     onlyTF <- onlyTF[c(-1,-2)]
#     #cat("files after init", filestomerge, "\n")
#     
#   } else {
#     
#     #cat("Files after reading init:", filestomerge, "\n")
#     addedfr <-read.table(filestomerge[1],quote="\"", comment.char="", 
#                          stringsAsFactors=FALSE)
#     cat("Merging additional TF", filestomerge[1], "\n")
#     names(addedfr)[1] <- "ID"
#     names(addedfr)[length(names(df2))] <- onlyTF[1]
#     #cat("toadd names", names(addedfr), "\n")
#     comb <- merge(comb, addedfr, all =T)
#     #print(head(comb,5))
#     filestomerge <- filestomerge[-1]
#     onlyTF <- onlyTF[-1]
#     
#   }
#   
# }
# comb[is.na(comb)] = 0
# write.table(comb, "allTFcounts.txt", sep = "\t", row.names = F, col.names = T)
# cat("Written to file!")
#le <- merge(`ZNF692(Zf).TFcountdf`, `ZSCAN22(Zf).TFcountdf`, all = T)

#for file in file name 

# ZSCAN22names <- ZSCAN22$ID
# ZNF692names <- ZNF692$ID
# 
# #remember union doesn't have them sorted
# both <- union(ZSCAN22names, ZNF692names) #should be 3660 names altog
# 
# bothsort <- sort(both)
# 
# leID <- le$ID
# inbothnotLE <- setdiff(bothsort, leID)
# nodup <- !duplicated(bothsort)



#####
##### sand-box shorter mergeFiles.R
####

## Take the TFcountdffiles and store the names to an object list

# files <- list.files(pattern = "^comb")
files <- list.files(pattern = "TFcountdf.txt")
#files <- files[1:5]
#files <- files[1:2]
#f1nm <- files[1]
#f1 <- read.delim(f1nm, stringsAsFactors = F, sep = " ")
#f1 <- read.delim(f1nm, stringsAsFactors = F, sep = " ")

#objlist <- lapply(files, read.delim)

#foxa1 <- read.delim(colnames(mergedframe)[9])

#objlist <- lapply(files, function(x) read.delim(x, stringsAsFactors = F, sep = " "))
#objlist <- lapply(files, function(x) read.delim(x, quote="\"", comment.char=""))



## Add the TF name as column header for each TFcountdffile
addcolnames <- function(TFcountdffile) {
  filefrm <- read.delim(TFcountdffile, stringsAsFactors = F, sep = " ", header = F)
  filenm <- TFcountdffile
  filenm <- gsub("-TFcountdf.txt", "", filenm)
  colnames(filefrm) <- c("genes", filenm)
  
  return(filefrm)
  
}

#newfile <- addcolnames("ZNF675(Zf)-TFcountdf.txt")

objlist <- lapply(files, function(x) addcolnames(x))

#Combine all the TFcountdffiles
mergedframe <- Reduce(function(x, y) merge(x, y, all=TRUE), objlist)
mergedframe[is.na(mergedframe)] = 0
write.table(mergedframe, "allTFcounts.txt", sep = "\t", row.names = F, col.names = T)

