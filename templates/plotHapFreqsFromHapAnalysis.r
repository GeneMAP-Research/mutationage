#!/usr/bin/env Rscript

getPackages <- function() 
{
    require(plyr)
    require(tools)
    #require(gsubfn)
}

getInputFiles <- function() 
{
    files <- Sys.glob("*.freq")
    hap <- read.csv2("haps.ids.txt", h=T, quote="")

    getInputFiles_output <- list()
    getInputFiles_output$files <- files
    getInputFiles_output$hap <- hap

    return(getInputFiles_output)
}

combineHapFiles <- function(files, hap)
{
    files <- getInputFiles_output$files
    hap <- getInputFiles_output$hap

    for (file in files) 
    {
        file_name <- file_path_sans_ext(file)
        file_name <- file_path_sans_ext(file_name)
        column_name <- c("Hap",file_name)
        hap_file <- read.csv2(file, h=T, sep="\t", quote="")
        hap_file <- data.frame(hap_file$Hap, hap_file$Freq)
        colnames(hap_file) <- column_name
        hap <- merge(hap, hap_file, by="Hap", all=T)
    }

    return(hap)
}

saveMergedHapFile <- function(hap)
{
    write.csv(hap, "test.haps.matrix.csv",  na="0", row.names=F, quote=F)
}

transposeMergedHapFile <- function() 
{
    mergedHap <- read.csv2("test.haps.matrix.csv", h=T, sep=",", quote="")
    mergedHap_rownames <- c()
    for (row in 1:nrow(mergedHap)) { mergedHap_rownames[row] <- paste0("Hap",row) }
    rownames(mergedHap) <- mergedHap_rownames
    hap_ids <- mergedHap$Hap
    mergedHap$Hap <- NULL
    transposedHap <- t(mergedHap)
    write.csv(transposedHap, "test.haps.transposed.matrix.csv",  na="0", row.names=T, quote=F)
}

getPlotInput <- function()
{
    transposedHap <- read.csv2("test.haps.transposed.matrix.csv", h=T, quote="", sep=",")
    rownames(transposedHap) <- transposedHap$X
    transposedHap$X <- NULL
#    orderedHap <- transposedHap[
#				order(transposedHap$Hap1, 
#				      transposedHap$Hap2, 
#				      transposedHap$Hap3, 
#				      transposedHap$Hap4, 
#				      transposedHap$Hap5, 
#				      transposedHap$Hap6, 
#				      transposedHap$Hap7, 
#				      transposedHap$Hap8, 
#				      transposedHap$Hap9, 
#				      transposedHap$Hap10)
#				,]

    orderedHap <- transposedHap[order(transposedHap$Hap2, transposedHap$Hap5, transposedHap$Hap1),]
    
    #colnames(orderedHap) <- hap_ids
    hapMatrix <- as.matrix(orderedHap)
    
    return(hapMatrix)
}

plotHapFreqs <- function(hapMatrix)
{
    pdf("test.pdf", colormodel="cmyk")
    par(fig=c(0,0.95,0.25,0.75), mar=c(4,3,4,2), bty="n", cex=0.7)
    barplot(t(hapMatrix), col=rainbow(ncol(hapMatrix)), border=0, space=0.0, las=3, cex.names=0.9, cex.axis=0.9)
    par(fig=c(0.87,1,0,1), mar=c(2,0,2,0), new=F, bty="n", cex=0.8)
    legend("center", legend=c(colnames(hapMatrix)), col=rainbow(ncol(hapMatrix)), horiz=F, pch=15, bty="n")
    dev.off()
}

getPackages()
getInputFiles_output <- getInputFiles()
hap <- combineHapFiles(getInputFiles_output)
saveMergedHapFile(hap)
transposeMergedHapFile()
hapMatrix <- getPlotInput()
plotHapFreqs(hapMatrix)

