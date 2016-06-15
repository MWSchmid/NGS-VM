#!/usr/bin/env Rscript

## arguments from commandline
myarg <- commandArgs()
argPos <- grep("--args", myarg, fixed = TRUE)
inputFile <- as.character(myarg[argPos+1])
outputFile <- as.character(myarg[argPos+2])
testSample <- as.character(myarg[argPos+3])
refSample <- as.character(myarg[argPos+4])

# a function to join overlapping fragments with the same sign for the average difference
f.join.overlapping.fragments.directional <- function(data, otherCols, summaryFunction, dirCol, gap = 0) {
  data <- data[with(data, order(chrom, start, end)),]
  data$regNum <- 1
  offset <- 1
  for (chr in unique(data$chrom)) {
    subData <- subset(data, chrom == chr)
    subData$regNum[1] <- offset
    if (nrow(subData) > 1) {
      temp <- (subData$start[2:nrow(subData)] - subData$end[1:(nrow(subData)-1)]) > gap
      dirCheck <- sign(subData[[dirCol]][2:nrow(subData)]) != sign(subData[[dirCol]][1:(nrow(subData)-1)])
      subData$regNum[2:nrow(subData)] <- offset + cumsum(temp | dirCheck)
    } 
    data[data$chrom == chr, ] <- subData
    offset <- max(data$regNum) + 1
  }
  out <- data.frame(
    chrom = aggregate(data$chrom, by=list(region=data$regNum), unique, simplify = TRUE)$x,
    start = aggregate(data$start, by=list(region=data$regNum), min, simplify = TRUE)$x,
    end = aggregate(data$end, by=list(region=data$regNum), max, simplify = TRUE)$x,
    count = aggregate(data$chrom, by=list(region=data$regNum), length, simplify = TRUE)$x,
    stringsAsFactors = FALSE
  )
  for (cn in otherCols) {
    out[[cn]] = aggregate(data[[cn]], by=list(region=data$regNum), summaryFunction, simplify = TRUE)$x
  }
  return(out)
}

data <- read.table(inputFile, sep = '\t', quote = '', col.names = c("chrom", "start", "end", "pVal", "aveDiff"), stringsAsFactors = FALSE)
cat("### Rscript: processChXPrunsTest.R:", "Loaded", nrow(data), "regions.\n")
data <- f.join.overlapping.fragments.directional(data, c("pVal", "aveDiff"), mean, "aveDiff")
cat("### Rscript: processChXPrunsTest.R:", nrow(data), "regions remain after merging overlapping regions.\n")
data$size <- data$end - data$start
cat("### Rscript: processChXPrunsTest.R: Distribution of the average difference:\n")
print(quantile(abs(data$aveDiff), seq(0,1,.1)))
cat("### Rscript: processChXPrunsTest.R: Distribution of the region size:\n")
print(quantile(data$size, seq(0,1,.1)))
cutoff <- quantile(abs(data$aveDiff), .9)
numOK <- with(data, sum((size > 5000) | (abs(aveDiff) > cutoff)))
cat("### Rscript: processChXPrunsTest.R:", numOK, "regions are > 5 kb or have an absolute average difference above 1.5.\n")
data$score <- data$aveDiff*100
data$name <- "toReplace"
data$name[data$aveDiff > 0] <- paste0(testSample, 1:sum(data$aveDiff > 0))
data$name[data$aveDiff < 0] <- paste0(refSample, 1:sum(data$aveDiff < 0))
write.table(data[,c("chrom", "start", "end", "name", "score")], outputFile, sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
cat("### Rscript: processChXPrunsTest.R: Done.\n")