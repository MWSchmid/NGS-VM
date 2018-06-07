#!/usr/bin/env Rscript

## arguments from commandline
myarg <- commandArgs()
argPos <- grep("--args", myarg, fixed = TRUE)
inputFile <- as.character(myarg[argPos+1])
outputFile <- as.character(myarg[argPos+2])

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

data <- read.table(inputFile, sep = '\t', quote = '', col.names = c("chrom", "start", "end", "leftBorder", "rightBorder", "posFraction"), stringsAsFactors = FALSE)
cat("### Rscript: processChXPedgeTest.R:", "Loaded", nrow(data), "regions.\n")
data <- f.join.overlapping.fragments.directional(data, c("leftBorder", "rightBorder", "posFraction"), mean, "posFraction")
cat("### Rscript: processChXPedgeTest.R:", nrow(data), "regions remain after merging overlapping regions.\n")
data$size <- data$end - data$start
cat("### Rscript: processChXPedgeTest.R: Distribution of the average positive fraction:\n")
print(quantile(abs(data$posFraction), seq(0,1,.1)))
cat("### Rscript: processChXPedgeTest.R: Distribution of the region size:\n")
print(quantile(data$size, seq(0,1,.1)))
data$score <- data$posFraction*abs(data$leftBorder)*abs(data$rightBorder)
data$score <- data$score/max(data$score)*100
data$name <- paste0("peak_", 1:nrow(data))
cat("### Rscript: processChXPedgeTest.R: Distribution of the compound score:\n")
print(quantile(data$score, seq(0,1,.1)))
write.table(data[,c("chrom", "start", "end", "name", "score")], outputFile, sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
cat("### Rscript: processChXPedgeTest.R: Done.\n")