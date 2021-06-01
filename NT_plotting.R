library(dplyr)
library(tidyverse)
library(RColorBrewer)

initialisespike <- function(ylabel, dataset, type) {
  if (type == 1) {
    limit = max(dataset)
    addition1 <- 160000
    addition2 <- 70000
  }
  else if (type == 2) {
    limit = max(dataset) + 0.1
    addition1 <- 0.1
    addition2 <- 0.05
  }
  else if (type == 3) {
    limit = max(dataset) + 1
    addition1 <- 2
    addition2 <- 1
  }
  plot(c(1, 3822), c(0, 0), ylim = c(0, (limit+addition1)), ylab = ylabel, 
       pch = ".", xlab = "genomic locus", xaxs="i", yaxs="i")
  rect(36, 0, 912, limit, 
       col = adjustcolor("red", alpha.f = 0.2), border = NA)
  rect(957, 0, 1623, limit, 
       col = adjustcolor("green", alpha.f = 0.2), border = NA)
  rect(1624, 0, 3624, limit, 
       col = adjustcolor("blue", alpha.f = 0.2), border = NA)
  lines(c(36, 1623), c(limit, limit), lwd = 5, col = "navy")
  lines(c(1624, 3624), c(limit, limit), lwd = 5, col = "red")
  text(200, (limit+addition2), "S1 Subunit")
  text(1750, (limit+addition2), "S2 Subunit")
}
initialiseRBD <- function(ylabel, dataset, type) {
  if (type == 1) {
    limit = max(dataset[957:1623]) + 10000
    addition1 <- 0
    addition2 <- 0
  }
  else if (type == 2) {
    limit = max(dataset[319:541])+0.05
    addition1 <- 0
    addition2 <- 0
  }
  else if (type == 3) {
    limit = max(dataset[319:541])+1
    addition1 <- 0
    addition2 <- 0
  }
  else if (type == 4) {
    limit = 100000
    addition1 <- 0
    addition2 <- 0
  }
  plot(c(1, 1273), c(0, 0), xlim = c(957,1623),ylim = c(0, (limit+addition1)), 
       ylab = ylabel, pch = ".", xlab = "genomic locus", xaxs="i", yaxs="i")
  rect(957, 0, 1623, limit, 
       col = adjustcolor("green", alpha.f = 0.2), border = NA)
  
}


counts <- as.numeric(readLines("nt_counts.txt"))

initialisespike("Number of Mutations", counts, 1)
lines(1:length(counts), counts, type = "l")

initialiseRBD("Number of Mutations", counts, 1)
lines(1:length(counts), counts, type = "l")

bases <- c("a","c","t","g","n","-")
initialiseRBD("Mutation Counts", counts, 4)
spikeseq <- paste0(readLines("spike_nt.txt"), collapse = "")
spikeseq <- unlist(str_split(spikeseq, ""))
colors <- c("red", "green", "blue", "orange", "purple", "pink")
for (pos in 957:1623) {
  x <- pos
  y1 <- 0
  for (nt in 1:(ncol(NTtable))) {
    color <- colors[nt]
    y2 <- y1 + NTtable[pos, nt]
    lines(c(x,x), c(y1, y2), col = color, lwd = 2, lend = 1)
    y1 <- y2 +1
  }
  text(x, (y1+1), bases[spikeseq[pos]], cex = 0.5)
}






