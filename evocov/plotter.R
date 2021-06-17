library(argparser)
library(dplyr)
library(tidyverse)
library(RColorBrewer)

initialisespikeNT <- function(ylabel, type) {
  if (type == 1) {
    ylimit = 1.1
    limit = 1
    addition2 <- 0.05
  }
  plot(c(1, 3822), c(0, 0), ylim = c(0, ylimit), ylab = ylabel, 
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
initialiseRBDNT <- function(ylabel) {
  plot(c(1, 1273), c(0, 0), xlim = c(957,1623),ylim = c(0, 1), 
       ylab = ylabel, pch = ".", xlab = "genomic locus", xaxs="i", yaxs="i")
  rect(957, 0, 1623, 1, 
       col = adjustcolor("green", alpha.f = 0.1), border = NA)
  text(1100, (0.9), "Receptor Binding Domain", cex = 1.5, col = "navy")
}
initialisespikeAA <- function(ylabel, dataset, type) {
  if (type == 1) {
    ylimit = 1.1
    limit = 1
    addition2 <- 0.05
  }
  plot(c(1, 1273), c(0, 0), ylim = c(0, (ylimit)), ylab = ylabel, pch = ".", xlab = "locus (by codon)", xaxs="i", yaxs="i")
  rect(13, 0, 304, limit, 
       col = adjustcolor("red", alpha.f = 0.2), border = NA)
  rect(319, 0, 541, limit, 
       col = adjustcolor("green", alpha.f = 0.2), border = NA)
  rect(543, 0, 1208, limit, 
       col = adjustcolor("blue", alpha.f = 0.2), border = NA)
  lines(c(12, 541), c(limit, limit), lwd = 5, col = "navy")
  lines(c(543, 1208), c(limit, limit), lwd = 5, col = "red")
  text(70, (limit+addition2), "S1 Subunit")
  text(600, (limit+addition2), "S2 Subunit")
}
initialiseRBDAA <- function(ylabel, dataset, type) {
  if (type == 1) {
    limit = max(dataset)
    addition1 <- 70000
    addition2 <- 30000
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
    limit = 10000
    addition1 <- 0
    addition2 <- 0
  }
  plot(c(1, 1273), c(0, 0), xlim = c(319,541),ylim = c(0, (limit+addition1)), 
       ylab = ylabel, pch = ".", xlab = "locus (by codon)", xaxs="i", yaxs="i")
  rect(13, 0, 304, limit, 
       col = adjustcolor("red", alpha.f = 0.2), border = NA)
  rect(319, 0, 541, limit, 
       col = adjustcolor("green", alpha.f = 0.2), border = NA)
  rect(543, 0, 1208, limit, 
       col = adjustcolor("blue", alpha.f = 0.2), border = NA)
}

args <- commandArgs(trailingOnly=TRUE)
nttable <- read.csv(args[1], header = F)
ntnons <- cbind(nttable[,1:4], nttable[,6])
ntcounts <- rowSums(nttable)
ntcountsnons <- rowSums(ntnons)

checkamplicons <- function(Ns){
  ampstarts <- c()
  ampstops <- c()
  inamplicon <- T
  counter <- 0
  for (i in Ns){
#    if (i>0.)
  }
}

pdf(file = "results.pdf", paper = "a4", width = 7, height = 10.5)
layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8), ncol = 1), heights = c(1,0.6,0.5,1,4,1, 4, 3))
par(mar=c(0,0,0,0))
plot(1, 1, col = "white", xaxt = "n", bty = "n", yaxt = "n")
text(1,1,"Evocov Epitope Selection Pipeline Results", cex = 2.5)
par(mar=c(0,0,0,0))
plot(1,1, col = "white", xaxt = "n", bty = "n", yaxt = "n")
text(1,1.1, cex = 1.5, paste0(c("Thank you for using the Evicov Evolution pipeline for SARS-CoV-2! The following \n results are based on ", 
                   args[2], " genome sequences from GISAID."), collapse = ""), font = 3)
par(mar=c(0,0,0,0))
plot(1, 1, col = "white", xaxt = "n", bty = "n", yaxt = "n")
text(1,0.8,"General Mutational Landscape by Nucleotide", cex = 2)
plot(1, 1, col = "white", xaxt = "n", bty = "n", yaxt = "n")
text(1,1, "The following plot is a general mutational landscape for the SARS-CoV-2 spike protein nucleotide sequence for the duration of the pandemic 
(i.e. using all available parsed sequences from GISAID).For all sequences, at each locus, if the nucleotide differed from the reference, the 
     count for that locus was increased by 1. The final counts for each locus were divided by the total number of sequences to give frequencies.")


par(mar=c(4.5, 4.5, 1, 1))
initialisespikeNT("Frequency Mutated", 1)
lines(1:length(ntcounts), ntcounts)
par(mar=c(0,0,0,0))
plot(1, 1, col = "white", xaxt = "n", bty = "n", yaxt = "n")
text(1,1, "For the above graph, N's are included as a difference from the reference, to allow instances of amplicon dropping to be detected.
     Specific suggestions of amplicon dropping in the data will be further elaborated on later in the variant-wise analysis. Below, the graph 
     is recreated without counting Ns, and a zoomed in view is shown of the receptor binding domain.")

par(mar=c(4.5, 4.5, 1, 1))
initialisespikeNT("Frequency Mutated", 1)
lines(1:length(ntcountsnons), ntcountsnons)

par(mar=c(4.5, 4.5, 1, 1))
initialiseRBDNT("Frequency Mutated")
lines(1:length(ntcountsnons), ntcountsnons)

