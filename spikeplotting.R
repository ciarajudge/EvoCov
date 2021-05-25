library(dplyr)
library(tidyverse)

initialisespike <- function(ylabel, dataset, type) {
  if (type == 1) {
    limit = max(dataset)
    addition1 <- 70000
    addition2 <- 30000
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
  plot(c(1, 1273), c(0, 0), ylim = c(0, (limit+addition1)), ylab = ylabel, pch = ".", xlab = "locus (by codon)", xaxs="i", yaxs="i")
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
initialiseRBD <- function(ylabel, dataset, type) {
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
    limit = max(dataset) + 1
    addition1 <- 2
    addition2 <- 1
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
addmutations <- function(b117, b1351, p1, b1427, b1429) {
  if (b117 == 1) {
    abline(v=69, col = "red")
    abline(v=70, col = "red")
    abline(v=144, col = "red")
    abline(v=484, col = "green")
    abline(v=494, col = "green")
    abline(v=501, col = "blue")
    abline(v=570, col = "blue")
    abline(v=614, col = "blue")
    abline(v=681, col = "blue")
    abline(v=716, col = "blue")
    abline(v=982, col = "blue")
    abline(v=1118, col = "blue")
    abline(v=1191, col = "green")
  }
  if (b1351 == 1) {
    abline(v=241, col = "red")
    abline(v=242, col = "red")
    abline(v=243, col = "red")
    abline(v=80, col = "blue")
    abline(v=215, col = "blue")
    abline(v=417, col = "blue")
    abline(v=484, col = "blue")
    abline(v=501, col = "blue")
    abline(v=614, col = "blue")
    abline(v=701, col = "blue")
  }
  if (p1 == 1) {
    abline(v=18, col = "blue")
    abline(v=20, col = "blue")
    abline(v=26, col = "blue")
    abline(v=138, col = "blue")
    abline(v=190, col = "blue")
    abline(v=417, col = "blue")
    abline(v=484, col = "blue")
    abline(v=501, col = "blue")
    abline(v=614, col = "blue")
    abline(v=655, col = "blue")
    abline(v=1027, col = "blue")
  }
  if (b1427 == 1) {
    abline(v=452, col = "blue")
    abline(v=614, col = "blue")
  }
  if (b1429 == 1) {
    abline(v=13, col = "blue")
    abline(v=152, col = "blue")
    abline(v=452, col = "blue")
    abline(v=614, col = "blue")
  }
}

counts <- as.numeric(readLines("fixedcounts.txt"))
plot(1:length(counts), counts, type = "l")

initialisespike("Number of Mutations", counts, 1)
lines(1:length(counts), counts, type = "l")

frequencies <- counts/1600000
initialisespike("Frequency mutated", frequencies, 2)
lines(1:length(counts), frequencies, type = "l")
addmutations(1,1,1,1,1)

initialiseRBD("Frequency mutated", frequencies, 2)
lines(1:length(counts), frequencies, type = "l")
addmutations(1,1,1,1,1)

uniquecounts <- as.numeric(readLines("uniquecounts.txt"))
initialisespike("Number of Mutation Events", uniquecounts, 3)
lines(1:length(counts), uniquecounts, type = "l")

