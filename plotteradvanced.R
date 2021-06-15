library(dplyr)
library(tidyverse)
library(progress)
library(gtools)
library(RColorBrewer)

MutationTable <- read.csv("Analysis/countstable_NT.csv", header = F)/902551
AAs <- c("A", "C", "T", "G", "N", "_")

AAs <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P",
         "S", "T", "W", "Y", "V", "X", "*", "_")
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
    limit = max(dataset[319:541])+1
    addition1 <- 0
    addition2 <- 0
  }
  else if (type == 4) {
    limit = 10000
    addition1 <- 0
    addition2 <- 0
  }
  plot(c(1, 1273), c(0, 0), xlim = c(500,520),ylim = c(0, 0.005), 
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


counts <- rowSums(MutationTable)

n <- length(AAs)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colors <- sample(col_vector, n, replace = F)

dev.new(width=12, height=4)
par(margin(1, 0, 0, 0))
initialiseRBD("Frequencies", counts, 2)
spikeseq <- paste0(readLines("Data/spike_NT.txt"), collapse = "")
spikeseq <- unlist(str_split(spikeseq, ""))
spikeseq <- spikeseq
spikeseq <- match(spikeseq, AAs)
for (pos in 2:nrow(MutationTable)) {
  x <- pos
  y1 <- 0
  for (AA in 1:(ncol(MutationTable)-1)) {
    color <- colors[AA]
    y2 <- y1 + MutationTable[pos, AA]
    lines(c(x,x), c(y1, y2), col = color, lwd = 6, lend = 1)
    y1 <- y2
  }
  if (counts[pos] > 0){
    text(x, (y1+0.0001), AAs[spikeseq[x]], cex = 0.7)
  }
}
legend("top", inset = c(0, -0.2),legend=AAs, pch = 15, col = colors,
       xpd = TRUE, horiz = T)



#for (pos in 1:ncol(MutationTable)) {
#  x <- pos + 318
#  y1 <- 0
#  y2 <- 2000
#  lines(c(x,x), c(y1, y2), col = colors[spikeseq[pos]], lwd = 6, lend = 1)
#}

#####Temporal#####
time20194 <- read.csv("Analysis/date/20194frequencies.csv", header = FALSE)
time20201 <- read.csv("Analysis/date/20201frequencies.csv", header = FALSE)
time20202 <- read.csv("Analysis/date/20202frequencies.csv", header = FALSE)
time20203 <- read.csv("Analysis/date/20203frequencies.csv", header = FALSE)
time20204 <- read.csv("Analysis/date/20204frequencies.csv", header = FALSE)
time20211 <- read.csv("Analysis/date/20211frequencies.csv", header = FALSE)
time20212 <- read.csv("Analysis/date/20212frequencies.csv", header = FALSE)

time20194counts <- rowSums(time20194)
time20201counts <- rowSums(time20201)
time20202counts <- rowSums(time20202)
time20203counts <- rowSums(time20203)
time20204counts <- rowSums(time20204)
time20211counts <- rowSums(time20211)
time20212counts <- rowSums(time20212)

par(mar=c(4, 4, 5, 4))
initialiseRBD("Frequency Mutated", time20212counts, 2)
points(1:length(time20194counts), time20194counts, pch = 12, col = "navy", cex = 1.5)
points(1:length(time20194counts), time20201counts, pch = 13, col = "orange", cex = 1.5)
points(1:length(time20194counts), time20202counts, pch = 13, col = "darkorchid", cex = 1.5)
points(1:length(time20194counts), time20203counts, pch = 13, col = "deeppink", cex = 1.5)
points(1:length(time20194counts), time20204counts, pch = 13, col = "navy", cex = 1.5)
points(1:length(time20194counts), time20211counts, pch = 14, col = "orange", cex = 1.5)
points(1:length(time20194counts), time20212counts, pch = 14, col = "darkorchid", cex = 1.5)

legend("top", inset = c(0, -0.2),legend=c("2019Q4", "2020Q1","2020Q2", "2020Q3","2020Q4","2021Q1", "2021Q2"), 
       pch = c(12, 13, 13, 13, 13, 14, 14), col = c("navy", "orange", "darkorchid","deeppink", "navy", "orange", "darkorchid"),xpd = TRUE, horiz = TRUE)


par(xpd = F)


##### Variant #####
B117 <- read.csv("Analysis/variant/B.1.1.7frequencies.csv", header = FALSE)
B1351 <- read.csv("Analysis/variant/B.1.351frequencies.csv", header = FALSE)
B1427 <- read.csv("Analysis/variant/B.1.427frequencies.csv", header = FALSE)
B1429 <- read.csv("Analysis/variant/B.1.429frequencies.csv", header = FALSE)
P1 <- read.csv("Analysis/variant/P.1frequencies.csv", header = FALSE)
B16172 <- read.csv("Analysis/variant/B.1.617.2frequencies.csv", header = FALSE)

B117counts <- rowSums(B117)
B1351counts <- rowSums(B1351)
B1427counts <- rowSums(B1427)
B1429counts <- rowSums(B1429)
B16172counts <- rowSums(B16172)
P1counts <- rowSums(P1)

initialisespike("Frequency Mutated", B117counts, 2)
points(1:length(P1counts), B117counts, col = "darkorchid")
points(1:length(P1counts), B1351counts, col = "navy")
points(1:length(P1counts), P1counts, col = "red")
points(1:length(P1counts), B1427counts, col = "cyan2")
points(1:length(P1counts), B1429counts, col = "darkolivegreen")
points(1:length(P1counts), B16172counts, col = "darkorange")
legend("top", inset = c(0, -0.2),legend=c("B.1.1.7", "B.1.351","P.1", "B.1.427","B.1.429","B.1.617.2"), 
       pch = 1, col = c("darkorchid", "navy", "red","cyan2", "darkolivegreen", "darkorange"),xpd = TRUE, horiz = TRUE)

initialisespike("Frequency Mutated", B117counts, 2)
points(1:length(P1counts), B117[,5], col = "darkorchid")
points(1:length(P1counts), B1351[,5], col = "navy")
points(1:length(P1counts), P1[,5], col = "red")
points(1:length(P1counts), B1427[,5], col = "darkslateblue")
points(1:length(P1counts), B1429[,5], col = "darkolivegreen")
points(1:length(P1counts), B16172[,5], col = "darkorange")
legend("top", inset = c(0, -0.2),legend=c("B.1.1.7", "B.1.351","P.1", "B.1.427","B.1.429","B.1.617.2"), 
       pch = 1, col = c("darkorchid", "navy", "red","darkslateblue", "darkolivegreen", "darkorange"),xpd = TRUE, horiz = TRUE)


