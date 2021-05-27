library(dplyr)
library(tidyverse)
library(progress)

MutationTable <- matrix(0, nrow = 22, ncol = 223)

AAs <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P",
         "S", "T", "W", "Y", "V", "X", "*")

admut <- function(pos, seq, MutationTable) {
  AAs <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P",
           "S", "T", "W", "Y", "V", "X", "*")
  column <- pos
  row <- match(seq[pos], AAs)
  if (is.na(row)) {
    row <- 22
  }
  MutationTable[row, column] <- MutationTable[row, column] + 1
  return(MutationTable)
}

tidyer <- function(s) {
  unlist(str_split(s, "', '"))[319:541]
}

sequences <- readLines("correctsequences.txt")
sequences1 <- lapply(sequences[1:100000], try(tidyer, TRUE))
pb <- progress_bar$new(total = length(sequences1))
pb$tick(0)
for (s in 1:length(sequences1)) {
  sequence <- sequences1[[s]]
  pb$tick()
  for (p in 1:length(sequence))
  MutationTable <- admut(p, sequence, MutationTable)
}

write.csv(MutationTable, "mutationtableprelim.csv")

n <- 22
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colors <- sample(col_vector, 22, replace = F)

initialiseRBD("Mutation Counts", counts, 4)
spikeseq <- paste0(readLines("spike.txt"), collapse = "")
spikeseq <- unlist(str_split(spikeseq, ""))
spikeseq <- spikeseq[319:541]
spikeseq <- match(spikeseq, AAs)
for (pos in 1:ncol(MutationTable)) {
  x <- pos + 318
  y1 <- 0
  for (AA in 1:(nrow(MutationTable)-1)) {
    if (AA == spikeseq[pos]) {
      next
    }
    color <- colors[AA]
    y2 <- y1 + MutationTable[AA, pos]
    lines(c(x,x), c(y1, y2), col = color, lwd = 6, lend = 1)
    y1 <- y2 +1
  }
  text(x, (y1+1), AAs[spikeseq[pos]], cex = 0.5)
}
legend(350, 10000, legend=AAs[1:21], pch = rep(15, 21), col = colors[1:21], horiz = T)

length(unique(colors))

for (pos in 1:ncol(MutationTable)) {
  x <- pos + 318
  y1 <- 0
  y2 <- 2000
  lines(c(x,x), c(y1, y2), col = colors[spikeseq[pos]], lwd = 6, lend = 1)

}