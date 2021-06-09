library(dplyr)
library(tidyverse)
library(progress)
library(gtools)

MutationTable <- matrix(0, nrow = 22, ncol = 223)

AAs <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P",
         "S", "T", "W", "Y", "V", "X", "*", "_")

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
sequences1 <- lapply(sequences[1:200000], try(tidyer, TRUE))
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






CorrelationTable <- matrix(0, nrow = 223, ncol = 223)
EitherOrTable <- matrix(0, nrow = 223, ncol = 223)
pb <- progress_bar$new(total = length(sequences1))
pb$tick(0)
baseline <- rep(0, length(spikeseq))
count <- 0
for (s in 1:length(sequences1)){
  sequence <- sequences1[[s]]
  if (any(is.na(sequence))) {next}
  count <- count + 1
  pb$tick()
  muts <- c()
  for (pos in 1:length(sequence)) {
    if (sequence[pos] =="X") {
      next
    }
    if (sequence[pos] != spikeseq[pos]){
      muts <- append(muts, pos)
      baseline[pos] <- baseline[pos]+1
    }
  }
  if (length(muts) > 1) {
    perms <- permutations(n=length(muts),r=2,v=muts,repeats.allowed=F)
    for (row in 1:nrow(perms)) {
      m1 <- perms[row, 1]
      m2 <- perms[row, 2]
      CorrelationTable[m1,m2] <- CorrelationTable[m1,m2] +1
    }
    for (mutation in muts) {
      EitherOrTable[,mutation] <- EitherOrTable[,mutation] + 1
      for (mutation2 in muts) {
        EitherOrTable[mutation2,mutation] <- EitherOrTable[mutation2,mutation] - 1
      }
    }
  }
  if (length(muts) == 1) {
    EitherOrTable[,muts[1]] <- EitherOrTable[,muts[1]] +1
    EitherOrTable[muts[1],muts[1]] <- 0
  }
}

phicoefficients <- matrix(0, nrow = 223, ncol = 223)
for (i in 1:nrow(CorrelationTable)) {
  for (j in 1:ncol(CorrelationTable)) {
    if (i==j) {next}
    i_on <- baseline[i]
    j_on <- baseline[j]
    i_off <- count - i_on
    j_off <- count - j_on
    both_on <- CorrelationTable[i,j]
    i_on_j_off <- EitherOrTable[j,i]
    j_on_i_off <- EitherOrTable[i,j]
    both_off <- count - (both_on + i_on_j_off +j_on_i_off)
    #phi <- ((both_on*both_off)-(i_on_j_off*i_on_j_off))/(sqrt(i_on*i_off*j_on*j_off))
    phi <- (((i_on + i_off)*both_on)-(i_on*j_on))/sqrt(i_on*j_on*(count-i_on)*(count-j_on))
    phicoefficients[i,j] <- phicoefficients[j,i] <- phi
  }
}

heatmapmaker <- function(matrix, names1, names2) {
  par(pty="s")
  library(RColorBrewer)
  dimnames(matrix) <- list(as.character(names1), as.character(names2))
  newdf <- as_tibble(matrix, rownames = NA)
  newdf <- cbind(firstletter = as.character(names1), newdf)
  newdf <- pivot_longer(newdf, !firstletter)
  newdf$firstletter <- as.character(newdf$firstletter)
  names(newdf)[1] <- "n"
  names(newdf)[2] <- "t"
  names(newdf)[3] <- "variable"
  column2 <- rep(names2, length(names1))
  newdf$t <- as.character(column2)
  newdf$t <- as.numeric(newdf$t)
  newdf$n <- as.numeric(newdf$n)
  uniquepi <- sort(unique(newdf$variable))
  newdf <- data.frame(cbind(newdf, rep(0, nrow(newdf))))
  colors <- c(rainbow(length(uniquepi)), "black")
  names(newdf)[4] <- "color"
  for (i in 1:nrow(newdf)){
    newdf[i, 4] <- match(newdf[i,3], uniquepi)
  }
  plot(newdf$t, newdf$n, col = colors[length(colors) + 1 - newdf$color], 
       pch = 15, cex = 1.6, ylab = "Locus (AA)", xlab = "Locus (AA)", xaxs = "i", yaxs = "i")
}

binaryphicoefficients <- phicoefficients > 0.5


heatmapmaker(CorrelationTable, 319:541, 319:541)
heatmapmaker(phicoefficients, 319:541, 319:541)
heatmapmaker(binaryphicoefficients, 319:541, 319:541)

write.csv(CorrelationTable, "prelimcorrelationanalysis.csv")

