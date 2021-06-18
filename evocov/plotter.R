library(argparser)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(stringr)

n <- 22
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
AAcolors <- sample(col_vector, 22, replace = F)

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
initialisespikeAA <- function(ylabel, type) {
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
initialiseRBDAA <- function(ylabel, ylimit, main = "") {
  plot(c(1, 1273), c(0, 0), xlim = c(319,541),ylim = c(0, ylimit), main = main,
         ylab = ylabel, pch = ".", xlab = "locus (by codon)", xaxs="i", yaxs="i")
  rect(319, 0, 541, ylimit, 
       col = adjustcolor("green", alpha.f = 0.2), border = NA)
}
stackedbarAA <- function(MutationTable, reference, counts) {
  AAs <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P",
           "S", "T", "W", "Y", "V", "X", "*", "_")
  spikeseq <- paste0(readLines(reference), collapse = "")
  spikeseq <- unlist(str_split(spikeseq, ""))
  spikeseq <- match(spikeseq, AAs)
  for (pos in 300:600) {
    x <- pos-1
    y1 <- 0
    for (AA in 1:20) {
      color <- AAcolors[AA]
      y2 <- y1 + MutationTable[pos, AA]
      lines(c(x,x), c(y1, y2), col = color, lwd = 3, lend = 1)
      y1 <- y2
    }
    if (counts[pos] > 0){
      text(x, (y1+0.01), AAs[spikeseq[x]], cex = 0.4)
    }
  }
  legend("top", inset = c(0, -0.2),legend=AAs[1:20], pch = 15, col = AAcolors[1:20],
         xpd = TRUE, horiz = T)
  

}
stackedbarAA2 <- function(MutationTable, reference, counts) {
  AAs <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P",
           "S", "T", "W", "Y", "V", "X", "*", "_")
  spikeseq <- paste0(readLines(reference), collapse = "")
  spikeseq <- unlist(str_split(spikeseq, ""))
  spikeseq <- match(spikeseq, AAs)
  for (pos in 300:600) {
    x <- pos-1
    y1 <- 0
    for (AA in 1:20) {
      color <- AAcolors[AA]
      y2 <- y1 + MutationTable[pos, AA]
      lines(c(x,x), c(y1, y2), col = color, lwd = 3, lend = 1)
      y1 <- y2
    }
    if (counts[pos] > 0){
      text(x, (y1+0.001), AAs[spikeseq[x]], cex = 0.5)
    }
  }
  
  
}
addfromfile <- function(filepath, type, dot, color) {
  table <- read.csv(filepath, header = T)
  counts <- rowSums(table)
  if (type == "l"){
    lines(1:length(counts), counts, col = color)
  }
  else if (type  == "p"){
    points(1:length(counts), counts, pch = dot, col = color)
  }
}

args <- commandArgs(trailingOnly=TRUE)
nttable <- read.csv(args[1], header = F)
ntnons <- cbind(nttable[,1:4], nttable[,6])
ntcounts <- rowSums(nttable)
ntcountsnons <- rowSums(ntnons)
aatable <- read.csv(args[3], header = F)
aacounts <- rowSums(aatable)


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
##Page 1
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

print("check")
##Page 2
layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8), ncol = 1), heights = c(1,2,5,5, 2, 2))
par(mar=c(0,0,0,0))
plot(1, 1, col = "white", xaxt = "n", bty = "n", yaxt = "n")
text(1,0.8,"Mutant Codon Frequency at Each Locus", cex = 2)
plot(1, 1, col = "white", xaxt = "n", bty = "n", yaxt = "n")
text(1,1, paste0(c("For each codon position in the RBD, a stacked bar is used to show the frequencies at that locus of codons that 
differ from the reference. Due to the heavily prevalent N501Y mutation (", AAcolors[19], ") the information in this plot is
difficult to view, so a zoomed in version is also provided."), collapse = ""))

par(mar=c(4.5, 4.5, 3, 1))
initialiseRBDAA("Frequency Mutated",  (max(aacounts[319:541])+0.05))
stackedbarAA(aatable, "Data/spike_AA.txt", aacounts)

par(mar=c(4.5, 4.5, 3, 1))
initialiseRBDAA("Frequency Mutated",  0.05)
stackedbarAA2(aatable, "Data/spike_AA.txt", aacounts)

par(mar=c(0,0,0,0))
plot(1, 1, col = "white", xaxt = "n", bty = "n", yaxt = "n")
text(1,1, "This information is also available for the genomic sequence rather than codon-wise, but is too wide to fit nicely into 
this PDF. For that reason, it has been saved separately in the 'Plots' folder. The information shown in these plots will be used
later in epitope scoring during the calculation of shannon entropies for each locus, and when the suggested epitope is selected,
a zoomed in version of this plot for the epitope in question will be presented once more.")
plot(1, 1, col = "white", xaxt = "n", bty = "n", yaxt = "n")
plot(1, 1, col = "white", xaxt = "n", bty = "n", yaxt = "n")
plot(1, 1, col = "white", xaxt = "n", bty = "n", yaxt = "n")

##Page 3
varfiles <- list.files("Analysis/variant")
details = file.info(list.files("Analysis/variant"))
details = details[with(details, order(as.POSIXct(ctime))), ]
varfiles = rownames(details)
variants <- c()
counts <- c()
for (var in varfiles) {
  v <- unlist(str_split(var, ".csv"))[1]
  c <- unlist (str_split(v, "_"))[2]
  v <- unlist (str_split(v, "_"))[1]
  variants <- append(variants, v)
  counts <- append(counts, c)
}
numvars <- length(variants)

layout(matrix(c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 7, 8, 9), ncol = 2, byrow = T), heights = c(1,0.8,4.2,1.5, 1, 4, 4))
par(mar=c(0,0,0,0))
plot(1, 1, col = "white", xaxt = "n", bty = "n", yaxt = "n")
text(1,0.8,"Variant Wise Analysis", cex = 2)
plot(1, 1, col = "white", xaxt = "n", bty = "n", yaxt = "n")
text(1,1, paste0(c("The following results are using the variants ", paste0(variants, collapse = ", "), 
                   " which 
                   may have been selected by you or chosen by default by the pipeline"), collapse = ""))

par(mar=c(4.5, 4.5, 3, 1))
initialisespikeAA("Frequency Mutated",  1)
Varcolors <- sample(brewer.pal(9,"Set1"), length(varfiles), replace = F)
for (x in 1:length(varfiles)) {
  addfromfile(paste0(c("Analysis/variant/",varfiles[x]), collapse = ""), "p", 19, Varcolors[x])
}
legend("top", inset = c(0, -0.2),legend=variants, pch = 19, col = Varcolors,
       xpd = TRUE, horiz = T)

par(mar=c(0,0,0,0))
plot(1, 1, col = "white", xaxt = "n", bty = "n", yaxt = "n")
text(1,1, "This plot can serve as a quality control check to ensure that the tabulation of mutations
has been carried out correctly. The mutations specific to a given variant should have a frequency of 
close to 1.0 in the variant in which it has been documented by the CDC, and closer to 0 in other variants.
For further clarity, the plots are separated out by variant below for the receptor binding domain (in the
interests of space). Ns have been included to provide evidence for instances of amplicon dropping specific
to some variants.")

par(mar=c(0,0,0,0))
plot(1, 1, col = "white", xaxt = "n", bty = "n", yaxt = "n")
text(1,1,"RBD Plots broken down by Variant", cex = 1.5, font =3)

graph <- 0

if ((graph != numvars)) {
  par(mar=c(4.5, 4.5, 3, 1))
  graph <- graph + 1
  initialiseRBDAA("Frequency Mutated", 1, paste0(c(variants[graph], ", based on ", counts[graph], " sequences"), collapse = ""))
  addfromfile(paste0(c("Analysis/variant/", varfiles[graph]), collapse = ""), 
              "l", 1, "black")
} else {
par(mar=c(0,0,0,0))
plot(1, 1, col = "white", xaxt = "n", bty = "n", yaxt = "n")
}

if ((graph != numvars)) {
  par(mar=c(4.5, 4.5, 3, 1))
  graph <- graph + 1
  initialiseRBDAA("Frequency Mutated", 1, paste0(c(variants[graph], ", based on ", counts[graph], " sequences"), collapse = ""))
  addfromfile(paste0(c("Analysis/variant/", varfiles[graph]), collapse = ""), 
              "l", 1, "black")
} else {
  par(mar=c(0,0,0,0))
  plot(1, 1, col = "white", xaxt = "n", bty = "n", yaxt = "n")
}

if ((graph != numvars)) {
  par(mar=c(4.5, 4.5, 3, 1))
  graph <- graph + 1
  initialiseRBDAA("Frequency Mutated", 1, paste0(c(variants[graph], ", based on ", counts[graph], " sequences"), collapse = ""))
  addfromfile(paste0(c("Analysis/variant/", varfiles[graph]), collapse = ""), 
              "l", 1, "black")
} else {
  par(mar=c(0,0,0,0))
  plot(1, 1, col = "white", xaxt = "n", bty = "n", yaxt = "n")
}

if ((graph != numvars)) {
  par(mar=c(4.5, 4.5, 3, 1))
  graph <- graph + 1
  initialiseRBDAA("Frequency Mutated", 1, paste0(c(variants[graph], ", based on ", counts[graph], " sequences"), collapse = ""))
  addfromfile(paste0(c("Analysis/variant/", varfiles[graph]), collapse = ""), 
              "l", 1, "black")
} else {
  par(mar=c(0,0,0,0))
  plot(1, 1, col = "white", xaxt = "n", bty = "n", yaxt = "n")
}

#Page 4
while (graph != numvars) {
  layout(matrix(c(1,1, 2, 3, 4, 5, 6, 7, 8, 9), ncol = 2, byrow = T), heights = c(1,5,5,5,5))
  par(mar=c(0,0,0,0))
  plot(1, 1, col = "white", xaxt = "n", bty = "n", yaxt = "n")
  text(1,0.8,"Variant Wise Analysis", cex = 2)
  while ((graph != numvars)) {
    par(mar=c(4.5, 4.5, 3, 1))
    graph <- graph + 1
    initialiseRBDAA("Frequency Mutated", 1, paste0(c(variants[graph], ", based on ", counts[graph], " sequences"), collapse = ""))
    addfromfile(paste0(c("Analysis/variant/", varfiles[graph]), collapse = ""), 
                "l", 1, "black")
  }
}

#Page 5
layout(matrix(c(1,2, 3, 4, 5), ncol = 1), heights = c(1,1,3,8))
par(mar=c(0,0,0,0))
plot(1, 1, col = "white", xaxt = "n", bty = "n", yaxt = "n")
text(1,0.8,"Epitope Scoring", cex = 2)
par(mar=c(0,0,0,0))
plot(1, 1, col = "white", xaxt = "n", bty = "n", yaxt = "n")
text(1,0.8,"How we carry out epitope scoring", cex = 1.5, font = 3)
par(mar=c(0,0,0,0))
plot(1, 1, col = "white", xaxt = "n", bty = "n", yaxt = "n")
text(1,0.8,"In assigning scores to epitopes, ", cex = 1)









