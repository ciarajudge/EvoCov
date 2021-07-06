library(tidyverse)

epitope1 <- read.csv("epitope1.csv", header = F)
epitope2 <- read.csv("epitope2.csv", header = F)
epitope1versions <- readLines("epitope1versions.txt")
epitope2versions <- readLines("epitope2versions.txt")
countries <- readLines("countries.txt")

dimnames(epitope1) <- list(epitope1versions,countries)
dimnames(epitope2) <- list(epitope2versions,countries)


