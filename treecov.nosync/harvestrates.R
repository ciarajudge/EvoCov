library(tidyverse)
library(argparser)

args <- commandArgs(trailingonly = TRUE)
if (args[1] == "iteration") {
if (file.exists("rates.csv")){
  master <- read.csv("rates.csv", header = T)
  rates <- read.table("rates", sep = "", skip = 9, header = F, fill = TRUE)
  sitewiserates <- rates[1:3822,4]
  combined <- cbind(master, sitewiserates)
  write.csv(combined, "rates.csv", row.names = FALSE)
} else {
  rates <- read.table("rates", sep = "", skip = 9, header = F, fill = TRUE)
  combined <- cbind(rates[1:3822,1], rates[1:3822,4])
  write.csv(combined, "rates.csv")
}
}
if (args[1]=="final"){
  rates <- read.csv("rates.csv")
  rates <- rates[,2:ncol(rates)]
  group1 <- sample(ncol(rates), round(0.5*ncol(rates)))
  group2 <- (1:ncol(rates))[! (1:ncol(rates)) %in% group1]
  rates1 <- rates[,group1]
  rates2 <- rates[,group2]
  
  l1 <- c()
  l2 <- c()
  total <- c()
  for (row in 1:nrow(rates)){
    l1[row] <- median(unlist(rates1[row,]))
    l2[row] <- median(unlist(rates2[row,]))
    total[row] <- median(unlist(rates[row,]))
  }
  plot(l1, l2)
  total <- as.character(total)
  fileConn <- file("evorates.txt") 
  writeLines(total, fileConn)    
  close(fileConn)
}




