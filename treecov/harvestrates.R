library(tidyverse)

if (file.exists("rates.csv")){
  master <- read.csv("rates.csv", header = F)
  rates <- read.table("rates", sep = "", skip = 9, header = F, fill = TRUE)
  sitewiserates <- rates[1:3822,4]
  combined <- cbind(master, sitewiserates)
  write.csv(combined, "rates.csv")
} else {
  rates <- read.table("rates", sep = "", skip = 9, header = F, fill = TRUE)
  combined <- cbind(rates[1:3822,1], rates[1:3822,4])
  write.csv(combined, "rates.csv")
}