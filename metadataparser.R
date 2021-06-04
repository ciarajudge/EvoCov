library(dplyr)
library(tidyverse)
library(progress)

meta <- read.table("metadata.tsv",header = T, fill = T, sep = "\t")
write.csv(meta, "metadata.csv")

accessions <- meta[,3]
write.table(accessions, "accessions.txt", sep = "\n", row.names = F, col.names = F)
locations <- meta[,5]
write.table(locations, "locations.txt", sep = "\n", row.names = F, col.names = F)


dates <- meta[,4]
write.table(dates, "dates1.txt", sep = "\n", row.names = F, col.names = F)
years <- c()
quarters <- c()
pb <- progress_bar$new(total = 1006876)
pb$tick(0)
for (d in dates) {
  pb$tick()
  year <- unlist(str_split(d, "-"))[1]
  month <- as.numeric(unlist(str_split(d, "-"))[2])

  if (is.na(month)) {
    quarter = 0
  }
  else if (any(month == c(01,02,03))){
    quarter = 1
  }
  else if (any(month == c(04,05,06))){
    quarter = 2
  }
  else if (any(month == c(07,08,09))){
    quarter = 3
  }
  else {
    quarter = 4
  }
  years <- append(years, year)
  quarters <- append(quarters, as.numeric(paste0(c(year, quarter), collapse = "")))
}

write.table(quarters, "quarters.txt", sep = "\n", row.names = F, col.names = F)
write.table(years, "years.txt", sep = "\n", row.names = F, col.names = F)

strains <- meta[,12]
write.table(strains, "strains.txt", sep = "\n", row.names = F, col.names = F)


