library(tidyverse)
library(hash)

epitope1 <- read.csv("Analysis/Epitopes/ByCountry/YYLQSYGFT.csv", header = F)
epitope2 <- read.csv("Analysis/Epitopes/ByCountry/WNRKRIS.csv", header = F)

countries <- readLines("Helpers/countries.txt")

epitope1versions <- epitope1[,1]
epitope2versions <- epitope2[,1]
epitope1 <- epitope1[,2:ncol(epitope1)]
epitope2 <- epitope2[,2:ncol(epitope2)]

dimnames(epitope1) <- list(epitope1versions,countries)
dimnames(epitope2) <- list(epitope2versions,countries)

WHOdata <- read.csv("/Users/ciarajudge/Desktop/WHO COVID-19 global table data July 9th 2021 at 3.51.09 PM.csv")
load("Helpers/countrydict.RData")
countrydata <- read.csv("epitope1wcountrydata.csv")

epitope1 <- rbind.data.frame(epitope1, colSums(epitope1))
epitope2 <- rbind.data.frame(epitope2, colSums(epitope2))
epitope1 <- rbind.data.frame(epitope1, rep(0, ncol(epitope1)))
epitope2 <- rbind.data.frame(epitope2, rep(0, ncol(epitope1)))

for (i in 2:nrow(WHOdata)){
  breaking <- F
  WHOcountry <- WHOdata[i, 1]
  GISAIDcountry <- suppressWarnings(dictionary[[WHOcountry]])
  ColumnNo <- match(GISAIDcountry, countries)
  if (is.null(GISAIDcountry)){
    next
  }
  epitope1[nrow(epitope1), ColumnNo] <- epitope1[nrow(epitope1), ColumnNo] + WHOdata[i,5]*5
  epitope2[nrow(epitope2), ColumnNo] <- epitope2[nrow(epitope2), ColumnNo] + WHOdata[i,5]*5
}

write.csv(epitope1, "epitope1wcountrydata.csv")
write.csv(epitope2, "epitope2wcountrydata.csv")


epitope1 <- epitope1[unlist(epitope1[nrow(epitope1), ]) != 0]
epitope2 <- epitope2[unlist(epitope2[nrow(epitope2), ]) != 0]

epitope1matrix <- matrix(0, nrow = nrow(epitope1)-2, ncol = ncol(epitope1)+1)
epitope2matrix <- matrix(0, nrow = nrow(epitope2)-2, ncol = ncol(epitope2)+1)


for (j in 1:(ncol(epitope1))){
  for (i in 1:(nrow(epitope1)-2)){
    if (epitope1[51,j]!=0){
    epitope1matrix[i,j] <- (epitope1[i,j]/epitope1[51,j])*epitope1[52,j]
    }
    else {
      epitope1matrix[i,j] <- 0
    }
  }
}

epitope1matrix[,ncol(epitope1matrix)] <- rowSums(epitope1matrix)
finalfreqsepi1 <- matrix(0, nrow = nrow(epitope1matrix)+1, ncol = 1)
finalfreqsepi1[,1] <- unlist(c(epitope1matrix[,nrow(epitope1matrix)]/sum(epitope1matrix[,nrow(epitope1matrix)]),0))
finalfreqsepi1[51,1] <- sum(finalfreqsepi1[,1])


for (j in 1:(ncol(epitope2))){
  for (i in 1:(nrow(epitope2)-2)){
    if (epitope2[54,j] != 0){
      epitope2matrix[i,j] <- (epitope2[i,j]/epitope2[54,j])*epitope2[55,j]
    }
    else{
      epitope2matrix[i,j] <- 0
    }
  }
}

epitope2matrix[,ncol(epitope2matrix)] <- rowSums(epitope2matrix)
finalfreqsepi2 <- matrix(0, nrow = nrow(epitope2matrix)+1, ncol = 1)
finalfreqsepi2[1:53,1] <- epitope2matrix[,ncol(epitope2matrix)]/sum(epitope2matrix[,ncol(epitope2matrix)])
finalfreqsepi2[54,1] <- sum(finalfreqsepi2[,1])

dimnames(finalfreqsepi1) <- list(unlist(c(epitope1versions, "total")),"Frequency")
dimnames(finalfreqsepi2) <- list(unlist(c(epitope2versions, "total")),"Frequency")


write.csv(epitope1matrix, "Epitope1NormFreqs.csv")
write.csv(epitope2matrix, "Epitope2NormFreqs.csv")

for (i in 2:nrow(WHOdata)){
  breaking <- F
  WHOcountry <- WHOdata[i, 1]
  GISAIDcountry <- match(WHOcountry, countries)
  while (is.na(GISAIDcountry)){
    WHOcountry <- readline(prompt = paste0(c(WHOcountry, "\n"),collapse = ""))
    GISAIDcountry <- match(WHOcountry, countries)
    if (WHOcountry==""){
      breaking <- T
      break
    }
  }
  if (breaking == T){
    next
  }
  epitope1[nrow(epitope1), GISAIDcountry] <- WHOdata[i,3]
  epitope2[nrow(epitope2), GISAIDcountry] <- WHOdata[i,3]
}


for (x in 1:length(keys(dictionary))){
  newtable[x, 1] <- keys(dictionary)[x]
  newtable[x,2] <- dictionary[[keys(dictionary)[x]]]
}

write.csv(newtable, "dictastable.csv")
