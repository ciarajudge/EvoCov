library(ape)

tree <- read.tree("global.tree")

for (i in 1:10){
  labelsample <- sample(tree$tip.label, 500)
  subsettedtree <- keep.tip(tree, labelsample)
  write.tree(subsettedtree, paste0(c("trees/sampledtree",as.character(i),".nhx"), collapse = ""))
  fileConn <- file(paste0(c("labels/subsetlabels",as.character(i),".txt"), collapse = ""))    
  writeLines(labelsample, fileConn)    
  close(fileConn)
}


