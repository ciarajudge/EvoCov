library(ape)

tree <- read.tree("global.tree")

labelsample <- sample(tree$tip.label, 100)
subsettedtree <- keep.tip(tree, labelsample)
write.tree(subsettedtree, "sampledtree.nhx")

fileConn <- file("subsetlabels.txt")    
writeLines(labelsample, fileConn)    
close(fileConn)
