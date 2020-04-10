library(ape)
library(MonoPhy)
library(phyloch)
library(phyloU)

setwd("/4_RAxML")
getwd() -> wd

###################################################
### Tribes
###################################################

setwd("/1_Aligns/3_singles_filtered")

read.csv("1 Markers_table_final.csv") -> meta
data.frame(Terminal=meta$spp, Clade=meta$Clade) -> clades
rownames(clades) <- clades$Terminal

setwd(wd)

###################################################
### Tree
###################################################

read.tree("RAxML_bipartitions.concatenated.rooted.tre") -> tree
ladderize(tree) -> tree

checkPhylo(tree, clades) -> c0
c0$data -> clades

Ntip(tree) == nrow(clades)

###################################################
### Check monophyly
###################################################

mono <- AssessMonophyly(tree, clades)

GetSummaryMonophyly(mono) -> mono.s
GetResultMonophyly(mono) -> mono.a
GetIntruderTips(mono) -> intruders
GetOutlierTips(mono) -> outliers
mono.s
mono.a

write.csv(mono.a, "monophyly_summary.csv")
write.csv(data.frame(unlist(outliers)), "monophyly_outliers.csv")

###################################################
### Node table and Plot
###################################################

GetAncNodes(mono) -> nodes
nodes$Clade -> nodes
as.numeric(as.character(nodes$MRCA)) -> nodes$MRCA
na.omit(nodes) -> nodes

cols <- rainbow(nrow(nodes))
sample(cols, length(cols), replace=F) -> cols

groups <- vector("list", length = nrow(nodes))
for (i in 1:nrow(nodes)) {
  descendants(tree, node=nodes$MRCA[i], labels=T) -> groups[[i]]
}

edge.color(tree, groups=groups, what="crown", col=cols) -> e.cols

plot(tree, show.tip.label = F, edge.color = e.cols)
