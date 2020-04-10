library(ape)
library(MonoPhy)
library(phyloch)
library(phyloU)

setwd("/6_RAxML_rogueless")
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

read.nexus("RAxML_bipartitions.concatenated.rooted.tre") -> tree
ladderize(tree) -> tree
Ntip(tree) == nrow(clades)

checkPhylo(tree, clades) -> c0
c0$data -> clades

###################################################
### Check monophyly
###################################################

mono <- AssessMonophyly(tree, clades)

GetSummaryMonophyly(mono)

GetResultMonophyly(mono) -> mono.sum
mono.sum

write.csv(mono.sum, "monophyly_summary.csv")

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
