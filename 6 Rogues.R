library(ape)
library(MonoPhy)
library(phyloch)
library(phyloU)
library(DECIPHER)
source("fastConc.R")
source("dstats.R")

setwd("/5_Rogues")
getwd() -> wd


###################################################
### Outliers
###################################################

read.csv("monophyly_outliers.csv", stringsAsFactors = F)[,2] -> outliers


###################################################
### Rogues
###################################################

file = "RogueNaRok_Maj_Con_concatenated_drop_1.txt"

read.table(file, header = T) -> rogues
rogues[-1,c(3:5)] -> rogues

as.character(rogues$taxon[which(rogues$rawImprovement > 1)]) -> rogues
sort(rogues) -> rogues

sort(unique(c(rogues,outliers))) -> drop

write.csv(drop, "Rogues_outliers_removed.csv")


###################################################
### Import fasta
###################################################

### Import fasta

setwd("/1_Aligns/5_singles_gene_trees")

read.dna("concatenated.phy") -> conc
read.csv("concatenated_map.csv", row.names=1) -> map

conc[-match(drop, rownames(conc)),] -> conc


###################################################
### Export new alignments
###################################################

setwd("/1_Aligns/6_singles_gene_trees_rogueless")

aligns <- vector("list", length = nrow(map))
names(aligns) <- rownames(map)
paste(rownames(map), ".fas", sep="") -> files

for (i in 1:length(aligns)) {
  conc[,map[i,1]:map[i,2]] -> a0
  delete.empty.cells(a0) -> a0
  write.dna(a0, file=files[i], "fasta")
  a0 -> aligns[[i]]
}


fastConc(aligns, fill.with.gaps = T, map=T) -> conc
conc$map -> map
conc$align -> conc


write.phy(conc, "concatenated.phy")
write.csv(map, "concatenated_map.csv")


###################################################
### Stats
###################################################

aligns$concatenated <- conc

lapply(aligns, dstats, missing.char = "n") -> stats.out
do.call(rbind, stats.out) -> stats.out
unlist(lapply(aligns, nrow)) -> n
stats.out[,c(1:5)] -> stats.out
data.frame(Terminals=n, stats.out) -> stats.out
colnames(stats.out)[2] <- "Aligned bp"

write.csv(stats.out, "DNA_stats_rogueless.csv")

pratchet(phyDat(conc)) -> tree
write.tree(tree, "concatenated_MP.tre")
