library(ape)
library(phyloch)
library(phyloU)
library(MonoPhy)

setwd("/11_treePL_secondary")
getwd() -> wd


###################################################
### Trees
###################################################

tree.file = "RAxML_bipartitions.concatenated.rooted.blfixed.tre"
boot.file = "RAxML_bootstrap.concatenated"

read.tree(tree.file) -> tree
read.tree(boot.file) -> boot

.compressTipLabel(boot) -> boot
attr(boot, "TipLabel") -> t
unlist(lapply(strsplit(t, "_"), FUN=function(x)(paste(x[2],x[3], sep="_")))) -> attr(boot, "TipLabel")

###################################################
### Tribes
###################################################

setwd("/1_Aligns/3_singles_filtered")

read.csv("1 Markers_table_final.csv", stringsAsFactors = F) -> meta
data.frame(Terminal=meta$spp, Clade=meta$Clade, stringsAsFactors = F) -> clades
rownames(clades) <- clades$Terminal

checkPhylo(tree, clades) -> c0
c0$data -> clades

setwd(wd)


###################################################
### Calibrations
###################################################

read.csv("secondary_calibrations.csv") -> dat.c

dat.c$Clades -> clades.c
calibrations <- vector("list", length = nrow(dat.c))
names(calibrations) <- toupper(clades.c)
calibrations.spp <- calibrations

for (i in 1:length(calibrations)) {
  as.numeric(dat.c[i,3:4]) -> calibrations[[i]]
  clades$Terminal[which(clades$Clade == clades.c[i])] -> spp0
  paste(spp0, collapse = " ") -> calibrations.spp[[i]]
}

### Add root and melas

melastomataceae <- c(56, 75.3)
all <- c(82.4, 98.8)

c2 <- list(melastomataceae, all)
names(c2) <- c("MELASTOMATACEAE", "ALL")

clades$Terminal[-which(clades$Clade == "CAP")] -> melastomataceae
clades$Terminal -> all
paste(melastomataceae, collapse = " ") -> melastomataceae
paste(all, collapse = " ") -> all

c2.spp <- list(melastomataceae, all)
names(c2.spp) <- c("MELASTOMATACEAE", "ALL")

c(calibrations, c2) -> calibrations
c(calibrations.spp, c2.spp) -> calibrations.spp


###################################################
### Root
###################################################

clades$Terminal[which(clades$Clade == "CAP")] -> outgroup

boot.r <- vector("list", length = length(boot))

for (i in 1:length(boot)) {
  boot[[i]] -> b0
  getMRCA(b0, outgroup) -> mrca
  root(b0, node=mrca, resolve.root = T) -> boot.r[[i]]
}


###################################################
### treePL configuration file
###################################################
# treefile = intree.tre
# smooth = 100
# numsites = 8502
# mrca = EUDICOT Papaver_nudicaule Lonicera_etrusca
# min = EUDICOT 123
# max = EUDICOT 140
# outfile = intree.dated.tre
# thorough
# prime
# randomcv 
# nthreads = 6
###################################################

options(scipen=999)

tree.file = "RAxML_bipartitions.concatenated.rooted.blfixed.tre"
tree$node.label <- NULL
write.tree(tree, tree.file)

.compressTipLabel(boot.r) -> boot.r

###################################################
### Estimate prime
###################################################

numsites = 14439
outfile = paste(sub(".tre", "", tree.file), "prime_treePL.tre", sep="_")
nthreads = 8

file = "treepl_prime_config.txt"
cat("treefile = ", tree.file, "\nnumsites = ", numsites,
    "\noutfile = ", outfile, "\nthorough \nprime \nnthreads = ", nthreads, sep="", fill=F, file=file)
for (i in 1:length(calibrations)) {
  cat("\nmrca =", names(calibrations)[i], calibrations.spp[[i]],  sep=" ", file=file, fill=F, append=T)
  cat("\nmin =", names(calibrations)[i], calibrations[[i]][1], sep=" ", file=file, append=T)
  cat("\nmax =", names(calibrations)[i], calibrations[[i]][2], sep=" ", file=file, append=T)
}
cat(" ", fill=T, file=file, append=T)


###################################################
### Estimate smooth parameter
###################################################

### From prime
### opt = 1
### optad = 1
### moredetailad
### optcvad = 1
### moredetailcvad


numsites = 14439
outfile = paste(sub(".tre", "", tree.file), "randomcv_treePL.tre", sep="_")
nthreads = 8
opt = 1
optad = 1
optcvad = 1
cvstart = 0.01
cvstop = 10
cvmultstep = 0.1

file = "treepl_randomcv_config.txt"
cat("treefile = ", tree.file, "\nnumsites = ", numsites, "\nopt = ", opt, "\noptad = ", optad, "\noptcvad = ", optcvad,
    "\noutfile = ", outfile, "\nmoredetailad \nmoredetailcvad \nthorough \nrandomcv \nnthreads = ", nthreads, "\ncvstart = ", cvstart,
    "\ncvstop = ", cvstop, "\ncvmultstep = ", cvmultstep, sep="", fill=F, file=file)
for (i in 1:length(calibrations)) {
  cat("\nmrca =", names(calibrations)[i], calibrations.spp[[i]],  sep=" ", file=file, fill=F, append=T)
  cat("\nmin =", names(calibrations)[i], calibrations[[i]][1], sep=" ", file=file, append=T)
  cat("\nmax =", names(calibrations)[i], calibrations[[i]][2], sep=" ", file=file, append=T)
}
cat(" ", fill=T, file=file, append=T)


###################################################
### treePL best tree
###################################################


### From randomcv
### smooth = 1 # max chisq

smooth = 1000
numsites = 14439
outfile = paste(sub(".tre", "", tree.file), "treePL.tre", sep="_")
nthreads = 6

file = "treepl_besttree_config.txt"
cat("treefile = ", tree.file, "\nsmooth = ", smooth, "\nnumsites = ", numsites,
    "\noutfile = ", outfile, "\nopt = ", opt, "\noptad = ", optad, "\noptcvad = ", optcvad,
    "\nmoredetailad \nmoredetailcvad \nthorough \nnthreads = ", nthreads, sep="", fill=F, file=file)
for (i in 1:length(calibrations)) {
  cat("\nmrca =", names(calibrations)[i], calibrations.spp[[i]],  sep=" ", file=file, fill=F, append=T)
  cat("\nmin =", names(calibrations)[i], calibrations[[i]][1], sep=" ", file=file, append=T)
  cat("\nmax =", names(calibrations)[i], calibrations[[i]][2], sep=" ", file=file, append=T)
}
cat(" ", fill=T, file=file, append=T)


###################################################
### Bootstraps
###################################################

setwd("F:/Melas_package/3_Phylogeny/2018-02/11_treePL_secondary/bootstrap")

for (i in 1:length(boot.r)) {
  tree.file = paste("boot_", i, ".tre", sep="")
  write.tree(boot.r[[i]], tree.file)
  outfile = paste(sub(".tre", "", tree.file), "treePL.tre", sep="_")
  file = paste("treepl_boot_", i, "_config.txt", sep="")
  cat("treefile = ", tree.file, "\nsmooth = ", smooth, "\nnumsites = ", numsites,
      "\noutfile = ", outfile, "\nopt = ", opt, "\noptad = ", optad, "\noptcvad = ", optcvad,
      "\nmoredetailad \nmoredetailcvad \nthorough \nnthreads = ", nthreads, sep="", fill=F, file=file)
  for (i in 1:length(calibrations)) {
    cat("\nmrca =", names(calibrations)[i], calibrations.spp[[i]],  sep=" ", file=file, fill=F, append=T)
    cat("\nmin =", names(calibrations)[i], calibrations[[i]][1], sep=" ", file=file, append=T)
    cat("\nmax =", names(calibrations)[i], calibrations[[i]][2], sep=" ", file=file, append=T)
  }
  cat(" ", fill=T, file=file, append=T)
}

### Master

list.files(pattern="config.txt") -> files
file = "run_treepl_boot.sh"
cat("#!/bin/bash", file=file, fill=T)
cat(paste("treePL", files), sep="\n", append=T, file=file)


###################################################
### Merge bootstraps
###################################################

setwd("F:/Melas_package/3_Phylogeny/2018-02/7_treePL/bootstrap")

list.files(pattern = "_treePL.tre") -> files
files[-grep(".r8s", files)] -> files
sapply(files, read.tree, simplify = F) -> boots
.compressTipLabel(boots) -> boots

setwd(wd)

write.nexus(boots, file = paste(boot.file, ".rooted_treePL.trees", sep=""))

outfile = paste(sub(".tre", "", tree.file), "treePL.tre", sep="_")

read.tree(outfile) -> out.tree
ladderize(out.tree, F) -> out.tree
write.nexus(out.tree, file="Target_tree.tre")
