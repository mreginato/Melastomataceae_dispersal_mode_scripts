library(ape)
library(MonoPhy)
library(phyloch)
library(phyloU)
library(DECIPHER)
library(phangorn)
source("fastConc.R")

setwd("/3_Gene_trees")
getwd() -> wd

###################################################
### Generate phylips
###################################################

### Import fasta

setwd("/1_Aligns/4_aliscored")
list.files(pattern=".fas") -> files
files[-grep(".txt", files)] -> files
sub(".fas", "", files) -> labs
sapply(files, read.dna, "fasta", simplify = F) -> dna
setwd(wd)

### Export for RAxML

for (i in 1:length(dna)) {
  write.phy(dna[[i]], sub(".fas", ".phy", files[i]))
}


###################################################
### Run RAxML
###################################################

### Ran on cipres

raxml.path = "raxmlHPC-PTHREADS-SSE3_8210.exe"
cores = 6

### ML only
### raxmlHPC-PTHREADS-SSE3_8210.exe -T 6 -f d -m GTRGAMMA -N 1 -O -p 462 
### -s F:/Melas_package/3_Phylogeny/2018-02/4_Gene trees/accD.phy -n accD.tre 
### -w F:/Melas_package/3_Phylogeny/2018-02/4_Gene trees/

### ML + rapid bootstrap
### raxmlHPC-PTHREADS-SSE3_8210.exe -T 2 -f a -x 139 -m GTRGAMMA -p 226 -N 100 
### -s F:/Melas_package/3_Phylogeny/2018-02/4_Gene_trees/atpB-rbcL.phy -n atpB-rbcL.tre 
### -O -w F:/Melas_package/3_Phylogeny/2018-02/4_Gene_trees/
  

for (i in 1:length(dna)) {
  call <- paste(raxml.path, " -T ", cores, " -f a -x 666 -m GTRGAMMA -N 100 -O -p 462 -s ", wd, "/", labs[i], ".phy -n ", labs[i], ".tre -w ", wd, sep="")
  system(call)
}


###################################################
### Checkk monophyly
###################################################

### Tribes

setwd("/1_Aligns/3_singles_filtered")

read.csv("1 Markers_table_final.csv") -> meta
data.frame(Terminal=meta$spp, Clade=meta$Clade) -> clades

setwd(wd)

###################################################
### Gene Trees
###################################################

list.files(pattern=".rooted.tre") -> files
sub("RAxML_bipartitions.", "", sub(".rooted.tre", "", files)) -> labs
sapply(files, read.tree, simplify = F) -> trees

### check labels

unlist(lapply(lapply(dna, rownames),length)) -> dna.l
unlist(lapply(trees, Ntip)) -> trees.l
dna.l == trees.l


###################################################
### Check monophyly
###################################################

rownames(clades) <- clades$Terminal

mono.a <- vector("list", length = length(labs))
names(mono.a) <- labs
intruders.l <- intruders <- outliers.l <- outliers <- mono.s <- mono.a

for (i in 1:length(labs)) {
  trees[[i]] -> t0
  unlist(lapply(strsplit(t0$tip.label, "-"), "[", 1)) -> t0$tip.label
  checkPhylo(t0, clades)$data -> d0
  mono <- AssessMonophyly(t0, d0)
  GetSummaryMonophyly(mono)$Clade -> mono.s[[i]]
  GetResultMonophyly(mono) -> mono.a[[i]]
  GetIntruderTips(mono) -> intruders[[i]]
  GetOutlierTips(mono) -> outliers[[i]]
  if (length(intruders[[i]]$Clade) > 0) {
    data.frame(marker=labs[i], type="intruder", spp=as.character(unlist(intruders[[i]]))) -> intruders.l[[i]]
  }
  if (length(outliers[[i]]$Clade) > 0) {
    data.frame(marker=labs[i], type="outlier", spp=as.character(unlist(outliers[[i]]))) -> outliers.l[[i]]
  }
}


do.call(cbind, mono.s) -> mono.s
do.call(rbind, unlist(mono.a, recursive=F)) -> mono.a

write.csv(mono.s, "monophyly_tab.csv")
write.csv(mono.a, "monophyly_summary.csv")


###################################################
### Export new alignments
###################################################

setwd("/1_Aligns/4_aliscored")
list.files(pattern=".fas") -> files
files[-grep(".txt", files)] -> files
sub(".fas", "", files) -> labs
sapply(files, read.dna, "fasta", simplify=F) -> dna
new.dna <- vector("list", length = length(files))
names(new.dna) <- labs

labs == names(outliers.l)

setwd("/1_Aligns/5_singles_gene_trees")

for (i in 1:length(files)) {
  dna[[i]] -> d0
  unlist(lapply(strsplit(rownames(d0), "-"), "[", 1)) -> rownames(d0)
  outliers.l[[i]] -> o0
  if (is.null(o0) == F) {
    as.character(o0$spp) -> o0
  }
  o0 -> r0
  if (is.null(r0) == F) {
    d0[-match(r0, rownames(d0)),] -> d0
  }
  d0 -> new.dna[[i]]
}


do.call(rbind, intruders.l) -> intruders.l
do.call(rbind, outliers.l) -> outliers.l
outliers.l -> rem
write.csv(rem, "monophyly_seqs_removed.csv", row.names = F)


### Align and Export new alignments

for (i in 1:length(files)) {
  new.dna[[i]] -> a0
  delete.empty.cells(a0) -> a0
  del.gaps(a0) -> a0
  write.dna(a0, file="temp.fas", "fasta")
  readDNAStringSet("temp.fas", format="fasta") -> a0
  unlink("temp.fas")
  AlignSeqs(a0, iterations = 50, refinements = 50, processors=6) -> a0
  AdjustAlignment(a0, processors=3) -> a0
  #write.phy(a0, files[i]) 
  writeXStringSet(a0, filepath = files[i], format="fasta")
  read.dna(files[i], "fasta") -> a0
  delete.empty.cells(a0) -> a0
  fillEndsWithN(a0) -> a0
  write.dna(a0, file=files[i], "fasta")
  a0 -> new.dna[[i]]
}


fastConc(new.dna, fill.with.gaps = T, map=T) -> conc
conc$map -> map
conc$align -> conc


write.phy(conc, "concatenated.phy")
write.csv(map, "concatenated_map.csv")
