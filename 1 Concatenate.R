library(ape)
library(phyloch)
source("dstats.R")
source("fastConc.R")

setwd("/1_Aligns/1_raw")
read.csv("Seqs_metadata.csv", stringsAsFactors = F) -> meta
unique(meta[,c(2,3,7)]) -> clades

in.dir <- "/1_Aligns/2_singles_aligned"
out.dir <- "/1_Aligns/3_singles_filtered"

###################################################
### Read singles
###################################################

setwd(in.dir)
list.files(pattern=".fas") -> files
sub(".fas", "", files) -> labels
sapply(files, FUN=function(x)(read.dna(x, "fasta"))) -> aligns
names(aligns) <- labels
setwd(out.dir)

###################################################
### Smart concatenation
###################################################
### Merge samples when possible; exclude duplicates 
###################################################

### Initial Markers Table

data.frame(table(unlist(lapply(aligns, rownames)))) -> all.labels
seqs.tab <- data.frame(spp=unlist(lapply(strsplit(as.character(all.labels[,1]), "-"), "[", 1)),
                       terminal=as.character(all.labels[,1]), markers=all.labels[,2])
seqs.tab[order(seqs.tab[,3], decreasing=T),] -> seqs.tab
as.character(unique(seqs.tab$spp)) -> spp
as.character(seqs.tab$terminal) -> terminals
markers.tab <- matrix(nrow=nrow(seqs.tab), ncol=length(aligns))
rownames(markers.tab) <- seqs.tab$terminal
colnames(markers.tab) <- names(aligns)
markers.tab[] <- 0
for (i in 1:nrow(markers.tab)) {
  terminals[i] -> t0
  unlist(lapply(aligns, FUN=function(x)(match(t0, rownames(x))))) -> x
  x[x > 0] <- 1
  which(is.na(x)) -> miss0
  if (length(miss0) > 0) {
    x[miss0] <- 0
  }
  x -> markers.tab[i,]
}
cbind(seqs.tab, markers.tab) -> markers.tab

clades$Clade[match(markers.tab$terminal, clades$Terminal)] -> markers.tab$Clade

write.csv(markers.tab, "0 Markers_table_initial.csv", row.names=F)

### Merge then exclude duplicates

markers.tab$merged <- ""
markers.tab$excluded <- ""

for (i in 1:length(spp)) {
  spp[i] -> sp0
  which(markers.tab$spp == sp0) -> x
  markers.tab[x,] -> tab0
  tab0[,-c(1:3)] -> tab1
  if (nrow(tab0) > 1) {
    merged <- vector()
    for (k in 2:(nrow(tab0))) {
      which(tab1[1,] == 0) -> miss0
      if (length(miss0) > 0) {
        miss0[which(tab1[k,miss0] == 1)] -> m0
        if (length(m0) > 0) {
          tab1[1,m0] <- 1
          c(merged,as.character(tab0[k,2])) -> merged
          for (j in 1:length(m0)) {
            aligns[[m0[j]]] -> a0
            rownames(a0)[match(as.character(tab0[k,2]), rownames(a0))] <- as.character(tab0[1,2])
            a0 -> aligns[[m0[j]]]
          }
        }
      }
    }
    if (length(merged) > 0) {
      paste(merged, collapse=", ") -> markers.tab$merged[x[1]]
    }
    paste(as.character(markers.tab$terminal[x[2:length(x)]]), collapse=", ") -> markers.tab$excluded[x[1]]
    markers.tab[-x[2:length(x)],] -> markers.tab
    as.character(tab0[2:nrow(tab0),2]) -> duplis
    for (w in 1:length(aligns)) {
      aligns[[w]] -> a0
      which(is.na(match(rownames(a0), duplis)) == F) -> rem
      if (length(rem) > 0) {
        a0[-rem,] -> a0
      }
      a0 -> aligns[[w]]
    }
  }
}

nrow(markers.tab) == length(unique(markers.tab$spp))

### Final markers table

markers.tab[order(markers.tab$spp),] -> markers.tab

markers.tab$merged -> merged
markers.tab$excluded -> excluded
markers.tab$spp -> spp.temp

data.frame(table(unlist(lapply(aligns, rownames)))) -> all.labels
seqs.tab <- data.frame(spp=unlist(lapply(strsplit(as.character(all.labels[,1]), "-"), "[", 1)),
                       terminal=as.character(all.labels[,1]), markers=all.labels[,2])
seqs.tab[order(seqs.tab$spp),] -> seqs.tab
as.character(unique(seqs.tab$spp)) -> spp
spp == spp.temp
as.character(seqs.tab$terminal) -> terminals
markers.tab <- matrix(nrow=nrow(seqs.tab), ncol=length(aligns))
rownames(markers.tab) <- seqs.tab$terminal
colnames(markers.tab) <- names(aligns)
markers.tab[] <- 0
for (i in 1:nrow(markers.tab)) {
  terminals[i] -> t0
  unlist(lapply(aligns, FUN=function(x)(match(t0, rownames(x))))) -> x
  x[x > 0] <- 1
  which(is.na(x)) -> miss0
  if (length(miss0) > 0) {
    x[miss0] <- 0
  }
  x -> markers.tab[i,]
}

clades$Clade[match(seqs.tab$terminal, clades$Terminal)] -> Clade
cbind(Clade,seqs.tab) -> seqs.tab

cbind(seqs.tab, markers.tab, merged, excluded) -> markers.tab

write.csv(markers.tab, "1 Markers_table_final.csv", row.names=F)

### Genbank table

genbank.tab <- markers.tab[,-c(3,18,19)]
genbank.tab[,4:ncol(genbank.tab)] <- NA

for (i in 1:nrow(genbank.tab)) {
  as.character(markers.tab$terminal[i]) -> p0
  as.character(markers.tab$merged[i]) -> p1
  unlist(strsplit(p1, ", ")) -> p1
  meta[which(meta$Terminal == p0),] -> m0
  genbank.tab[i,match(m0$Marker, colnames(genbank.tab))] <- m0$Accession
  if (is.na(match("", p1))) {
    for (k in 1:length(p1)) {
      colnames(genbank.tab)[which(is.na(genbank.tab[i,]))] -> keep
      meta[which(meta$Terminal == p1[k]),] -> m0
      m0[which(is.na(match(m0$Marker, keep)) == F),] -> m0
      if (nrow(m0) > 0) {
        genbank.tab[i,match(m0$Marker, colnames(genbank.tab))] <- m0$Accession
      }
    }
  }
}

sub("_", " ", genbank.tab$spp) -> genbank.tab$spp
as.matrix(genbank.tab) -> genbank.tab
genbank.tab[is.na(genbank.tab[])] <- "-"
data.frame(genbank.tab) -> genbank.tab
write.csv(genbank.tab, "1 Markers_table_final_genbank.csv", row.names=F)


### Export new alignments

paste(labels, ".fas", sep="") -> files
for (i in 1:length(aligns)) {
  aligns[[i]] -> a0
  delete.empty.cells(a0) -> a0
  del.gaps(a0) -> a0
  write.dna(a0, file="temp.fas", "fasta")
  readDNAStringSet("temp.fas", format="fasta") -> a0
  unlink("temp.fas")
  AlignSeqs(a0, iterations = 50, refinements = 50, processors=6) -> a0
  AdjustAlignment(a0, processors=6) -> a0
  writeXStringSet(a0, filepath = files[i], format="fasta")
  read.dna(files[i], "fasta") -> a0
  fillEndsWithN(a0) -> a0
  delete.empty.cells(a0) -> a0
  unlist(lapply(strsplit(rownames(a0), "-"), "[", 1)) -> rownames(a0)
  write.dna(a0, file=files[i], "fasta")
  write.phy(a0, file=paste(labels[i], ".phy", sep=""))
  a0 -> aligns[[i]]
}


fastConc(aligns, fill.with.gaps = T, map = T) -> conc
conc$map -> map
conc$align -> conc

rownames(conc) -> spp
spp[which(duplicated(spp))]

write.phy(conc, "concatenated.phy")
write.csv(map, "concatenated_map.csv")

jpeg("conc.jpg")
image(conc)
dev.off()

### Stats

aligns$concatenated <- conc

lapply(aligns, dstats, missing.char = "n") -> stats.out
do.call(rbind, stats.out) -> stats.out
unlist(lapply(aligns, nrow)) -> n
stats.out[,c(1:5)] -> stats.out
data.frame(Terminals=n, stats.out) -> stats.out
colnames(stats.out)[2] <- "Aligned bp"

write.csv(stats.out, "DNA_stats.csv")

