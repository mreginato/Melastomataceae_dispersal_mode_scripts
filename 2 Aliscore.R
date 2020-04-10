library(ape)
library(phyloch)
library(phangorn)
source("dstats.R")
source("fastConc.R")

setwd("/1_Aligns/4_aliscored")
getwd() -> wd


###################################################
### Generate bash script
###################################################

aligns.dir = "/1_Aligns/3_singles_filtered"
setwd(aligns.dir)
list.files(pattern = '.fas') -> files
sub(".fas", ".tre", files) -> files.t
sapply(files, read.dna, "fasta") -> aligns
setwd(wd)

file = "run_aliscore.sh"
cat("#!/bin/bash", fill=T, file=file)
sub("F:", "/mnt/f", aligns.dir) -> aligns.dir
sub("F:", "/mnt/f", wd) -> trees.dir
cat(paste("cp ", aligns.dir, "/", files, " ~/aliscore/", files, sep=""), sep="\n", file=file, append=T)

### Aliscore

for (i in 1:length(files)) {
  cat("\nperl Aliscore.02.2.pl -i ", files[i], " -r", sep="", file=file, append=T)
}
cat("\n", file=file, append=T)


###################################################
### Generate guide trees
###################################################

for (i in 1:length(aligns)) {
  pratchet(phyDat(aligns[[i]])) -> t0
  write.tree(t0, files.t[i])
}

###################################################
### Import aliscore output
###################################################

aligns -> old.aligns
list.files(pattern="List_random.txt") -> txt.files

for (i in 1:length(files)) {
  aligns[[i]] -> a0
  id <- scan(txt.files[i], what = "char", quiet = TRUE)
  id <- as.numeric(id)
  if (length(id) > 0) {
    a0[,-id] -> a0
    delete.empty.cells(a0) -> a0
  }
  write.dna(a0, files[i], "fasta")
  aligns[[i]] <- a0
}

unlist(lapply(aligns, ncol))
unlist(lapply(old.aligns, ncol))
unlist(lapply(aligns, ncol))/unlist(lapply(old.aligns, ncol)) -> removed.chars.n
write.csv(removed.chars.n, "Kept_chars_percent.csv")

fastConc(aligns, fill.with.gaps = T, map = T) -> conc
conc$map -> map
conc$align -> conc

jpeg("conc.jpg")
image(conc)
dev.off()

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

write.csv(stats.out, "DNA_stats_aliscored.csv")

