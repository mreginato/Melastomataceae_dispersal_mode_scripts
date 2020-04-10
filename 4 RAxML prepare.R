library(ape)
library(phyloch)

setwd("/4_RAxML")
getwd() -> wd

##################################################
### Import concatenated and map
##################################################

setwd("/1_Aligns/5_singles_gene_trees")

read.dna("concatenated.phy") -> conc
read.csv("concatenated_map.csv", row.names=1) -> map


##################################################
### Best scheme from partitionFinder
##################################################

s1 <- c("accD", "psbK", "trnS")
s2 <- c("atpB-rbcL", "ndhF", "trnL")
s3 <- c("atpF", "rpl16")
s4 <- c("psbA")
s5 <- c("rbcL")
s6 <- c("ETS")
s7 <- c("ITS")
s8 <- c("waxy")

part <- list(s1,s2,s3,s4,s5,s6,s7,s8)


##################################################
### Generate new concatenated 
##################################################

setwd(wd)

pos <- vector(length=length(part))
conc[,1] -> concatenated

for (i in 1:length(part)) {
  match(part[[i]], rownames(map)) -> c1
  for (k in 1:length(c1)) {
    c1[k] -> c0
    conc[,map[c0,1]:map[c0,2]] -> a0
    cbind(concatenated, a0) -> concatenated
  }
}
concatenated[,-1] -> concatenated

ncol(conc) == ncol(concatenated)

write.phy(concatenated, "concatenated.phy")


##################################################
### Generate part file
##################################################

pos <- matrix(nrow=length(part), ncol=2)
e0 = 0

map$end - (map$start-1) -> map$length

for (i in 1:length(part)) {
  match(part[[i]], rownames(map)) -> c1
  sum(map$length[c1]) -> l0
  e0++1 -> start
  e0++l0 -> end
  end -> e0
  start -> pos[i,1]
  end -> pos[i,2]
}

paste("DNA, gene", 1:length(part), "=", pos[,1], "-", pos[,2], sep="") -> part
cat(part, sep="\n", file="concatenated.part")
