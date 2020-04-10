library(ape)
library(DECIPHER)
source("checkNamesMelas.R")
source("trim_align_edges.R")

setwd("/1_Aligns/1_raw")
getwd() -> wd

###################################################
### Import sequences and metadata from Geneious 
###################################################

file = "All.csv"
keep.cols = c("Marker..Marker", "Organism", "Clade..Clade", "Sequence", "isolate", "specimen_voucher", "Accession")

read.csv(file, stringsAsFactors = F) -> all
subset(all, select=keep.cols) -> all
colnames(all) <- c("Marker", "Species", "Clade", "Sequence", "Isolate", "Voucher", "Accession")


###################################################
### Get rid of "cf.", "aff." and double spaces 
### between genus and epithet
###################################################

get.rid.indets = T

if (get.rid.indets) {
  all$Species -> spp
  spp[agrep("cf.", spp, max.distance = 0)]
  spp[agrep("aff.", spp, max.distance = 0)]
  spp[agrep("  ", spp, max.distance = 0)]
  sub(" cf.", "", spp, fixed=T) -> spp
  sub(" aff.", "", spp, fixed=T) -> spp
  spp -> all$Species
  
  ### Get rid of indets 
  
  all$Species -> spp
  spp[agrep("sp.", spp, max.distance=0)]
  all[-agrep("sp.", spp, max.distance=0),] -> all
} else {
  all$Species -> spp
  gsub("\\.", "", spp) -> all$Species
}


###################################################
### check Taxonomy
###################################################

read.csv("Melnet_data_2018_Feb_21.csv") -> melnet


all[order(all$Species),] -> all
all$Species -> spp

name.check <- checkNamesMelas(spp, author=NULL, database=melnet) 
unique(name.check[which(is.na(name.check$Confidence)),1])

c <- name.check[which(name.check$Confidence < 0.9),]
name.check[which(name.check$Confidence < 0.9),c(1,3)]

manual = T

if (manual) {
  c[,c(1,3)]
  ### dont change
  dc <- c("Tibouchina papyrus", "Miconia polychaete", "Miconia navifolia", "Leandra purpureovillosa")
  name.check[which(is.na(match(name.check$Old_spp, dc)) == F),1] -> name.check[which(is.na(match(name.check$Old_spp, dc)) == F),3]
}

### Missing

name.check[which(is.na(name.check$New_spp)),1] -> name.check[which(is.na(name.check$New_spp)),3]

tail(spp) == tail(name.check$Old_spp)
name.check$New_spp -> all$Species


###################################################
### Check consistency voucher/spp/isolate
###################################################

unique(all$Voucher) -> vouchers
vouchers.c <- data.frame(vouchers=vouchers, check.spp=NA, check.iso=NA)

for (i in 1:length(vouchers)) {
  all[which(all$Voucher == vouchers[i]),] -> a0
  length(unique(a0$Species)) -> l.spp
  length(unique(a0$Isolate)) -> l.iso
  if (l.spp == 1) {
    vouchers.c[i,2] <- "Ok"
  }
  if (l.iso == 1) {
    vouchers.c[i,3] <- "Ok"
  }
}

sort(vouchers[which(is.na(vouchers.c$check.spp))])
sort(vouchers[which(is.na(vouchers.c$check.iso))])

unique(all$Isolate) -> isolates
isolates.c <- data.frame(isolates=isolates, check.spp=NA, check.iso=NA)

for (i in 1:length(isolates)) {
  all[which(all$Isolate == isolates[i]),] -> a0
  length(unique(a0$Species)) -> l.spp
  length(unique(a0$Isolate)) -> l.iso
  if (l.spp == 1) {
    isolates.c[i,2] <- "Ok"
  }
  if (l.iso == 1) {
    isolates.c[i,3] <- "Ok"
  }
}

sort(isolates[which(is.na(isolates.c$check.spp))])
sort(isolates[which(is.na(isolates.c$check.iso))])


### check same T-number different species

all$Isolate -> t.number
all$Species -> spp
unique(t.number) -> Ts
check.t.spp <- vector()
for (i in 1:length(Ts)) {
  Ts[i] -> t0
  spp[which(t.number == t0)] -> spp0
  c(check.t.spp, length(unique(spp0))) -> check.t.spp
}

Ts[which(check.t.spp > 1)]

### check same T-number different voucher

all$Isolate -> t.number
all$Voucher -> voucher
unique(t.number) -> Ts
check.t.vou <- vector()
for (i in 1:length(Ts)) {
  Ts[i] -> t0
  voucher[which(t.number == t0)] -> spp0
  c(check.t.vou, length(unique(spp0))) -> check.t.vou
}

Ts[which(check.t.vou > 1)]


### check same voucher different species

unique(voucher) -> Vs
check.v.spp <- vector()
for (i in 1:length(Vs)) {
  Vs[i] -> v0
  spp[which(voucher == v0)] -> spp0
  c(check.v.spp, length(unique(spp0))) -> check.v.spp
}

Vs[which(check.v.spp > 1)]


###################################################
### Sort seqs by number of markers
###################################################

all$Voucher[which(is.na(all$Voucher))] <- all$Isolate[which(is.na(all$Voucher))]
gsub("\\.", "", all$Voucher) -> all$Voucher
gsub(" ", "", all$Voucher) -> all$Voucher
gsub("_", "", all$Voucher) -> all$Voucher
gsub("-", "", all$Voucher) -> all$Voucher

paste(gsub(" ", "_", all$Species), all$Voucher, sep="-") -> all$Terminal


###################################################
### Remove markers with less than x seqs
###################################################

x = 100

all[-which(duplicated(all[,c(1,2)])),] -> less

table(less$Marker) -> markers.t
sort(markers.t)
names(which(markers.t < x)) -> rem.markers
#unique(c("psbA", rem.markers)) -> rem.markers

all[-which(is.na(match(all$Marker, rem.markers)) == F),] -> all


###################################################
### Remove duplicated seqs by marker
###################################################

paste(all$Marker, all$Terminal) -> x
sort(x[which(duplicated(x))])
which(duplicated(x)) -> y

if (length(y) > 0) {
  all[-y,] -> all
}


###################################################
### Generate alignments
###################################################

unique(all$Marker) -> markers

for (i in 1:length(markers)) {
  markers[i] -> m0
  all[which(all$Marker == m0),] -> d0
  writeLines(paste(paste(">",d0$Terminal, sep=""), d0$Sequence, sep="\n"), sep="\n", con=paste(m0, ".fas", sep=""))
}


###################################################
### Export New metadata
###################################################

subset(all, select=c("Terminal", "Species", "Clade")) -> meta1
meta1[-which(duplicated(meta1)),] -> meta1

markers.table <- matrix(ncol=length(markers), nrow=nrow(meta1))
colnames(markers.table) <- markers
rownames(markers.table) <- meta1$Terminal
markers.table[] <- ""

for (i in 1:length(markers)) {
  markers[i] -> m0
  all[which(all$Marker == m0),] -> d0
  markers.table[match(d0$Terminal, rownames(markers.table)),i] <- d0$Accession
}

data.frame(meta1,markers.table) -> meta2
write.csv(meta2, "Markers_matrix.csv", row.names=F)

all[,-4] -> all
write.csv(all, "Seqs_metadata.csv", row.names=F)


###################################################
### Align
###################################################


setwd(wd)

in.dir <- "/1_Aligns/1_raw"
out.dir <- "/1_Aligns/2_singles_aligned"

list.files(pattern=".fas") -> files

method = "decipher"

if (method == "muscle") {
  maxit = 2 # 16 default
  for (i in 1:length(files)) {
    call <- paste("muscle -maxiters ", maxit, " -in ", in.dir, "/", files[i], " -out ", out.dir, "/", files[i], sep="")
    system(call)
  }
}

if (method == "decipher") {
  for (i in 1:length(files)) {
    readDNAStringSet(files[i], format="fasta") -> a0
    OrientNucleotides(a0) -> a0
    AlignSeqs(a0) -> a0
    writeXStringSet(a0, filepath = paste(out.dir, files[i], sep="/") , format="fasta")
  }
}

trim = T

if (trim) {
  setwd(out.dir)
  for (i in 1:length(files)) {
    read.dna(files[i], "fasta") -> a0
    trim_align_edges(a0, min.missing = 0.75) -> a0
    write.dna(a0, file=files[i], "fasta")
  }
  setwd(wd)
}

