checkNamesMelas <- function(species, author=NULL, database=NULL, max.distance=list(insertions=0.01,deletions=0.01,substitutions=0.01)) {
  if (is.null(database)) {
    #melnet <- data(xxx)
  } else {
    database -> melnet
  }
  if (is.null(author) == F) {
    paste(species, author) -> old.spp.full
    check.full = T
  } else {
    old.spp.full <- vector(length=length(species))
    old.spp.full[] <- ""
    check.full = F
  }
  species -> old.spp
  t.check <- as.data.frame(matrix(ncol=5, nrow=length(old.spp)))
  t.check[,1] <- old.spp
  t.check[,2] <- old.spp.full
  t.check[,5] <- NA
  colnames(t.check) <- c("Old_spp", "Old_spp_full", "New_spp", "New_spp_full", "Confidence")
  
  ### Match 1 - full spp.
  if (check.full) {
    match(old.spp.full, melnet$Full) -> m1
    t.check[,4] <- as.character(melnet$Accepted[m1])
    t.check[which(is.na(m1) == F),5] <- 1
  }
  
  ### Match 2 - short spp.
  which(is.na(t.check$New_spp_full)) -> m0
  if (length(m0) > 0) {
    t.check[m0,] -> t.missing
    t.check[-m0,] -> t.check
    t.missing$Old_spp -> old.spp
    match(old.spp, melnet$Species) -> m2
    t.missing[,4] <- as.character(melnet$Accepted[m2])
    t.missing[which(is.na(m2)==F),5] <- 0.9
    rbind(t.check,t.missing) -> t.check
  }

  ### Match 3 - agrep full spp.
  
  which(is.na(t.check$Confidence)) -> m0
  if (length(m0) > 0 && check.full == T) {
    t.check[m0,] -> t.missing
    t.check[-m0,] -> t.check
    t.missing$Old_spp_full -> old.spp.full
    agrep0 <- cbind(old.spp.full, NA)
    for (i in 1:nrow(agrep0)) {
      agrep(agrep0[i,1], melnet$Full, max.distance = max.distance) -> g0
      if (length(g0) == 1) {
        as.character(melnet$Accepted[g0]) -> agrep0[i,2]
      }
      cat("\r",i)
    }
    t.missing[,4] <- agrep0[,2]
    t.missing[which(is.na(agrep0[,2])==F),5] <- 0.8
    rbind(t.check,t.missing) -> t.check
  }
  
  ### Match 4 - agrep short spp.
  which(is.na(t.check$Confidence)) -> m0
  if (length(m0) > 0) {
    t.check[m0,] -> t.missing
    t.check[-m0,] -> t.check
    t.missing$Old_spp -> old.spp
    agrep0 <- cbind(old.spp, NA)
    for (i in 1:nrow(agrep0)) {
      agrep(agrep0[i,1], melnet$Species, max.distance = max.distance) -> g0
      if (length(g0) == 1) {
        as.character(melnet$Accepted[g0]) -> agrep0[i,2]
      }
      cat("\r",i)
    }
    t.missing[,4] <- agrep0[,2]
    t.missing[which(is.na(agrep0[,2])==F),5] <- 0.7
    rbind(t.check,t.missing) -> t.check
    tail(t.check)
  }
  spp.less <- unlist(lapply(strsplit(t.check$New_spp_full, " "), FUN=function(x)(paste(x[1],x[2]))))
  spp.less[which(spp.less == "NA NA")] <- NA
  t.check$New_spp <- spp.less
  t.check[order(t.check$Old_spp),] -> t.check
  return(t.check)
}