dstats <- function(x, include.amb=T, missing.char=c("n", "-", "?")) {
  x -> a0
  ncol(a0) -> total0
  if (include.amb) {
    tapply(a0, INDEX=col(a0), base.freq, freq=T, all=T) -> freq0
    unlist(lapply(freq0, "[", "-")) -> gap0
    length(which(gap0 > 0)) -> gap.n
    unlist(lapply(freq0, FUN=function(x)(sum(x[c(5:14)])))) -> het0
    length(which(het0 > 0)) -> het.n
    lapply(freq0, "[", -c(15:17)) -> freq0
  } else {
    tapply(a0, INDEX=col(a0), base.freq, freq=T, all=F) -> freq0
    unlist(lapply(freq0, "[", "-")) -> gap0
    length(which(gap0 > 0)) -> gap.n
    het.n <- NA
  }
  lapply(freq0, FUN=function(x)(length(which(x > 0)))) -> v0
  lapply(freq0, FUN=function(x)(length(which(x > 1)))) -> p0
  length(which(v0 > 1)) -> var0
  length(which(p0 > 1)) -> pis0
  freqs <- base.freq(x, all = T)
  missing.data <- sum(freqs[missing.char] * 100)
  missing.data <- round(missing.data, 1)
  c(Total=total0, Variable=var0, PIS=pis0, Missing.data=missing.data, Sites_with_Gaps=gap.n, Sites_with_Heterozygotes=het.n) -> out
  return(out)
}


