trim_align_edges <- function(dna, min.missing=0.5, quiet=F, edges.only=T) {
  dna -> a0
  tapply(a0, INDEX=col(a0), base.freq, freq=T, all=T) -> freq0
  lapply(freq0, "[", c(15:17)) -> freq0
  unlist(lapply(freq0, sum))/nrow(a0) -> freq0
  if (min.missing == 1) {
    freq0 == min.missing -> x
  } else {
    freq0 > min.missing -> x
  }
  if (edges.only) {
    trim.r <- trim.l <- vector()
    if (length(which(x == F)) > 0) {
      if (x[1]) {
        c(1:which(x == F)[1]) -> trim.l
      }
      if (x[length(x)]) {
        which(x == F) -> y
        c(y[length(y)]:length(x)) -> trim.r
      }
      c(trim.l, trim.r) -> trim
    } else {
      c(1:length(x)) -> trim
    }
  } else {
    which(x) -> trim
  }
  if (length(trim) > 0) {
    a0[,-trim] -> a0
    if (quiet == F) {
      cat(length(trim), "sites were trimmed")
    }
    if (ncol(a0) > 0) {
      tapply(a0, INDEX=row(a0), base.freq, freq=T, all=T) -> freq0
      lapply(freq0, "[", c(15:17)) -> freq0
      unlist(lapply(freq0, sum))/ncol(a0) -> freq0
      which(freq0 == 1) -> drop
    } else {
      c(1:nrow(a0)) -> drop
    }
    if (length(drop) > 0) {
      a0[-drop,] -> a0
      if (quiet == F) {
        cat("\n", length(drop), "sequence(s) were dropped")
      }
    }
  } else {
    if (quiet == F) {
      cat("No trimming was necessary")
    }
  }
  return(a0)
}
