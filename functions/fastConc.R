fastConc <- function (x, check.names = TRUE, fill.with.gaps = FALSE, quiet = FALSE, map=F) {
  obj <- x
  n <- length(obj)
  if (n == 1) 
    return(obj[[1]])
  for (i in 1:n) if (!is.matrix(obj[[1]])) 
    stop("the 'cbind' method for \"DNAbin\" accepts only matrices")
  NR <- unlist(lapply(obj, nrow))
  for (i in 1:n) class(obj[[i]]) <- NULL
  if (check.names) {
    nms <- unlist(lapply(obj, rownames))
    if (fill.with.gaps) {
      NC <- unlist(lapply(obj, ncol))
      nms <- unique(nms)
      ans <- matrix(as.raw(240), length(nms), sum(NC))
      rownames(ans) <- nms
      from <- 1
      for (i in 1:n) {
        to <- from + NC[i] - 1
        tmp <- rownames(obj[[i]])
        nmsi <- tmp[tmp %in% nms]
        ans[nmsi, from:to] <- obj[[i]][nmsi, , drop = FALSE]
        from <- to + 1
      }
    }
    else {
      tab <- table(nms)
      ubi <- tab == n
      nms <- names(tab)[which(ubi)]
      ans <- obj[[1]][nms, , drop = FALSE]
      for (i in 2:n) ans <- cbind(ans, obj[[i]][nms, , 
        drop = FALSE])
      if (!quiet && !all(ubi)) 
        warning("some rows were dropped.")
    }
  }
  else {
    if (length(unique(NR)) > 1) 
      stop("matrices do not have the same number of rows.")
    ans <- matrix(unlist(obj), NR)
    rownames(ans) <- rownames(obj[[1]])
  }
  class(ans) <- "DNAbin"
  if (map) {
    unlist(lapply(x, ncol)) -> ncols
    matrix(ncol=2, nrow=length(ncols)) -> map
    rownames(map) <- names(x)
    colnames(map) <- c("start", "end")
    map[1,1] <- 1
    map[1,2] <- ncols[1]
    for (i in 2:length(ncols)) {
      map[i,1] <- map[i-1,2]+1
      map[i,2] <- map[i-1,2]+ncols[i]  
    }
    list(align=ans, map=map) -> ans
  }
  return(ans)
}

