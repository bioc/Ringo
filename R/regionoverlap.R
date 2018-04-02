# Function to compute overlaps between (enriched) genomic regions in one list and those in another list
regionOverlap <- function(xdf, ydf, chrColumn="chr",startColumn="start",endColumn="end", mem.limit=1e8) {
  stopifnot(is.data.frame(xdf), is.data.frame(ydf),
            all(c(chrColumn, startColumn, endColumn) %in% names(xdf)),
            all(c(chrColumn, startColumn, endColumn) %in% names(ydf)))
  if (!all(xdf[[startColumn]]<=xdf[[endColumn]]))
    stop("Some regions in",deparse(substitute(xdf)),"have end postions that are smaller than their start positions.\n")
  if (!all(ydf[[startColumn]]<=ydf[[endColumn]]))
    stop(paste("Some regions in",deparse(substitute(ydf)),"have end postions that are smaller than their start positions.\n"))
  ### avoid creation of too big matrices here:
  if (log2(nrow(xdf))+log2(nrow(ydf))<log2(mem.limit)){
    res <- .Call(ringoRegionOverlap, as.character(xdf[[chrColumn]]), as.integer(as.character(xdf[[startColumn]])), as.integer(as.character(xdf[[endColumn]])), as.character(ydf[[chrColumn]]), as.integer(as.character(ydf[[startColumn]])), as.integer(as.character(ydf[[endColumn]])))
    ## convert the result to a dgCMatrix (package matrix)
    res <- as(res, "dgCMatrix")
  } else {
    ## iterative procedure:
    if (log2(nrow(xdf))> log2(mem.limit))
      stop("Matrix ", deparse(substitute(xdf)), " has too many rows. Work with subsets of this matrix or increase parameter 'mem.limit'.\n") ## maybe implement a workaround that as well later
    maxYRows <- floor(2^(log2(mem.limit)-log2(nrow(xdf))))
    yRowSeq  <- unique(as.integer(c(seq(maxYRows, nrow(ydf), by=maxYRows), nrow(ydf))))
    partRess <- vector("list", length(yRowSeq))
    firstYCol <- as.integer(1)
    for (j in seq(length(yRowSeq))){
      lastYCol <-  yRowSeq[j]
      partYdf <- ydf[firstYCol:lastYCol,,drop=FALSE]
      thisRes <- .Call(ringoRegionOverlap, as.character(xdf[[chrColumn]]), as.integer(as.character(xdf[[startColumn]])), as.integer(as.character(xdf[[endColumn]])), as.character(partYdf[[chrColumn]]), as.integer(as.character(partYdf[[startColumn]])), as.integer(as.character(partYdf[[endColumn]])))
      thisRes <- as(thisRes, "dgCMatrix")
      partRess[[j]] <- thisRes
      firstYCol <-  lastYCol+1
    }# for (j in seq(length(yRowSeq)))
    res <- do.call("cbind", partRess)
  }
  return(res)
}#regionOverlap
