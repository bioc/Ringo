
sliding.quantile <- function(positions, scores, half.width, prob=0.5, return.counts=TRUE) {
  stopifnot(!is.unsorted(positions), length(positions) == length(scores), half.width >= 0, prob >= 0, prob <= 1)
  res <- .Call(ringoSlidingQuantile, as.integer(positions), as.numeric(scores), as.integer(half.width), as.numeric(prob))
  if (return.counts){
    colnames(res) <- c("quantile","count")
    rownames(res) <- positions
  } else {
    res <- res[,1, drop=TRUE]
    names(res) <- positions
  }
  return(res)
}#sliding.quantile

computeRunningMedians <- function(xSet, probeAnno, modColumn="Cy5", allChr,
   winHalfSize=400, min.probes=5, quant=0.5, combineReplicates=FALSE,
   nameSuffix=".sm", checkUnique=TRUE, uniqueCodes=c(0), verbose=TRUE)
{
  ## 0. check arguments:
  stopifnot(inherits(xSet,"ExpressionSet"),
            inherits(probeAnno,"probeAnno"), validObject(probeAnno),
            is.numeric(quant), (quant>=0)&(quant<=1), length(quant)==1 )
  if (!modColumn %in% varLabels(xSet))
    stop(paste("There is no column ",deparse(substitute(modColumn)),
               " defined in the phenoData of ",
               deparse(substitute(xSet)),".\n", sep=""))
  ## which chromosomes to do the smoothing for:
  if (missing(allChr)) allChr <- chromosomeNames(probeAnno)
  stopifnot(is.character(allChr))
  # initialize result matrix:
  if (combineReplicates)
    grouping <- factor(pData(xSet)[[modColumn]])
  else
    grouping <- factor(sampleNames(xSet), levels=sampleNames(xSet))

  newExprs <- matrix(NA, nrow=nrow(exprs(xSet)), ncol=nlevels(grouping))
  rownames(newExprs) <- allProbes <- featureNames(xSet)
  
  for (chr in allChr){
    if (verbose) cat("\nChromosome",chr, "...\n")
    chrsta <- probeAnno[paste(chr,"start",sep=".")]
    chrend <- probeAnno[paste(chr,"end",sep=".")]
    chrmid <- round((chrsta+chrend)/2)
    chridx <- probeAnno[paste(chr,"index",sep=".")]
    if (checkUnique){
      chruni <- probeAnno[paste(chr,"unique",sep=".")]
      stopifnot(length(chruni)==length(chridx))
      chridx <- chridx[chruni %in% uniqueCodes]
      chrmid <- chrmid[chruni %in% uniqueCodes]
      if (length(chrmid)==0){
        warning(paste("No reporters with unique hits on chromosome",chr,".\n"))
        next}
    } #  if (checkUnique)
    if (is.character(chridx))
        chridx <- match(chridx, allProbes)
    ## check if any NAs or non-matching numbers are in chridx:
    areProbes <- is.element(chridx, seq(nrow(xSet)))
    if (!all(areProbes)){
      warning(paste(sum(!areProbes),"probes on chromosome",chr,
                    "are listed in",deparse(substitute(probeAnno)),
                    ", but do not correspond to features of",
                    deparse(substitute(xSet)),"."))
      chridx <- chridx[areProbes]
      chrmid <- chrmid[areProbes]
    }
    ### iterate over sample groups:
    for (i in 1:nlevels(grouping)){
      modSamples   <- which(grouping == levels(grouping)[i])
      if (verbose) cat(sampleNames(xSet)[modSamples],"... ")
      combined.dat <- as.vector(t(exprs(xSet)[chridx,modSamples,drop=FALSE]))
      # as.vector(t(X)) leads to columns (samples) being appended one value by one value into long vector
      combined.pos <- rep(chrmid, each=length(modSamples))
      slidingRes <- sliding.quantile(positions=combined.pos, scores=combined.dat, half.width=winHalfSize, prob=quant, return.counts=TRUE)
      slidingRes <- slidingRes[seq(1, nrow(slidingRes)+1-length(modSamples), by=length(modSamples)),,drop=FALSE]
      chrrm <- slidingRes[,"quantile"]
      slidingRes[,"count"] <- slidingRes[,"count"]/length(modSamples)
      areBelow <- slidingRes[,"count"] < min.probes
      if (any(areBelow)) chrrm[areBelow] <- NA
      stopifnot(length(chrrm) == length(chrmid))
      newExprs[chridx,i] <- chrrm
    }#for (i in 1:length(all.mods))
  } #for (chr in allChr)
  if (verbose) cat("\nConstruction result ExpressionSet...")
  # cat construct ExpressionSet of results:
  sample.labels <- paste(as.character(levels(grouping)), nameSuffix, sep="")
  if (!combineReplicates)
      newPD <- phenoData(xSet)[order(grouping),]
  else
      newPD <- new("AnnotatedDataFrame", data=data.frame(label=sample.labels, row.names=sample.labels), varMetadata=data.frame("varLabel"=c("label"),row.names=c("label")))
  newEset <- new('ExpressionSet',exprs=newExprs,  phenoData = newPD)
      # experimentData = [MIAME], annotation = [character], ...
  featureNames(newEset) <- featureNames(xSet)
  featureData(newEset)  <- featureData(xSet)
  sampleNames(newEset)  <- sample.labels
  if (verbose) cat("Done.\n")
  return(newEset)
}#computeRunningMedians
