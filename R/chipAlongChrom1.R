## function to plot ChIP-chip intensities along chromosome, inspired by
##  tilingArray's plotAlongChrom

chipAlongChrom1 <- function (eSet, chrom, probeAnno, xlim, ylim=NULL, samples=NULL, paletteName="Dark2", colPal=NULL, byStrand = FALSE, ylabel="fold change [log]", rugCol="#000010", itype="r", ipch=20,icex=1, ilwd=3, ilty=1, useGFF=TRUE, gff=NULL, featCol="darkblue", zero.line=TRUE, putLegend=TRUE, add=FALSE, maxInterDistance=200, coord=NULL, verbose=TRUE, ...)
{
  # 0. check arguments  
  stopifnot(inherits(eSet,"ExpressionSet"), inherits(probeAnno, "probeAnno"), validObject(probeAnno))
  eSetProbeNames <- featureNames(eSet)
  if (is.null(eSetProbeNames))
    stop("Could not determine probe identifiers from expression set.\nCheck 'featureNames' or 'geneNames' of expression set.\n")
  if (is.null(samples)) samples <- 1:ncol(exprs(eSet))
  if (!is.null(coord) && missing(xlim)){
      stopifnot(is.numeric(coord), length(coord)==2)
      xlim <- coord
  }
  thisCall <- match.call()
  # 1. get intensities in selected region
  if (verbose) cat("Getting probe intensities in selected regions..,\n")

  # get probes on chromosome:
  sta = probeAnno[paste(chrom, "start", sep=".")]
  end = probeAnno[paste(chrom, "end",   sep=".")]
  mid <- round((sta+end)/2)
  ind = probeAnno[paste(chrom, "index", sep=".")]
  uni = probeAnno[paste(chrom, "unique",   sep=".")]
  names(mid) <- ind
  nProbesOnChr <- length(mid)
  if (!missing(xlim)){
    stopifnot(length(xlim)==2, xlim[1]>0, xlim[2]>xlim[1])
    areProbesInLimits <- (mid>=xlim[1])&(mid<=xlim[2])
    usedProbes <- mid[areProbesInLimits]
    usedProbesCol <- as.numeric(uni[areProbesInLimits]!=0)+1
  } else {usedProbes <- mid}

  if (length(usedProbes) < 1)
    stop("No reporter-mapped positions in specified region!\n")
  nSamples <- length(samples)
  usedProbesIdx <- match(names(usedProbes),eSetProbeNames)
  if (any(is.na(usedProbesIdx)))
    warning(paste("The identifiers of", sum(is.na(usedProbesIdx)),
       "reporters in the region to plot are not found as",
       "'featureNames' of", deparse(substitute(eSet)),"\n"))
  
  chromExprs <- exprs(eSet)[usedProbesIdx, samples, drop=FALSE]
  if (all(is.na(as.vector(chromExprs))))
    warning("Only NA values in specified region.\n")
  
  # 2. select colors for plotting
  if (verbose) cat("Preparing color scheme...\n")
  if (is.null(colPal)){
    if (nSamples > 9)
      colPal <- sample(colors(), 9)
    else 
      colPal <- suppressWarnings(brewer.pal(nSamples, paletteName))
  }
  
  # 3. do plotting
  if (is.null(ylim))
    ylim <- range(exprs(eSet), na.rm=TRUE)
  
  absProbes <- abs(unlist(usedProbes))
  rangeX   <- xlim  # previously: range(absProbes)
  rangeX <- rangeX + round(diff(rangeX)*c(-0.05,0.05))

  if (itype %in% c("r","u")){
    interProbeDistances <- diff(absProbes)
    closeProbeClusters <- Ringo:::clusters(interProbeDistances <= maxInterDistance)
  }# if (itype %in% c("r","u"))
  if (verbose) cat("Plotting intensities...\n")
  xaxisSide <- ifelse(useGFF, 3, 1) # on which side of plot to draw genome coordinates, top or bottom
  if (!add) {# should a new plot be generated? (DEFAULT:yes)
    plot(chromExprs[,1], xlim=rangeX, ylim=ylim,
         xlab=NA, xaxt="n", ylab=ylabel, type="n", frame.plot=FALSE, ...)
    xaxisSide <- ifelse(useGFF, 3, 1) # on which side of plot to draw genome coordinates, top or bottom
    axis(side=xaxisSide)  # add x-axis and axis label next line
    mtext(paste("Chromosome",chrom,"Coordinate [bp]"), side=xaxisSide, line=2.5, font=2)
    if (zero.line) abline(h=0, lty=2)
  }#if (!add)
  
  if (length(usedProbes) > 0) {
    for (i in 1:nSamples){
      if (itype %in% c("r","u")){
        if (nrow(closeProbeClusters)>0){
          for (j in 1:nrow(closeProbeClusters)){
            clusterPos <- closeProbeClusters[j,1]+(0:closeProbeClusters[j,2])
            lines(x=absProbes[clusterPos], y=chromExprs[clusterPos,i], col=colPal[i], lwd=ilwd, lty=ilty, type=switch(itype,"r"="l","u"="c"), cex=icex)
          }#for (j in 1:nrow(closeProbeClusters))
        }#if (nrow(closeProbeClusters)>0)
        points(absProbes, y=chromExprs[,i], col=c(colPal[i],ifelse(itype=="p","grey",colPal[i]))[usedProbesCol], lwd=ilwd, type="p", pch=ipch, cex=icex)
      } else {
        points(x=absProbes, y=chromExprs[,i], col=c(colPal[i],ifelse(itype=="p","grey",colPal[i]))[usedProbesCol], lwd=ilwd, lty=ilty, type=itype, pch=ipch, cex=icex)
      }# if (itype %in% c("r","u"))
    }#for (i in 1:nSamples)
    # if we are plotting points, we can colour non-unique probe signals grey,
    #  if not, we have to take care that we don't color everything grey
    if (!add){
      rug(absProbes, side=xaxisSide, col=rugCol, lwd=2)
      if (any(usedProbesCol==2))
        rug(absProbes[usedProbesCol==2], side=xaxisSide, col="grey", lwd=2)
    }#if (!add)
  } #if (length(usedProbes) < 1)
  
  if (putLegend){
    if (is.character(samples))
      samples <- match(samples,sampleNames(eSet))
    legend(x=ifelse(add,"bottomleft","topleft"), legend=sampleNames(eSet)[samples], fill=colPal, bty="n")
  }#if (putLegend)
  
  # 4. annotate genomic features as well
  if (all(useGFF,!is.null(gff),!add)){
    stopifnot(is.data.frame(gff), all(c("name","chr","strand","start","end")%in%names(gff)))
    if (verbose) cat("Obtain genomic features...\n")
    if (! "symbol" %in% names(gff)){ gff$symbol <- gff$name}
    else { gff$symbol[gff$symbol==""] <- gff$name[gff$symbol==""]}
    mtline <- 0 ; mcex <- 1
    areOnChrom <- (gff$chr==chrom)
    genestrand <- ifelse(gff$strand[areOnChrom]%in%c("+",1),1,-1)
    genestarts <- ifelse(genestrand>0,gff$start[areOnChrom],
                         gff$end[areOnChrom])
    geneends <- ifelse(genestrand>0,gff$end[areOnChrom],
                       gff$start[areOnChrom])
    extendBeyond <-  (abs(gff$start[areOnChrom])<=rangeX[1] & abs(gff$end[areOnChrom])>=rangeX[2])
    areIn <- ((abs(genestarts)>=rangeX[1] & abs(genestarts)<=rangeX[2])|
              (abs(geneends)>=rangeX[1] & abs(geneends)<=rangeX[2]) |
              extendBeyond)
    genestarts <- genestarts[areIn]
    geneends  <- geneends[areIn]
    genestrand <- genestrand[areIn]
    transdirection <- ifelse(genestrand>0,">","<")
    chrgene <- gff$symbol[areOnChrom][areIn]
    
    if (length(genestarts)>0) {
      symbolSpacing <- round(diff(rangeX)/100)
      for (i in 1:length(genestarts)){
        mtext(transdirection[i], at=seq(abs(genestarts[i]),abs(geneends[i]), by=genestrand[i]*symbolSpacing), line=mtline, col=featCol, cex=mcex, side=1)
        sym.pos <- abs(genestarts[i])
        sym.label <- chrgene[i]
        sym.font <- ifelse((sym.pos<rangeX[1])|(sym.pos>rangeX[2]),3,1)
        sym.adj <- 0.5
        if (sym.pos < rangeX[1]){
          sym.pos   <- rangeX[1]
          sym.label <- paste("->",sym.label)
          sym.adj <- 0
        }
        if (sym.pos > rangeX[2]){
          sym.pos   <- rangeX[2]
          sym.label <- paste(sym.label,"<-")
          sym.adj <- 1
        }
        mtext(sym.label, at=sym.pos, line=mtline+1, adj=sym.adj, col=featCol, side=1, cex=mcex, font=sym.font)
      }#for
        
    }#if (length(chrstarts)<0)
  }#if (useGFF & ("gff" %in% names(chromLocObj)))  
  invisible(chromExprs)  
}# chipAlongChrom1
