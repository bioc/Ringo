
### new version of chipAlongChrom using lattice graphics
##   adopted from the tilingArray package

# plotting along chromosome (extra bar/s, gff optional,
#  legend at bottom optional)
# also need ability to change the color scheme for panel


### also provide an S4 method as a wrapper
setMethod("plot",signature=c("ExpressionSet","probeAnno"),function(x, y,...){
  chipAlongChrom(eSet=x, probeAnno=y, ...)
})

chipAlongChrom <- function(eSet, probeAnno, chrom, xlim, ylim,
                           samples=NULL, paletteName="Set2",
                           colPal=NULL, ylab="Fold change [log]",
                           ipch=16, ilwd=3, ilty=1, icex=3, gff = NULL,
                           featureExclude=c("chromosome", "nucleotide_match",
                             "insertion"), zeroLine=TRUE, sampleLegend=TRUE,
                           sampleLegendPos="topleft",
                           featureLegend=FALSE, maxInterDistance=200,
                           coord=NULL, highlight, main, ...)
{
  require("grid")
  ########################################################
  ## check arguments
  ##########################################################
  stopifnot(inherits(eSet,"ExpressionSet"),
            inherits(probeAnno, "probeAnno"),
            validObject(probeAnno))
  eSetProbeNames <- featureNames(eSet)
  if (is.null(samples)) samples <- 1:ncol(exprs(eSet))
  if (!is.null(coord) && missing(xlim)){
    stopifnot(is.numeric(coord), length(coord)==2)
    xlim <- coord
  }
  thisCall <- match.call()
  sampleLegendPos <- match.arg(sampleLegendPos,
     c("topleft","topright","bottomleft", "bottomright"))
  
  ################################################
  ## set up the viewports of the plot layout.
  #########################################################
  n <- length(samples)
  VP <- c("title"=0.4, "expr+"=10, "gff+"=1, "coord"=1, "gff-"=1, "legend"=1)
  
  if(!featureLegend)
    VP = VP[-which(names(VP)=="legend")]
  if(missing(gff))
     VP = VP[-which(names(VP)=="gff+" | names(VP)=="gff-")]
  defaultColors = c("+" = "#00441b", "-" = "#081d58", "duplicated" = "grey",
    "cp" = "#555555", "ci" = "#777777", "highlight" = "red",
    "threshold" = "black","rmp"= "#010101") # rmp: reporter-match position
  if (!is.null(colPal))
    colPal <- rep(colPal, length.out=ncol(eSet))
  else
    colPal <- brewer.pal(8, paletteName)[1:ncol(eSet)]

  ## plot margin
  pushViewport(viewport(width=0.85, height=0.95)) ## plot margin
  pushViewport(viewport(layout=grid.layout(length(VP), 1, heights=VP)))

  strand <- "+"
   ## interpret input data and probe mapping
  ## eSet and probeAnno
  y <- exprs(eSet)
  stopifnot(is.matrix(y))
  
  ## get probe annotation on that chromosome from probeAnno
  sta   <- probeAnno[paste(chrom, "start", sep=".")]
  end   <- probeAnno[paste(chrom, "end",   sep=".")]
  mid   <- round((sta+end)/2)
  index <- probeAnno[paste(chrom, "index", sep=".")]
  uni   <- probeAnno[paste(chrom, "unique",   sep=".")]
  names(mid) <- index
  nProbesOnChr <- length(mid)
  if (!missing(xlim)){
    stopifnot(length(xlim)==2, xlim[1]>0, xlim[2]>xlim[1])
    areProbesInLimits <- (mid>=xlim[1])&(mid<=xlim[2])
    usedProbes <- mid[areProbesInLimits]
    usedProbesCol <- as.numeric(uni[areProbesInLimits]!=0)+1
    uni <- uni[areProbesInLimits]
  } else {usedProbes <- mid}
  if (length(usedProbes) < 1)
    stop("No reporter-mapped positions in specified region!\n")
  nSamples <- length(samples)
  usedProbesIdx <- match(names(usedProbes),eSetProbeNames)
  if (any(is.na(usedProbesIdx)))
    warning(paste("The identifiers of", sum(is.na(usedProbesIdx)),
                  "reporters in the region to plot are not found as",
                  "'featureNames' of", deparse(substitute(eSet)),"\n"))

  dat <- list(x = usedProbes,
              y = y[usedProbesIdx, samples, drop=FALSE],
              flag = uni)
  stopifnot(is.numeric(dat$flag))
  lengthChr <- end[length(end)]
    
  ## At this point, no matter what the input to the function was,
  ##    we have the list 'dat' with elements
  ## x: x coordinate
  ## y: y coordinate
  ## flag
  
  ## if no region is specified, plot the whole chromosome
  if(missing(coord) && missing(xlim))
    xlim <- coord <- c(1, lengthChr)

  ## set up the y-axis limits
  ylimdata <- range(as.vector(dat[["y"]][dat[["x"]]>=xlim[1] & dat[["x"]]<=xlim[2],]), na.rm=TRUE)
  if (missing(ylim)) ylim <- ylimdata + c(-0.15, 0.15)*diff(ylimdata)
  
  ## plot the data
  vpr <- which(names(VP)==sprintf("expr%s", strand))
  
  ##############################################################
  ### set up viewport for expression plot ######################
  ##############################################################
  pushViewport(dataViewport(xData=xlim, yData=ylim, extension=0, clip="off",
                            layout.pos.col=1, layout.pos.row=vpr))
  grid.yaxis(gp=gpar(cex=0.8),...)
  grid.text(ylab, x=-0.075, y=0.6, just=c("right", "center"), default.units="npc", rot=90, gp=gpar(cex=0.9, font=2),...)
  pushViewport(dataViewport(xData=xlim, yData=ylim, extension=0, clip="on",
                            layout.pos.col=1, layout.pos.row=vpr))
  
  for(j in seq(1:n)) {
    datj <- dat
    datj$y <- dat$y[,j]
    ## plot the data
    plotOneChIPSample(datj, xlim=xlim, chr=chrom, strand="+", vpr=vpr,
                      sampleColor=colPal[j],
                      maxInterDistance=maxInterDistance,
                      ipch=ipch, ilwd=ilwd, icex=icex,...)
  }# for(j in seq(1:n))

  ###-----------------------------------------------------------------
  ### sample legend
  ###-----------------------------------------------------------------
  if (sampleLegend){
    # determine positions of sample legend
    sLPos <- switch(sampleLegendPos,
                    "topleft"=list(x=unit(0, "npc"), y = unit(1, "npc"),
                      just=c("left","top")),
                    "bottomleft"=list(x=unit(0, "npc"), y = unit(0, "npc"),
                      just=c("left","bottom")),
                    "topright"=list(x=unit(1, "npc"), y = unit(1, "npc"),
                      just=c("right","top")),
                    "bottomright"=list(x=unit(1, "npc"), y = unit(0, "npc"),
                      just=c("right","bottom")))
    
    draw.key(key=list(rectangles=list(col=colPal),
               text=list(sampleNames(eSet)[samples], cex=0.8)),
             vp=viewport(x = sLPos$x, y = sLPos$y,
               width=max(strwidth(sampleNames(eSet)[samples],
                 units="figure"))*1.2+0.05,
               height=strheight("abc", units="figure")*1.2*n,
               just=sLPos$just, clip="off"), draw=TRUE)
  } # if (sampleLegend)
    
  ### end plotting the data
  popViewport(2)

  ############################################################################
  ### plot annotated genome features (supplied in the gff) 
  ############################################################################
  if(!missing(gff))
    for (strand in c("+","-"))
      plotFeatures(gff=gff, chr=chrom, xlim=xlim, strand=strand, 
                   featureExclude=featureExclude, featureColorScheme=1,
                   vpr=which(names(VP)==sprintf("gff%s", strand)),...)

  #############################################################################
  ######### plot  chromosomal coordinate axis #################################
  #############################################################################
  pushViewport(dataViewport(xData=xlim, yscale=c(-0.4,0.8), extension=0, 
                            layout.pos.col=1, layout.pos.row=which(names(VP)=="coord")))
  grid.lines(xlim, c(0,0), default.units = "native")
  tck <- alongChromTicks(xlim)
  grid.text(label=formatC(tck, format="d"), x=tck, y=0.2, just=c("centre", "bottom"), gp = gpar(cex=.6), default.units="native")
  grid.segments(x0 = tck, x1 = tck, y0 = -0.17, y1 = 0.17,  default.units = "native")

  #############################################################################
  ####### plot probe positions ################################################
  #############################################################################
  rmpCols <- ifelse(dat$flag==0, defaultColors["rmp"],
                    defaultColors["duplicated"])
  stopifnot(length(rmpCols)==length(usedProbes))
  grid.segments(x0 = usedProbes, x1 = usedProbes, y0 = -0.01, y1 = 0.1,
                default.units="native", gp=gpar(lwd=2, col=rmpCols))
  
  if(!missing(highlight)){
    ## this part was modified to draw arrows for transcripts rather than bars
    mt = (match(highlight$strand, c("-", "+"))-1.5)*2
    if (!is.null(highlight$coord))
      co <- highlight$coord
    else
      co <- highlight$x
    if(is.na(mt) || !is.numeric(co))
      stop("Invalid parameter 'highlight'.")
    strand.num <- ifelse(highlight$strand=="-",-1,1)
    grid.segments(x0=co, x1=co+(500*strand.num), y0=c(0.4,0.4)*mt,
                  y1=c(0.4,0.4)*mt, default.units = "native", arrow=arrow(),
                  gp=gpar(col="violetred4", lwd=4))
  }
  popViewport()
  
  #### TITLE
  pushViewport(viewport(layout.pos.col=1,
                        layout.pos.row=which(names(VP)=="title")))
  grid.text(label=paste("Chromosome", chrom, "coordinate [bp]"), x=0.5, y=1,
            just="centre", gp=gpar(cex=1, font=2))
  if(!missing(main))
    grid.text(label=main, x=0.05, y=1, just="centre", gp=gpar(cex=1, font=2))
  popViewport()
  
  #### LEGEND
  if(featureLegend)
    plotAlongChromLegend(which(names(VP)=="legend"),
         featureColorScheme=1, featureExclude=featureExclude)
  popViewport(2)
  invisible(dat)
  
}#chipAlongChrom2


## ------------------------------------------------------------
## auxiliary function: plot Features
## ------------------------------------------------------------
plotFeatures <- function(gff, chr, xlim, strand, vpr, featureColorScheme=1,
                         featureExclude=c("chromosome", "nucleotide_match",
                           "insertion"), featureNoLabel=c("uORF", "CDS"),...)
{
  #### check arguments:
  stopifnot(is.data.frame(gff),
            all(c("name","chr","strand","start","end")%in%names(gff)),
            length(strand)==1, strand %in% c("+", "-"))
  
  if (! "symbol" %in% names(gff)){ gff$symbol <- gff$name}
  else { gff$symbol[gff$symbol==""] <- gff$name[gff$symbol==""]}
  stopifnot(length(strand)==1, strand %in% c("+", "-"))
  
  translateStrand <- function(x){
    if (x==-1) return("-")
    if (x==1) return("+")
    return(NA)}
  stopifnot(all(gff[,"start"] <= gff[, "end"]))
  if (is.numeric(gff$strand))
    gff$strand <- sapply(gff$strand, translateStrand)
   
  pushViewport(dataViewport(xData=xlim, yscale=c(-1.2,1.2),  extension=0,
                            clip="on", layout.pos.col=1, layout.pos.row=vpr))

  sel <- which(gff[, "chr"]     == chr &
               gff[, "strand"]  == strand &
               gff[, "start"]   <= xlim[2] &
               gff[, "end"]     >= xlim[1])
  
  ## for label, use "symbol" if available, otherwise "name"
  featName = gff[sel, "symbol"]
  
  ## split by feature type (e.g. CDS, ncRNA)
  feature  = as.character(gff[sel, "feature"])
  featsp = split(seq(along=sel), feature)
  
  ## There are now five different cases, and we need to deal with them:
  ## - ignorable features, given by featureExclude
  ## - genes: a horizontal line + name
  ## - introns: a caret
  ## - CDS: a box + no name
  ## - all others: a colored box + name

  ## in this vector we save those features for which we want to have names
  whnames = integer(0)

  ## 1. drop the ignorable ones
  featsp = featsp[ ! (names(featsp) %in% featureExclude) ]
  
  ## 3.introns
  wh = ("intron" == names(featsp))
  if(any(wh)) {
    i = featsp[["intron"]]
    s = sel[i]
    mid = (gff$start[s]+gff$end[s])/2
    wid = (gff$end[s]-gff$start[s])/2 
    for(z in c(-1,1))
      grid.segments(x0 = mid,
                    x1 = mid+z*wid,
                    y0 = 1.20*c("+"=1, "-"=-1)[strand],  ## istrand is 1 or 2
                    y1 = 0.95*c("+"=1, "-"=-1)[strand],
                    default.units = "native",
                    gp = gpar(col="black"))
     featsp = featsp[!wh]
  } ## if
  
  ## 4. colors for boxes
  ## check that we know how deal with all features
  featCols = featureColors(featureColorScheme)

  whm = names(featsp) %in% rownames(featCols)
  if(!all(whm))
    warning("Don't know how to handle feature of type(s) '",
            paste(names(featsp)[!whm], collapse=", "), "' in gff.", sep="")

  sfeatsp  = featsp[rownames(featCols)]
  ll       = listLen(sfeatsp)
  
  if(any(ll>0)) {
    i  = unlist(sfeatsp)
    gp = gpar(col = rep(featCols$col,  ll),
                 fill = rep(featCols$fill, ll))
    s  = sel[i]
    grid.rect(x     = gff$start[s],
              y     = 0,
              width = gff$end[s]-gff$start[s],
              height= 2,
              default.units = "native",
              just  = c("left", "center"),
              gp    = gp)
    whnames = c(whnames, unlist(sfeatsp[!(names(sfeatsp) %in% featureNoLabel)]))
    ## additional potentially useful values for featureNoLabel: "binding_site", "TF_binding_site"
  }

  ## labels
  if( !all(tolower(featureNoLabel)=="all") && (length(whnames)>0)) {

    ## this is a bit of a hack to abbreviate the labels of
    ##  "binding site" features:
    bindingRegexpr = "binding.?site.*$"
    isBindingSite = (regexpr(bindingRegexpr, featName[whnames]) > 0)
    if(any(isBindingSite)) {
      ## replace long labels
      featName[whnames] = gsub(bindingRegexpr, "bs", featName[whnames])
    }

    ## remove duplicated names that are not binding sites
    whnames = whnames[isBindingSite | !duplicated(featName[whnames])]

    txtcex = 0.6
    txtdy  = 0.7
    s      = sel[whnames]
    txtx   = (pmax(min(xlim), gff$start[s])+ pmin(max(xlim), gff$end[s]))/2
    txty   = numeric(length(s))
    ord    = order(txtx)
    whnames = whnames[ord]
    s      = s[ord]
    txtx   = txtx[ord]
    
    strw   = convertWidth(stringWidth(featName[whnames]), "native", valueOnly=TRUE)*txtcex
    rightB = txtx[1] + 0.5*strw[1]
    doText = rep(TRUE, length(whnames))

    # adjust text labels to be still readable in feature-dense areas:
    if(length(whnames) >1) {
      for(k in 2:length(whnames)) {
        leftB = txtx[k] - 0.5*strw[k]
        if(leftB > rightB) { # all texts not overlapping next to each other?
          rightB = txtx[k] + 0.5*strw[k]
        } else { # any overlaps?
          if(!any(txty[k-(1:2)]==txtdy)) {#  2 previous labels not moved up?
            txty[k]= txtdy                #   then this one 
          } else {                        #  else try move down:
            if(!any(txty[k-(1:2)]== -txtdy)) { 
              txty[k]= -txtdy             #  if 2 previous ones weren't
            } else {
              doText[k] = FALSE           #  otherwise don't put the label
            }
          }
        } ##  else
      } ## for
    }
    grid.text(label = featName[whnames][doText],
              x = txtx[doText], y = txty[doText], gp=gpar(cex=txtcex), 
              default.units = "native")
  } ## if
  popViewport()
} ## plotFeatures

##------------------------------------------------------------
##
##------------------------------------------------------------
alongChromTicks <- function(x){
  rx = range(x)
  lz = log((rx[2]-rx[1])/3, 10)
  fl = floor(lz)
  if( lz-fl > log(5, 10))
    fl = fl +  log(5, 10)
  tw = round(10^fl)
  i0 = ceiling(rx[1]/tw)
  i1 = floor(rx[2]/tw)
  seq(i0, i1)*tw
}# alongChromTicks

##------------------------------------------------------------
## featureColors
## note that features are drawn in the order in which they appear
## here, this can be used to let important features overdraw less
## important ones (e.g. tRNA is more specific than ncRNA)
## to test, say tilingArray:::plotAlongChromLegend()
##------------------------------------------------------------
featureColors <- function(scheme=1){
  defaultColors = c(
    "centromere"               = "#FFEDA0",    ## orange
    "telomere"                 = "#FFEDA0",    ## orange
    "novel_pseudogene"         = "#f0f0f0",    ## light gray
    "novel_retrotransposed"    = "#f0f0f0",    ## light gray
    "novel_coding"             = "#e0f1f2",    ## lighter blue
    "novel_RNA"                = "#b6e97a",    ## green 
    "novel_tRNA"               = "#b6e97a",    ## green 
    "novel_scRNA"              = "#9C9BC1",    ## purple
    "novel_snRNA"              = "#9C7BC1",    ## purple
    "novel_rRNA"               = "#ffbe71",    ## meat
    "novel_snoRNA"             = "#8F6A68",    ## red brown
    "novel_miRNA"              = "#DC76DC",    ## light red-violet
    "pseudogene"               = "#e0e0e0",    ## light gray
    "uORF"                     = "#FED976" ,   ## orange
    "nc_primary_transcript"    = "#a0a0a0",    ## grey
    "repeat_family"            = "#CC6666",    ## light red
    "repeat_region"            = "#e31a1c",    ## bright red
    "retrotransposed"          = "#f1b6da",    ## pink
    "transposable_element"     = "#f1b6da",    ## pink
    "transposable_element_gene"= "#f1b6da",
    "ARS"         = "#CC9966",             ## light brown
    "insertion"   = "#FFEDA0",             ## orange
    "CDS_dubious" = "#e0f1f2",             ## lighter blue
    "gene"        = "#addfff",             ## light blue
    "CDS"         = "#addfff",             ## light blue
    "coding"      = "#5E88B0",             ## light blue
    "exon"        = "#5E88B0",             ## aquamarine
    "transcript"  = "#5E88B0",             ## aquamarine
    "ncRNA"       = "#a6d96a",             ## green 
    "tRNA"        = "#a6d96a",             ## green
    "snRNA"       = "#8C6BB1",             ## purple
    "rRNA"        = "#fdae61",             ## meat
    "snoRNA"      = "#7F5A58",             ## red brown
    "miRNA"       = "#cc66cc",             ## light red-violet
    "nucleosome_binding_motif" = "#C9C299",## lemon chiffon
    "TF_binding_site" = "#C9C299"          ## lemon chiffon
    )
  darkenborder <- as.logical(c(rep(1, length(defaultColors)-2),0, 0))
  stopifnot(length(darkenborder)==length(defaultColors))
  
  fill = switch(scheme,
    default  = defaultColors,
    unicolor = ifelse(is.na(defaultColors), NA,  "#addfff"),  ## light blue
    stop("Encountered error when filling in colors."))
  
  ## calculate hex string for a color that is a little bit darker than the
  ## hex string in the argument
  darken <- function(x, factor=0.5) {
    wh = which(!is.na(x))
    hex = sapply(x[wh], substring, first=c(2,4,6), last=c(3,5,7))
    hex = apply(hex, 2, function(h) as.integer(factor*as.integer(paste("0x", h, sep=""))))
    res = rep(as.character(NA), length(x))
    res[wh] = apply(hex, 2, function(h) sprintf("#%02x%02x%02x", h[1], h[2], h[3]))
    return(res)
  }
  border <- ifelse(darkenborder, darken(fill), fill)
  res <- data.frame(fill=I(fill), col =I(border))
  rownames(res) <- names(defaultColors) 
  return(res)
} # function featureColors

##------------------------------------------------------------
## legend
##------------------------------------------------------------
plotAlongChromLegend <- function(vpr, nr=2, featureColorScheme=1,
    featureExclude=c("chromosome", "nucleotide_match", "insertion"),
    mainLegend, cexLegend=0.35, cexMain=1)
{
  endVP = FALSE
  # when this function is called on its own 
  # set up a viewport
  if(missing(vpr)) { 
     endVP=TRUE      
     vpr = newVP(main=mainLegend, dataPanelHeight=1, cexMain=cexMain) # newVP sets up a new viewport
  }
  formatRow = function(featColsOneRow, row) {
    ## print(featColsOneRow)
    strWid   = convertWidth(stringWidth(rownames(featColsOneRow)), "npc", valueOnly=TRUE)
    n        = length(strWid)
    inbetWid = 0.2*min(strWid)
    totWid   = sum(strWid)+(n-1)*inbetWid
    x        = c(0, cumsum(strWid[-n])) + (0:(n-1))*inbetWid 
    y        = numeric(length(x))

    x      = x/totWid
    strWid = strWid/totWid
    grid.rect(x = x, width = strWid, 
              y = unit(row, "native"), height = unit(1, "native")- unit(1, "mm"), 
              just  = c("left", "center"), default.units="npc",
              gp    = do.call(gpar, featColsOneRow))
    
    grid.text(label = rownames(featColsOneRow),
              x = unit(x + strWid/2, "native"), y = unit(row, "native"),
              just  = c("center", "center"), gp=gpar(cex=cexLegend))
  } 

  featCols = featureColors(featureColorScheme)
  featCols = featCols[ !(rownames(featCols) %in% featureExclude), ]

  pushViewport(viewport(layout.pos.col=1, layout.pos.row=vpr, yscale=c(0.5, nr+0.5)))

  i = 1:nrow(featCols)
  for(r in 1:nr)
    formatRow(featCols[ceiling(i/nrow(featCols)*nr-1e-10)==r, ], row=nr-r+1)
  
  popViewport()

  if(endVP)
     popViewport(2)
}

# this function sets up a new viewport.  It is used by plotAlongChromLegend, 
# plotSegmentationHeatmap and plotOneChIPSample when they are called as 
# stand-alone functions (ie when vpr is not specified)
newVP <- function(main, cexMain=1, dataPanelHeight=1, vpHeight=0.7, titleOffSet=0) {
  if(!missing(main)) {
    vpr = c("title"=0.1, "data"=dataPanelHeight)
    pushViewport(viewport(width=0.85, height=vpHeight)) ## plot margin
    pushViewport(viewport(layout=grid.layout(length(vpr), 1, heights=vpr)))  
    pushViewport(viewport(layout.pos.col=1, layout.pos.row=which(names(vpr)=="title")))
    grid.text(label=main, x=0.5, y=1.1+titleOffSet, just="centre", gp=gpar(cex=cexMain))  
    popViewport()
    vpr = which(names(vpr)=="data")
  } else {
    vpr = c("data"=dataPanelHeight)
    pushViewport(viewport(width=0.85, height=vpHeight)) ## plot margin
    pushViewport(viewport(layout=grid.layout(length(vpr), 1, heights=vpr)))
  }
  return(vpr)
}#newVP


##-------------------------------------------------------------
##   plot one ChIP-chip sample with lines and dots
##-------------------------------------------------------------
plotOneChIPSample <- function(dat, xlim, ylim, ylab, 
  chr=1, strand="+", vpr, sampleColor=NULL, zeroLine=TRUE,
  main, pointSize=unit(1.0, "mm"), maxInterDistance=200, itype="r",
  ilwd=3, ilty=1, icex=4, ipch=16, cexAxisLabel=1, cexAxis=1,...)
{
  endVP = FALSE
  if(missing(vpr)) {
     endVP=TRUE
     vpr = newVP(main=main, dataPanelHeight=1, vpHeight=0.95, titleOffSet=-0.9)
  }

  if(is.matrix(dat$y))
    dat$y = rowMeans(dat$y) ##  if >1 samples, take mean over samples
  stopifnot(length(dat$y)==length(dat$x), length(dat$flag)==length(dat$x))
  
  xorg  = dat$x
  if(missing(xlim)) {
    xlim=range(dat$x, na.rm=TRUE)
  } else {
    sel  = (dat$x>=xlim[1])&(dat$x<=xlim[2])
    dat$x = dat$x[sel]
    dat$y = dat$y[sel]
    dat$flag = dat$flag[sel]
  }
  
  if(missing(ylim))
    ylim = quantile(dat$y, c(0,1), na.rm=TRUE)

  ## the expression data. use two viewports for different clipping behavior
  if (missing(vpr)){
    if(!missing(ylab)) {
      pushViewport(dataViewport(xData=xlim, yData=ylim, extension=0, clip="off",
                                layout.pos.col=1, layout.pos.row=vpr))
      grid.yaxis(gp=gpar(cex=cexAxis),...)
      grid.text(ylab, x=-0.075, y=0.5, rot=90, gp=gpar(cex=cexAxisLabel),...)
      pushViewport(dataViewport(xData=xlim, yData=ylim, extension=0, clip="on",
                                layout.pos.col=1, layout.pos.row=vpr))
    } else {
      pushViewport(dataViewport(xData=xlim, yData=ylim, extension=0, clip="off",
                                layout.pos.col=1, layout.pos.row=vpr))
      grid.yaxis(gp=gpar(cex=cexAxis))
      pushViewport(dataViewport(xData=xlim, yData=ylim, extension=0, clip="on",
                                layout.pos.col=1, layout.pos.row=vpr))
    }
  }
  
  defaultColors <- c("+" = "#081d58", "-" = "#081d58", "duplicated" = "grey",
                    "cp" = "#555555", "ci" = "#777777", "highlight" = "red",
                     "threshold" = "grey")
  if (is.null(sampleColor))
    sampleColor <- "#081d58"

  ord  <- sort(which(dat$flag==0))# , which(dat$flag==0))
  colFlag <- ifelse(dat$flag==0, sampleColor, defaultColors["duplicated"])

  ### draw zero line
  if (zeroLine)
    grid.lines(y=unit(0, "native"), gp=gpar(col="#000000", lty=2, lwd=2))

  ### probe positions
  absProbes <- abs(unlist(dat$x[ord]))
  rangeX   <- xlim  # previously: range(absProbes)
  rangeX <- rangeX + round(diff(rangeX)*c(-0.05,0.05))
  interProbeDistances <- diff(absProbes)
  closeProbeClusters <- Ringo:::clusters(interProbeDistances <= maxInterDistance)
  
  if (length(absProbes) > 0) {
    if (itype %in% c("r","u")){
      if (nrow(closeProbeClusters)>0){
        areInClusters <- vector("logical", length(absProbes))
        ## this next for-loop is ugly, there must be a nicer way
        for (j in 1:nrow(closeProbeClusters))
          areInClusters[closeProbeClusters[j,1]+(0:closeProbeClusters[j,2])] <- TRUE
        grid.polyline(x=absProbes[areInClusters], y=dat$y[ord][areInClusters], default.units="native", id=rep(seq(nrow(closeProbeClusters)), closeProbeClusters[,2]+1), gp=gpar(col=sampleColor, lwd=ilwd, lty=ilty))
      }#if (nrow(closeProbeClusters)>0)
      grid.points(dat$x, y=dat$y,  pch=ipch, size=pointSize*icex, gp=gpar(col=colFlag))
    } else {
      grid.points(dat$x, y=dat$y,  pch=ipch, size=pointSize, gp=gpar(col=colFlag))
    }# if (itype %in% c("r","u"))
  }
  #popViewport(2)

  if(endVP)
     popViewport(2)

} ## plotOneChIPSample
