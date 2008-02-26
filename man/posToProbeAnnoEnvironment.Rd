\name{posToProbeAnnoEnvironment}
\alias{posToProbeAnnoEnvironment}
\alias{createProbeAnnoEnvironment}
\title{Function for creating a probeAnno environment}
\description{
  This function allows the user to create a probeAnno environment that
  holds the mapping between probes on the array and their genomic match
  position(s). As input, the function takes either
  a.) one of NimbleGen's POS file or a similar file that holds the
  mapping of probes to the genome. OR
  b.) a \code{data.frame} holding this information
}
\usage{
posToProbeAnnoEnvironment(pos, chrNameColumn = "CHROMOSOME",
   probeColumn = "PROBE_ID", chrPositionColumn = "POSITION",
   lengthColumn = "LENGTH", verbose = TRUE, ...)
}
\arguments{
  \item{pos}{either a file-name that specifies the path to the POS or
    other mapping file OR a data.frame holding the mapping}
  \item{chrNameColumn}{name of the column in the file or
    \code{data.frame} that holds the chromosome name of the match}
  \item{probeColumn}{name of the column that holds the matching probe's
    unique identifier}
  \item{chrPositionColumn}{name of the column that holds the match
    genomic position/coordinate on the chromosome}
  \item{lengthColumn}{name of the column that holds the length of the 
    match position, in case of perfect match should correspond to the
    sequence length of the probe}
  \item{verbose}{logical; should progress be written to STDOUT?}
  \item{\dots}{further arguments passed on to function \code{read.delim}
    if \code{pos} is a file path.}
}
\details{
  The default column names correspond to the column names in a NimbleGen
  POS file.

  For custom mappings, using the tools Exonerate, BLAT or MUMmer, the
  scripts directory of this package holds Perl scripts to generate such
  a POS file from the respective output files.
}
\value{
  The results is an environment that holds four vectors for each
  chromosome/strand, namely, say for chromosome 1:
  \item{1.start}{genomic start coordinates of probe matches on
    chromosome 1}
  \item{1.end}{genomic start coordinates of probe matches on
    chromosome 1}
  \item{1.index}{identifier of probes matching at these coordinates}
  \item{1.unique}{vector of the same length as the three before; encoding
    how many matches the corresponding probe has in the given file or
    \code{data.frame}. An entry of '0' indicates that the probe matching
    at this position has only this one match.}
}
\note{Users do not have to use this function but can generate such an
  environment any way they like. We recommend to check the result using
  the function \code{validProbeAnno}.}
\author{Joern Toedling \email{toedling@ebi.ac.uk}}
\seealso{\code{\link{validProbeAnno}},\code{\link{read.delim}}}
\examples{
 exPos <- read.delim(file.path(system.file("exData",package="Ringo"),
                     "MOD_2003-12-05_SUZ12_1in2.pos"),
                     header=TRUE,as.is=TRUE)
 str(exPos)
 exProbeAnno <- posToProbeAnnoEnvironment(exPos)
 ## is equivalent to 
 exProbeAnno2 <-
posToProbeAnnoEnvironment(file.path(system.file("exData",package="Ringo"),
                          "MOD_2003-12-05_SUZ12_1in2.pos"))
 ls(exProbeAnno)
}
\keyword{manip}