\name{corPlot}
\alias{corrPlot}
\alias{corPlot}
\title{Function to plot correlation of different samples}
\description{
  This function can be used to visualise the (rank) correlation in
  expression data between different samples or sample groups.
}
\usage{
corPlot(eset, samples = NULL, grouping = NULL, ref = NULL,
   useSmoothScatter = TRUE, ...)
}
\arguments{
  \item{eset}{object of class \code{ExpressionSet} holding the array
    data, or a numeric matrix instead}
  \item{samples}{which samples' expression shall be correlated to each
    other; either a numeric vector of sample numbers in the
    \code{ExpressionSet} or a character vector that must be contained in
    the \code{sampleNames} of the \code{ExpressionSet}, default
    \code{NULL} means take all samples in the \code{ExpressionSet}}
  \item{grouping}{an optional factor vector defining if the correlation
    should be assessed between groups of samples, rather than individual
  samples. If two or more samples are assigned into the same group,
  the mean over these samples' expression values is taken before
  computing correlation. Default NULL means assess correlation between
  individual samples only.}
  \item{ref}{reference than only applies if argument \code{grouping} is
    given; see \code{\link[stats]{relevel}}}
  \item{useSmoothScatter}{logical; should the function
    \code{\link[graphics]{smoothScatter}} be used?}
  \item{\dots}{additional arguments, not used yet}
}
\value{
  The function only returns NULL (invisible). 
  The function is called for its side-effect
  producing the pairs plot.
}
\author{Joern Toedling}
\seealso{\code{\link[Biobase:class.ExpressionSet]{ExpressionSet}},
  \code{\link[stats]{relevel}},  
  \code{\link[graphics]{pairs}},
  \code{\link[graphics]{smoothScatter}}}
\examples{
   data(sample.ExpressionSet)
   if (interactive()) 
   corPlot(sample.ExpressionSet,
    grouping=paste(sample.ExpressionSet$sex,
    sample.ExpressionSet$type, sep="."))
}
\keyword{hplot}% at least one, from doc/KEYWORDS
