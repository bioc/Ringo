\name{nonzero-methods}
\docType{methods}
\alias{nonzero-methods}
\alias{nonzero,dgCMatrix-method}
\alias{nonzero,matrix.csr-method}
\alias{nonzero,matrix-method}
\alias{nonzero}
\title{Methods for Function nonzero}
\description{
  Auxiliary functions to retrieve the indices of non-zero elements in
  sparse matrices.
}
\section{Methods}{
\describe{
  \item{x = "dgCMatrix"}{returns the indices of non-zero elements in
    matrices of class \code{dgCMatrix}}
  \item{x = "matrix.csr"}{returns the indices of non-zero elements in
    matrices of class \code{matrix.csr}}
  \item{x = "matrix"}{returns the indices of non-zero elements in
    matrices of base class \code{matrix}; equivalent to
    \code{which(x != 0, arr.ind=TRUE)}}
}}
\value{
  A two-column matrix. Each row gives the row and column index of a
  non-zero element in the supplied matrix \code{x}.
}
\note{
  Originally we used the \code{matrix.csr} class from SparseM, but we
  have switched to the class \code{dgCMatrix} from package Matrix, as
  that package is part of the R distribution bundle now.

  The idea is to have a function similar to
  \code{which(x != 0, arr.ind=TRUE)} if \code{x} is a matrix.
}
\seealso{\code{\link[Matrix]{dgCMatrix-class}}}
\examples{
  (A <- matrix(c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,-34),
               nrow=5, byrow=TRUE))
  str(A.dgc <- as(A, "dgCMatrix"))
  nonzero(A.dgc)
  A2.dgc <- cbind(A.dgc, A.dgc)
  as.matrix(A2.dgc)
  nonzero(A2.dgc)
}
\keyword{methods}
