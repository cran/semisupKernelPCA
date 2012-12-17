\name{computeStandardKernel}
\alias{computeStandardKernel}
\title{
computeStandardKernel
}
\description{
Compute a kernel matrix from a data set.
}
\usage{
computeStandardKernel(data, pgaussian=TRUE)
}
\arguments{
  \item{data}{
n x d matrix of d-dimensional row-elements.
}
  \item{pgaussian}{
if TRUE, use the p-gaussian kernel function (see Francois et al. 2005). Else use the classic gaussian kernel, ie p-gaussian function with p=2 and sigma set to the maximal pairwise distance between elements of the data set.
}
}
\value{
n x n kernel matrix.
}
\author{
Pierrick Bruneau
}
\references{
Francois, D., Wertz, V. and Verleysen, M. (2005) About the locality of kernels in high-dimensional spaces. \emph{International Symposium on Applied Stochastic
Models and Data Analysis}. 238-245.
}
\seealso{
computeCompositeKernel
}
\examples{
data(iris)
irisdat <- as.matrix(iris[,1:4])
kernel <- computeStandardKernel(irisdat)
}

