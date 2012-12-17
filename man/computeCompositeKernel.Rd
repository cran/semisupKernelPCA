\name{computeCompositeKernel}
\alias{computeCompositeKernel}
\title{
computeCompositeKernel
}
\description{
compute a kernel matrix that combines a standard kernel with a user-specified partial labelling. Labellings are handled according to "method" parameter.
}
\usage{
computeCompositeKernel(data, alpha=3, method="neighbors", pgaussian=TRUE)
}
\arguments{
  \item{data}{
n x (d+1) matrix of d-dimensional row-elements. Last matrix column is a label specification in 1..n (NA if element is unlabelled)
}
  \item{alpha}{
parametrizes the kernel matrix transformation (see Bruneau and Otjacques, 2012)
}
  \item{method}{
transformation method, being either "simple" or "neighbors" (see Bruneau and Otjacques, 2012)
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
Bruneau, P. and Otjacques, B. (2012) Including semi-supervision in a kernel matrix, with a view to interactive visual clustering. Tech Report hal-00751407, CRP Gabriel Lippmann.
Francois, D., Wertz, V. and Verleysen, M. (2005) About the locality of kernels in high-dimensional spaces. \emph{International Symposium on Applied Stochastic Models and Data Analysis}. 238-245.
}
\seealso{
computeStandardKernel
computeKernelProjection
}
\examples{
data(iris)
irisdat <- as.matrix(cbind(iris[,1:4], rep(NA, 150)))
irisdat[sample(1:150,4),5] <- 1:4
kernel <- computeCompositeKernel(irisdat)
}



