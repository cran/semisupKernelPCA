\name{computeKernelProjection}
\alias{computeKernelProjection}
\title{
computeKernelProjection
}
\description{
Perform the semi-supervised Kernel PCA projection (see Bishop, 2006) directly from a data set, using the semi-supervised kernel transformation method in (Bruneau and Otjacques, 2012). Somehow combines computeCompositeKernel and computeProjectionFromKernel.
}
\usage{
computeKernelProjection(data, dims=2, alpha=3, method="neighbors", pgaussian=TRUE)
}
\arguments{
  \item{data}{
n x (d+1) matrix of d-dimensional row-elements. Last matrix column is a label specification in 1..n (NA if element is unlabelled)
}
\item{dims}{
number of dimensions to which the projection is restricted.
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
n x dims matrix of row-elements.
}
\author{
Pierrick Bruneau
}
\references{
Bishop, C. M. (2006) \emph{Pattern recognition and machine learning}. Springer.
Bruneau, P. and Otjacques, B. (2012) Including semi-supervision in a kernel matrix, with a view to interactive visual clustering. Tech Report hal-00751407, CRP Gabriel Lippmann.
Francois, D., Wertz, V. and Verleysen, M. (2005) About the locality of kernels in high-dimensional spaces. \emph{International Symposium on Applied Stochastic Models and Data Analysis}. 238-245.
}
\seealso{
computeCompositeKernel
computeProjectionFromKernel
}
\examples{
data(iris)
irisdat <- as.matrix(cbind(iris[,1:4], rep(NA, 150)))
irisdat[sample(1:150,4),5] <- 1:4
proj <- computeKernelProjection(irisdat)
}

