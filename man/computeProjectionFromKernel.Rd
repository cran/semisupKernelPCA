\name{computeProjectionFromKernel}
\alias{computeProjectionFromKernel}
\title{
computeProjectionFromKernel
}
\description{
Compute the kernel PCA projection from an input kernel matrix.
}
\usage{
computeProjectionFromKernel(kernel, dims=2, eigentype=c("basic", "irlba"))
}
\arguments{
  \item{kernel}{
n x n kernel matrix
}
  \item{dims}{
number of output dimensions for the projection. Cannot exceed n.
}
	\item{eigentype}{
	Indicates the eigendecomposition routine that should be used, either the standard ("basic"), or an optimization for the extraction of a few major eigenpairs ("irlba")
}
}
\value{
n x dims matrix of the projected data.
}
\author{
Pierrick Bruneau
}
\references{
Bishop, C. M. (2006) \emph{Pattern recognition and machine learning}. Springer.
}
\seealso{
computeStandardKernel
}
\examples{
data(iris)
irisdat <- as.matrix(iris[,1:4])
kernel <- computeKernel(irisdat)
proj <- computeProjectionFromKernel(kernel)
}


