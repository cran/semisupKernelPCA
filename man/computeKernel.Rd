\name{computeKernel}
\alias{computeKernel}
\title{
computeKernel
}
\description{
Compute a kernel matrix from a data set.
}
\usage{
computeKernel(data, type=c("gaussian", "pgaussian"))
}
\arguments{
  \item{data}{
n x d matrix of d-dimensional row-elements.
}
  \item{type}{
if "pgaussian", use the p-gaussian kernel function (see Francois et al. 2005). If "gaussian", use the classic gaussian kernel, ie p-gaussian function with p=2 and sigma set to the maximal pairwise distance between elements of the data set.
}
}
\value{
n x n kernel matrix.
}
\author{
Pierrick Bruneau
}
\references{
  Andrew Y. Ng, Michael I. Jordan, Yair Weiss\cr
  \emph{On Spectral Clustering: Analysis and an Algorithm}\cr
  Neural Information Processing Symposium 2001\cr
  \url{http://www.nips.cc/NIPS2001/papers/psgz/AA35.ps.gz}
Francois, D., Wertz, V. and Verleysen, M. \cr
\emph{About the locality of kernels in high-dimensional spaces} 
\emph{International Symposium on Applied Stochastic Models and Data Analysis 2005}, pp 238-245.
}
\seealso{
computeProjectionFromKernel
}
\examples{
data(iris)
irisdat <- as.matrix(iris[,1:4])
kernel <- computeKernel(irisdat, "pgaussian")
}

