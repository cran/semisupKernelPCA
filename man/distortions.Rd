\name{distortions}
\alias{distortions}
\title{
distortions
}
\description{
the same objects are represented by two different datasets: this function computes the pairwise distances in both datasets, and the distortions induced by these distance matrices (see Aupetit, 2007)
}
\usage{
distortions(origdata, projdata)
}
\arguments{
  \item{origdata}{
n x d original dataset
}
  \item{projdata}{
n x q projected/transformed dataset
}
}
\value{
list object, with components:
- compress : vector of n compressive distortion values
- stretch : vector of n stretching distortion values
}
\author{
Pierrick Bruneau
}
\references{
Aupetit, M. (2007) Visualizing distortions and recovering topology in continuous projection techniques. \emph{Neurocomputing}. Pages 1304-1330.
}
\seealso{
computeKernel
computeProjectionFromKernel
}
\examples{
data(iris)
irisdat <- as.matrix(cbind(iris[,1:4], rep(NA, 150)))
irisdat[sample(1:150,4),5] <- 1:4
proj <- computeProjectionFromKernel(computeKernel(irisdat))
dists <- distortions(irisdat[,1:4], proj)
}

