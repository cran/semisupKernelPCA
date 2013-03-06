\name{transformKernel}
\alias{transformKernel}
\alias{transformKernel,matrix-method}
\title{Transforms a Kernel matrix as specified by a linkinfo object}
\description{
Transforms a Kernel matrix as specified by a linkinfo object
}
\usage{
\S4method{transformKernel}{matrix}(kern, linkinfo, type=c("none", "simple", "extended"), 
		linkfun=function(x) x^(1/6), notlinkfun=function(x) {1 - (1-x)^(1/6)})
}

\arguments{
 \item{kern}{the kernel matrix to transform}
 \item{linkinfo}{linkinfo object that specifies the similarities to transform}
 \item{type}{among c("none", "simple", "extended"), see reference for explanation}
 \item{linkfun}{function applied to similarities that match the linkinfo$link relation}
 \item{notlinkfun}{function applied to similarities that match the linkinfo$notlink relation}
}

\references{
  Pierrick Bruneau and Benoit Otjacques\cr
  \emph{An Interactive, Example-Based, Visual Clustering System}\cr
  Tech Report, CRP Gabriel Lippmann, 2013\cr
  \url{http://hal.archives-ouvertes.fr/hal-00797367}
}

\author{Pierrick Bruneau\cr \email{pbruneau@gmail.com}}

\seealso{
  \code{\link{linkinfo}} 
}
\examples{

data(iris)

krn <- computeKernel(iris[,1:4])
lk <- linkinfo(krn)
lk <- add(lk, 1, 5, "link")
lk <- add(lk, 3, 55, "notlink")
krn <- transformKernel(krn, lk, "extended")
}

