\name{linkinfo}
\alias{linkinfo}
\alias{linkinfo,matrix-method}
\title{Linking Information for Kernel Transformation}
\description{
Returns an empty object (ie no link specified) from the linkinfo class.
}
\usage{
\S4method{linkinfo}{matrix}(mat)
}

\arguments{
 \item{mat}{the kernel matrix to which we intend to associate some linking info. Determines the size of the object.}   
}

\references{
  Pierrick Bruneau and Benoit Otjacques\cr
  \emph{An Interactive, Example-Based, Visual Clustering System}\cr
  Tech Report, CRP Gabriel Lippmann, 2013\cr
  \url{http://hal.archives-ouvertes.fr/hal-00797367}
}

\author{Pierrick Bruneau\cr \email{pbruneau@gmail.com}}

\examples{

data(iris)

krn <- computeKernel(iris[,1:4])
lk <- linkinfo(krn)
}

