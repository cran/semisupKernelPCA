\name{check.linkinfo}
\alias{check.linkinfo}
\alias{check,linkinfo-method}
\title{Checks the consistency of the object}


\description{Checks the consistency of the object through temporary application of the transitive/composite closures}

\usage{
\S4method{check}{linkinfo}(object, ...)
}

\arguments{
  \item{object}{a linkinfo class object}
  \item{...}{unused additional arguments, for compliance with generic mechanism}
}

\value{
	TRUE if the relations are consistent, FALSE else.
}

\references{
  Pierrick Bruneau and Benoit Otjacques\cr
  \emph{An Interactive, Example-Based, Visual Clustering System}\cr
  Tech Report, CRP Gabriel Lippmann, 2013\cr
  \url{http://hal.archives-ouvertes.fr/hal-00797367}
}

\author{Pierrick Bruneau\cr \email{pbruneau@gmail.com}}

\seealso{
  \code{\link{transform.linkinfo}} 
}


\examples{
data(iris)

krn <- computeKernel(iris[,1:4])
lk <- linkinfo(krn)
lk <- add(lk, 1, 5, "link")
lk <- add(lk, 3, 55, "notlink")
check(lk)
}

