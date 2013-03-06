\name{transform.linkinfo}
\alias{transform.linkinfo}
\alias{transform,linkinfo-method}
\title{Transform the linkinfo object to their closures}


\description{The parameterized linkinfo object is augmented by its transitive/composite closures}

\usage{
\S4method{transform}{linkinfo}(`_data`, ...)
}

\arguments{
  \item{_data}{a linkinfo class object}
  \item{...}{unused additional arguments, for compliance with generic mechanism}
}

\value{
	Augmented linkinfo object
}

\references{
  Pierrick Bruneau and Benoit Otjacques\cr
  \emph{An Interactive, Example-Based, Visual Clustering System}\cr
  Tech Report, CRP Gabriel Lippmann, 2013\cr
  \url{http://hal.archives-ouvertes.fr/hal-00797367}
}

\author{Pierrick Bruneau\cr \email{pbruneau@gmail.com}}

\seealso{
  \code{\link{check.linkinfo}} 
}


\examples{
data(iris)

krn <- computeKernel(iris[,1:4])
lk <- linkinfo(krn)
lk <- add(lk, 1, 5, "link")
lk <- add(lk, 3, 55, "notlink")
lk <- transform(lk)
}

