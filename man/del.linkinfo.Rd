\name{del.linkinfo}
\alias{del.linkinfo}
\alias{del,linkinfo-method}
\title{Deletes links from the object}


\description{removes one or several couples from both relations recorded in the object.}

\usage{
\S4method{del}{linkinfo}(object, inda, indb, ...)
}

\arguments{
  \item{object}{a linkinfo class object}
  \item{inda}{indexes in \code{1:length(linkinfo)} of the 1st elements in the couples to be removed}
  \item{indb}{indexes in \code{1:length(linkinfo)} of the 2nd elements in the couples to be removed}
  \item{...}{unused additional arguments, for compliance with generic mechanism}
}

\value{
	updated linkinfo object.
}

\author{Pierrick Bruneau\cr \email{pbruneau@gmail.com}}

\examples{
data(iris)

krn <- computeKernel(iris[,1:4])
lk <- linkinfo(krn)
lk <- add(lk, 1, 5, "link")
lk <- add(lk, 3, 55, "notlink")
lk <- del(lk, 1, 5)
}

