\name{add.linkinfo}
\alias{add.linkinfo}
\alias{add,linkinfo-method}
\title{Add links to the object}


\description{adds one or several couples to either link or not-link relations recorded in the object.}

\usage{
\S4method{add}{linkinfo}(object, inda, indb, type=c("link", "notlink"), ...)
}

\arguments{
  \item{object}{a linkinfo class object}
  \item{inda}{indexes in \code{1:length(linkinfo)} of the 1st elements in the couples to be added}
  \item{indb}{indexes in \code{1:length(linkinfo)} of the 2nd elements in the couples to be added}
  \item{type}{the relation where the couples should be added, in \code{c("link", "notlink")}}
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
}

