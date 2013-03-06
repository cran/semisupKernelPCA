\name{del}
\alias{del}
\title{Deletes from an object}


\description{Generic method to remove elements from an object.}

\usage{
del(object, ...)
}

\arguments{
  \item{object}{object from which elements should be removed}
  \item{...}{to allow further arguments by implementations}
}

\value{
	updated object.
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

