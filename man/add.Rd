\name{add}
\alias{add}
\title{Generic add method}


\description{Generic method to add elements to an object.}

\usage{
add(object, ...)
}

\arguments{
  \item{object}{object where elements should be added}
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
}

