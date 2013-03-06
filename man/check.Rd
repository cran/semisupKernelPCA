\name{check}
\alias{check}
\title{Generic method to check an object}


\description{Checks an object, returning TRUE or FALSE.}

\usage{
check(object, ...)
}

\arguments{
  \item{object}{object to be checked}
  \item{...}{to allow further arguments by implementations}
}

\value{
	TRUE if checking procedure successful, FALSE else.
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
check(lk)
}

