\name{length.linkinfo}
\alias{length.linkinfo}
\alias{length,linkinfo-method}
\title{returns dimension of parameterized linkinfo object}


\description{returns dimension of parameterized linkinfo object}

\usage{
\S4method{length}{linkinfo}(object)
}

\arguments{
  \item{object}{a linkinfo class object, as created by \code{linkinfo}}
  \item{...}{unused additional arguments, for compliance with generic mechanism}
}

\value{
	dimension of the object
}

\author{Pierrick Bruneau\cr \email{pbruneau@gmail.com}}

\examples{
data(iris)

krn <- computeKernel(iris[,1:4])
lk <- linkinfo(krn)

length(lk)
}

