\name{linkinfo-class}
\docType{class}
\alias{linkinfo-class}

\title{Class "linkinfo"}
\description{Linking Information Class for Kernel Transformation}
\section{Objects from the Class}{
Objects can be created by calls to the function \code{linkinfo}. 
}

\section{Slots}{
  \describe{
    \item{\code{link}:}{Binary matrix containing the link-type contraints, \cr
    ie items of the respective kernel matrix should be linked if TRUE.}
    \item{\code{notlink}:}{Binary matrix containing the not-link-type contraints}
  }
}

\section{Methods}{
  \describe{
    \item{length}{\code{signature(object = "linkinfo")}: returns the
      number of elements associated to the object}
    \item{add}{\code{signature(object = "linkinfo")}: adds link or not-link instances in the relations}
    \item{del}{\code{signature(object = "linkinfo")}: removes instances from the relations}
	\item{check}{\code{signature(object = "linkinfo")}: checks the consistency of the relations}
	\item{transform}{\code{signature(object = "linkinfo")}: compute the transitive/composite closures of the relations (see reference)}
  }
}

\author{Pierrick Bruneau\cr \email{pbruneau@gmail.com}}


\references{
  Pierrick Bruneau and Benoit Otjacques\cr
  \emph{An Interactive, Example-Based, Visual Clustering System}\cr
  Tech Report, CRP Gabriel Lippmann, 2013\cr
  \url{http://hal.archives-ouvertes.fr/hal-00797367}
}

\seealso{
  \code{\link{linkinfo}} 
}
\examples{

data(iris)

krn <- computeKernel(iris[,1:4])
lk <- linkinfo(krn)

length(lk)
lk <- add(lk, 1, 5, "link")
check(lk)
}
\keyword{classes}

