\docType{methods}
\name{altgr}
\alias{altgr}
\alias{altgr,ANY-method}
\title{altgr: Open a Shiny Application for ...}

\description{
  This opens a shiny visualization application in the browser
  based on ...
}
\usage{
altgr(object, ...)
}
\arguments{
  \item{object}{data object to display}
  \item{...}{additional arguments passed to methods; currently unused.}
}
\value{
  Any ...
}
\seealso{
  \url{http://bioconductor.org/packages/2.13/bioc/html/interactiveDisplay.html}
}
\examples{

if(interactive()) {

## Open an browser application for the purpose of manually biclustering an  
## ExpressionSet object and obtaining a GO summary for a specific bicluster.

data(mmgr)
altgr(mmgr)

}

}
\keyword{methods}

