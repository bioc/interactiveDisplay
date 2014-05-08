\docType{methods}
\name{simplenet}
\alias{simplenet}
\alias{simplenet,ANY-method}
\title{simplenet: Open a Shiny Application for ...}

\description{
  This opens a shiny visualization application in the browser
  based on ...
}
\usage{
simplenet(object, ...)
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

simplenet(mtcars)

}

}
\keyword{methods}

