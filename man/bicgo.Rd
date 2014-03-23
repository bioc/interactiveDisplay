\docType{methods}
\name{bicgo}
\alias{bicgo}
\alias{bicgo,ANY-method}
\title{bicgo: Open a Shiny Application for manual/interactive biclustering and GO exploration}

\description{
  This opens a shiny visualization application in the browser
  based on the submitted ExpressionSet object.
}
\usage{
bicgo(object, ...)
}
\arguments{
  \item{object}{data object to display}
  \item{...}{additional arguments passed to methods; currently unused.}
}
\value{
  Any ExpressionSet object.
}
\seealso{
  \url{http://bioconductor.org/packages/2.13/bioc/html/interactiveDisplay.html}
}
\examples{

if(interactive()) {

## Open an browser application for the purpose of manually biclustering an  
## ExpressionSet object and obtaining a GO summary for a specific bicluster.

data(expr)
bicgo(expr)

}

}
\keyword{methods}

