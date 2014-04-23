\docType{methods}
\name{gridtweak}
\alias{gridtweak}
\alias{gridtweak,ANY-method}
\title{gridtweak: Open a Shiny Application for the purpose of tweaking grid plots}

\description{
  This opens a shiny visualization application in the browser.
}
\usage{
gridtweak(...)
}
\arguments{
  \item{...}{additional arguments passed to methods; currently unused.}
}
\value{
  Any grid based plot.  For example: a plot produced with lattice, ggplot2 or 
  biobase libraries.
}
\seealso{
  \url{http://bioconductor.org/packages/2.13/bioc/html/interactiveDisplay.html}
}
\examples{

if(interactive()) {

## Send a grid based plot to a browser as a Javascript interactive SVG

gridtweak()

}

}
\keyword{methods}

