\docType{methods}
\name{gridsvgjs}
\alias{gridsvgjs}
\alias{gridsvgjs,ANY-method}
\title{gridsvgjs: Open a Shiny Application for a Grid Plot}

\description{
  This opens a shiny visualization application in the browser
  based on the submitted plot.
}
\usage{
gridsvgjs(object, ...)
}
\arguments{
  \item{object}{data object to display}
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

library(ggplot2)
data(mtcars)
qp <- qplot(mpg, data=mtcars, geom="density", fill=factor(cyl), alpha=I(.4))
gridsvgjs(qp)

}

}
\keyword{methods}

