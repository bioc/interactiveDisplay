\docType{methods}
\name{gridsvgjs}
\alias{gridsvgjs}
\alias{gridsvgjs,ANY-method}
\title{gridsvgjs: open shiny application for grid plot}

\description{
  This opens a shiny visualation application in the browser
  based on the submitted plot.
}
\usage{
gridsvgjs(object)
}
\arguments{
  \item{object}{a grid based plot to display in browser}
}
\value{
  Any grid based plot.  For example: a plot produced with lattice, ggplot2 or biobase libraries.
}
\seealso{
  \url{http://bioconductor.org/packages/2.11/bioc/html/interactiveDisplay.html}
}
\examples{

if(interactive()) {

## Send grid based plot to browser as a Javascript interactive SVG

library(ggplot2)
data(mtcars)
qp <- qplot(mpg, data=mtcars, geom="density", fill=factor(cyl), alpha=I(.4))
gridsvgjs(qp)

}

}
\keyword{methods}

