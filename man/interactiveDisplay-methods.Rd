\docType{methods}
\name{display}
\alias{display}
\alias{display,ANY-method}
\alias{display,missing-method}
\alias{display,ExpressionSet-method}
\alias{display,GRanges-method}
\alias{display,GRangesList-method}
\alias{display,RangedSummarizedExperiment-method}
\alias{display,MRexperiment-method}


\title{display: Open a Shiny application for a Bioconductor object}
\description{
  This opens a shiny visualization application in the browser
  based on the submitted object.
}
\usage{
display(object, ...)
}
\arguments{
  \item{object}{data object to display}
  \item{...}{additional arguments passed to methods; currently unused.}
}
\value{
  Usually some variation of the initial input object, but it may be
  altered by the display widget (subset for example).
}
\seealso{
  \url{http://bioconductor.org/packages/2.13/bioc/html/interactiveDisplay.html}
}
\author{Shawn Balcome and Marc Carlson}
\examples{

if(interactive()) {

## draw a RangedSummarizedExperiment object
data(se)
display(se)

## draw a GRanges object
data(mmgr)
display(mmgr)

## some display methods allow subsetting.
## To take advantage, just use an assignment operator like this:
mmgr2 <- display(mmgr)


## draw a GRangesList object
data(mmgrl)
display(mmgrl)

## draw an ExpressionSet object
data(expr)
display(expr)

## draw an MRexperiment object (placeholder!!!)
data(mr)
display(mr)

}

}

\keyword{methods}

