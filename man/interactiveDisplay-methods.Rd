\docType{methods}
\name{display}
\alias{display}
\alias{display,ANY-method}
\alias{display,missing-method}
\alias{display,ExpressionSet-method}
\alias{display,GRanges-method}
\alias{display,GRangesList-method}
\alias{display,SummarizedExperiment-method}
\alias{display,AnnotationHub-method}

\title{display: open shiny application for object}
\arguments{
  \item{object}{data object to display}
}
\value{
  Usually some variation of the initial input object, but it may be
  altered by the display widget (subset for example).
}
\description{
  This opens a shiny visualation application in the browser
  based on the submitted object.
}
\examples{

if(interactive()) {

## draw an SummarizedExperiment object
data(se)
display(se)

## draw a GRanges object
data(mmgr)
display(mmgr)

## draw a GRangesList object
data(mmgrl)
display(mmgrl)

## draw an ExpressionSet object
data(expr)
display(expr)

}

}
\seealso{
  \url{http://bioconductor.org/packages/2.11/bioc/html/interactiveDisplay.html}
}
\keyword{methods}

