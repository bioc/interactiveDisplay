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

\title{display: Open a Shiny application for a Bioconductor object}
\description{
  This opens a shiny visualation application in the browser
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

\keyword{methods}

