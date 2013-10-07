\name{expr}
\alias{expr}
\title{An Example ExpressionSet object}
\description{
The expression data are real but anonymized. The data are from an experiment that used Affymetrix U95v2 chips. The data were processed by dChip and then exported to R for analysis.
The data illustrate ExpressionSet-class, with assayData containing the required matrix element exprs and an additional matrix se.exprs. se.exprs has the same dimensions as exprs.
The phenoData and standard error estimates (se.exprs) are made up. The information in the "description" slot is fake.
}
\details{
The data for 26 cases, labeled A to Z and 500 genes. Each case has three covariates: sex (male/female); type (case/control); and score (testing score).
}
\examples{
data(expr)
}
\keyword{datasets}
