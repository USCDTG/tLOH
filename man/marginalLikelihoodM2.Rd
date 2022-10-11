\name{marginalLikelihoodM2}
\alias{marginalLikelihoodM2}
\title{
Marginal M2 Calculation
}
\description{
Calculation of the marginal likelihood of Model 2, HET
}

\arguments{
  \item{x}{
Dataframe output by tLOHDataImport
}
  \item{a}{
Alpha value
}
  \item{b}{
Beta value
}

}
\details{
The reference and total counts should come from a .csv output by the spatial LOH pre-processing pipeline.
The recommended values for both Alpha2 and Beta2 is 500.
}
\value{
The value returned from marginalLikelihoodM1 is numeric
}
\author{
Michelle Webb
}

\examples{
test <- data.frame(REF=c(10,2,3,4,5,10),TOTAL=c(20,20,20,20,20,20))
apply(test, MARGIN = 1, FUN = marginalLikelihoodM1, a = 500, b = 500)
}
