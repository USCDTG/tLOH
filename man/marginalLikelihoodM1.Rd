\name{marginalLikelihoodM1}
\alias{marginalLikelihoodM1}
\title{
Marginal M1 Calculation
}
\description{
Calculation of the marginal likelihood of Model 1, LOH
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
The recommended values for both Alpha1 and Beta1 is 1.25.
}
\value{
The value returned from marginalLikelihoodM1 is numeric
}
\author{
Michelle Webb
}

\examples{
test <- data.frame(REF=c(10,2,3,4,5,10),TOTAL=c(20,20,20,20,20,20))
apply(test, MARGIN = 1, FUN = marginalLikelihoodM1, a = 1.25, b = 1.25)
}
