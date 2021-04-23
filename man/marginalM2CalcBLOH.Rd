\name{marginalM2CalcBLOH}
\alias{marginalM2CalcBLOH}

\title{
Marginal M2 Calculation
}
\description{
Calculation of the marginal for Model 2
}
\usage{
marginalM2CalcBLOH(x, a, b)
}

\arguments{
  \item{x}{
Counts for the reference allele
}
  \item{a}{
Alpha value
}
  \item{b}{
Beta value
}
}

\author{
Michelle Webb
}

\examples{
test <- data.frame(REF=c(10,2,3,4,5,10),TOTAL=c(20,20,20,20,20,20))
apply(test, MARGIN = 1, FUN = marginalM2CalcBLOH, a = 10,b = 10)

}
