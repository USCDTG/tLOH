\name{marginalM2CalcBHET}
\alias{marginalM2CalcBHET}

\title{
Calculation of marginal M2 het
}
\description{
Calculation of marginal M2 het
}
\usage{
marginalM2CalcBHET(x, a, b)
}

\arguments{
  \item{x}{
Number of counts for the reference allele
}
  \item{a}{
Alpha value
}
  \item{b}{
Beta value
}
}

\value{
The value returned is an integer
}

\author{
Michelle Webb
}

\examples{
save <- data.frame(REF=c(10,2,3,4,5,10),TOTAL=c(20,20,20,20,20,20))
apply(save, MARGIN = 1, FUN = marginalM2CalcBHET, a = 10,b = 10)
}

\keyword{ marginalM2Calc }

