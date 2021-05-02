\name{marginalM1Calc}
\alias{marginalM1Calc}
\title{
Calculate marginal of Model 1
}
\description{
This function takes the number of counts for a reference allele as x,
    and the number of total allele counts as y.
}

\arguments{
  \item{x}{
Number of counts for the reference allele.
}
  \item{y}{
Number of counts total at this SNP position.
}
}
\details{
The reference and total counts should come from a .csv output by the spatial LOH pre-processing pipeline.
}
\value{
The value returned from marginalM1Calc is numeric
}
\author{
Michelle Webb
}

\examples{
marginalM1Calc(10, 0.5)
}
