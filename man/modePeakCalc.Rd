\name{modePeakCalc}
\alias{modePeakCalc}
\title{
Calculation of mode peak
}
\description{
This function takes a set of numbers and outputs a mode peak value. To be used 
in a larger function that will be updated.
}

\arguments{
  \item{x}{
List of allele fraction values
}
}
\details{
List of values should be the allele fractions of SNPs with the top 25 percent 
of counts in a region. If only one value is input, that value is returned.
}
\value{
The value returned is numeric
}
\author{
Michelle Webb
}

\examples{
test <- c(1,2,3,4,5)
modePeakCalc(test)
}
