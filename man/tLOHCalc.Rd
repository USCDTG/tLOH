\name{tLOHCalc}
\alias{tLOHCalc}
\title{
Assement of evidence for spatial LOH
}
\description{
Takes the output from the tLOH pre-processing pipeline and calculates
Bayes factors for allele fractions at each SNP position.
}
\usage{
tLOHCalc(input, sample)
}

\arguments{
  \item{input}{
Input directory
}
  \item{sample}{
Sample Name
}
}

\value{
Output is a dataframe
}

\author{
Michelle Webb
}

\examples{
\dontrun{tLOHCalc(input, sample)}
}
