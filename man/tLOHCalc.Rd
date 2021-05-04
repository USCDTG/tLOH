\name{tLOHCalc}
\alias{tLOHCalc}
\title{
Assement of evidence for LOH in clusters from spatial transcriptomics data
}
\description{
Calculates Bayes factors for allele fractions at each SNP position. Uses 
dataframe output by tLOHDataImport
}
\usage{
tLOHCalc(forCalcDF)
}

\arguments{
  \item{forCalcDF}{
  Input dataframe generated from the tLOHDataImport function
}
}

\value{
Output is a dataframe with values that can be visualized with 
alleleFrequencyPlot() or aggregateCHRPlot()
}

\author{
Michelle Webb
}

\examples{
data('dataExample')
df <- tLOHCalc(dataExample)
head(df)
}
