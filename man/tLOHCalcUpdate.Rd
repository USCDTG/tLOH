\name{tLOHCalcUpdate}
\alias{tLOHCalcUpdate}
\title{
Assement of evidence for LOH in clusters from spatial transcriptomics allele count data
}
\description{
Calculates Bayes factors for allele fractions at each SNP position. Uses 
dataframe output by tLOHDataImport.
}
\usage{
tLOHCalcUpdate(forCalcDF, alpha1, beta1,alpha2, beta2, countThreshold)
}

\arguments{
  \item{forCalcDF}{
  Input dataframe generated from the tLOHDataImport function
}
  \item{alpha1}{
  Model 1 alpha value
}
  \item{beta1}{
  Model 1 beta value
}
  \item{alpha2}{
  Model 2 alpha value
}
  \item{beta2}{
  Model 2 beta value
}
  \item{countThreshold}{
  Threshold for minimum number of read counts
}
}

\value{
Output is a dataframe with Bayes Factor values
}

\author{
Michelle Webb
}

\examples{
data('humanGBMsampleAC')
df <- tLOHCalcUpdate(humanGBMsampleAC, 1.25,1.25,500,500,4)
head(df)
}
