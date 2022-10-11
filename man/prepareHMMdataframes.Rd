\name{prepareHMMdataframes}
\alias{prepareHMMdataframes}
\title{
Prepare dataframes for HMM analysis
}
\description{
Split output from tLOHCalc or tLOHCalcUpdate into a list of cluster and chromosome separated dataframes. Applies an ordered quantile normalization on the bayes factor K values in each dataset.
}
\usage{
prepareHMMdataframes(importedData)
}

\arguments{
  \item{importedData}{
  Input dataframe generated from the tLOHCalc or tLOHCalcUpdate function
}
}

\value{
Output is a list of dataframes separated by chromosome and cluster
}

\author{
Michelle Webb
}

\examples{
data('humanGBMsampleAC')
df <- tLOHCalcUpdate(humanGBMsampleAC,1.25,1.25,500,500,4)
output <- prepareHMMdataframes(df)
}
