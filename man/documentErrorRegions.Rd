\name{documentErrorRegions}
\alias{documentErrorRegions}
\title{Identification of Errors
}
\description{
Generates dataframes equal to the length of original data containing NA if errors were found during the run HMM process. Important for final overview of results
}
\usage{
documentErrorRegions(a,b)
}

\arguments{
  \item{a}{
  List of dataframes from prepareHMMdataframes
}
  \item{b}{
  List of HMM state determination dataframes from the run_HMM series of functions
}
}

\value{
Output is a list of dataframes to be used by the regionAnalysis function
}

\author{
Michelle Webb
}

\examples{
estimatedStates <- data.frame(state = c(1,1,1), S1 = c(1,1,1), S2 = c(1,1,1))
sampleDataFrame <- data.frame(data = c(1,1,1))
list1 <- list(sampleDataFrame,sampleDataFrame)
list2 <- list(estimatedStates,estimatedStates)
documentErrorRegions(list1,list2)
}
