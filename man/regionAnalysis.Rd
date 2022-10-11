\name{regionAnalysis}
\alias{regionAnalysis}
\title{Summary of HMM Regions
}
\description{
Generates summary metrics for HMM regions
}
\usage{
regionAnalysis(originalDF,dataframeList)
}

\arguments{
  \item{originalDF}{
  Original imported dataframe from tLOHCalcUpdated
}
  \item{dataframeList}{
  List of HMM state determination dataframes from the run_HMM series of functions
}
}

\value{
Output is a dataframe containing region metrics and data for each HMM segment
}

\author{
Michelle Webb
}

\examples{
data('humanGBMsampleAC')
data('initialStartProbabilities')
df <- tLOHCalcUpdate(humanGBMsampleAC,1.25,1.25,500,500,4)
dataframeList <- prepareHMMdataframes(df)
trProbs <- cbind(c(0.8999,0.1001),c(0.1001,0.8999))
dataframeList2 <- runHMM_1(dataframeList, initialStartProbabilities, trProbs)
dataframeList3 <- runHMM_2(dataframeList2)
output <- runHMM_3(dataframeList3)
final <- regionAnalysis(dataframeList,output)
}
