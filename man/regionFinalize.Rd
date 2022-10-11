\name{regionFinalize}
\alias{regionFinalize}
\title{Summary of HMM Regions
}
\description{
Final metrics and summary for regions
}
\usage{
regionFinalize(finalList1)
}

\arguments{
  \item{finalList1}{
  List of dataframes output by the regionAnalysis function
}
}

\value{
Output is a table containing all calculations from the bayes factor and HMM analysis
}

\author{
Michelle Webb
}

\examples{
\dontrun{
data('humanGBMsampleAC')
data('initialStartProbabilities')
df <- tLOHCalcUpdate(humanGBMsampleAC,1.25,1.25,500,500,4)
dataframeList <- prepareHMMdataframes(df)
trProbs <- cbind(c(0.8999,0.1001),c(0.1001,0.8999))
dataframeList2 <- runHMM_1(dataframeList, initialStartProbabilities, trProbs)
dataframeList3 <- runHMM_2(dataframeList2)
output <- runHMM_3(dataframeList3)
intermediate <- regionAnalysis(dataframeList,output)
final <- regionFinalize(intermediate)
}}
