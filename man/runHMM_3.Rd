\name{runHMM_3}
\alias{runHMM_3}
\title{Step 3 of HMM process
}
\description{
Applies the depmixS4 method posterior on normalized K values.
}
\usage{
runHMM_3(dataframeList)
}

\arguments{
  \item{dataframeList}{
  List of depmixS4 depmix.fitted class output generated from the runHMM_2 step
}
}

\value{
Output is a list of depmixS4 posterior state classifications for each input dataframe
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
}
