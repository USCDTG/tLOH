\name{runHMM_2}
\alias{runHMM_2}
\title{Step 2 of HMM process
}
\description{
Applies the depmixS4 method fit on .
}
\usage{
runHMM_2(dataframeList)
}

\arguments{
  \item{dataframeList}{
  List of depmixS4 depmix class objects generated from the runHMM_1 step
}
}

\value{
Output is a list of depmixS4 depmix.fitted class output for each input dataframe
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
output <- runHMM_2(dataframeList2)
}
