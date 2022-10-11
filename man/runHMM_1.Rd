\name{runHMM_1}
\alias{runHMM_1}
\title{Step 1 of HMM process
}
\description{
Applies the depmixS4 method depmix on normalized K values.
}
\usage{
runHMM_1(dataframeList, initProbs, trProbs)
}

\arguments{
  \item{dataframeList}{
  List of dataframes separated by cluster and chromosome from the prepareHMMdataframes function.
}
  \item{initProbs}{
  Dataframe containing 22 rows and two columns, initProb1 and initProb2. Each row represents a chromosome in sequential order, with initProb1 being the probability of state1 and initProb2 being the probability of state2.
}
  \item{trProbs}{
  Matrix of transition start probabilities for HMM
}
}

\value{
Output is a list of depmixS4 depmix class objects for each input dataframe
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
output <- runHMM_1(dataframeList, initialStartProbabilities, trProbs)
}
