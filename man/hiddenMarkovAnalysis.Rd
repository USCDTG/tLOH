\name{hiddenMarkovAnalysis}
\alias{hiddenMarkovAnalysis}
\title{Run Multi-Step HMM Analysis on tLOHCalcUpdate output
}
\description{
Applies the depmixS4 method of HMM Analysis on tLOHCalcUpdate output to obtain segments, noted by the output 'state' column
}
\usage{
hiddenMarkovAnalysis(df, initProbs, trProbs)
}

\arguments{
  \item{df}{
  A dataframe output by tLOHCalcUpdate
}
  \item{initProbs}{
  Dataframe containing 22 rows and two columns, initProb1 and initProb2. Each row represents a chromosome in sequential order, with initProb1 being the probability of state1 and initProb2 being the probability of state2.
}
  \item{trProbs}{
  Matrix of transition start probabilities for HMM
}
}

\value{
Output is a dataframe containing HMM analysis output and tLOHCalcUpdate output summary
}

\author{
Michelle Webb
}

\examples{
data('humanGBMsampleAC')
data('initialStartProbabilities')
df <- tLOHCalcUpdate(humanGBMsampleAC,1.25,1.25,500,500,4)
trProbs <- cbind(c(0.8999,0.1001),c(0.1001,0.8999))
output <- hiddenMarkovAnalysis(df,initialStartProbabilities,trProbs)
}
