\name{summarizeRegions1}
\alias{summarizeRegions1}
\title{Step 1 of region summary
}
\description{
Function used by regionFinalize to group segments
}
\usage{
summarizeRegions1(x)
}

\arguments{
  \item{x}{
  A list of dataframes containing calculated values and JMM state determinations
}
}

\value{
Output is a list of dataframes
}

\author{
Michelle Webb
}

\examples{
\donttest{data('humanGBMsampleAC')
data('initialStartProbabilities')
df <- tLOHCalcUpdate(humanGBMsampleAC,1.25,1.25,500,500,4)
dataframeList <- prepareHMMdataframes(df)
trProbs <- cbind(c(0.8999,0.1001),c(0.1001,0.8999))
dataframeList2 <- runHMM_1(dataframeList, initialStartProbabilities, trProbs)
dataframeList3 <- runHMM_2(dataframeList2)
output <- runHMM_3(dataframeList3)
intermediate <- regionAnalysis(dataframeList,output)
finalList1 <- purrr::map(intermediate, 
~dplyr::mutate(.x, state = as.character(state)))
sampleValues <- as.data.frame(purrr::reduce(finalList1,full_join))
sampleData <- summarizeRegions1(finalList1)}
}
