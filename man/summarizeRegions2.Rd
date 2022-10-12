\name{summarizeRegions2}
\alias{summarizeRegions2}
\title{Step 2 of region summary
}
\description{
Function used by regionFinalize to identify segment start and end positions
}
\usage{
summarizeRegions2(finalTable)
}

\arguments{
  \item{finalTable}{
  A dataframe containing metrics from the Bayes Factor and HMM analysis
}
}


\value{
Output is a dataframe
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
finalList1 <- purrr::map(intermediate, 
~dplyr::mutate(.x, state = as.character(state)))
sampleValues <- as.data.frame(purrr::reduce(finalList1,full_join))
sampleData <- summarizeRegions1(finalList1)
sampleData$lengthOfInterval <- sampleData$intervalEnd - sampleData$intervalStart
sampleDF <- summarizeRegions2(sampleValues)}
}
