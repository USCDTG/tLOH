\name{splitByChromosome}
\alias{splitByChromosome}
\title{
Split dataframe into individual chromosome dataframes
}
\description{
Creates individual chromosome dataframes
}
\usage{
splitByChromosome(listOfDataframes,numberOfDataframes)
}

\arguments{
  \item{listOfDataframes}{
  Input dataframe generated from the tLOHDataImport function
}
  \item{numberOfDataframes}{
  Number of dataframes in list
}
}

\value{
Output is a list of dataframe separated by chromosome
}

\author{
Michelle Webb
}

\examples{
data('humanGBMsampleAC')
df <- tLOHCalcUpdate(humanGBMsampleAC,1.25,1.25,500,500,4)
output <- splitByChromosome(list(df),1)
}
