\name{alleleFrequencyPlot}
\alias{alleleFrequencyPlot}
\title{
Visualization of data output from the tLOHCalc function
}
\description{
Creates a plot with panels for each cluster. The x-axis is chromosome, y-axis is
allele frequency. Point color is Log10(1/K) where K is a Bayes factor
}
\arguments{
  \item{df}{
  An input dataframe with merged cluster data
}
  \item{sample}{
  Name of sample for plot title
  }
}

\value{
Output is a plot of allele frequency for each cluster. Can be assigned to 
object and visualized individually. For each panel, the y axis has a min of 0 
and max of 1
}

\author{
Michelle Webb
}

\examples{
data('dataExample')
df <- tLOHCalc(dataExample)
alleleFrequencyPlot(df,"Example")
}