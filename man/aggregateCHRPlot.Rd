\name{aggregateCHRPlot}
\alias{aggregateCHRPlot}
\title{
Visualization of data output from the tLOHCalc function, aggregated per 
chromosome
}
\description{
Output is a plot of the sum of Log10(1/K) values (K is a Bayes factor) per 
chromosome for each cluster. The dotted line at y=3 represents threshold for 
substantial evidence toward Model 2
}
\arguments{
  \item{df}{
  An input dataframe with merged cluster data output by tLOHCalc
}
  \item{sample}{
  Name of sample for plot title
  }
}

\value{
Output is a plot where the y axis is sum of Log10(1/K) values (K is a 
Bayes factor) per chromosome and the x axis is chromosome
}

\author{
Michelle Webb
}