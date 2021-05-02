\name{tLOHCalc}
\alias{tLOHCalc}
\title{
Assement of evidence for LOH in clusters from spatial transcriptomics data
}
\description{
Calculates Bayes factors for allele fractions at each SNP position. Uses 
dataframe output by either tLOHDataImport or tLOHimportCollapsedVCF
}
\usage{
tLOHCalc(forCalcDF)
}

\arguments{
  \item{forCalcDF}{
  input dataframe that can be generated through VCF import or VariantAnnotation 
  CollapsedVCF import
}
}

\value{
Output is a dataframe with values that can be visualized with 
alleleFrequencyPlot() or aggregateCHRPlot()
}

\author{
Michelle Webb
}

