\name{tLOHDataImport}
\alias{tLOHDataImport}
\title{
 Import VCF for tLOHCalc
}
\description{
Import a VCF with per-cluster allele count information at heterozygous SNP 
positions for the tLOHCalc calculation function.
}
\arguments{
  \item{vcf}{
  An input VCF file. Spatial transcriptomics clusters make up the sample 
  columns. AD and DP fields are required. Each SNP should be annotated with
  dbSNP rsIDs.
}
}

\value{
Output is a dataframe with required fields for tLOHCalc
}

\author{
Michelle Webb
}

