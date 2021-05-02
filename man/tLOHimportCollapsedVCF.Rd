\name{tLOHimportCollapsedVCF}
\alias{tLOHimportCollapsedVCF}
\title{
 Import VariantAnnotation CollapsedVCF for tLOHCalc
}
\description{
Import a CollapsedVCF with per-cluster allele count information at heterozygous 
SNP positions for the tLOHCalc calculation function.
}
\arguments{
  \item{variantAnnotationVCF}{
  An input CollapsedVCF object.
}
}

\value{
Output from variantAnnotationVCF is a dataframe with required fields for 
tLOHCalc
}

\author{
Michelle Webb
}
