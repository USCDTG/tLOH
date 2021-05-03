#' Imported dataset processed with tLOHImportData.
#' @name dataExample
#' @description A dataset containing the allele count information for 9 spatial
#'transcriptomics clusters
#'
#' @format A data frame with 34601 rows and 7 variables:
#' \describe{
#'   \item{rsID}{dbSNP rs identifier}
#'   \item{CLUSTER}{cluster number}
#'   \item{TOTAL}{total number of counts}
#'   \item{REF}{counts for the reference allele}
#'   \item{ALT}{counts for the alternative allele}
#'   \item{CHR}{chromosome number}
#'   \item{POS}{genomic position}
#' }
#' @docType data
#' @keywords datasets
#' @source  Lab data repository
#' @examples 
#' data("dataExample")
"dataExample"