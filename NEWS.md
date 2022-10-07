## tLOH ChangeLog

Changes in version 1.5.1
+ Updated version numbering to match Bioconductor standard
+ Updated DESCRIPTION file
+ Updated README file

Changes in version 0.99.10 (2021-05-05)
+ Replaced DelayedArray in DESCRIPTION imports field with MatrixGenerics,
    updated NAMESPACE
+ Updated R dependency to R (>= 4.0)
+ Updated tLOHDataImport and tLOHCalc functions

Changes in version 0.99.9 (2021-05-4)
+ Added humanGBMsampleAC.Rd file to man/
+ Added RoxygenNote line to DESCRIPTION

Changes in version 0.99.8 (2021-05-4)
+ Added import(DelayedArray) to NAMESPACE, was previously already in 
    DESCRIPTION imports field
+ Changed matrixStats rowRanges to DelayedArray rowRanges in functions.R
+ Removed import of matrixStats rowRanges from NAMESPACE
+ Adjusted formatting in humanGBMsampleAC field of R/data.R documentation

Changes in version 0.99.7 (2021-05-4)
+ Renamed the data/dataExample.rda to data/humanGBMsampleAC.rda, which is a set 
    of pre-processed allele counts from a human glioblastoma sample analyzed 
    with the 10X Genomics Visium Spatial Gene Expression platform

Changes in version 0.99.6 (2021-05-4)
+ Testing dependency fix
+ Fixed typos

Changes in version 0.99.5 (2021-05-3)
+ Added back in vignette
+ Added running examples to documentation
+ Modified required input format to be only a VCF. Future option will include             pre-processed CollapsedVCF object as well

Changes in version 0.99.4 (2021-05-1)
+ Modified required input format from directory of pre-processed .csv files to 
    a VCF or VariantAnnotation CollapsedVCF object containing columns for each 
    cluster
+ Split the plotLOH function into two shorter functions: alleleFrequencyPlot() 
    and aggregateCHRPlot
+ Removed filtering within tLOHCalc. Detailed information on filtering for next
    update.
+ Modified NAMESPACE file
+ Added VariantAnnotation and GenomicRanges to README acknowledgments and
    DESCRIPTION file
+ Removed inst/extdata/sampleData.tar.gz and replaced with Example.vcf
+ Temporarily removed running examples in documentation - will update in 
    next patch
+ Temporarily removed vignette - will update in next patch

Changes in version 0.99.3 (2021-04-27)
+ Updated plotLOH.Rd to not run the code example

Changes in version 0.99.2 (2021-04-27)
+ Updated VignetteEngine to knitr::rmarkdown from knitr::knitr in tLOH.Rmd
+ Updated lines in tLOH.Rmd

Changes in version 0.99.1 (2021-04-27)
+ Removed LazyData: true from DESCRIPTION file
+ Updated Description field in the DESCRIPTION file
+ Consolidated Author and Maintainer fields to Authors@R field in DESCRIPTION 
    file
+ Removed tracking of system files .DS_Store and tLOH.Rproj by adding to 
    .gitignore
+ Removed data/exampleDF.rda and added df.rda to inst/extdata
+ Updated vignette to reflect .rda location and name change
+ Updated vignette text to be more informative and contain link to 
    pre-processing steps
+ Added documentation on contents of inst.extdata/df.rda in the vignette
+ Reduced size of sample images and tar.gz in inst/extdata

Changes in version 0.99.0 (2021-04-22)    
+ Initial Code Edits            
+ Updated README, CITATION, and vignette