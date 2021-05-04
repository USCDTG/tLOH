# tLOH
### v0.99.6
Assessment of evidence for loss of heterozygosity in spatial transcriptomics pre-processed data using Bayes factor calculations.

## About
Loss of heterozygosity (LOH) refers to a genomic event where a chromosomal region on one allele is lost. Evidence for this event can be calculated by examining allele frequency (AF), or the ratio of alleles, at known heterozygous positions. AF at a known heterozygous (single nucleotide polymorphism) SNP should be around 0.50, where half of the sequencing reads at that position align to the reference allele and half to the alternative. As AF approaches 0 or 1, this represents an uneven distribution of counts to either reference or alternative. Bayesian statistics can be used to assess the ratio of likelihood of a heterozygous or LOH event at each SNP.         

The functions included in this package allow for calculation of a Bayes factor at each SNP provided. A 10X Genomics Visium spatial transcriptomics BAM must be pre-processed to obtain a VCF with sample columns for each cluster (graph or k-means). A detailed pipeline will be provided in future release, though steps are as follows: Separate a spatial transcriptomics BAM into per-cluster BAM files, filter reads for heterozygous or likely heterozygous SNP positions, obtain allele counts, and store data in VCF format with columns for each cluster. Required fields are DP (read depth) and AD (counts for reference and alternative alleles). Data import format is a VCF. Output from this R package includes a dataframe with Bayes factor calculations for all clusters at all sites. There are two separate plotting function in the package to visualize allele fraction and aggregated Bayes factors per chromosome.

![alt text](https://github.com/USCDTG/tLOH/blob/main/inst/extdata/bayesFactor.png)

M1 and M2 are independent events                

Pr(M1|D) - Probability of Model 1 given data            
Pr(M2|D) - Probability of Model 2 given data             
Pr(M1) - Probability of Model 1                 
Pr(M2) - Probabiliy of Model 2                
              
For this tool, the Pr(M1) is set at 0.5 for a heterozygous event. Alpha and Beta in the beta distribution are set at 10 and 10, respectively.         

## Installation
After downloading the R package .zip file from GitHub, convert to tar.gz and run the following commands:

```
R
> library('devtools')
> devtools::install_local('/path/to/tLOH.tar.gz')
> library('tLOH')
```
or

```
R
> install.packages('/path/to/tLOH.tar.gz')
> library('tLOH')
```

## Usage
An example VCF is available in the /inst/extdata directory as Example.vcf. tLOH can read in a VCF with the function tLOHDataImport().

####  Data Import
```
myDF <- tLOHDataImport('/full/path/to/Example/vcf')    
```

#### tLOH Calculation of Bayes factors
```
tLOHOutput <- tLOHCalc(myDF)
```

#### Visualization
```
alleleFrequencyPlot(tLOHOutput,'SampleNameForPlotTitle')
```
![alt text](https://github.com/USCDTG/tLOH/blob/main/inst/extdata/Example_alleleFrequencyPlot.png)              

```
aggregateCHRPlot(tLOHOutput,'SampleNameForPlotTitle')
```
![alt text](https://github.com/USCDTG/tLOH/blob/main/inst/extdata/Example_aggregateCHRPlot.png)

Dotted line represents stringent threshold for substantial evidence toward Model 2.

## Notes
This version is optimized for human data aligned to GRCh38. The HLA region on chromosome 6 is omitted from this analysis (chr6:28510120-33500500), but will be analyzed in further versions. SNP positions with total allele counts above 2000 were not included, but will be considered in future release. Additional visualizations and SNP annotation/filtering guidelines and are planned for future patch.

## Prerequisites
- R (>= 3.5.0)
- scales    
- tidyverse
- ggplot2
- data.table
- purrr
- dplyr
- VariantAnnotation
- GenomicRanges


## Contact
Michelle Webb  
michelgw@usc.edu

## Acknowledgments
**10X Visium Spatial Gene Expression** https://www.10xgenomics.com/products/spatial-gene-expression              
**R:** R Core Team (2019). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.     
**scales:** Hadley Wickham and Dana Seidel (2020). scales: Scale Functions for Visualization. R package version 1.1.1. https://CRAN.R-project.org/package=scales                 
**tidyverse:** Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686
**ggplot2:** H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York 2016.         
**data.table:** Matt Dowle and Arun Srinivasan (2020). data.table: Extension of \`data.frame\`. R package version 1.13.0. https://CRAN.R-project.org/package=data.table          
**purrr:** Lionel Henry and Hadley Wickham (2020). purrr: Functional Programming Tools. R package version 0.3.4. https://CRAN.R-project.org/package=purrr               
**dplyr:** Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2020). dplyr: A Grammar of Data Manipulation. R package version 1.0.0. https://CRAN.R-project.org/package=dplyr           
**VariantAnnotation:** Obenchain V, Lawrence M, Carey V, Gogarten S, Shannon P, Morgan M (2014).
“VariantAnnotation: a Bioconductor package for exploration and annotation of
genetic variants.” _Bioinformatics_, *30*(14), 2076-2078. doi:
10.1093/bioinformatics/btu168 (URL:
https://doi.org/10.1093/bioinformatics/btu168).   
**GenomicRanges:** Lawrence M, Huber W, Pag\`es H, Aboyoun P, Carlson M, et al. (2013) Software
  for Computing and Annotating Genomic Ranges. PLoS Comput Biol 9(8): e1003118.
  doi:10.1371/journal.pcbi.1003118             
**Bayes factors** Jeffreys, Harold (1998) [1961]. The Theory of Probability(3rd ed.). 
Oxford, England. p. 432. ISBN 9780191589676.            
