# tLOH v0.99.0

Assessment of evidence for loss of heterozygosity in spatial transcriptomics pre-processed data using Bayes factor calculations.

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
## Information
The functions included in this package allow for calculation of a Bayes factor using allele fraction data at each single nucleotide polymorphic site that is provided. Visium spatial data must be pre-processed using the instructions at https://github.com/USCDTG/spatialAlleleCountPipeline to obtain an individual sample directory of per-cluster .csv's that will be the input for tLOH. The output of this R package is a sample .csv with Bayes factor calculations for all clusters at all sites. The naming convention is 'sample\_tLOH_output.csv'.There is a separate plotting function, plotSpatialLOH(), in the package to visualize allele fraction and Bayes factors per chromosome.

![alt text](https://github.com/USCDTG/tLOH/blob/main/inst/extdata/bayesFactor.png)

M1 and M2 are independent events                

Pr(M1|D) - Probability of Model 1 given data            
Pr(M2|D) - Probability of Model 2 given data             
Pr(M1) - Probability of Model 1                 
Pr(M2) - Probabiliy of Model 2                
              
For this tool, the Pr(M1) is set at 0.5 for a heterozygous event. Alpha and Beta in the beta distribution are set at 10 and 10, respectively.         
 
## To Run
The input directory specified for tLOHCalc must be sample-specific, and contain files with the naming convention 'sample_cluster\[clusterNumber\]\_alleleCounts.csv'. There should be as many input .csv's as there are calculated clusters.

Each input .csv should contain the following columns:

"CHR" - chromosome number            
"CHR\_withLabel" - chomosome containing the 'chr' prefix            
"POS" - SNP position        
"POS\_oneBefore" - one position before
"rsID" - SNP rsID            
"vcfQUAL" - QUAL field from VCF           
"vcf_refNucleotide" - reference nucleotide from the VCF          
"vcf\_altNucleotide" - alternative nucleotide from the VCF       
"INFO" - VCF Information field          
"vcfAF" - VCF allele fraction            
"genotype" - genotype from the VCF, such as '0/1'        
"A" - coverage of A           
"C" - coverage of C       
"G" - coverage of G         
"T" - coverage of T        
"REF" - reference counts           
"ALT" - alternative counts         

```
myDF <- tLOHCalc('inputDirectory','sampleName')
listOfPlots <- plotSpatialLOH(myDF,'sampleName')

## Visualize
listOfPlots[[1]]
```
![alt text](https://github.com/USCDTG/tLOH/blob/main/inst/extdata/Example_alleleFractionPlot.png)              

```
listOfPlots[[2]]
```
![alt text](https://github.com/USCDTG/tLOH/blob/main/inst/extdata/Example_linePlot.png)  

```
listOfPlots[[3]]
```
![alt text](https://github.com/USCDTG/tLOH/blob/main/inst/extdata/Example_columnPlot.png)

## Notes
This version is optimized for human data aligned to GRCh38. The HLA region on chromosome 6 is omitted from this analysis (chr6:28510120-33500500), but will be analyzed in further versions. SNP positions with total allele counts above 2000 were not included, but will be considered in future release.

## Dependencies
- R (>= 3.5.0)
- scales    
- tidyverse
- ggplot2
- data.table
- purrr
- dplyr


## Contact
Michelle Webb  
michelgw@usc.edu

## References
**R:** R Core Team (2019). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.     
**scales:** Hadley Wickham and Dana Seidel (2020). scales: Scale Functions for Visualization. R package version 1.1.1. https://CRAN.R-project.org/package=scales                 
**tidyverse:** Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686
**ggplot2:** H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York 2016.         
**data.table:** Matt Dowle and Arun Srinivasan (2020). data.table: Extension of \`data.frame\`. R package version 1.13.0. https://CRAN.R-project.org/package=data.table          
**purrr:** Lionel Henry and Hadley Wickham (2020). purrr: Functional Programming Tools. R package version 0.3.4. https://CRAN.R-project.org/package=purrr               
**dplyr:** Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2020). dplyr: A Grammar of Data Manipulation. R package version 1.0.0. https://CRAN.R-project.org/package=dplyr               
**Bayes factors** Jeffreys, Harold (1998) [1961]. The Theory of Probability(3rd ed.). 
Oxford, England. p. 432. ISBN 9780191589676.            
        