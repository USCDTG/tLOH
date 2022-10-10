marginalM1Calc <- function(x,y){
    successes <- x["REF"]
    trials <- x["TOTAL"]
    dbinom(x = successes, size = trials, prob = y)
}

marginalM2CalcBLOH <- function(x,a,b){
    successes <- x["REF"]
    trials <- x["TOTAL"]
    p <- function(theta){dbinom(x = successes, size = trials, prob = theta) *
            dbeta(x = abs(0.5 - theta), shape1 = a, shape2 = b)}
    integrate(f = p, lower = 0, upper = 1)[[1]]
}


marginalM2CalcBHET <- function(x,a,b){
    successes <- x["REF"]
    trials <- x["TOTAL"]
    p <- function(theta){dbinom(x = successes, size = trials, prob = theta) *
            dbeta(x = theta, shape1 = a, shape2 = b)}
    integrate(f = p, lower = 0, upper = 1)[[1]]
}

removeOutlierFromCalc <- function(dataframe, cols, rows, newValue = NA) {
    if (any(rows)) {
        data.table::set(dataframe, rows, cols, newValue)
    }
}

tLOHDataImport <- function(vcf){
    inputVCF <- readVcf(vcf)
    sampleClusters <- seq(1,dim(inputVCF)[[2]],by=1)
    intermediateDFList <- lapply(sampleClusters, function(x){
        counts <- t(data.frame(geno(inputVCF)$AD[,x]))
        depth <- suppressWarnings(data.frame(data.table::melt(as.data.table(
            as.data.frame(geno(inputVCF)$DP), keep.rownames = "rsID"), 
            variable.name = "CLUSTER", value.name = "TOTAL")))
        names(depth) <- c("rsID", "CLUSTER", "TOTAL")
        intermediate <- data.frame(rsID = rownames(counts), 
                                   REF = counts[,1], 
                                   ALT = counts[,2],
                                   CLUSTER = x)
        intermediate <- merge(depth,intermediate,by=c("rsID", "CLUSTER"))
        intermediate <- intermediate[intermediate$TOTAL > 0,]
        rownames(intermediate) <- NULL
        return(intermediate)})
    preMergeData <- reduce(intermediateDFList,full_join)
    gr <- rowRanges(inputVCF)
    positionList <- data.frame(CHR = seqnames(gr),
                               POS = start(gr), 
                               rsID = names(gr))
    importedDF <- merge(preMergeData, positionList, by = c("rsID"))
    importedDF <- removeOutlierFromCalc(importedDF,"TOTAL",
                                        which(importedDF$TOTAL > 2000),NA)
    return(importedDF)
}

tLOHCalc <- function(forCalcDF){
    try({
        forCalcDF <- forCalcDF[complete.cases(forCalcDF), ]
        marginalM1 <- apply(forCalcDF[,c("REF","TOTAL")], 
                            MARGIN = 1,
                            FUN = marginalM1Calc, y = 0.5)
        marginalM2_HET <-  apply(forCalcDF[,c("REF","TOTAL")],
                                 MARGIN = 1,
                                 FUN = marginalM2CalcBHET, a = 10, b = 10)
        marginalM2_LOH <- apply(forCalcDF[,c("REF","TOTAL")],
                                MARGIN = 1,
                                FUN = marginalM2CalcBLOH, a = 10, b = 10)
        forCalcDF$`p(D|het)` <- marginalM2_HET
        forCalcDF$`p(D|loh)` <- marginalM2_LOH
        added <- forCalcDF$`p(D|het)` + forCalcDF$`p(D|loh)`
        forCalcDF$`p(het|D)` <- forCalcDF$`p(D|het)` / added
        forCalcDF$`p(loh|D)` <- forCalcDF$`p(D|loh)` / added
        forCalcDF$bayesFactors <- forCalcDF$`p(D|het)` / forCalcDF$`p(D|loh)`
        forCalcDF$inverseBayes <- 1 / forCalcDF$bayesFactors
        forCalcDF$LogBayesFactors <- log(forCalcDF$bayesFactors)
        forCalcDF$LogInverseBayes <- log(forCalcDF$inverseBayes)
        forCalcDF$Log10BayesFactors <- log10(forCalcDF$bayesFactors)
        forCalcDF$Log10InverseBayes <- log10(forCalcDF$inverseBayes)
        forCalcDF$AF <- forCalcDF$ALT / forCalcDF$TOTAL},
        silent=TRUE)
    forCalcDF$CLUSTER <- as.numeric(forCalcDF$CLUSTER)
    forCalcDF$`CLUSTER_AF` <- forCalcDF$CLUSTER + forCalcDF$AF
    forCalcDF <- forCalcDF[!(forCalcDF$CHR == 6 
                                           & forCalcDF$POS > 28510120 
                                           & forCalcDF$POS < 33500500),]
    forCalcDF$CHR_F <- factor(gsub("chr","",forCalcDF$CHR), 
                                     levels=c('1','2','3','4','5','6','7',
                                              '8','9','10','11','12','13',
                                              '14','15','16','17','18',
                                              '19','20','21','22','23',
                                              '24','25'))
    forCalcDF <- forCalcDF[complete.cases(forCalcDF), ]
    return(forCalcDF)
} 

alleleFrequencyPlot <- function(df,sample){
    toPlot <- df
    uniqueClusters <- seq(1.5,1+length(unique(toPlot$CLUSTER)), by = 1)
    topLimit <- tail(uniqueClusters, n=1) + 0.7
    ylines <- seq(1,max(unique(toPlot$CLUSTER)), by = 1)
    labels <- sprintf("Cluster %s", seq(1,length(unique(toPlot$CLUSTER)),
                                        by = 1))
    p1 <- ggplot2::ggplot(toPlot, aes(x = POS, y = `CLUSTER_AF`,
                                      size = Log10InverseBayes)) +
        ggtitle(sample) +
        geom_point(size = 0.50, aes(color = Log10InverseBayes), alpha = 0.75) +
        scale_color_gradient2(low = "darkblue",
                              mid = "transparent",
                              high = "red",
                              midpoint = 0, limits=c(-3, 3),
                              oob = scales::squish) +
        scale_y_continuous(limits = c(0.95,topLimit), expand = c(0,0),
                           breaks = uniqueClusters, labels = labels) +
        facet_grid(~CHR_F, scales = 'free_x', space = 'free_x', switch = 'x') +
        labs(color = expression(paste("Log10( ", frac("1","K"),")"))) +
        xlab("Chromosome") +
        ylab("Allele Fraction") +
        geom_hline(yintercept=ylines, color = "black", size=0.2, alpha = 0.75) +
        theme(panel.spacing = unit(0, "lines"),
              panel.background = element_rect(fill = NA, color = "grey"),
              legend.key = element_rect(colour = NA, fill = NA),
              plot.background = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              panel.grid = element_blank(),
              axis.title.x = element_text(face = "bold"),
              axis.title.y = element_text(face = "bold"),
              axis.text.y = element_text(face = "bold", size = 10))
    return(p1)
}

aggregateCHRPlot <- function(df,sample){
    intermediate <- data.table::setDT(df)
    toPlot2 <- intermediate[, sum(Log10InverseBayes),by=list(CLUSTER,`CHR_F`)]
    names(toPlot2) <- c("CLUSTER","CHR_F","SumLog10InverseBF")
    p2 <- ggplot2::ggplot(toPlot2, aes(x = as.factor(CLUSTER),
                                       y = SumLog10InverseBF,
                                       fill = as.factor(CLUSTER),
                                       color = as.factor(CLUSTER))) +
        geom_bar(stat = "identity") +
        facet_grid(~`CHR_F`, scales = 'free_x', space = 'free_x', 
                   switch = 'x') +
        ggtitle(sample) +
        xlab("Chromosome") +
        ylab("Sum of Log10(1/BF)") +
        geom_hline(yintercept=3, linetype="dashed",
                   color = "black", size=0.50) +
        scale_y_continuous(n.breaks = 15) +
        labs(fill = "CLUSTER", color = "CLUSTER") +
        theme(panel.spacing = unit(0.0001, "lines"),
              panel.background = element_rect(fill = NA, color = NA),
              legend.key = element_rect(colour = NA, fill = NA),
              plot.background = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              panel.grid = element_blank(),
              axis.title.x = element_text(face = "bold"),
              axis.title.y = element_text(face = "bold"),
              axis.text.y = element_text(face = "bold", size = 10))
    return(p2)
}

modePeakCalc <- function(x) {
    if(length(x) > 1){
        den <- density(x, kernel=c("gaussian"))
        ( den$x[den$y==max(den$y)] )     
    } else if(length(x) == 1){
        den <- x
    }
}

