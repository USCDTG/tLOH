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

tLOHCalc <- function(input,sample){
    setupDF <- NULL
    AC <- list.files(input, pattern=".csv")
    setupDF <- lapply(AC, function(x){
        setwd(input)
        message(paste0("- Running analysis for " ,basename(x)))
        ac <- read.csv(x, header=TRUE)
        names(ac) <- c("CHR","CHR_withLabel","POS","POS_oneBefore",
                       "rsID","vcfQUAL","vcf_refNucleotide","vcf_altNucleotide",
                       "INFO","vcfAF","genotype","A","C","G","T","REF","ALT")
        ac$TOTAL <- ac$ALT + ac$REF
        removeOutlierFromCalc(ac,"TOTAL",which(ac$TOTAL > 2000),NA)
        ac$CHR <- gsub("chr","",ac$CHR)
        ac$CHR <- suppressWarnings(as.numeric(ac$CHR))
        ac <- dplyr::filter(ac, grepl('(0, 1)', genotype) &
                                grepl('UTR|synonymous', INFO))
        subset <- ac[c("CHR","POS","REF","ALT","TOTAL")]
        subset <- subset[complete.cases(subset), ]
        subset <- subset[which(subset$TOTAL > 10),]

        try({marginalM1 <- apply(subset, MARGIN = 1,
                                 FUN = marginalM1Calc, y = 0.5)
        marginalM2_HET <-  apply(subset, MARGIN = 1,
                                 FUN = marginalM2CalcBHET, a = 10, b = 10)
        marginalM2_LOH <- apply(subset, MARGIN = 1,
                                FUN = marginalM2CalcBLOH, a = 10, b = 10)
        subset$`p(D|het)` <- marginalM2_HET
        subset$`p(D|loh)` <- marginalM2_LOH
        added <- subset$`p(D|het)` + subset$`p(D|loh)`
        subset$`p(het|D)` <- subset$`p(D|het)` / added
        subset$`p(loh|D)` <- subset$`p(D|loh)` / added
        subset$bayesFactors <- subset$`p(D|het)` / subset$`p(D|loh)`
        subset$inverseBayes <- 1 / subset$bayesFactors
        subset$LogBayesFactors <- log(subset$bayesFactors)
        subset$LogInverseBayes <- log(subset$inverseBayes)
        subset$Log10BayesFactors <- log10(subset$bayesFactors)
        subset$Log10InverseBayes <- log10(subset$inverseBayes)
        subset$AF <- subset$ALT / subset$TOTAL}, silent=TRUE)
        subset$sample <- sample
        subset$Cluster <- unlist(strsplit(basename(x),split="_"))[2]
        subset$Cluster <- gsub("cluster","",subset$Cluster)
        subset$`Cluster_AF` <- as.integer(subset$Cluster) + subset$AF
        subset$CHR_F <- factor(subset$CHR, levels=c('1','2','3','4','5','6','7',
                                                    '8','9','10','11','12','13',
                                                    '14','15','16','17','18',
                                                    '19','20','21','22','23',
                                                    '24'))
        subset <- subset[!(subset$CHR == 6 & subset$POS > 28510120
                           & subset$POS < 33500500),]
        return(subset)
    })
    finalDF <- reduce(setupDF,full_join)
    return(finalDF)
}

plotLOH <- function(df,sample){
    toPlot <- df
    uniqueClusters <- seq(1.5,1+length(unique(toPlot$Cluster)), by = 1)
    topLimit <- tail(uniqueClusters, n=1) + 0.7
    ylines <- seq(1,max(unique(toPlot$Cluster)), by = 1)
    labels <- sprintf("Cluster %s", seq(1,length(unique(toPlot$Cluster)),
                                        by = 1))
    p1 <- ggplot2::ggplot(toPlot, aes(x = POS, y = `Cluster_AF`,
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

    intermediate <- data.table::setDT(toPlot)
    toPlot2 <- intermediate[, sum(Log10InverseBayes),by=list(Cluster,CHR)]
    names(toPlot2) <- c("Cluster","CHR","SumLog10InverseBF")
    p2 <- ggplot2::ggplot(toPlot2, aes(x = as.factor(Cluster),
                                       y = SumLog10InverseBF,
                                       fill = as.factor(Cluster),
                                       color = as.factor(Cluster))) +
        geom_bar(stat = "identity") +
        facet_grid(~CHR, scales = 'free_x', space = 'free_x', switch = 'x') +
        ggtitle(sample) +
        xlab("Chromosome") +
        ylab("Sum of Log10(1/BF)") +
        geom_hline(yintercept=3, linetype="dashed",
                   color = "black", size=0.50) +
        scale_y_continuous(n.breaks = 15) +
        labs(fill = "Cluster", color = "Cluster") +
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
    listOfPlots <- list(p1,p2)
    return(listOfPlots)
}

