marginalLikelihoodM1 <- function(x,a,b){
    e=0.01
    successes <- x["REF"]
    trials <- x["TOTAL"]
    p <- function(theta){dbinom(x = successes, size = trials, prob = theta) *
            (((1-e)*dbeta(x = abs(0.5 - theta), 
                          shape1 = a, 
                          shape2 = b))+((e)*dbeta(x = theta, 
                                                  shape1 = 1, shape2 = 1)))
    }
    integrate(f = p, lower = 0, upper = 1)[[1]]
}

marginalLikelihoodM2 <- function(x,a,b){
    e=0.01
    successes <- x["REF"]
    trials <- x["TOTAL"]
    p <- function(theta){dbinom(x = successes, size = trials, prob = theta) *
            (((1-e)*dbeta(x = abs(theta), 
                          shape1 = a, 
                          shape2 = b))+((e)*dbeta(x = theta, 
                                                  shape1 = 1, shape2 = 1)))
    }
    integrate(f = p, lower = 0, upper = 1)[[1]]
}

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


tLOHCalc <- function(forCalcDF,alpha1,beta1,alpha2,beta2,countThreshold){
    try({
        marginalM1_LOH <- apply(forCalcDF[,c("REF","TOTAL")],
                                MARGIN = 1,
                                FUN = purrr::possibly(marginalLikelihoodM1,NA), 
                                a = alpha1, b = beta1)
        
        marginalM2_HET <-  apply(forCalcDF[,c("REF","TOTAL")],
                                 MARGIN = 1,
                                 FUN = purrr::possibly(marginalLikelihoodM2,NA), 
                                 a = alpha2, b = beta2)
        forCalcDF$`p(D|loh)` <- marginalM1_LOH
        forCalcDF$`p(D|het)` <- marginalM2_HET
        added <- forCalcDF$`p(D|het)` + forCalcDF$`p(D|loh)`
        forCalcDF$`p(het|D)` <- forCalcDF$`p(D|het)` / added
        forCalcDF$`p(loh|D)` <- forCalcDF$`p(D|loh)` / added
        forCalcDF$bayesFactors <-  forCalcDF$`p(D|loh)` / forCalcDF$`p(D|het)`
        forCalcDF$inverseBayes <- 1 / forCalcDF$bayesFactors
        forCalcDF$LogBayesFactors <- log(forCalcDF$bayesFactors)
        forCalcDF$LogInverseBayes <- log(forCalcDF$inverseBayes)
        forCalcDF$Log10BayesFactors <- log10(forCalcDF$bayesFactors)
        forCalcDF$Log10InverseBayes <- log10(forCalcDF$inverseBayes)
        forCalcDF$AF <- forCalcDF$ALT / forCalcDF$TOTAL
    })
    forCalcDF$CLUSTER <- as.numeric(forCalcDF$CLUSTER)
    forCalcDF$`CLUSTER_AF` <- forCalcDF$CLUSTER + forCalcDF$AF
    forCalcDF <- forCalcDF[!(forCalcDF$CHR == 'chr6' 
                             & forCalcDF$POS > 28510120 
                             & forCalcDF$POS < 33500500),]
    forCalcDF <- forCalcDF[!(forCalcDF$CHR == 6 
                             & forCalcDF$POS > 28510120 
                             & forCalcDF$POS < 33500500),]
    forCalcDF$CHR <- gsub("chr","",forCalcDF$CHR)
    forCalcDF$CHR_F <- factor(gsub("chr","",forCalcDF$CHR), 
                              levels=c('1','2','3','4','5','6','7',
                                       '8','9','10','11','12','13',
                                       '14','15','16','17','18',
                                       '19','20','21','22','23',
                                       '24','25'))
    forCalcDF <- forCalcDF[complete.cases(forCalcDF), ]
    forCalcDF <- forCalcDF[which(forCalcDF$TOTAL >= countThreshold),]
    forCalcDF$alpha <- alpha1
    forCalcDF$beta <- beta1
    forCalcDF$alpha2 <- alpha2
    forCalcDF$beta2 <- beta2
    forCalcDF <- dplyr::arrange(forCalcDF,CLUSTER,CHR,POS)
    return(forCalcDF)
}


tLOHCalcUpdate <- function(forCalcDF,alpha1,beta1,alpha2,beta2,countThreshold){
    try({
        marginalM1_LOH <- apply(forCalcDF[,c("REF","TOTAL")],
                                MARGIN = 1,
                                FUN = purrr::possibly(marginalLikelihoodM1,NA), 
                                a = alpha1, b = beta1)
        
        marginalM2_HET <-  apply(forCalcDF[,c("REF","TOTAL")],
                                 MARGIN = 1,
                                 FUN = purrr::possibly(marginalLikelihoodM2,NA), 
                                 a = alpha2, b = beta2)
        forCalcDF$`p(D|loh)` <- marginalM1_LOH
        forCalcDF$`p(D|het)` <- marginalM2_HET
        added <- forCalcDF$`p(D|het)` + forCalcDF$`p(D|loh)`
        forCalcDF$`p(het|D)` <- forCalcDF$`p(D|het)` / added
        forCalcDF$`p(loh|D)` <- forCalcDF$`p(D|loh)` / added
        forCalcDF$bayesFactors <-  forCalcDF$`p(D|loh)` / forCalcDF$`p(D|het)`
        forCalcDF$inverseBayes <- 1 / forCalcDF$bayesFactors
        forCalcDF$LogBayesFactors <- log(forCalcDF$bayesFactors)
        forCalcDF$LogInverseBayes <- log(forCalcDF$inverseBayes)
        forCalcDF$Log10BayesFactors <- log10(forCalcDF$bayesFactors)
        forCalcDF$Log10InverseBayes <- log10(forCalcDF$inverseBayes)
        forCalcDF$AF <- forCalcDF$ALT / forCalcDF$TOTAL
    })
    forCalcDF$CLUSTER <- as.numeric(forCalcDF$CLUSTER)
    forCalcDF$`CLUSTER_AF` <- forCalcDF$CLUSTER + forCalcDF$AF
    forCalcDF <- forCalcDF[!(forCalcDF$CHR == 'chr6' 
                             & forCalcDF$POS > 28510120 
                             & forCalcDF$POS < 33500500),]
    forCalcDF <- forCalcDF[!(forCalcDF$CHR == 6 
                             & forCalcDF$POS > 28510120 
                             & forCalcDF$POS < 33500500),]
    forCalcDF$CHR <- gsub("chr","",forCalcDF$CHR)
    forCalcDF$CHR_F <- factor(gsub("chr","",forCalcDF$CHR), 
                              levels=c('1','2','3','4','5','6','7',
                                       '8','9','10','11','12','13',
                                       '14','15','16','17','18',
                                       '19','20','21','22','23',
                                       '24','25'))
    forCalcDF <- forCalcDF[complete.cases(forCalcDF), ]
    forCalcDF <- forCalcDF[which(forCalcDF$TOTAL >= countThreshold),]
    forCalcDF$alpha <- alpha1
    forCalcDF$beta <- beta1
    forCalcDF$alpha2 <- alpha2
    forCalcDF$beta2 <- beta2
    forCalcDF <- dplyr::arrange(forCalcDF,CLUSTER,CHR,POS)
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

splitByChromosome <- function(listOfDataframes,numberOfDataframes){
    output <- list()
    numberOfChromosomes <- seq(1,22,1)
    output <- lapply(numberOfDataframes, 
                     function(x) lapply(numberOfChromosomes, 
                                        function(s) 
                                            dplyr::filter(listOfDataframes[[x]], 
                                                          CHR == s)))
    return(output)
}


prepareHMMdataframes <- function(importedData){
    numberOfClusters <- unique(importedData$CLUSTER)
    numberOfChromosomes <- seq(1,22,1)
    df_by_cluster <- lapply(numberOfClusters, 
                            function(x) 
                                importedData[which(importedData$CLUSTER == 
                                                       as.numeric(x)),])
    df_by_cluster_by_chr <- splitByChromosome(df_by_cluster,numberOfClusters)
    nObs <- seq(1,length(numberOfClusters) * length(numberOfChromosomes),1)
    toAnalyze <- purrr::flatten(df_by_cluster_by_chr)
    orderedNorm <- list()
    output <- list()
    
    for(i in nObs){
        if(nrow(toAnalyze[[i]]) == 0){
            toAnalyze[[i]]$orderNorm <- NULL
            output[[i]] <- toAnalyze[[i]]
        } else{
            output[[i]] <- cbind(toAnalyze[[i]],
                                 orderNorm = bestNormalize::orderNorm(
                                     toAnalyze[[i]]$bayesFactors, 
                                     warn = FALSE)$x.t)
        }
    }
    return(output)
}

runHMM_1 <- function(dataframeList, initProbs, trProbs){
    initProbs <- initProbs
    trX <- trProbs
    output <- list()
    for(i in 1:length(dataframeList)){
        CHRvalue <- unique(dataframeList[[i]]$CHR)
        initProb1 <- initProbs[CHRvalue,1] #notLOH
        initProb2 <- initProbs[CHRvalue,2] #LOH
        tryCatch(output[[i]] <- depmixS4::depmix(orderNorm~1,
                                                 data = dataframeList[[i]],
                                                 family = gaussian(),
                                                 trstart = trX,
                                                 instart = c(initProb1,
                                                             initProb2),
                                                 nstates = 2),
                 error=function(cond) {
                     message(cond)
                     return(NA)
                 })
    }
    return(output)
}

runHMM_2 <- function(dataframeList){
    output <- list()
    for(i in 1:length(dataframeList)){
        tryCatch(output[[i]] <- depmixS4::fit(dataframeList[[i]],
                                              emc=depmixS4::
                                                  em.control(rand=FALSE),
                                              type = 'viterbi'),
                 error=function(cond) {
                     message(cond)
                     return(output[[i]] <- 'NA')
                 })
    }
    return(output)
}

runHMM_3 <- function(dataframeList){
    output <- list()
    for(i in 1:length(dataframeList)){
        tryCatch(output[[i]] <- depmixS4::
                     posterior(dataframeList[[i]],
                               type = 'viterbi'),
                 error=function(cond) {
                     message(cond)
                     return(NA)
                 })
    }
    return(output)
}


regionAnalysis <- function(originalDF,dataframeList){
    finalList <- documentErrorRegions(originalDF,dataframeList)
    for(i in 1:length(finalList)){
        df <- finalList[[i]] |> dplyr::arrange(CLUSTER,CHR,POS)
        rlel <- rle(df$state)$lengths
        df$group <- unlist(lapply(1:length(rlel), 
                                  function(i) rep(i, rlel[i])))
        df$Log10BayesFactorsEdited <- df$Log10BayesFactors
        intermediateSum <- df |>
            dplyr::group_by(group) |>
            naniar::replace_with_na_at(.vars = 'Log10BayesFactorsEdited',
                               condition = ~(.x > -0.5 & .x < 0.5)) |>
           dplyr::summarize(sequentialSum = sum(Log10BayesFactorsEdited, 
                                          na.rm = TRUE),
                      sequentialSumOriginalBF = sum(Log10BayesFactors),
                      medianBF = median(Log10BayesFactors),
                      meanBF = mean(Log10BayesFactors))
        toSave <- merge(df,intermediateSum,by="group")
        toSave$state <- as.character(toSave$state)
        finalList[[i]] <- toSave
    }
    return(finalList)
} 

documentErrorRegions <- function(a,b){
    intermediate <- a
    est.states <- b
    for(i in 1:length(est.states)){
        if(is.null(est.states[[i]]) == TRUE){
            lengthForNewDF <- nrow(intermediate[[i]])
            est.states[[i]] <- data.frame('state' = rep('errored',
                                                        lengthForNewDF),
                                          'S1' = rep(0,lengthForNewDF),
                                          'S2' = rep(0,lengthForNewDF))
        } else{
            test <- 'false'
        }
    }
    emptyList <- list()
    nDF <- seq(1,as.numeric(length(intermediate)),1)
    loop = nDF
    loop = loop[-length(loop)]
    
    loop = seq(1,length(est.states),1)
    for(i in loop){
        if(is.null(est.states[[i]]) == TRUE){
            emptyList[[i]] <- data.frame(NULL)
        } else if(nrow(est.states[[i]]) == 0){
            emptyList[[i]] <- data.frame(NULL)
        } else if(nrow(est.states[[i]] > 1)){
            emptyList[[i]] <- cbind(intermediate[[i]],est.states[[i]])
        }
    }
    finalList <- emptyList[sapply(emptyList, function(x) dim(x)[1]) > 0]
}

summarizeRegions1 <- function(x){
    finalList <- x
    ListOfSegments <- list()
    for(i in 1:length(finalList)){
        df1 <- data.frame(state = paste0('state',finalList[[i]]$state), 
                          start = finalList[[i]]$POS)
        segments <- df1[order( df1$start ),] |> 
            dplyr::group_by(state, rleid = with(rle(state), 
                                         rep(seq_along(lengths), 
                                             lengths))) |> 
                                         dplyr::summarize(
                                             intervalStart = min(start),
                                             intervalEnd  = max(start)) |> 
                                         dplyr::arrange(rleid)
        segments$CHR <- unique(finalList[[i]]$CHR)
        segments$CLUSTER <- unique(finalList[[i]]$CLUSTER)
        ListOfSegments[[i]] <- segments
    }
    output <- purrr::reduce(ListOfSegments,full_join)
    return(output)
}

summarizeRegions2 <- function(finalTable){
    sampleValues <- finalTable
    collapseTable <- sampleValues[c('group','CHR','POS',
                                    'sequentialSumOriginalBF',
                                    'medianBF','meanBF','CLUSTER')]
    collapsed <- collapseTable |>
        dplyr::group_by_at(dplyr::vars(CLUSTER,group,CHR,
                                       sequentialSumOriginalBF,
                         medianBF,meanBF)) |>
        dplyr::summarize_all(paste, collapse=",")
    collapsed$intervalStart <- unlist(lapply(
        stringr::str_split(collapsed$POS,','),
        function(x) min(as.numeric(x))))
    collapsed$intervalEnd <- unlist(lapply(
        stringr::str_split(collapsed$POS,','),
        function(x) max(as.numeric(x))))
    collapsed$POS <- NULL
    collapsed$lengthOfInterval <- collapsed$intervalEnd - 
        collapsed$intervalStart
    collapsed$lengthOfInterval[collapsed$lengthOfInterval == 0] <- 1
    sampleDF <- collapsed
    sampleDF$CHR <- as.numeric(sampleDF$CHR)
    return(sampleDF)
}

regionFinalize <- function(finalList1){
    finalList1 <- purrr::map(finalList1, 
                             ~dplyr::mutate(.x, state = as.character(state)))
    sampleValues <- as.data.frame(purrr::reduce(finalList1,full_join))
    sampleData <- summarizeRegions1(finalList1)
    sampleData$lengthOfInterval <- sampleData$intervalEnd -
        sampleData$intervalStart
    sampleDF <- summarizeRegions2(sampleValues)
    output <- merge(sampleValues,sampleDF,
                    by=c('group','sequentialSumOriginalBF','medianBF',
                         'meanBF','CHR','CLUSTER'))
    output$CHR <- as.numeric(output$CHR)
    data2 <- output |>
        dplyr::arrange(CLUSTER,CHR,POS) |>
        dplyr::group_by(CLUSTER,CHR,medianBF) |>
        dplyr::summarise(nSNPs_in_Region = dplyr::n(), dplyr::across()) |>
        dplyr::arrange(CLUSTER,CHR,POS) |>
        dplyr::group_by(CLUSTER,CHR,POS,
                        segmentNumber = with(rle(state),
                                             rep(seq_along(lengths), 
                                                 lengths))) |>
        dplyr::arrange(CLUSTER,CHR,POS)
    
    valuesToTest <- data2 |> dplyr::group_by(segmentNumber) |>
        dplyr::filter(quantile(TOTAL, 0.75)<TOTAL) |>
        dplyr::select(AF)
    valuesToTest$DeltaAlleleFraction <- abs(valuesToTest$AF-0.5)
    dt <- data.table::data.table(valuesToTest)
    dt2 <- dt[, modePeakCalc(DeltaAlleleFraction), by= segmentNumber]
    names(dt2)[2] <- 'modePeak'
    dt2 <- data.frame(merge(dt2,dt,by='segmentNumber'))
    a <- merge(dt2,data2,by=c('segmentNumber','AF'),all=TRUE)
    a$modePeak <- NULL
    lookup <- unique(dt2[c('segmentNumber','modePeak')])
    final <- merge(a,lookup,by=c('segmentNumber'))
    final$group  <- NULL
    outputFinal <- final
    outputFinal <- dplyr::arrange(outputFinal,CLUSTER,CHR,POS)
    return(outputFinal)
}

hiddenMarkovAnalysis <- function(df, initProbs, trProbs){
    list1 <- prepareHMMdataframes(df)
    step1 <- runHMM_1(list1, initProbs, trProbs)
    step2 <- runHMM_2(step1)
    step3 <- runHMM_3(step2)
    step4 <- regionAnalysis(list1,step3)
    done <- regionFinalize(step4)
    done <- done[c('CLUSTER','rsID','CHR','POS','REF','ALT','TOTAL','AF',
                   'bayesFactors','Log10BayesFactors','p(D|loh)','p(D|het)',
                   'p(het|D)','p(loh|D)','orderNorm','state','S1','S2',
                   'modePeak','medianBF','meanBF','sequentialSum',
                   'sequentialSumOriginalBF','intervalStart','intervalEnd',
                   'lengthOfInterval','nSNPs_in_Region')] |> 
        dplyr::arrange(CLUSTER,CHR,POS) |> 
        dplyr::filter(lengthOfInterval > 1000)
    names(done) <- c(
        'cluster','rsID','chromosome','position','referenceCounts',
        'alternativeCounts','totalCounts','alleleFraction','bayesFactors',
        'log10BayesFactors','p(D|loh)','p(D|het)','p(het|D)','p(loh|D)',
        'orderedNormalization','state','probabilityOfS1','probabilityOfS2',
        'segment_modePeak','segment_medianBF','segment_meanBF',
        'segment_sequentialSumAdjusted','segment_sequentialSum',
        'segment_intervalStart',
        'segment_intervalEnd','segment_lengthOfInterval','segment_nSNPs')
    return(done)
}


