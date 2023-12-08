#' helper functions from
#' https://github.com/koenvandenberge/bulkATACGC/blob/master/methods/gcqn_validated.R
FQnorm <- function(counts, type="mean"){
    rk <- apply(counts,2,rank,ties.method='min')
    counts.sort <- apply(counts,2,sort)
    if(type=="mean"){
        # refdist <- apply(counts.sort,1,mean)
        refdist <- base::rowMeans(counts.sort)
    } else if(type=="median"){
        #refdist <- apply(counts.sort,1,median)
        refdist <- matrixStats::rowMedians(counts.sort)
    }
    norm <- apply(rk,2,function(r){ refdist[r] })
    rownames(norm) <- rownames(counts)
    return(norm)
}

gcqn <- function(counts, gcGroups, summary='mean', round=TRUE){
    gcBinNormCounts <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts), 
        dimnames=list(rownames(counts),colnames(counts)))
    for(ii in 1:nlevels(gcGroups)){
        id <- which(gcGroups==levels(gcGroups)[ii])
        if(length(id) == 1){
            normCountBin <- counts[id,]
            if(round) normCountBin <- round(normCountBin)
            gcBinNormCounts[id,] <- normCountBin
            next
        }
        countBin <- counts[id,,drop=FALSE]
        normCountBin <- FQnorm(countBin, type=summary)
        if(round) normCountBin <- round(normCountBin)
        normCountBin[normCountBin<0] <- 0
        gcBinNormCounts[id,] <- normCountBin
    }
    return(gcBinNormCounts)
}

gcqn_qsmooth <- function(counts, gcGroups, bio, round=TRUE){
    gcBinNormCounts <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts), 
        dimnames=list(rownames(counts),colnames(counts)))
    for(ii in 1:nlevels(gcGroups)){
        id <- which(gcGroups==levels(gcGroups)[ii])
        countBin <- counts[id,]
        qs <- qsmooth::qsmooth(countBin, group_factor=bio)
        normCountBin <- qs@qsmoothData
        if(round) normCountBin <- round(normCountBin)
        normCountBin[normCountBin<0L] <- 0L
        gcBinNormCounts[id,] <- normCountBin
    }
    return(gcBinNormCounts)
}

#' wrapper function by Jiayi Wang

GCQuantileNorm <- function(se, genome=NULL, g=20, 
                           summary='median', round=TRUE) {
    
    stopifnot(!is.null(genome) || !is.null(rowData(se)$bias))
    # get GC content
    gr <- rowRanges(se)
    if(!is.null(rowData(se)$bias)){
        gcContent <- rowData(se)$bias
    }else{
        peakSeqs <- getSeq(x = genome, gr)
        gcContent <- letterFrequency(peakSeqs, "GC", as.prob = TRUE)[,1]
    }
    gcGroups <- Hmisc::cut2(gcContent, g = g)
    
    # run QC quantile normalization
    countsGCQN <- gcqn(counts(se), gcGroups, summary = summary, round = round)
}


GCSmoothQuantile <- function(se, genome=NULL, g=20, bio) {
    library(qsmooth)
    # get GC content
    gr <- rowRanges(se)
    se <- as(se,"SummarizedExperiment")
    if(!is.null(rowData(se)$bias)){
        gcContent <- rowData(se)$bias
    }else{
        peakSeqs <- getSeq(x = genome, gr)
        gcContent <- letterFrequency(peakSeqs, "GC", as.prob = TRUE)[,1]
    }
    gcGroups <- Hmisc::cut2(gcContent, g = g)
    
    # run QC smooth quantile normalization
    # colData(se)$group_id <- factor(substr(colnames(se), 1, 
    #     nchar(colnames(se)) - 5))
    
    # countsGCSQ <- gcqn_qsmooth(counts(se), gcGroups,
    #        bio=droplevels(colData(se)$group_id))
    countsGCSQ <- gcqn_qsmooth(se, gcGroups,
                               bio=droplevels(as.factor(colData(se)[,bio])))
}

