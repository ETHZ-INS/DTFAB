# First flavor of monaLisa: MOtif aNAlysis with Lisa

runmonaLisa <- function(DAR, motifs, peaks, genome, nBins=11, minAbsLfc=0.3,
                        background = c("otherBins","zeroBin"), minMatchScore=10){
  
  # At first a normalization should be done
  
  # Compute FC which --> later appended to the GR obj.
  
  peak_FC <- jitter(DAR$logFC, 0.001) # To avoid non-unique breaks
  
  values(peaks) <- DataFrame(peak_FC)

  fc2 <- peak_FC[which(abs(peak_FC)>=minAbsLfc)]
  nElements <- ceiling(length(fc2)/nBins)
  
  
  
  ptm <- proc.time()
  bins <- monaLisa::bin(x = peaks$peak_FC, 
                        binmode = "equalN", 
                        nElements = nElements,
                        minAbsX = minAbsLfc)
  
  # BinDensity <- plotBinDensity(peaks$peak_FC,
  #               bins)

  DARseqs <- getSeq(genome, peaks)
  if(!is.null(names(peaks))) names(DARseqs) <- names(peaks)
  
  # Checking Biases via visualization
  
  # GC_fraction <- plotBinDiagnostics(seqs = DARseqs, 
  #                    bins = bins, 
  #                    aspect = "GCfrac")
  # 
  # Dinucl_fraction <- plotBinDiagnostics(seqs = DARseqs, 
  #                    bins = bins, 
  #                    aspect = "dinucfreq")
  
  # Motif enrichment analysis
  
  se <- calcBinnedMotifEnrR(seqs = DARseqs, 
                            bins = bins, 
                            pwmL = motifs, min.score=minMatchScore,
                            BPPARAM = BiocParallel::MulticoreParam(8),
                            background = match.arg(background))
  
  # Calculate bin-level p-values based on Simes method with code from 
  # https://github.com/markrobinsonuzh/DAMEfinder/blob/master/R/simes_pval.R
  simes <- function(pval){ 
    min((length(pval)*pval[order(pval)])/seq(from=1, to=length(pval), by=1))
  }
  
  ML <- se
  zerobin <- which(colData(ML)$bin.nochange)
  # assays(ML)$signedSigns <- t(t(sign(assays(ML)$log2enr))*
  #                               rep(c(-1,0,1),c(zerobin-1L, 1, ncol(ML)-zerobin)))
  ML <- ML[,-zerobin]
  MLp <- 10^-assays(ML)$negLog10P
  # MLsimes2 <- sapply(seq_len(nrow(ML)), FUN=function(i){
  #   tryCatch({
  #     ss <- assays(ML)$signedSigns[i,]
  #     strongestSign <- ss[which.min(MLp[i,])]
  #     simes(MLp[i,which(ss==strongestSign)])
  #   }, error=function(e) return(NA_real_))
  # })
  # names(MLsimes2) <- row.names(ML)
  
  MLsimes <- apply(MLp, 1, simes)
  
  # Correct the p-values using FDR correction 
  MLdf <- data.frame(p=MLsimes, padj=p.adjust(MLsimes))
  MLdf <- MLdf[order(MLdf$p),]
  MLdf$rank <- seq_along(row.names(MLdf))
  
  # calculate correlation across bins
  cors <- cor(t(assays(ML)$log2enr), seq_len(ncol(ML)), method="spearman")[,1]
  names(cors) <- row.names(ML)
  MLdf$binSpearman <- cors[row.names(MLdf)]
  
  runtime <- proc.time()-ptm
  
  return(list(res=se, runtime=runtime, df=MLdf))
  
}



# legacy, for testing pval aggregation variants
getMLglobalP <- function(ml, method=c("simes","cauchy","weighted.cauchy","geomean")){
  method <- match.arg(method)
  simes <- function(pval){ 
    min((length(pval)*pval[order(pval)])/seq(from=1, to=length(pval), by=1))
  }
  MLp <- 10^-assays(ml)$negLog10P[,-which(colData(ml)$bin.nochange)]
  i <- seq_len(ncol(MLp))
  il <- split(i,cut(i,2))
  apply(MLp, 1, FUN=function(x){
    if(min(x[il[[1]]],na.rm=TRUE)>min(x[il[[2]]],na.rm=TRUE)){
      x <- x[il[[2]]]
      weights <- seq_along(x)
    }else{
      x <- x[il[[1]]]
      weights <- rev(seq_along(x))
    }
    w <- which(!is.na(x))
    if(length(w)==0) return(NA)
    switch(method,
      simes=simes(x[w]),
      cauchy=ACAT::ACAT(x[w]),
      weighted.cauchy=ACAT::ACAT(x[w], weights=weights[w]),
      geomean=10^mean(log10(x)))
  })
}

