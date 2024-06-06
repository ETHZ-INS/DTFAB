#@' author Emanuel Sonder

sampleSwitch <- function(total, size, seed=42){ 
  set.seed(seed)
  data.table::setDTthreads(2)
  if(total >= size)
  {
    s <- sample(total, size=size, replace=FALSE)
  }
  else
  {
    s <- sample(total, size=size, replace=TRUE)
  }
  
  return(s)
}

.importFrags <- function(bamPath, which=NULL, fracSub=1, annotationStyle="NCBI")
{
  data.table::setDTthreads(2)
  if(is.null(which))
  {
    # currently only set for human
    which <- GRanges(Rle(1:22), IRanges(start=rep(1,22), width=536870912))
  }
  
  if(endsWith(bamPath, ".bam"))
  {
    param <- ScanBamParam(what=c('pos', 'qwidth', 'isize'),
                          which=which, 
                          flag=scanBamFlag(isUnmappedQuery=FALSE))
  
   readPairs <- readGAlignmentPairs(bamPath, param=param)
  
   # get fragment coordinates from read pairs
   frags <- GRanges(seqnames(GenomicAlignments::first(readPairs)), 
                    IRanges(start=pmin(GenomicAlignments::start(GenomicAlignments::first(readPairs)), 
                                       GenomicAlignments::start(GenomicAlignments::second(readPairs))), 
                            end=pmax(GenomicAlignments::end(GenomicAlignments::first(readPairs)), 
                                     GenomicAlignments::end(GenomicAlignments::second(readPairs)))))
   frags <- keepStandardChromosomes(frags, pruning.mode="coarse")
  
   frags <- granges(frags, use.mcols=TRUE)
  }
  else if(endsWith(bamPath, ".bed")){
    frags <- fread(bamPath)
    frags <- makeGRangesFromDataFrame(as.data.frame(frags)) 
    frags <- keepStandardChromosomes(frags, pruning.mode="coarse")
    frags <- subsetByOverlaps(frags, which)
    readPairs <- NULL
  }
  else if(endsWith(bamPath, ".rds")){
    frags <- readRDS(bamPath)
    frags <- keepStandardChromosomes(frags, pruning.mode="coarse")
    frags <- subsetByOverlaps(frags, which)
    readPairs <- NULL
  }
  
  if(fracSub<1)
  {
    n <- length(frags)
    nSub <- ceiling(n*fracSub)
    
    frags <- frags[sample(1:n, nSub, replace=FALSE)]
    readPairs <- readPairs[sample(1:n, nSub, replace=FALSE)]
  }
  
  # change annotation style
  seqlevelsStyle(frags) <- annotationStyle

  return(list(fragments=frags,
              readPairs=readPairs))
}

.importPeaks <- function(peaksPath, 
                         which=NULL, 
                         annotationStyle="NCBI",
                         colNames=c("chr","start", "end", 
                                    "name", "score", "strand",
                                    "signalValue", "pValue", "qValue", "peak"))
{ 
  data.table::setDTthreads(2)
  peaks <- fread(peaksPath, col.names=colNames)

  peaks <- makeGRangesFromDataFrame(peaks, keep.extra.columns=TRUE)
  seqlevelsStyle(peaks) <- annotationStyle
  peaks <- keepStandardChromosomes(peaks, pruning.mode="coarse")
  
  peaks <- subsetByOverlaps(peaks, which)
  
  return(peaks)
}


# estimate the log fold change based on the enrichment over the input
.estLfc <- function(enrichment, 
                    lfcDistDt,
                    enrCol="enr",
                    lfcCol="lfc",
                    seed=42)
{
  data.table::setDTthreads(2)
  set.seed(seed)
  
  #TODO: Change
  enrPeakDt <- as.data.table(enrichment)
  lfcDistDt <- as.data.table(lfcDistDt)
  lfcDistDt <- copy(lfcDistDt)

  # empirical cdf of observed enrichments
  lfcCDF <- ecdf(lfcDistDt[[enrCol]])
  
  # quantile normalization of enrichment over input of ChIP-peaks
  normEnr <- preprocessCore::normalize.quantiles.use.target(matrix(enrPeakDt[[enrCol]], ncol=1),
                                            target=lfcDistDt[[enrCol]], copy=FALSE)
  enrSampQuant <- lfcCDF(normEnr)
  
  # match corresponding log foldchanges
  setorderv(lfcDistDt, enrCol)
  ind <- round(enrSampQuant*nrow(lfcDistDt))
  ind[ind==0] <- 1
  enrPeakDt$lfc <- lfcDistDt$lfc[ind]
  enrPeakDt[,peak_id:=1:nrow(enrPeakDt)]
  
  return(enrPeakDt)
}

# biasFileDir: deepTools output
.varyGCBias <- function(frags, 
                        biasFileDir, 
                        frac=1, 
                        minGC=0, 
                        maxGC=1,
                        annotationStyle="NCBI",
                        genome=BSgenome.Hsapiens.UCSC.hg38,
                        seed=42)
{
  set.seed(seed)
  data.table::setDTthreads(2)
  
  # Import GC bias table (output deeptools)
  gcBiasTable <- fread(biasFileDir, col.names=c("bin", "prob"), header=FALSE, sep="\t")
  gcBiasTable[,gc_bin:=seq(from=0,to=1, by=1/nrow(gcBiasTable))[1:nrow(gcBiasTable)]]
  
  # subset gc rates
  gcBiasTableSub <- subset(gcBiasTable, gc_bin>=minGC & gc_bin<=maxGC)
  
  gcBiasTableSub[,gc_bin_start:=gc_bin]
  gcBiasTableSub[,gc_bin_end:=data.table::shift(gc_bin_start, 1, type="lead")]
  gcBiasTableSub[nrow(gcBiasTableSub),]$gc_bin_end <- maxGC
  
  # get gc content per frag
  frags <- makeGRangesFromDataFrame(as.data.frame(frags))
  seqlevelsStyle(frags) <- "UCSC"
  frags$GC_content <- suppressWarnings(Repitools::gcContentCalc(frags, 
                                                                organism=genome))
  seqlevelsStyle(frags) <- annotationStyle
  
  # annotate each frag with gc bin
  fragsDt <- as.data.table(frags) # convert to data.table for simplicity
  fragsDt$min_GC_content <- fragsDt$GC_content
  fragsDt$max_GC_content <- fragsDt$GC_content  
  
  setkey(gcBiasTableSub, gc_bin_start, gc_bin_end)
  fragsDt <- foverlaps(fragsDt, gcBiasTableSub,
                          by.x=c("min_GC_content", "max_GC_content"),
                          by.y=c("gc_bin_start", "gc_bin_end"))
  
  # Calculate probability for sampling for each gc bin
  nSample <- round(nrow(fragsDt)*frac)
  fragsDt[,sample_size:=floor(prob*nrow(fragsDt)), by=bin]
  
  # sample from gc bins
  fragsSubDt <- fragsDt[,.SD[sampleSwitch(.N, unique(sample_size), seed=seed)], 
                        by=gc_bin]

  return(fragsSubDt)
}

# fragsDt needs to contain width
.varyFragDist <- function(fragsDt, 
                          frac=0.8,
                          nClust=4, 
                          fitGMM=FALSE,
                          referenceData=NULL,
                          prob=c(0.5, 0.4, 0.08, 0.02),
                          cuts=NULL,
                          estimateProb=FALSE,
                          seed=42)
{ 
  set.seed(seed)
  data.table::setDTthreads(2)
  fragsDt[,width:=end-start]
  if(is.null(cuts)) cuts <- c(0, 120, 300, 500, max(fragsDt$width))
  else cuts <- unique(c(0, cuts, max(max(fragsDt$width), max(cuts))))
  
  if(fitGMM)
  {
    # Estimate fragment length distributions
    if(is.null(referenceData))
    {
      fit = mclust::Mclust(fragsDt$width, G=nClust, model="V")
    }
    else
    {
      fit = mclust::Mclust(referenceData$width, G=nClust, model="V")
    }
    
    if(estimateProb)
    {
      prob <- fit$parameters$pro
    }
    
    # Annotate each fragment with the most likely cluster
    fragsTable$cluster <- fit$classification
  }
  else
  {
    fragsDt[,cluster:=as.numeric(cut(get("width"), breaks=cuts, include.lowest=TRUE))]
  }
  
  # Number of reads to sample for each type
  nSample <- round(nrow(fragsDt)*frac)
  nFragsType <- round(nSample*prob)
  
  fragsSubDts <- lapply(1:nClust, function(i){
    data.table::setDTthreads(2)
    sampledFrags <- fragsDt[cluster==i, ][sampleSwitch(.N, nFragsType[i], seed=seed),]
  })
  
  # sample according to probability of each fragment type
  typeDt <- data.table(cluster=unique(fragsDt$cluster))
  setorder(typeDt, cluster)
  typeDt$sample_size <- nFragsType
  fragsDt <- merge(fragsDt, typeDt, by=c("cluster"))
  
  fragsSubDt <- fragsDt[,.SD[sampleSwitch(.N, unique(sample_size), seed=seed)], 
                        by=cluster]
  
  return(fragsSubDt)
}

.varEffectSize <- function(fragDt, 
                           peaks,
                           logFCs=NULL,
                           effectStrength=0.5,
                           noiseLevel=0.1,
                           maxReadPerPeak=300,
                           seed=42){
  data.table::setDTthreads(2)
  set.seed(seed)
  
  if(is.null(logFCs))
  {
    logFCs <- rep(0, nrow(peakDt))
  }
  
  # convert to data.table for simplicity
  peakDt <- as.data.table(peaks)
  peakDt$peak_id <- 1:nrow(peakDt)
  peakDt$logFC <- logFCs
  
  fragDt[,frag_id:=paste(seqnames, start, end, sep="_")]
  
  # get frags within peaks
  setkey(peakDt, seqnames, start, end)
  fragInPeakDt <- foverlaps(fragDt, 
                            peakDt,
                            type="any",
                            nomatch=NULL,
                            by.x=c("seqnames", "start", "end"),
                            by.y=c("seqnames", "start", "end"),
                            minoverlap=1,
                            mult="all")
  
  # Get frags outside of peaks
  fragOutPeakDt <- fragDt[!unique(fragInPeakDt, by="frag_id"), on=.(frag_id)]
  
  # Calculate the number of fragments to sample per peak
  fragInPeakDt[,fc:=2^abs(logFC)]
  if(effectStrength>0)
  {
    fragInPeakDt[, n_frags:=min(.N, ceiling(.N/(data.table::first(fc)*effectStrength))), 
                   by=c("peak_id")]
    #fragInPeakDt[, n_frags:=max(.N, floor(.N*data.table::first(fc)*effectStrength)), 
    #             by=c("peak_id")]
  }
  else
  {
    fragInPeakDt[, n_frags:=.N, by=c("peak_id")]
  }
  
  # Only sample from positive logFcs run for multiple samples
  #if(nrow(fragInPeakDt[logFC>0,])>0 & effectStrength>0)
  if(nrow(fragInPeakDt[logFC<0,])>0 & effectStrength>0)
  {
   # fragInPeakSubDt <- fragInPeakDt[logFC>0,][,.SD[
  #    sampleSwitch(.N, min(maxReadPerPeak, data.table::first(n_frags)), seed=seed)], by=peak_id]
   # fragInPeakSubDt <- rbind(fragInPeakSubDt, fragInPeakDt[logFC<=0,])
    
    fragInPeakSubDt <- fragInPeakDt[logFC<0,][,.SD[
      sampleSwitch(.N, min(maxReadPerPeak, data.table::first(n_frags)), seed=seed)], by=peak_id]
    fragInPeakSubDt <- rbind(fragInPeakSubDt, fragInPeakDt[logFC>=0,])
  }
  else
  {
    fragInPeakSubDt <- fragInPeakDt
  }
  
  fragInPeakSubDt <- fragInPeakSubDt[,c("seqnames", "i.start", "i.end"), with=FALSE]
  colnames(fragInPeakSubDt) <- c("seqnames", "start", "end")
  
  # add back fragments outside of peaks
  fragSubDt <- rbind(fragInPeakSubDt, 
                     fragOutPeakDt[,  c("seqnames", "start", "end"), with=FALSE])
  
  return(fragSubDt)
}

#' Function to vary different ATAC-seq properties
#'
#' Function to vary different Atac-seq related properties by subsampling fragments.
#' It is possible to vary GC-Bias, Fragment size distribution and effect strenght (of an experimental condition)
#' 
#' @param bamPath path to the .Bam file to subsample fragments from
#' @param bedPath path to the .Bed file with the peaks
#' @param sampleName name of the .Bam file to be analyized (subsampled .Bam gets saved under this name)
#' @param design a vector of 1 & -1 indicating to which of the two groups a
#' sample belongs to
#' @param logFCs vector with log fold changes for each peak. Needs to be of the same length 
#' as the peaks in the bed file. 
#' @param which GRanges object stating which regions should be selected
#' @param fracSub proportion of reads to subsample in each step 
#' (three steps in total, GC-bias based, fragment size based and effect strength based subsampling)
#' @param fragFragsOutside proportion of fragments outside of peaks
#' @param effectStrength Factor to which the fragments of differentially accessible peaks should be
#' sub- / upsampled. 
#' @param nFragTypes Number of fragment types (nucleosome-free, mono-, di-, tri-nucleosome, etc. containing fragments)
#' @param fitGMM TRUE/FALSE if GMM should be fitted for fragment length distribution or fixed thresholds should be used.
#' @param prob vector of probabilities of fragment types, needs same length as nFragTypes
#' @param estimateProb TRUE/FALSE if probabilities of fragment types should be estimated from data
#' @param minOverlap minimal overlap to consider a fragment being within a peak
#' @param minGC minimal gc content of fragments to be considered (to be removed)
#' @param maxGC maximal gc content of fragments to be considered (to be removed)
#' @param annotationStyle Either NCBI or UCSC
#' @param genome BSgenome object to be used. 
#' @return data.table with subsampled fragment coordinates & .bam file of these fragments saved on disk.
varyAtacSignal <- function(bamPath, 
                           peaks, 
                           atacPeaks,
                           sampleName,
                           chIPlogFCs,
                           atacLogFCs,
                           depth,
                           biasFileDir=NULL,
                           which=NULL, 
                           fracSub=0.8,
                           effectStrength=0.5,
                           nFragTypes=4,
                           fitGMM=FALSE,
                           prob=c(0.6, 0.3, 0.07, 0.03),
                           estimateProb=FALSE,
                           minOverlap=1,
                           minGC=0,
                           maxGC=1,
                           simGCBias=TRUE,
                           simFLD=TRUE,
                           simEffectStrength=TRUE,
                           varyAtacPeaks=TRUE,
                           maxReadPerPeak=300,
                           matchCounts=FALSE,
                           annotationStyle="NCBI",
                           genome=BSgenome.Hsapiens.UCSC.hg38,
                           seed=42)
{
  set.seed(seed)
  data.table::setDTthreads(2)
  
  # Import fragments
  bamData <- .importFrags(bamPath, which, fracSub, annotationStyle)
  frags <- bamData$fragments
  readPairs <- bamData$readPairs
  fragRanges <- frags
  fragsDt <- as.data.table(frags) 
  
  fragsDt <- fragsDt[,.SD[sample(.N, depth)]]
  
  if(simEffectStrength)
  {
    print("vary chIP peaks")
    # Vary Effect size of ChIP-peaks
    fragSubDt <- .varEffectSize(fragsDt, peaks, 
                                effectStrength, 
                                logFCs=chIPlogFCs,
                                maxReadPerPeak=maxReadPerPeak,
                                seed=seed)
    
    # like this correspondance to logFCs is lost
    if(varyAtacPeaks)
    {
      print("Vary atac peaks")
      atacPeakDt <- cbind(as.data.table(atacPeaks), 
                          data.table(logFCs=atacLogFCs))
      atacSolePeakDt <- atacPeakDt[!overlapsAny(atacPeaks, peaks),]
      
      if(matchCounts)
      {
        atacSolePeakRanges <- makeGRangesFromDataFrame(as.data.frame(atacSolePeakDt))
        atacCounts <- countOverlaps(atacSolePeakRanges, fragRanges)
        chIPCounts <- countOverlaps(peaks, fragRanges)
        iqrChIP <- iqr(chIPCounts)
        medianChIP <- median(chIPCounts)
      
        idAtacPeaks <- which((atacCounts<=(medianChIP+iqrChIP)) & 
                             (atacCounts>=(medianChIP-iqrChIP)))
        atacSolePeakDt <- atacSolePeakDt[idAtacPeaks, ]
      
        # match width
        #atacSolePeakDt[,med:=floor((end-start)/2+start)]
        #medWidth <- median(width(peaks))
        #atacSolePeakDt[,start:=med-floor(medWidth/2)]
        #atacSolePeakDt[,end:=med+floor(medWidth/2)]
      }
      
      nPeaksOverlaps <- length(peaks[overlapsAny(peaks, atacPeaks)])
      atacSubSolePeakDt <-  atacSolePeakDt[sample(1:nrow(atacSolePeakDt), 
                                                  min(nrow(atacSolePeakDt), 
                                                      nPeaksOverlaps)),]
      atacSubSolePeakRanges <- makeGRangesFromDataFrame(as.data.frame(atacSubSolePeakDt))
      fragSubDt <- .varEffectSize(fragSubDt, atacSubSolePeakRanges, 
                                  effectStrength, 
                                  logFCs=atacSubSolePeakDt$logFCs,
                                  maxReadPerPeak=maxReadPerPeak,
                                  seed=seed)
    }
  }
  else fragSubDt <- fragsDt
  
  # Vary GC Bias
  if(simGCBias)
  {
    if(is.null(biasFileDir))
    {
      stop("Bias-file directory has to be provided if gc bias should be simulated (simGCBias=TRUE)")  
    }
    
    if(!is.na(biasFileDir))
    {
      print("Bias file is specified")
      fragSubDt <- .varyGCBias(fragSubDt, biasFileDir, fracSub,
                           minGC, maxGC, annotationStyle, genome, seed=seed)
    }
    else
    {
      print("Bias file is not specified, observed bias is used")
      fragSubDt <- as.data.table(fragSubDt) 
    }
  }
  else
  {
    fragSubDt <- as.data.table(fragSubDt) 
  }
    
  if(simFLD)
  {
    # Vary Frag Dist
    fragSubDt <- .varyFragDist(fragSubDt, fracSub, nClust=nFragTypes,
                               fitGMM=fitGMM, prob=prob, 
                               estimateProb=estimateProb, seed=seed)
  }
  
  fragSubDt[,frag_id:=1:nrow(fragSubDt)]
  
  # Convert to GRanges and get GC content
  fragSubDt <- makeGRangesFromDataFrame(fragSubDt)
  seqlevelsStyle(fragSubDt) <- "UCSC"
  fragSubDt$gc_content <- Repitools::gcContentCalc(fragSubDt, organism=genome)
  
  return(fragSubDt)
}

#' Function to simulate a two condition ATAC seq experiment
#'
#' Single samples are varied by group specific parameters in a two group setting.
#' Fragments are sub-samples based on a GC-Bias, Fragment size distribution and
#' effect strength which can be varied (see varyAtacSignal). 
#' 
#' @param bamPaths a vector of the paths of Bam files to subsample from. 
#' @param chIPPeakDir path to the .Bed file with the chIP peaks
#' @param atacPeakDir path to the .bed file with the ATAC peaks
#' @param sampleNames a vector with the names of Bam files analyzed
#' @param design a vector of 1 & -1 indicating to which of the two groups a
#' sample belongs to
#' @param paramsGroup1 data.table with parameters for group 1 (see function: varyAtacSignal)
#' @param paramsGroup2 data.table with parameters for group 2 (see function: varyAtacSignal)
#' @param logFCs vector with log fold changes for each peak. Needs to be of the same length 
#' as the peaks in the bed file. 
#' @param which GRanges object stating which regions should be selected
#' @param effectStrength Factor to which the fragments of differentially accessible peaks should be
#' sub- / upsampled. 
#' @param nFragTypes Number of fragment types (nucleosome-free, mono-, di-, tri-nucleosome, etc. containing fragments)
#' @param minOverlap minimal overlap to consider a fragment being within a peak
#' @param minGC minimal gc content of fragments to be considered (to be removed)
#' @param maxGC maximal gc content of fragments to be considered (to be removed)
#' @param annotationStyle Either NCBI or UCSC
#' @param genome BSgenome object to be used. 
#' @return data.table with subsampled fragment coordinates & .bam file of these fragments saved on disk.
simAtacData <- function(bamPaths, 
                        chIPPeakDir, 
                        atacPeakDir,
                        sampleNames,
                        gcBiases, 
                        design,
                        paramsGroup1,
                        paramsGroup2, 
                        lfcDist,
                        which=NULL,
                        effectStrength=2,
                        nFragTypes=4,
                        minOverlap=1,
                        fracSub=1,
                        minGC=0,
                        maxGC=1,
                        simGCBias=TRUE,
                        simFLD=TRUE,
                        simEffectStrength=TRUE,
                        varyAtacPeaks=TRUE,
                        maxReadPerPeak=300,
                        matchCounts=FALSE,
                        annotationStyle="NCBI",
                        colNamesChIPPeaks=c("chr","start", "end", 
                                            "name", "score", "strand",
                                            "signalValue", "pValue", "qValue", "peak"),
                        colNamesAtacPeaks=c("chr","start", "end", 
                                            "name", "score", "strand",
                                            "signalValue", "pValue", "qValue", "peak"),
                        enrColChIPName="signalValue",
                        enrColAtacName="signalValue",
                        lfcCol="lfc",
                        equalLib=TRUE,
                        seed=42,
                        genome=BSgenome.Hsapiens.UCSC.hg38,
                        BPPARAM=SerialParam()){
  
  set.seed(seed)
  seqlevelsStyle(which) <- annotationStyle
  
  # get minimal library size
  depths <- lapply(bamPaths, function(bamPath){
  bamData <- .importFrags(bamPath, which, fracSub, annotationStyle)
  frags <- bamData$fragments
  length(frags)
  })
  
  minDepth <- min(unlist(depths))
  
  if(simEffectStrength){
    
  # import peaks 
  chIPPeaks <- .importPeaks(chIPPeakDir, which, 
                            annotationStyle, 
                            colNames=colNamesChIPPeaks)

  atacPeaks <- .importPeaks(atacPeakDir, which, 
                            annotationStyle,
                            colNames=colNamesAtacPeaks)
  
  # estimate ChIP-peak logFCs
  colnames(lfcDist) <- c(enrColChIPName, lfcCol)
  peakDt <- .estLfc(chIPPeaks, lfcDist,  enrCol=enrColChIPName, 
                    lfcCol=lfcCol, seed=seed)
  chIPlogFCs <- peakDt$lfc

  # estimate ATAC-peak logFCs
  setnames(lfcDist, enrColChIPName, enrColAtacName)
  atacPeakDt <- .estLfc(atacPeaks, lfcDist,  enrCol=enrColAtacName, 
                        lfcCol=lfcCol, seed=seed)
  atacLogFCs <- atacPeakDt$lfc
  }
  else{
    chIPPeaks=NULL
    chIPlogFCs=NULL
    
    atacPeaks=NULL
    atacLogFCs=NULL
  }
  
  # Positive samples
  posSamples <- bamPaths[which(design==1)]
  posSampleNames <- sampleNames[which(design==1)]

  if(simGCBias)
  {
    posGcBiases <- gcBiases[which(design==1)]
  }
  else
  {posGcBiases <- NULL}
  
  simSamplesPos <- bplapply(1:length(posSamples), function(i){
    
    data.table::setDTthreads(2)
    simData <- varyAtacSignal(bamPath=posSamples[i], 
                              peaks=chIPPeaks,
                              atacPeaks=atacPeaks,
                              biasFileDir=posGcBiases[i],
                              sampleName=posSampleNames[i],
                              chIPlogFCs=chIPlogFCs,
                              atacLogFCs=atacLogFCs*-1,
                              depth=minDepth,
                              effectStrength=effectStrength,
                              nFragTypes=nFragTypes,
                              fracSub=fracSub,
                              prob=c(paramsGroup1$prob_nf[1], 
                                     paramsGroup1$prob_mono[1], 
                                     paramsGroup1$prob_di[1], 
                                     paramsGroup1$prob_tri[1]),
                              estimateProb=paramsGroup1$estimateProb,
                              which=which,
                              minOverlap=minOverlap,
                              simGCBias=simGCBias,
                              simFLD=simFLD, 
                              simEffectStrength=simEffectStrength,
                              varyAtacPeaks=varyAtacPeaks,
                              maxReadPerPeak=maxReadPerPeak,
                              matchCounts=matchCounts,
                              minGC=minGC,
                              maxGC=maxGC,
                              annotationStyle=annotationStyle,
                              genome=genome,
                              seed=seed)
    simData$sample <- posSampleNames[i]
    simData$group <- paramsGroup1$name[1]
    simData <- as.data.table(simData)
    
    return(simData)
  }, 
  BPPARAM=BPPARAM)
  
  simSamplesPos <- rbindlist(simSamplesPos)
  
  # Negative samples
  negSamples <- bamPaths[which(design==-1)]
  negSampleNames <- sampleNames[which(design==-1)]
  
  if(simGCBias)
  {
    negGcBiases <- gcBiases[which(design==-1)]
  }
  else
  {negGcBiases <- NULL}
  
  simSamplesNeg <- bplapply(1:length(negSamples), function(i){
    data.table::setDTthreads(2)
    
    simData <- varyAtacSignal(bamPath=negSamples[i], 
                              peaks=chIPPeaks,
                              atacPeaks=atacPeaks,
                              biasFileDir=negGcBiases[i],
                              sampleName=negSampleNames[i],
                              chIPlogFCs=chIPlogFCs*-1,
                              atacLogFCs=atacLogFCs,
                              depth=minDepth,
                              effectStrength=effectStrength,
                              nFragTypes=nFragTypes,
                              fracSub=fracSub,
                              prob=c(paramsGroup2$prob_nf[1], 
                                     paramsGroup2$prob_mono[1], 
                                     paramsGroup2$prob_di[1], 
                                     paramsGroup2$prob_tri[1]),
                              estimateProb=paramsGroup2$estimateProb,
                              which=which,
                              minOverlap=minOverlap,
                              simGCBias=simGCBias,
                              simFLD=simFLD,
                              simEffectStrength=simEffectStrength,
                              varyAtacPeaks=varyAtacPeaks,
                              maxReadPerPeak=maxReadPerPeak,
                              matchCounts=matchCounts,
                              minGC=minGC,
                              maxGC=maxGC,
                              annotationStyle=annotationStyle,
                              genome=genome,
                              seed=seed)
    
    simData$sample <- negSampleNames[i]
    simData$group <- paramsGroup2$name[1]
    simData <- as.data.table(simData)
    
    return(simData)
  }, 
  BPPARAM=BPPARAM)
  
  simSamplesNeg <- rbindlist(simSamplesNeg)
  simSamples <- rbind(simSamplesPos, simSamplesNeg)
  
  # even sequencing depth
  if(equalLib)
  {
    minDepth <- min(simSamples[,.N,by=sample]$N)
    simSamples <- simSamples[,.SD[sample(.N, minDepth)],by = sample]
  }
  # adding dummy cols 
  #simSamples$strand <- "*"
  #simSamples$name <- "."
  #simSample$score <- 0L
    
  return(list(sim=simSamples, lfc=chIPlogFCs, depths=minDepth))
}


saveToBed <- function(simData, 
                      outDir,
                      folderName,
                      saveMerged=TRUE){
  
    if(!dir.exists(file.path(outDir,folderName))) dir.create(file.path(outDir,folderName), recursive=TRUE)
    simDataPath <-  file.path(outDir, folderName, "simfrags.rds")
    saveRDS(simData, simDataPath)
    
    # retrieve the simulated fragments
    simFrags <- simData[[1]]

    # bed file format: chrom, start, end, name
      bedSim <- simFrags[,c("seqnames", "start", "end", "sample", "gc_content", "strand"), with=FALSE]
      bedSim$strand <- "."
      bedSim$score <- 0L
      bedSim$name <- "."
      
      if(saveMerged){
      bedSimMerged <- copy(bedSim)
      bedSimMerged$sample <- "merged"
      colnames(bedSimMerged) <- c("chrom", "chromStart", "chromEnd", "sample", "gc_content", "name", "score", "strand")

      # save merged bed file
      mergedFilePath <- file.path(outDir, folderName, "mergedSamples.bed")
      setorder(bedSimMerged, chrom, chromStart, chromEnd)
      write.table(bedSimMerged[,c("chrom", "chromStart", "chromEnd", "name", "score", "strand"), with=FALSE],
                  mergedFilePath,
                  quote=FALSE, col.names=FALSE,
                  sep="\t",
                 row.names=FALSE)}

    # save .bed file by sample
    bedSim <- split(bedSim, by="sample")
    lapply(names(bedSim),function(sample){
      setorder(bedSim[[sample]], seqnames, start, end)
      singleFilePaths <- file.path(outDir, folderName, paste0(sample, ".bed"))
      write.table(bedSim[[sample]][,c("seqnames", "start", "end", "name", "score", "strand"), with=FALSE],
                  singleFilePaths,
                  quote=FALSE, col.names=FALSE,
                  sep="\t",
                  row.names=FALSE)})
}