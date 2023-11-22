# chromVAR

register(MulticoreParam(8))

runchromVAR <- function(counts, 
                        genome, 
                        motifs,
                        design, doNorm=FALSE, TMM = FALSE){
  # Adding GC-bias
  ptm <- proc.time()
  
  counts <- addGCBias(counts, 
            genome)
  
  # for downstream analyses...
  counts$meanGC <- sapply(as.data.frame(as.matrix(assay(counts))), FUN=function(x){
    weighted.mean(rowData(counts)$bias, w=x)
  })
  gcBins <- c(0,0.4,0.45,0.48,0.51,0.54,0.57,0.6,0.65,0.7,1)
  gcBins <- cut(rowData(counts)$bias, breaks=gcBins, labels=FALSE, include.lowest=TRUE)
  counts$GCcor <- sapply(as.data.frame(as.matrix(assay(counts))), FUN=function(x){
    x <- split(x,gcBins)
    cor(sapply(x,FUN=mean),seq_along(x))
  })
  
  # Extract non-overlapping peaks
  
  if(is.null(row.names(counts))) row.names(counts) <- as.character(granges(counts))
  counts_filtered <- filterPeaks(counts, 
                                 non_overlapping = TRUE)
  

  # check if match matrix already given, otherwise
  #  Finding which peaks contain which motifs
  if(is(motifs, "sparseMatrix")){
    motif_matrix <- motifs[row.names(counts_filtered),]
  }else{
    motif_matrix <- matchMotifs(motifs, 
                            counts_filtered, 
                            genome)
  }
  
  # Compute deviations
  
  dev <- computeDeviations(object = counts_filtered, 
                           annotations = motif_matrix,
                           expectation = computeExpectations(counts_filtered),
                           background_peaks=getBackgroundPeaks(counts_filtered, 
                                                               niterations=2000))
  
  # Compute variability
  
#  variability <- computeVariability(dev)
  # variability_plot <- plotVariability(variability, use_plotly = FALSE)
  
  # variability_plot
  
  # Limma analysis
  
    # First, the normalization of z-scores
  if(doNorm){
    assays(dev)$norm <- scale(assays(dev)$z)
    devMat <- assays(dev)$norm
  }else{
    devMat <- assays(dev)$z
  }
  
  # Estimate the fold changes and standard errors by fitting a linear model for each gene
  
  fit <- lmFit(devMat, 
               design)
  
  # Apply empirical Bayes smoothing to the standard errors
  
  fit <- eBayes(fit)
  
  # Show statistics for the top 10 genes
  
  topTFs <- topTable(fit, number = Inf)
  
  # List Readouts
  # topTFs
  
  runtime <- proc.time()-ptm
  
  return(list(res=topTFs, runtime=runtime, obj1=dev, obj2=counts))
}
