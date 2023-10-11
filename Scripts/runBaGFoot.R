#' getModelBasedActivity
#' 
#' @param x A named vector of bam files.
#' @param pmoi A GRanges of motif instances, with a 'motif_id' column.
#' @param wSize Positive integer scaler of the window size around motifs. (200 is bagfoot's)
#' @param paired Logical; whether or not the data is paired.
#' @param peaks A GRanges of peaks (analyses will be restricted to that). If 
#'   not given, signal around all motif instances in `pmoi` is used.
#' @param shift How to shift the reads (default ATAC shift)
#' @param minFragLength Minimum fragment length to consider
#' @param maxFragLength Maximum fragment length to consider
#' @param BPPARAM A BiocParallel param object for multithreading
#'
#' @return A summarized experiment containing the motif scores.
#' @author Pierre-Luc Germain
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom BiocParallel SerialParam bplapply
#' @import S4Vectors IRanges GenomicRanges
#' @export
getModelBasedActivity <- function(x, pmoi, paired, peaks=NULL, wSize=200,
                                  shift=c(4L,-5L), minFragLength=30L,
                                  maxFragLength=2000,
                                  BPPARAM=BiocParallel::SerialParam(progress=TRUE)){
  if(is.null(peaks)){
    peaks <- GenomicRanges::reduce(pmoi)
  }else{
    pmoi <- pmoi[overlapsAny(pmoi, peaks)]
  }
  onlyReg <- GenomicRanges::reduce(resize(peaks, width=width(peaks)+(2*wSize), fix="center"))
  if(is.null(names(x))) names(x) <- make.names(basename(x), unique=TRUE)
  message("Reading in the data...")
  x <- pbapply::pblapply(x, FUN=function(x){
    if(grepl("\\.bam$", x, ignore.case=TRUE)){
      co <- epiwraps::bam2bw(x, paired=paired, output_bw=NA, only=onlyReg,
                             binWidth=1L, type=ifelse(paired,"ends","start"), shift=shift,
                             minFragLength=minFragLength, maxFragLength=maxFragLength,
                             scaling=FALSE, verbose=FALSE)
    }else{
      co <- epiwraps::frag2bw(x, paired=paired, output_bw=NA, only=onlyReg,
                              binWidth=1L, type=ifelse(paired,"ends","start"), shift=shift,
                              minFragLength=minFragLength, maxFragLength=maxFragLength,
                              scaling=FALSE, verbose=FALSE)
    }
    co[seqlevels(pmoi)]
  })
  pmoi <- split(pmoi, pmoi$motif_id)
  message("Computing motifs signals...")
  res <- bplapply(pmoi, BPPARAM=BPPARAM, FUN=function(m){
    motifSize <- width(m)[[1]]
    m <- resize(m, width=2*wSize, fix="center")
    vs <- lapply(x, FUN=function(x){
      v <- Views(x, m)
      suppressWarnings(as.integer(colSums(.views2Matrix(v, padVal=0L))))
    })
    rm(m)
    global <- colSums(do.call(rbind,vs))
    # to avoid seq bias, use the median coverage as weight for inside the footprint
    moRange <- round(wSize+c(-1,1)*motifSize/2)
    global[moRange] <- round(median(global[moRange]))
    bgP <- mean(global[c(1:2,length(global)-1,length(global))])
    global <- (1L+global)/(1L+bgP)
    list(scores=sapply(vs, FUN=function(x){
      sum(x*global)
    }), global=global)
  })
  global <- do.call(rbind, lapply(res, FUN=function(x) t(x$global)))
  res <- do.call(rbind, lapply(res, FUN=function(x) x$scores))
  se <- SummarizedExperiment(list(scores=res))
  S4Vectors::metadata(se)$globalProfiles <- global
  rowData(se)$nsites <- lengths(pmoi)
  rowData(se)$motifSize <- sapply(pmoi,FUN=function(x) width(x)[1])
  se
}

#' getBagFootLike
#' 
#' Get motif*sample insertion counts within (footpring) and flanking each motif 
#'
#' @param x A named vector of bam files.
#' @param pmoi A GRanges of motif instances, with a 'motif_id' column.
#' @param fSize Positive integer scalar of the flanking window size around motifs. (200 is bagfoot's)
#' @param fdSize Positive integer scalar of the flanking used to calculate footprint depth.
#' @param paired Logical; whether or not the data is paired.
#' @param peaks A GRanges of peaks (analyses will be restricted to that)
#' @param shift How to shift the reads (default ATAC shift)
#' @param BPPARAM A BiocParallel param object for multithreading
#'
#' @return A summarized experiment containing the motif scores.
#' @author Pierre-Luc Germain
#' @importFrom SummarizedExperiment SummarizedExperiment rowData<-
#' @importFrom BiocParallel SerialParam bplapply
#' @import S4Vectors IRanges GenomicRanges
#' @export
getBagFootLike <- function(x, 
                           pmoi, 
                           paired, 
                           peaks=NULL, 
                           fSize=200L, 
                           fdSize=20L, 
                           shift=c(4L,-5L), 
                           BPPARAM=BiocParallel::SerialParam(progress=TRUE)){
  if(is.null(peaks)) peaks <- GenomicRanges::reduce(pmoi)
  onlyReg <- GenomicRanges::reduce(resize(peaks, width=width(peaks)+(2*fSize), fix="center"))
  if(is.null(names(x))) names(x) <- make.names(basename(x), unique=TRUE)
  pmoi <- split(pmoi, pmoi$motif_id)
  # eventually replace with a bamChrChunkApply version to reduce mem usage
  res <- lapply(x, FUN=function(x){
    message("Reading ", x)
    if(grepl("\\.bam$", x, ignore.case=TRUE)){
      x <- epiwraps::bam2bw(x, paired=paired, output_bw=NA, only=onlyReg,
                            binWidth=1L, type=ifelse(paired,"ends","start"), shift=shift,
                            scaling=FALSE, verbose=FALSE)
    }else{
      x <- epiwraps::frag2bw(x, paired=paired, output_bw=NA, only=onlyReg,
                             binWidth=1L, type=ifelse(paired,"ends","start"), shift=shift,
                             scaling=FALSE, verbose=FALSE)
    }
    x <- x[seqlevels(pmoi[[1]])]
    message("Computing motif profiles...")
    bplapply(pmoi, BPPARAM=BPPARAM, FUN=function(m){
      motifSize <- width(m)[[1]]
      moRange <- Reduce(":",round(fSize+c(-1,1)*motifSize/2))
      fdRange <- setdiff(Reduce(":",round(fSize+c(-1,1)*fdSize)),
                         moRange)
      m <- resize(m, width=2*fSize, fix="center")
      prof <- as.integer(colSums(.views2Matrix(Views(x, m), padVal=0L)))
      rm(m)
      moScore <- sum(prof[moRange])
      flScore <- sum(prof[fdRange])
      fdScore <- (1+flScore)/length(fdRange)-(1+moScore)/motifSize
      prof[moRange] <- 0L
      list(flankingArea=sum(prof), footprint=moScore, immediateFlanking=flScore,
           footprintDepth=fdScore/pmax(flScore, fdScore))
    })
  })
  names(n) <- n <- names(res[[1]][[1]])
  se <- SummarizedExperiment(lapply(n, FUN=function(n){
    x <- do.call(cbind, lapply(res, FUN=function(x){
      as.matrix(sapply(x, FUN=function(x) x[[n]]))
    }))
    colnames(x) <- names(res)
    x
  }))
  rowData(se)$nsites <- lengths(pmoi)
  rowData(se)$motifSize <- sapply(pmoi,FUN=function(x) width(x)[1])
  se
}


# converts a RleViews or RleViewsList with views of the same width to a matrix,
# setting out-of-bounds regions to `padVal`
.views2Matrix <- function(v, padVal=NA_integer_){
  if(!is(v, "RleViewsList")) v <- RleViewsList(v)
  ws <- width(v)[[1]][[1]]
  stopifnot(all(unlist(width(v))==ws))
  x <- Reduce(c, lapply(v[lengths(v)>0], padVal=padVal, FUN=.view2paddedAL))
  matrix(unlist(x), byrow=TRUE, ncol=ws)
}

# converts a RleViews to an AtomicList, setting out-of-bounds regions to padVal
.view2paddedAL <- function(v, padVal=NA_integer_, forceRetAL=TRUE){
  v2 <- trim(v)
  if(isInt <- is.integer(runValue(v2[[1]]))){
    stopifnot(is.integer(padVal))
  }else{
    stopifnot(is.numeric(padVal))
  }
  toAL <- function(v){
    if(isInt) return(IntegerList(v))
    NumericList(v)
  }
  if(all(width(v2)==width(v))){
    if(forceRetAL) return(toAL(v))
    return(v)
  }
  if(any(w <- width(v2)==0)){
    v <- v[which(!w)]
    v2 <- v2[which(!w)]
    warning(sum(w), " views were excluded as completely out of range.")
  }
  # figure out how much is trimmed on either side
  pleft <- start(v2)-start(v)
  pright <- end(v)-end(v2)
  # concatenate the list elements with their padding
  n <- seq_along(v2)
  v <- splitAsList(c( rep(padVal, sum(pleft)),
                      unlist(toAL(v2), use.names=FALSE),
                      rep(padVal, sum(pright))),
                   c(rep(n, pleft), rep(n, width(v2)), rep(n, pright)))
  names(v) <- names(v2)
  v
}

runBaGFoot <- function(readlist,
                       pmoi,
                       paired_arg){
  ptm <- proc.time()
  se.bf <- getBagFootLike(readlist,
                 pmoi,
                 paired_arg,
                 BPPARAM=MulticoreParam(6, progress=TRUE))
  runtime <- proc.time()-ptm
  groups <- rep(LETTERS[1:2], each=(length(readlist)/2))
  design <- model.matrix(~groups)
  # get common normalization factors
  dds <- calcNormFactors(DGEList(assays(se.bf)$footprint+assays(se.bf)$flanking))
  # differential analysis for the footprint depth
  fit <- eBayes(lmFit(assays(se.bf)$footprintDepth,design))
  resfo <- as.data.frame(topTable(fit, number=Inf))
  # differential analysis for the flanking
  dds$counts <- assays(se.bf)$flanking
  fit <- eBayes(lmFit(voom(dds, design),design))
  resfl <- as.data.frame(topTable(fit, number=Inf))
  # combine the two sets of results:
  m <- merge(resfo, resfl, by="row.names", suffix=c(".fo",".fl"))
  m$combined.P <- apply(m[,grep("P.Value",colnames(m))],1,FUN=aggregation::fisher)
  # w <- which(sign(m$logFC.fo)==sign(m$logFC.fl))
  # m$combined.P[w] <- 10^apply(log10(as.matrix(m[w,grep("P.Value",colnames(m))])),1,FUN=mean)
  m$FDR <- p.adjust(m$combined.P)
  row.names(m) <- m$Row.names
  head(m[order(m$combined.P),])
  return(list(m, runtime, se.bf))
}

runMBA <- function(readlist,
                   pmoi,
                   paired_arg){
  ptm <- proc.time()
  se.mba <- getModelBasedActivity(readlist, 
                                  pmoi, 
                                  paired_arg, 
                                  BPPARAM=MulticoreParam(6, progress=TRUE))
  runtime <- proc.time()-ptm
  dds <- calcNormFactors(DGEList(assay(se.mba)))
  groups <- rep(LETTERS[1:2],each=(length(readlist)/2))
  design <- model.matrix(~groups)
  fit <- eBayes(lmFit(voom(dds, design), design))
  topTFs <- topTable(fit, number = Inf)
  return(list(res=topTFs, runtime=runtime, obj1=se.mba))
}
