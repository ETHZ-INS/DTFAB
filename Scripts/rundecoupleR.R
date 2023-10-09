rundecoupleR <- function(counts,
                         genome,
                         motifs,
                         DAR,
                         decoupleR_modes=c("mlm", "ulm", "norm_wsum")){

  # Generate required matrix of logFCs
  
  DARmat <- as.matrix(DAR$logFC)
  row.names(DARmat) <- rownames(DAR)
  
  # Generate required network
  
  if(is(motifs, "GRanges")){
    o <- findOverlaps(granges(counts), motifs)
    network <- data.frame(source=as.character(motifs[to(o)]$motif_id),
                          target=row.names(counts)[from(o)])
    network <- network[!duplicated(network),]
  }else{
    matchingMotifs <- matchMotifs(motifs, counts, genome)
    row.names(matchingMotifs) <- as.character(granges(matchingMotifs))
    matchingMtx <- as(assay(matchingMotifs), "TsparseMatrix")
    network <- data.frame(source=factor(colnames(matchingMtx), colnames(matchingMtx))[matchingMtx@j+1L], 
                          target=factor(row.names(matchingMtx), row.names(matchingMtx))[matchingMtx@i+1L])
  }
  network$mor <- 1L
  network$target <- sub(":\\+$", "", network$target)
  
  # remove colinearity
  network <- .removeColinearity(network)
  
  ptm <- proc.time()
  if(is.null(colnames(DARmat))) colnames(DARmat) <- "A"
  dcplR <- decouple(DARmat, network, .source='source', .target='target',
                    minsize=5, include_time=TRUE, consensus_score=TRUE,
                    statistics=decoupleR_modes)
  runtime <- proc.time()-ptm
  runtime2 <- sapply(split(dcplR$statistic_time, dcplR$statistic), unique)
  runtime2["consensus"] <- sum(runtime2)
  
  dcplR <- dcplR[dcplR$condition=="A",]
  res <- lapply(split(as.data.frame(dcplR[,c("source","score","p_value")]),
                      dcplR$statistic), FUN=function(x){
    x <- x[order(x$p_value, abs(x$score)),]
    colnames(x)[3] <- "p"
    row.names(x) <- x$source
    x$source <- NULL
    x$rank <- seq_len(nrow(x))
    x$padj <- p.adjust(x$p)
    x
  })
  
  return(list(raw=list(res=dcplR, runtime=runtime, runtime2=runtime2),
              res=res))
}


.removeColinearity <- function(reg, minsize=5, truth=c(), cor.thres=0.96){
  # remove colinearity
  m <- reshape2::dcast(reg, target~source, value.var = "mor", fill = 0)
  row.names(m) <- m[,1]
  m <- as.matrix(m[,-1])
  m <- m[row.names(m),]
  m <- m[,colSums(m!=0)>=minsize]
  cc <- cor(m)
  cc[upper.tri(cc, diag=TRUE)] <- NA_real_
  removed <- vector("character")
  while(nrow(w <- which(cc>cor.thres, arr.ind=TRUE))>0){
    # if the true TF is colinear with something else, keep the true one
    isTrueReg <- intersect(row.names(w), truth)
    removed <- c(removed, setdiff(row.names(w),truth))
    if(length(isTrueReg)>0) removed <- c(removed, colnames(cc)[w[isTrueReg,2]])
    keep <- setdiff(row.names(cc),row.names(w))
    cc <- cc[keep,][,keep]
  }
  message(paste0(ncol(cc),"/",length(unique(reg$source))," regulons kept."))
  if(length(removed)>0)
    message("The following factors were removed due to collinearity with other factors:\n",
            paste(removed, collapse=", "))
  reg[reg$source %in% row.names(cc),]
}