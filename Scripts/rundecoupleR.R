rundecoupleR <- function(counts, genome, motifs, DAR,
                         norm=c("TMM","GCSmoothQuantile"),
                         mode=c("logFC", "sample"),
                         decoupleR_modes=c("mlm", "ulm", "norm_wsum")){

  norm <- match.arg(norm)
  mode <- match.arg(mode)
  
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
    network <- data.frame(source=factor(colnames(matchingMtx), 
                                        colnames(matchingMtx))[matchingMtx@j+1L], 
                          target=factor(row.names(matchingMtx), 
                                        row.names(matchingMtx))[matchingMtx@i+1L])
  }
  network$mor <- 1L
  network$target <- sub(":\\+$", "", network$target)
  
  # remove colinearity
  network <- .removeColinearity(network)
  
  ptm <- proc.time()
  if(mode=="sample"){
    res <- .decoupleR.sample(counts, network, decoupleR_modes, norm=norm)
  }else{
    # Generate required matrix of logFCs
    DARmat <- as.matrix(DAR$logFC)
    row.names(DARmat) <- rownames(DAR)
    res <- .decoupleR.logFC(DARmat, network, decoupleR_modes)
  }
  runtime2 <- sapply(split(res$raw$res$statistic_time,
                           res$raw$res$statistic), unique)
  runtime2["consensus"] <- sum(runtime2)
  
  res$raw$runtime <- proc.time()-ptm
  res$raw$runtime2=runtime2
  return(res)
}

.decoupleR.logFC <- function(DARmat, network, decoupleR_modes){
  if(is.null(colnames(DARmat))) colnames(DARmat) <- "A"
  dcplR <- decouple(DARmat, network, .source='source', .target='target',
                    minsize=5, include_time=TRUE, consensus_score=TRUE,
                    statistics=decoupleR_modes)
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

  return(list(raw=list(res=dcplR), res=res))
}

.decoupleR.sample <- function(se, network, decoupleR_modes, norm="TMM"){
  se$group <- rep(LETTERS[1:2],each=ncol(se)/2)
  if(norm=="TMM"){
    y <- log1p(cpm(calcNormFactors(DGEList(assay(se)))))
  }else{
    y <- GCSmoothQuantile(se, bio="group")
  }
  dcplR <- decouple(y, network, .source='source', .target='target',
                    minsize=5, include_time=TRUE, consensus_score=TRUE,
                    statistics=decoupleR_modes)
  res <- lapply(split(dcplR, dcplR$statistic), FUN=function(dc){
    m <- reshape2::dcast(dc, source~condition, value.var="score")
    row.names(m) <- m[,1]
    m <- as.matrix(m[,-1])
    m <- m[apply(m,1,FUN=sd)>0,]
    fit <- eBayes(lmFit(m, model.matrix(~se$group)))
    res <- topTable(fit, number=Inf)
    data.frame(row.names=row.names(res), logFC=res$logFC, t=res$t,
               p=res$P.Value, padj=res$adj.P.Val, rank=seq_len(nrow(res)))
  })
  return(list(raw=list(res=dcplR), res=res))
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
