ulm.ebayes <- function(signal, matchMatrix, use.intercept=TRUE, minSize=5){
    if(is.data.frame(signal)) signal <- setNames(signal$logFC, row.names(signal))
    signal <- signal[names(signal) %in% row.names(matchMatrix)]
    matchMatrix <- matchMatrix[names(signal),]
    if(is(matchMatrix, "sparseMatrix")) matchMatrix <- as(matchMatrix, "dgCMatrix")
    matchMatrix <- matchMatrix[,which(colSums(matchMatrix!=0L)>=minSize)]
    if(use.intercept){
      modelmtx <- model.matrix(~signal)
    }else{
      modelmtx <- model.matrix(~0+signal)
    }
    lmfit <- lmFit(as.matrix(t(matchMatrix)), modelmtx)
    ulmres <- as.data.frame(topTable(eBayes(lmfit), coef="signal", number = Inf)[,c(1,4,5)])
    colnames(ulmres) <- c("score", "p", "padj")
    ulmres
}

ulmGC.ebayes <- function(DAR, matchMatrix, use.intercept=TRUE, minSize=5){
    if(is.null(DAR$bias)) stop("GC bias data not found!")
    signal <- as.numeric(DAR$logFC)
    names(signal) <- rownames(DAR)
    GC <- DAR$bias
    signal <- signal[names(signal) %in% row.names(matchMatrix)]
    matchMatrix <- matchMatrix[names(signal),]
    matchMatrix <- matchMatrix[,which(colSums(matchMatrix!=0L)>=minSize)]
    if(use.intercept){
      modelmtx <- model.matrix(~GC+signal)
    }else{
      modelmtx <- model.matrix(~0+GC+signal)
    }
    lmfit <- lmFit(as.matrix(t(matchMatrix)), modelmtx)
    ulmres <- as.data.frame(topTable(eBayes(lmfit), coef="signal", number = Inf)[,c(1,4,5)])
    colnames(ulmres) <- c("score", "p", "padj")
    ulmres
}


runulm <- function(DARmat,
                   matchMtx, GC=FALSE){

  ptm <- proc.time()
  if(GC){
    ulm <- ulmGC.ebayes(DARmat, matchMtx)
  }else{
    ulm <- ulm.ebayes(DARmat, matchMtx)
  }
  runtime <- proc.time()-ptm

  return(list(res=ulm, runtime=runtime))
}


