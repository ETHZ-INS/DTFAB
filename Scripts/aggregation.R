runAggregations <- function(resList, methods=c("fisher","simes","stouffer","npPvalAgg")){
  motifs <- unique(unlist(lapply(resList, row.names)))
  pmat <- sapply(resList, FUN=function(x) x[motifs,"p"])
  row.names(pmat) <- motifs
  simes <- function(pval){ 
    pval <- pval[!is.na(pval)]
    min((length(pval)*pval[order(pval)])/seq(from = 1, to = length(pval), by = 1))
  }
  lapply(setNames(methods,methods), FUN=function(meth){
    if(meth=="npPvalAgg"){
      y <- npPvalAgg(pmat)
      row.names(y) <- y$TF
      return(y[order(y$rank),c("p","padj","rank")])
    }
    ag <- apply(pmat, 1, function(x){
      x <- x[!is.na(x)]
      if(length(x)==1) return(x)
      if(length(x)==0) return(NA_real_)
      switch(meth, simes=simes(x), fisher=metap::sumlog(x)$p,
             stouffer=metap::sumlog(x)$p, min=min(x))
    })
    res <- data.frame(row.names=row.names(pmat), p=ag,
                      padj=p.adjust(ag, method="fdr"))
    res <- res[order(res$p),]
    res$rank <- rank(res$p)
    res
  })
}

#' npPvalAgg
#'
#' @param x A matrix of unadjusted p-values (methods as columns, TFs as rows)
#' @param niter The number of iterations (increase to get lower possible p-values)
#' @param trans The transformation to be applied on the p-values, either rank,
#'   scale, or sqrtrank
#' @param fn the aggregation function to be applied on the transformed matrix
#'
#' @return A data.frame of p-values and adjusted p-values
npPvalAgg <- function(x, niter=2000, trans=c("sqrtrank","rank","logitscale"), fn=rowMeans){
  trans <- match.arg(trans)
  transfn <- switch(trans,
                    rank=rank,
                    logitscale=function(x) scale(gtools::logit(x)),
                    sqrtrank=function(x) sqrt(rank(x))
  )
  x <- as.data.frame(lapply(as.data.frame(x), FUN=transfn), row.names=row.names(x))
  ag <- fn(as.matrix(x), na.rm=TRUE)
  if(trans %in% c("rank","sqrtrank")){
    xnull <- as.data.frame(lapply(x, FUN=function(j){
      x <- sample.int(length(j), length(j)*niter, replace=TRUE)
      if(trans=="sqrtrank") x <- sqrt(x)
      x
    }))
  }else{
    xnull <- as.data.frame(lapply(x, FUN=function(j){
      as.numeric(j)[sample.int(length(j),niter*length(j),replace=TRUE)]
    }))
  }
  xnull <- fn(as.matrix(xnull), na.rm=TRUE)
  pag <- sapply(ag,FUN=function(x) sum(xnull<=x, na.rm=TRUE))/length(xnull)
  out <- data.frame(TF=row.names(x), padj=p.adjust(pag, method = "fdr"),p=pag)
  out <- out[order(out$p),]
  out$rank <- seq_along(row.names(out))
  out
}
