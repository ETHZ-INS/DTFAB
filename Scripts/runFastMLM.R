#' fastMLM
#'
#' @param accmat A matrix of peak (rows) accessibility per sample/cell (columns)
#' @param annotation A matrix of motif (columns) matches per peak (rows)
#' @param useIntercept Logical; whether to use an intercept in the mode
#' @param poisson Logical; whether to use poisson regression (assumes `accmat`
#'   is a count matrix).
#' @param BPPARAM BiocParallel param for multithreading.
#'
#' @return A matrix (of `ncol(annotation)` rows and `ncol(accmat)` columns) with
#'   the activity scores (model t-values) of each motif in each sample.
#'   
#' @details 
#' Regresses each column of `accmat` on `annotation`, and uses the coefficients'
#' t-values as activity scores.
#' 
#' @export
fastMLM <- function(accmat, annotation, useIntercept=TRUE, poisson=FALSE,
                    BPPARAM=BiocParallel::SerialParam()){
  stopifnot(is.matrix(accmat) && is.matrix(annotation))
  if(useIntercept) annotation <- cbind(rep(1L,nrow(annotation)), annotation)
  res <- bplapply(seq_len(ncol(accmat)), BPPARAM=BPPARAM, FUN=function(i){
    if(!isTRUE(poisson)){
      mod <- RcppArmadillo::fastLmPure(annotation, accmat[,i])
      tvals <- mod$coefficients/mod$stderr
    }else{
      if(FALSE && require("Rfast", quietly=TRUE, include.only="glm_poisson")){
        mod <- glm_poisson(a, y[,1], full=TRUE)$info
      }else{
        mod <- glm(accmat[,i]~0+annotation, family="poisson")
        mod <- coef(summary(mod))
      }
      tvals <- mod[,1]/mod[,2]
    }
    tvals
  })
  res <- matrix(unlist(res), nrow=ncol(annotation))
  row.names(res) <- colnames(annotation)
  if(useIntercept) res <- res[-1,,drop=FALSE]
  colnames(res) <- colnames(accmat)
  res
}

runFastMLM <- function(counts, momat, design=NULL, intercept=TRUE,
                       BPPARAM=BiocParallel::MulticoreParam(4)){
  if(is.null(design)){
    design <- rep(LETTERS[1:2],each=ncol(counts)/2)
  }
  counts$condition <- design
  ptm <- proc.time()
  e <- GCSmoothQuantile(counts, bio="condition")
  e <- fastMLM(e, as.matrix(momat), useIntercept=intercept, BPPARAM=BPPARAM)
  fit <- eBayes(lmFit(e, model.matrix(~design)))
  res <- topTable(fit, number=Inf)
  res <- data.frame(row.names=row.names(res), logFC=res$logFC, t=res$t,
                    p=res$P.Value, padj=res$adj.P.Val)
  res <- res[order(res$p, -abs(res$logFC)),]
  res$rank <- seq_len(nrow(res))
  runtime <- proc.time()-ptm
  raw <- list(runtime=runtime, runtime2=runtime, obj1=e)
  return(list(raw=raw, res=res))
}
