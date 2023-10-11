# taken from enrichMiR, with slight modifications
regreg <- function(signal, matchMatrix, use.intercept=TRUE, minSize=0, binary=NULL, alpha=0){
  signal <- signal[names(signal) %in% row.names(matchMatrix)]
  matchMatrix <- matchMatrix[names(signal),]
  matchMatrix <- matchMatrix[,which(colSums(matchMatrix!=0L)>=minSize)]
  allTFs <- colnames(matchMatrix)
  suppressPackageStartupMessages({
    library(glmnet)
  })
  res1 <- ulm.ebayes(signal, matchMatrix, use.intercept=use.intercept, minSize=minSize)
  sig <- head(row.names(res1), max(c(20,sum(res1$p_value<0.1))))
  matchMatrix <- matchMatrix[,sig]
  
  
  # regularized regression with cross-validation
  if(isLogistic <- is.logical(signal)){
    fits <- cv.glmnet(matchMatrix, signal, standardize=FALSE, alpha=alpha, 
                      family="binomial", lower.limits=0, intercept=use.intercept)
  }else{
    fits <- cv.glmnet(matchMatrix, signal, standardize=FALSE, alpha=alpha, 
                      family="gaussian", intercept=use.intercept)
  }
  
  l <- fits$lambda.1se
  if(fits$nzero[fits$lambda==l]<2){
    w1 <- which(fits$lambda==fits$lambda.min)
    w2 <- which(fits$lambda==fits$lambda.1se)
    l <- fits$lambda[floor((w1-w2)/2+w2)]
  }
  if(fits$nzero[fits$lambda==l]<2) l <- fits$lambda.min
  co <- coef(fits, s=l)
  co[co[,1]!=0,,drop=FALSE]
  
  beta <- co[,1]
  names(beta) <- row.names(co)
  
  signald <- data.frame( y=signal, as.matrix(matchMatrix[,intersect(colnames(matchMatrix),names(beta)),drop=FALSE]) )
  
  # new fit to get significance estimates
  if(use.intercept){
    form <- y~.
  }else{
    form <- y~0+.
  }
  res <- tryCatch({
    mod <- lm( form, data=signald )
    # we extract the coefficients and p-values, and reorganize the output:
    res <- coef(summary(mod))
    res[order(res[,4]),,drop=FALSE]
  }, error=function(e){
    warning(e)
    data.frame(row.names=names(beta), beta=beta, z=rep(NA_real_, length(beta)), pvalue=1)
  })
  colnames(res) <- c("score","stderr","t","p_value")
  res <- res[grep("^\\(Intercept\\)$|FALSE$", row.names(res), invert=TRUE),c(1,4),drop=FALSE]
  row.names(res) <- gsub("TRUE","",row.names(res))
  
  missingTFs <- setdiff(allTFs, row.names(res))
  missingTFs <- data.frame(row.names=missingTFs, score=rep(0,length(missingTFs)),
                           p_value=rep(1,length(missingTFs)))
  res <- rbind(res,missingTFs)
  res <- res[row.names(res)!="median",]
  # res <- DataFrame(res)
  # we adjust using all features as number of comparisons
  if(nrow(res)>0){
    res$FDR <- p.adjust(res$p_value, n=max(nrow(res),length(allTFs),na.rm=TRUE))
  }
  res
}

runregreg <- function(DARmat,
                      matchMtx, alpha=1){  
  
  ptm <- proc.time()
  regres <- regreg(DARmat, matchMtx, alpha=alpha)
  runtime <- proc.time()-ptm
  
  return(list(res=regres, runtime=runtime))
}
