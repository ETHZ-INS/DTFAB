# viper


runviper <- function(counts_control, 
                     counts_perturbed,
                     DAR,
                     peaks,
                     regulons,
                     mode=c("msviper", "viper"),
                     method=c("ttest", "none"),
                     binarizeRegulon=FALSE){
  
  if(binarizeRegulon){
    regulons <- lapply(regulons, FUN=function(x){
      x$likelihood <- pmax(1L,x$likelihood)
      x
    })
  }
  
  cond <- c(rep("ctrl", ncol(counts_control)), 
            rep("mut", ncol(counts_perturbed)))
  se <- cbind(counts_control, counts_perturbed)
  
  mode <- match.arg(mode, c("msviper", "viper"))
  if(mode=="viper")
  {
    # vst transform (important that we don't scale later in the viper method)
    transCounts <- DESeq2::vst(as.matrix(assays(se)$counts))
    rownames(transCounts) <- as.character(rowRanges(se))
    ptm <- proc.time()
    vpres <- viper(transCounts, regulons, method=method, verbose=TRUE) # method="none"
    runtime <- proc.time()-ptm
    return(list(vpres, runtime))
  }
  else if(mode=="msviper")
  {
    signature <- DAR$logFC # not sure whats the best choice for signature: lets discuss
    names(signature) <- rownames(DAR)
    
    ptm <- proc.time()
    mres <- msviper(signature, regulons, nullmodel=NULL, verbose=TRUE)
    runtime <- proc.time()-ptm
    return(list(res=mres, runtime=runtime))
  }
}