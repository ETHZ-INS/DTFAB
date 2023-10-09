# viper

# dATestedgeR by E. Sonder

dATestedgeR <- function(counts_control, 
                        counts_perturbed)
{
  # literature: https://www.nature.com/articles/s41598-020-66998-4#data-availability
  # code from: https://github.com/Zhang-lab/ATACseq_benchmarking/blob/master/simulated_ATAC_tests_edgeR.R
  
  cond <- c(rep("ctrl", ncol(counts_control)), 
            rep("mut", ncol(counts_perturbed)))
  se <- cbind(counts_control, counts_perturbed)
  
  modelMat <- model.matrix(~cond)
  dge <- DGEList(counts=assays(se)$counts)
  dge <- calcNormFactors(dge)
  
  dge <- estimateGLMCommonDisp(dge, modelMat)
  dge <- estimateGLMTrendedDisp(dge, modelMat) 
  dge <- estimateGLMTagwiseDisp(dge, modelMat)
  fit <- glmFit(dge, modelMat)
  fit <- glmLRT(fit, coef="condmut")
  res <- topTags(fit, n=nrow(dge), sort="none")$table
  if(!is.null(rowData(counts_control)$bias)) res$bias <- rowData(counts_control)$bias
  rownames(res) <- as.character(rowRanges(se))
  
  return(res)
}


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