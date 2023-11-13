# dATestedgeR by E. Sonder

dATestedgeR <- function(counts_control, counts_perturbed,
                        norm.method=c("TMM","GCQuantileNorm","GCSmoothQuantile")){
  # literature: https://www.nature.com/articles/s41598-020-66998-4#data-availability
  # code from: https://github.com/Zhang-lab/ATACseq_benchmarking/blob/master/simulated_ATAC_tests_edgeR.R
  
  norm.method <- match.arg(norm.method)

  cond <- c(rep("ctrl", ncol(counts_control)), 
            rep("mut", ncol(counts_perturbed)))
  se <- cbind(counts_control, counts_perturbed)
  if(!inherits(se,"SummarizedExperiment")) se <- SummarizedExperiment(list(counts=se))
  se$group <- cond

  modelMat <- model.matrix(~cond)
  
  dge <- switch(norm.method,
    TMM=calcNormFactors(DGEList(counts=assays(se)$counts)),
    GCQuantileNorm=DGEList(GCQuantileNorm(se)),
    GCSmoothQuantile=DGEList(GCSmoothQuantile(se, bio="group"))
  )
  
  dge <- estimateGLMCommonDisp(dge, modelMat)
  dge <- estimateGLMTrendedDisp(dge, modelMat) 
  dge <- estimateGLMTagwiseDisp(dge, modelMat)
  fit <- glmFit(dge, modelMat)
  fit <- glmLRT(fit, coef="condmut")
  res <- topTags(fit, n=nrow(dge), sort="none")$table
  if(!is.null(rowData(counts_control)$bias)) res$bias <- rowData(counts_control)$bias
  
  return(res)
}

