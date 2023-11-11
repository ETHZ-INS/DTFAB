runATACseqTFEA <- function(bams, pmoi){
  ns <- length(bams)
  pmoi$motif_alt_id <- pmoi$pvalue <- pmoi$qvalue <- pmoi$matched_sequence <- NULL
  colnames(mcols(pmoi)) <- gsub("motif_id","motif",colnames(mcols(pmoi)))
  ptm <- proc.time()
  res <- TFEA(bams[(ns/2)+1:(ns/2)], bams[1:(ns/2)], bindingSites=pmoi)
  res2 <- res$resultsTable
  colnames(res2)[4:5] <- c("p","padj")
  set.seed(123)
  res2 <- res2[order(res2$p, -abs(res2$normalizedEnrichmentScore),
                     sample.int(nrow(res2))),]
  res2$rank <- seq_len(nrow(res2))
  runtime <- proc.time()-ptm
  
  return(list(res=res2, runtime=runtime))
}

