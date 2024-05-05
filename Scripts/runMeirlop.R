runMeirlop <- function(path){
  ptm <- proc.time()
  system(paste0("meirlop --jobs 10 --fa ", path, "/others/scored.fasta ", path,
                "/others/motifs.jaspar ", path, "/raw/meirlop.out"))
  e <- read.delim(paste0(path, "/raw/meirlop.out/lr_results.tsv"), header=TRUE, row.names=1)
  row.names(e) <- gsub(" .*","",row.names(e))
  e <- e[order(e$pval, -abs(e$coef)),]
  res <- data.frame(row.names=row.names(e), logFC=e$coef, p=e$pval, padj=e$padj)
  res$rank <- seq_len(nrow(res))
  runtime <- proc.time()-ptm
  raw <- list(runtime=runtime, runtime2=runtime, obj1=e)
  return(list(raw=raw, res=res))
}
