# fGSEA

rungsea <- function(DAR,
                    peakset,
                    peaks)
{
  peak_FC <- jitter(DAR$logFC, 0.001) # To avoid non-unique breaks
  values(peaks) <- DataFrame(peak_FC)
  peakdf <- as.data.frame(peaks)
  peakdf$cat <- str_c(peakdf$seqnames, sep = ":", peakdf$start)
  peakdf$cat <- str_c(peakdf$cat, sep = "-", peakdf$end)
  peak_ranks <- setNames(peakdf$peak_FC, peakdf$cat)
  peak_ranks
  
  ptm <- proc.time()
  
  gseares <- fgsea::fgseaMultilevel(peakset, 
                                    peak_ranks,
                                    eps = 0,
                                    nPermSimple = 10000)
  gseares <- gseares[order(gseares$padj),]
  
  runtime <- proc.time()-ptm
  
  return(list(res=gseares, runtime=runtime))
}