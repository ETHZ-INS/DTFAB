# Second flavor of monaLisa: Regression Based Approach for Motif Selection

runStabSel <- function(DAR,
                         motifs,
                         peaks,
                         genome){
  ptm <- proc.time()
  
  # Adding peak logFCs to the peaks
  motifs <- do.call(TFBSTools::PWMatrixList,mapply(mo=motifs, n=names(motifs), FUN=function(mo,n){ mo@name<-n; mo}))
  peak_FC <- DAR$logFC
  
  values(peaks) <- DataFrame(peak_FC)
  
  # Obtaining the sequences for differentially accessible regions
  
  DARseqs <- getSeq(genome, peaks)
  if(!is.null(names(peaks))) names(DARseqs) <- names(peaks)
  
  # Finding motif hits within the differentially accessible regions.
  
  suppressWarnings({
    hits <- findMotifHits(query = motifs,
                        subject = DARseqs,
                        min.score = 10.0,
                        BPPARAM = MulticoreParam(8))
  })

  TFBSmatrix <- unclass(table(factor(seqnames(hits), levels = seqlevels(hits)),
                            factor(hits$pwmname, levels = names(motifs))))

  zero_TF <- colSums(TFBSmatrix) == 0
  sum(zero_TF)

  TFBSmatrix <- TFBSmatrix[, !zero_TF]

  fMono <- oligonucleotideFrequency(DARseqs, width = 1L, as.prob = TRUE)
  fDi <- oligonucleotideFrequency(DARseqs, width = 2L, as.prob = TRUE)
  fracGC <- fMono[, "G"] + fMono[, "C"]
  oeCpG <- (fDi[, "CG"] + 0.01) / (fMono[, "G"] * fMono[, "C"] + 0.01)

  # add GC and oeCpG to predictor matrix
  TFBSmatrix <- cbind(fracGC, oeCpG, TFBSmatrix)
  TFBSmatrix[1:6, 1:6]

  set.seed(123)
  se <- randLassoStabSel(x = TFBSmatrix, y = peaks$peak_FC,
                       cutoff = 0.9)

  Stab_plot <- plotStabilityPaths(se)
  
  runtime <- proc.time()-ptm
  return(list(res=se, runtime=runtime, obj1=Stab_plot))
}