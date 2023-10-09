# First flavor of monaLisa: MOtif aNAlysis with Lisa

runmonaLisa <- function(DAR,
                        motifs,
                        peaks,
                        genome){
  
  # At first a normalization should be done
  
  # Compute FC which --> later appended to the GR obj.
  
  peak_FC <- jitter(DAR$logFC, 0.001) # To avoid non-unique breaks
  
  values(peaks) <- DataFrame(peak_FC)

  fc2 <- peak_FC[which(abs(peak_FC)>=0.3)]
  nElements <- ceiling(length(fc2)/11)
  
  
  
  ptm <- proc.time()
  bins <- monaLisa::bin(x = peaks$peak_FC, 
                        binmode = "equalN", 
                        nElements = nElements,
                        minAbsX = 0.3)
  
  BinDensity <- plotBinDensity(peaks$peak_FC,
                bins)

  DARseqs <- getSeq(genome, 
                    peaks) # This is how they got sequences. I will also try it with our sequences from PLs motif prep.
  
  # Checking Biases via visualization
  
  # GC_fraction <- plotBinDiagnostics(seqs = DARseqs, 
  #                    bins = bins, 
  #                    aspect = "GCfrac")
  # 
  # Dinucl_fraction <- plotBinDiagnostics(seqs = DARseqs, 
  #                    bins = bins, 
  #                    aspect = "dinucfreq")
  
  # Motif enrichment analysis
  
  se <- calcBinnedMotifEnrR(seqs = DARseqs, 
                            bins = bins, 
                            pwmL = motifs, 
                            BPPARAM = BiocParallel::MulticoreParam(8),
                            background = c("zeroBin"))
  
  runtime <- proc.time()-ptm
  return(list(res=se, runtime=runtime))
  
  # This is commented out since it is currently not relevant.
  #
  # sel <- apply(assay(se, "negLog10Padj"), 
  #              1, 
  #              function(x) max(abs(x),
  #                              0, 
  #                              na.rm = TRUE)) > 1.3
  # 
  # seSel <- se[sel, ]
  # 
  # enrich_plot <- plotMotifHeatmaps(x = seSel, 
  #                   which.plots = c("log2enr", 
  #                                   "negLog10Padj"), 
  #                   width = 2.0, 
  #                   cluster = TRUE, 
  #                   maxEnr = 2, 
  #                   maxSig = 10, 
  #                   show_motif_GC = TRUE)
  # SimMatSel <- motifSimilarity(rowData(seSel)$motif.pfm)
  # 
  # # Create hclust object, similarity defined by 1 - Pearson correlation
  # 
  # hcl <- hclust(as.dist(1 - SimMatSel), 
  #               method = "average")
  # hcl_plot <- plotMotifHeatmaps(x = seSel, 
  #                   which.plots = c("log2enr", 
  #                                   "negLog10Padj"), 
  #                   width = 1.8, 
  #                   cluster = hcl, 
  #                   maxEnr = 2, 
  #                   maxSig = 10,
  #                   show_dendrogram = TRUE,
  #                   show_seqlogo = TRUE,
  #                   width.seqlogo = 1.2)
  #  
  # # Selecting readout
  # 
  # return(list("BinDensity", 
  #             "GC_fraction", 
  #             "Dinucl_fraction", 
  #             "enrich_plot", 
  #             "hcl_plot"))
}

