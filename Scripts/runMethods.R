#' List of (real) benchmark datasets
#'
#' @param onlyPE Logical; whether to include only paired-end datasets (default F)
#'
#' @return A named list of datasets, where the names also indicate the 
#'   corresponding subfolder. Each dataset is then also a list containing, 
#'   minimally, the `truth` (character vector of true perturbed TF(s)) and
#'   `species` (either h or m) slots. If the fragments are in bed/bed.gz format
#'   (instead of bam), also specify `readType="bed"`.
getDatasets <- function(onlyPE=FALSE){
  datasets <- list(
    BANP=list(truth="BANP", species="m", type="dTag"),
    ESR1=list(truth=c("ESR1","ESR2"), species="h", type="ligand"),
    GATA1=list(truth=c("GATA1"), species="h", readType="bed", type="deletion"),
    GATA2=list(truth=c("GATA2"), species="h", readType="bed", type="deletion"),
    RUNX1=list(truth=c("RUNX1"), species="h", readType="bed", type="deletion"),
    RUNX2=list(truth=c("RUNX2"), species="h", readType="bed", type="deletion"),
    KLF1=list(truth=c("KLF1"), species="h", readType="bed", type="deletion"),
    MYC=list(truth=c("MYC","MAX"), species="h", readType="bed", type="deletion"),
    NR1H3=list(truth=c("NR1H3"), species="m", type="ligand"),
    NR1H4=list(truth=c("NR1H4","RXRA","RXRB"), species="m", type="ligand"),
    NR3C1=list(truth=c("NR3C1","GCR"), species="h", type="ligand", paired=FALSE)
  )
  if(onlyPE) datasets$NR3C1 <- NULL
  return(datasets)
}

#' Gets the names of the methods to run
#'
#' @param onlyTop Logical, whether to get only the top methods
getMethods <- function(onlyTop=FALSE){
  if(onlyTop)
    return(c("chromVAR", "minaLisa.others", "monaLisa.zero", "StabSel", 
             "msVIPER", "MBA"))
  return(c( "chromVAR", "monaLisa", "StabSel", "GSEA", "decoupleR", 
            "VIPER", "VIPERb", "msVIPER", "msVIPERb",
            "ulm", "ulmB", "ulmGC", "regreg", "regregR",
            "BaGFoot", "MBA", "ATACseqTFEA" ))
}

#' Run all methods on all datasets
#'
#' @param datasets A named list of datasets, as produced by `getDatasets()`
#' @param ... Passed to `runMethods`
#'
#' @return Nothing (results saved to disk)
runAll <- function(datasets=getDatasets(), methods=getMethods(), ...){
  wd <- getwd()
  for(dn in names(datasets)){
    ds <- datasets[[dn]]
    if(!is.null(ds$folder)) dn <- ds$folder
    if(dir.exists(dn)){
      print(dn)
      runMethods(ds, dn, methods=methods, ...)
    }else{
      warning("Could not find dataset ",dn)
    }
    setwd(wd)
  }
}

# multi-threaded version of runAll()
runAllMt <- function(datasets=getDatasets(), nthreads=3, ...){
  wd <- getwd()
  library(BiocParallel)
  bplapply(datasets, BPPARAM=MulticoreParam(nthreads), FUN=function(x){
    dn <- head(x$truth,1)
    if(!dir.exists(dn)){
      warning("Could not find dataset ",tf)
      return(0)
    }
    runMethods(x, dn, ...)
    return(1)
  })
}

#' Run all methods on one dataset
#'
#' @param dataset A dataset object (one of the list elements as produced by
#'   `getDatasets()`)
#' @param folder The folder where the dataset data is (and results will be saved).
#'   Defaults to the current folder.
#' @param scriptsFolder The path to the Scripts folder, relative to `folder`
#' @param methods The methods to run (default all). Use `getMethods(onlyTop=TRUE)`
#'   to run only the top methods.
#' @param decoupleR_modes The decoupleR methods (beside consensus) to run.
#' @param rndSeed The random seed
#' @param peakWidth The peak width to enforce (default 300, use 150 for NF peaks)
#' @param forceRerun Whether to regenerate peak counts and motif instances
#' @param outSubfolder Subfolder in which to save the results.
#' @param readlist Vector of paths to aligned files, otherwise will be detected 
#'   from folder.
#' @param DA.norm Normalization method to use for peak differential accessibility
#' 
#' @return A list of dataframes for each method (also saves them to disk)
runMethods <- function(dataset, folder=".", scriptsFolder="../../Scripts",
                       methods=getMethods(), rndSeed=1997, peakWidth=300, 
                       decoupleR_modes=c("mlm", "ulm", "udt", "wsum"),
                       forceRerun=FALSE, outSubfolder="runATAC_results",
                       readlist=NULL, DA.norm="TMM"){
  
  methods <- match.arg(methods, several.ok = TRUE)
  
  suppressPackageStartupMessages({
    library(viper)
    library(SummarizedExperiment)
    library(Matrix)
    library(GenomicRanges)
    library(edgeR)
    library(limma)
    library(monaLisa)
    library(fgsea)
    library(stringr)
    # for variations to the chromVAR normalization, install:
    # devtools::install_github("Jiayi-Wang-Joey/chromVAR")
    library(chromVAR)
    library(motifmatchr)
    library(BiocParallel)
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(BSgenome.Mmusculus.UCSC.mm10)
    library(Rsamtools)
    library(dplyr)
    library(MotifDb)
    library(universalmotif)
    library(memes)
    library(GenomicAlignments)
    library(data.table)
    library(decoupleR)
    library(AUCell)
    library(TFBSTools)
    library(ggplot2)
    library(epiwraps)
    library(ATACseqTFEA)
  })
  register(MulticoreParam(8))
  
  setwd(folder)
  
  mypath <- function(x, folder=".") paste0(folder,"/",x)
  if(is.null(dataset$readType)) dataset$readType <- "bam"
  if(is.null(dataset$paired)) dataset$paired <- TRUE
  if(is.null(dataset$peakFile)) dataset$peakFile <- head(
    list.files(mypath("peaks"), pattern="narrowPeak$|narrowPeak\\.gz$|bed\\.gz$|bed$", full=TRUE),1)
  
  # Obtaining and re-sizing peaks
  peaks <- sort(rtracklayer::import(dataset$peakFile))
  if(!is.na(peakWidth))
    peaks <- resize(peaks, width=peakWidth, fix="center")
  peaks <- keepStandardChromosomes(peaks, pruning.mode="coarse")
  
  # Paths to the method wrappers:
  for(f in c(list.files(scriptsFolder, pattern="run.+R"), "getpmoi.R",
             "differentialAccessibility.R", "QCnorm.R")){
    source(file.path(scriptsFolder, f))
  }
  
  dir.create(mypath(outSubfolder), showWarnings=FALSE)
  for(f in c("others","raw","with_pvalues","scores_only")){
    dir.create(mypath(f,outSubfolder), showWarnings=FALSE)
  }
  genome <- if(dataset$species=="h"){
    BSgenome.Hsapiens.UCSC.hg38
  } else {
    BSgenome.Mmusculus.UCSC.mm10
  }
  spec <- ifelse(dataset$species=="h", "Hsapiens", "Mmusculus")
  
  if(!is.null(dataset$seqStyle)){
    seqStyle <- dataset$seqStyle
  }else{
    if(dataset$readType=="bam"){
      b1 <- head(list.files(pattern="bam$", "seq_files", full=TRUE),1)
      seqStyle <- seqlevelsStyle(Rsamtools::BamFile(b1))
    }else{
      seqStyle <- "UCSC"
    }
    if(length(seqStyle)>1)
      seqStyle <- ifelse(any(grepl("Ensembl|ensembl",seqStyle)), "ensembl", "UCSC")
  }
  if(any(seqStyle=="NCBI")) seqStyle <- "UCSC"
  seqlevelsStyle(peaks) <- seqStyle
  seqlevelsStyle(genome) <- seqStyle  

  # use pmoi if available
  if(file.exists(pmoiPath <- mypath("others/pmoi.rds",outSubfolder)) && 
     !forceRerun){
    pmoi <- readRDS(pmoiPath)
  }else{
    pmoi <- getpmoi(genome=genome,
                    peaks=peaks,
                    spec=spec,
                    seqStyle=seqStyle, srcFolder=scriptsFolder)
    saveRDS(pmoi, pmoiPath)
  }
  if(is.null(readlist))
    readlist <- list.files(mypath("seq_files"), 
                           pattern=paste0(dataset$readType,"$"), full=TRUE)
  if(is.null(dataset$design))
    dataset$design <-
    c(rep(-1, length(readlist)/2), rep(1, length(readlist)/2))
  design <- dataset$design
  
  # Obtaining read counts if not already available
  if(file.exists(cntsPath <- mypath("others/countmatrix.rds",outSubfolder)) &&
     !forceRerun){
    counts <- readRDS(cntsPath)
  }else{
    counts <- chromVAR::getCounts(readlist, peaks,
                                  paired = dataset$paired,
                                  format = dataset$readType)
    counts <- chromVAR::addGCBias(counts, genome=genome)
    row.names(counts) <- as.character(granges(counts))
    saveRDS(counts, cntsPath)
  }
  
  
  #############################################
  ## FROM HERE ON TAKEN FROM runATAC.R (written by Felix Gerbaldo) with small modifs
  #############################################
  
  seqlevelsStyle(pmoi) <- seqStyle
  
  pmoi2 <- pmoi # needed unmodified for BaGFoot
  
  # Compute regulons and GSEA genesets from pmoi object
  
  pmoi$motif_id <- factor(pmoi$motif_id)
  
  # divide motif scores by the max for that motif
  
  msmax <- max(splitAsList(pmoi$score, pmoi$motif_id))
  pmoi$score2 <- pmoi$score/msmax[pmoi$motif_id]
  
  # overlap peaks with motif instances
  
  o <- findOverlaps(peaks, pmoi, ignore.strand=TRUE)
  
  # build a matrix of scores (e.g. for alternative use with chromVAR)
  
  m <- sparseMatrix(from(o), as.integer(pmoi$motif_id[to(o)]), x=pmoi$score2[to(o)])
  mb <- sparseMatrix(from(o), as.integer(pmoi$motif_id[to(o)]))
  row.names(mb) <- row.names(m) <- as.character(granges(peaks))
  colnames(mb) <- colnames(m) <- levels(pmoi$motif_id)
  
  # build a regulon (for use with viper)
  # here we need to name the peaks (I use the coordinates) and make sure the names
  #  corresponding to those of the differential accessibility signature later on
  
  ll <- split( data.frame(peak=as.character(granges(peaks))[from(o)],
                          score=pmoi$score2[to(o)]),
               pmoi$motif_id[to(o)] )
  
  regulons <- lapply(ll, FUN=function(x){
    # count a peak only once if it has two instances of the same motif
    x <- x[order(-x$score),]
    x <- x[!duplicated(x$score),]
    # reformat to regulon
    list( tfmode=setNames(rep(1L,nrow(x)), x$peak),
          likelihood=setNames(x$score, x$peak) )
  })
  
  # Changing the names in the regulons to the conventional ones
  
  if (spec=="Hsapiens") {
    motifnames <- fread(file.path(scriptsFolder, "HOCOMOCOv11_core_annotation_HUMAN_mono.tsv"))
  } else if (spec=="Mmusculus") {
    motifnames <- fread(file.path(scriptsFolder, "HOCOMOCOv11_core_annotation_MOUSE_mono.tsv"))
  }
  
  # These will be turned into the motifs used for chromVAR and monaLisa later but first, they serve as donor of "wrong" TF names
  
  motifs <- getNonRedundantMotifs(format = "PWMatrix", 
                                  species = spec)
  
  # Listing TF names we don't want (wrong)
  
  motiflist <- c()
  
  for (i in seq_along(motifs)){
    for (j in seq_along(motifnames$Model)) {
      if (motifs[[i]]@ID == motifnames$Model[[j]])
        motiflist$wrong <- append(motiflist$wrong,motifs[[i]]@name)
    }
  }
  
  # Adding a list entry with the TF names we want (correct)
  
  for (i in seq_along(motifs)){
    for (j in seq_along(motifnames$Model)) {
      if (motifs[[i]]@ID == motifnames$Model[[j]])
        motiflist$correct <- append(motiflist$correct,motifnames$`Transcription factor`[[j]])
    }
  }
  
  # Replacing the undesired names in the regulons with the desired ones
  
  for (i in seq_along(regulons)){
    for (j in seq_along(motiflist$wrong)){
      if (names(regulons)[[i]] == motiflist$wrong[[j]])
        names(regulons)[[i]] <- motiflist$correct[[j]]
    }
  }
  
  saveRDS(regulons, mypath("others/regulons.rds",outSubfolder))
  
  # get normal lists of peaks per motif (for GSEA)
  
  genesets <- split( as.character(granges(peaks))[from(o)], pmoi$motif_id[to(o)] )
  
  # Replacing the undesired names in the genesets with the desired ones
  
  for (i in seq_along(genesets)){
    for (j in seq_along(motiflist$wrong)){
      if (names(genesets)[[i]] == motiflist$wrong[[j]])
        names(genesets)[[i]] <- motiflist$correct[[j]]
    }
  }
  
  saveRDS(genesets, mypath("others/genesets.rds",outSubfolder))
  
  # Here, the names in the motifs for chromVAR and monaLisa are changed to the correct ones, as well
  
  for (i in seq_along(motifs)){
    for (j in seq_along(motifnames$Model)) {
      if (motifs[[i]]@ID == motifnames$Model[[j]])
        names(motifs)[[i]] <- motifnames$`Transcription factor`[[j]]
    }
  }
  
  BANP_motif <- readRDS(file.path(scriptsFolder, "BANP.PFMatrix.rds"))
  BANP_motifPW <- TFBSTools::toPWM(BANP_motif) 
  motifs$BANP <- BANP_motifPW
  if(spec=="Mmusculus")
    motifs$NR1H3 <- readRDS(file.path(scriptsFolder, "NR1H3.PWMatrix.rds"))
  
  saveRDS(motifs, mypath("others/motifs.rds",outSubfolder))
  
  # Compute differentially accessible regions required to run monaLisa, StabSel, fGSEA, VIPER, and msVIPER 
  
  if (any(grepl("mona|Stab|GSEA|VIPER|decoupleR|ulm|regreg",methods))){
    
    set.seed(rndSeed)
    npos <- sum(dataset$design == 1)
    nneg <- sum(dataset$design == -1)
    counts_control <- counts[, colnames(counts)[1:npos]]
    counts_perturbed <- counts[, colnames(counts)[(npos+1):(npos+nneg)]]
    DAR <- dATestedgeR(counts_control, 
                       counts_perturbed, norm.method=DA.norm)
    
    
    
    # Generate required matrix of logFCs
    
    DARmat <- as.numeric(DAR$logFC)
    names(DARmat) <- rownames(DAR)
  }
  
  # Start running methods
  
  readouts <- list()
  
  # Run ATACseqTFEA
  
  if("ATACseqTFEA" %in% methods){
    set.seed(rndSeed)
    res <- runATACseqTFEA(readlist, pmoi)
    saveRDS(res, mypath("raw/ATACseqTFEA_raw.rds",outSubfolder))
    saveRDS(res$res, mypath("with_pvalues/ATACseqTFEA.rds",outSubfolder))
    readouts$ATACseqTFEA <- res$res
  }
  
  
  # Run chromVAR without normalization
  
  if ("chromVAR" %in% methods){
    set.seed(rndSeed)
    CV <- runchromVAR(counts, 
                           genome, 
                           m, 
                           design)
    
    saveRDS(CV, mypath("raw/CV_raw.rds",outSubfolder))
    
    # Select desired information from chromVAR readout
    
    CVsel <- CV[[1]][order(CV[[1]]$adj.P.Val),]
    CVdf <- data.frame(row.names = rownames(CVsel),
                       logFC = CVsel[, "logFC"], 
                       padj = CVsel[, "adj.P.Val"], 
                       p = CVsel[, "P.Value"])
    CVdf$rank = seq_along(row.names(CVdf))
    
    saveRDS(CVdf, mypath("with_pvalues/CV.rds",outSubfolder))
    readouts$CV <- CVdf
    
    dev <- CV[[3]]
    dev$condition <- factor(design)
    dd <- chromVAR::differentialDeviations(dev, "condition")
    dd <- dd[order(dd$p_value),]
    colnames(dd) <- c("p","padj")
    dd$rank <- seq_len(nrow(dd))
    saveRDS(dd, file.path(outSubfolder, "with_pvalues", "CVoriginal.rds"))
    
  }
  
  # Run chromVAR with normalization
  
  if ("chromVAR" %in% methods){
    # re-use stock chromVAR object:
    dev <- CV[[3]]
    assays(dev)$norm <- scale(assays(dev)$z)
    devMat <- assays(dev)$norm
    fit <- eBayes(lmFit(devMat, design))
    topTFs <- topTable(fit, number = Inf)
    CV <- list(topTFs, CV[[2]], dev)
    saveRDS(CV, mypath("raw/CVnorm_raw.rds",outSubfolder))
    
    # Select desired information from chromVAR readout
    
    CVsel <- CV[[1]][order(CV[[1]]$adj.P.Val),]
    CVdf <- data.frame(row.names = rownames(CVsel),
                       logFC = CVsel[, "logFC"], 
                       padj = CVsel[, "adj.P.Val"], 
                       p = CVsel[, "P.Value"])
    CVdf$rank = seq_along(row.names(CVdf))
    
    saveRDS(CVdf, mypath("with_pvalues/CVnorm.rds",outSubfolder))
    readouts$CVnorm <- CVdf
  }
  
  if ("decoupleR" %in% methods){
    
    set.seed(rndSeed)
    res <- rundecoupleR( counts,
                         genome,
                         motifs=pmoi,#motifs,
                         DAR,
                         decoupleR_modes = decoupleR_modes)
    saveRDS(res$raw, "./runATAC_results/raw/decoupleR_raw.rds")
    for(f in names(res$res))
      saveRDS(res$res[[f]],
              paste0("./runATAC_results/with_pvalues/decoupleR",f,".rds"))
  }
  
  # Run monaLisa 
  
  if ("monaLisa" %in% methods || "minaLisa.zero" %in% methods){
    
    ML <- runmonaLisa(DAR, 
                      motifs, 
                      peaks, 
                      genome, background="zeroBin")
    saveRDS(ML, "./runATAC_results/raw/MLzero_raw.rds")
    saveRDS(ML$df, "./runATAC_results/with_pvalues/MLzero.rds")
    readouts$MLzero <- ML$df
    
  }
  if ("monaLisa" %in% methods || "minaLisa.others" %in% methods){

    set.seed(rndSeed)
    ML <- runmonaLisa(DAR, 
                      motifs, 
                      peaks, 
                      genome)
    
    saveRDS(ML, "./runATAC_results/raw/ML_raw.rds")
    saveRDS(ML$df, "./runATAC_results/with_pvalues/ML.rds")
    readouts$ML <- ML$df
  }
  
  if ("monaLisa" %in% methods){
    set.seed(rndSeed)
    ML <- runmonaLisa(DAR, 
                      motifs, 
                      peaks, 
                      genome, minAbsLfc=0.15)
    saveRDS(ML, "./runATAC_results/raw/MLlower_raw.rds")
    saveRDS(ML$df, "./runATAC_results/with_pvalues/MLlower.rds")

    set.seed(rndSeed)
    ML <- runmonaLisa(DAR, 
                      motifs, 
                      peaks, 
                      genome, nBins=7)
    saveRDS(ML, "./runATAC_results/raw/MLfewerBins_raw.rds")
    saveRDS(ML$df, "./runATAC_results/with_pvalues/MLfewerBins.rds")

    # use correlation across bins
    MLdf <- readouts$ML
    MLdf <- MLdf[order(abs(MLdf$binSpearman)*-log10(MLdf$p), decreasing=TRUE),]
    MLdf$rank <- seq_along(row.names(MLdf))
    saveRDS(MLdf, "./runATAC_results/with_pvalues/MLsp.rds")
    readouts$MLsp <- MLdf
    
  }
  
  # Run StabSel
  
  if ("StabSel" %in% methods){
    set.seed(rndSeed)
    MLStabSel <- runStabSel(DAR, 
                            motifs, 
                            peaks, 
                            genome)
    
    saveRDS(MLStabSel, "./runATAC_results/raw/MLStabSel_raw.rds")
    
    MLStabSeldf <- data.frame(TF = colnames(MLStabSel[[1]]),
                              sel_Prob = MLStabSel[[1]]$selProb)
    MLStabSeldf <- MLStabSeldf[order(MLStabSeldf$sel_Prob, decreasing = TRUE),]         
    MLStabSeldf <- data.frame(row.names = MLStabSeldf$TF,
                              sel_Prob = MLStabSeldf$sel_Prob)
    MLStabSeldf <- subset(MLStabSeldf, !(row.names(MLStabSeldf) %in% c("fracGC", "oeCpG")))
    MLStabSeldf$rank <- seq_along(rownames(MLStabSeldf))
    p <- (5-MLStabSeldf$sel_Prob)/100
    w <- which(MLStabSeldf$sel_Prob<0.9)
    p[w] <- (10-MLStabSeldf$sel_Prob[w])/100
    MLStabSeldf$p <- MLStabSeldf$padj <- p
    saveRDS(MLStabSeldf, "./runATAC_results/with_pvalues/MLStabSel.rds")
    readouts$MLStabSel <- MLStabSeldf
  }
  
  # Run fGSEA
  
  if ("GSEA" %in% methods){
    
    set.seed(rndSeed)
    GSEA <- rungsea(DAR, 
                    genesets, 
                    peaks)
    
    saveRDS(GSEA, "./runATAC_results/raw/GSEA_raw.rds")
    
    # Select desired information from fgsea readout
    
    GSEA <- GSEA[[1]][order(GSEA[[1]]$pval, -abs(GSEA[[1]]$NES)),]
    GSEAdf <- data.frame(row.names = GSEA$pathway,
                         NES = GSEA$NES,
                         p = GSEA$pval,
                         padj = GSEA$padj)
    
    GSEAdf$rank <- seq_along(row.names(GSEAdf))
    
    saveRDS(GSEAdf, "./runATAC_results/with_pvalues/GSEA.rds")
    readouts$GSEA <- GSEAdf
  }
  
  # Run VIPER
  if ("VIPER" %in% methods){
    
    set.seed(rndSeed)
    VIPER <- runviper(counts_control, 
                      counts_perturbed, 
                      DAR, 
                      peaks, 
                      regulons, 
                      mode = "viper", 
                      method = "ttest")
    
    saveRDS(VIPER, "./runATAC_results/raw/VIPER_raw.rds")
    
    VIPERtTest <- rowTtest(VIPER[[1]][,1:npos], VIPER[[1]][,(npos+1):(npos+nneg)])
    
    # Select desired information from VIPER readout
    
    VIPERdf <- data.frame(VIPERtTest)
    VIPERdf <- VIPERdf[order(VIPERdf$p.value),]
    VIPERdf$padj <- p.adjust(VIPERdf$p.value, method="fdr")
    VIPERdf <- data.frame(row.names = rownames(VIPERdf),
                          statistic = VIPERdf$statistic,
                          padj = VIPERdf$padj,
                          p = VIPERdf$p.value)
    VIPERdf$rank <- seq_along(row.names(VIPERdf))
    
    saveRDS(VIPERdf, "./runATAC_results/with_pvalues/VIPER.rds")
    readouts$VIPER <- VIPERdf
  }
  
  # Run VIPER
  if ("VIPERb" %in% methods){
    
    set.seed(rndSeed)
    VIPERb <- runviper(counts_control, 
                      counts_perturbed, 
                      DAR, 
                      peaks, 
                      regulons, 
                      mode = "viper", 
                      method = "ttest", binarizeRegulon=TRUE)
    
    saveRDS(VIPERb, "./runATAC_results/raw/VIPERb_raw.rds")
    
    VIPERtTest <- rowTtest(VIPERb[[1]][,1:npos], VIPERb[[1]][,(npos+1):(npos+nneg)])
    
    # Select desired information from VIPER readout
    
    VIPERbdf <- data.frame(VIPERtTest)
    VIPERbdf <- VIPERbdf[order(VIPERbdf$p.value),]
    VIPERbdf$padj <- p.adjust(VIPERbdf$p.value, method="fdr")
    VIPERbdf <- data.frame(row.names = rownames(VIPERbdf),
                          statistic = VIPERbdf$statistic,
                          padj = VIPERbdf$padj,
                          p = VIPERbdf$p.value)
    VIPERbdf$rank <- seq_along(row.names(VIPERbdf))
    
    saveRDS(VIPERbdf, "./runATAC_results/with_pvalues/VIPERb.rds")
    readouts$VIPERb <- VIPERbdf
  }
  
  # Run msVIPER
  
  if ("msVIPER" %in% methods){
    set.seed(rndSeed)
    msVIPER <- runviper(counts_control, 
                        counts_perturbed,
                        DAR,
                        peaks,
                        regulons,
                        mode = "msviper",
                        method = "ttest")
    
    saveRDS(msVIPER, "./runATAC_results/raw/msVIPER_raw.rds")
    
    # Select desired information from msVIPER readout
    
    msVIPERdf <- data.frame(NES = msVIPER[[1]]$es$nes,
                            padj = p.adjust(msVIPER[[1]]$es$p.value, method="fdr"),
                            p = msVIPER[[1]]$es$p.value)
    
    msVIPERdf <- msVIPERdf[order(msVIPERdf$p),]
    msVIPERdf$rank <- seq_along(row.names(msVIPERdf))
    
    saveRDS(msVIPERdf, "./runATAC_results/with_pvalues/msVIPER.rds")
    readouts$msVIPER <- msVIPERdf
  }
  
  
  if ("msVIPERb" %in% methods){
    set.seed(rndSeed)
    msVIPERb <- runviper(counts_control, 
                        counts_perturbed,
                        DAR,
                        peaks,
                        regulons,
                        mode = "msviper",
                        method = "ttest", binarizeRegulon=TRUE)
    
    saveRDS(msVIPERb, "./runATAC_results/raw/msVIPERb_raw.rds")
    
    # Select desired information from msVIPERb readout
    
    msVIPERbdf <- data.frame(NES = msVIPERb[[1]]$es$nes,
                            padj = p.adjust(msVIPERb[[1]]$es$p.value, method="fdr"),
                            p = msVIPERb[[1]]$es$p.value)
    
    msVIPERbdf <- msVIPERbdf[order(msVIPERbdf$p),]
    msVIPERbdf$rank <- seq_along(row.names(msVIPERbdf))
    
    saveRDS(msVIPERbdf, "./runATAC_results/with_pvalues/msVIPERb.rds")
    readouts$msVIPERb <- msVIPERbdf
  }
  
  # Run ulm
  
  if ("ulm" %in% methods){
    set.seed(rndSeed)
    ulm <- runulm(DARmat = DARmat,
                  matchMtx = m)
    
    saveRDS(ulm, "./runATAC_results/raw/ulm_raw.rds")
    
    ulmdf <- data.frame(row.names = rownames(ulm[[1]]),
                        score = ulm[[1]]$score,
                        padj = ulm[[1]]$padj,
                        p = ulm[[1]]$p)
    ulmdf <- ulmdf[order(ulmdf$p, -abs(ulmdf$score)),]
    ulmdf$rank <- seq_along(row.names(ulmdf))
    
    
    saveRDS(ulmdf, "./runATAC_results/with_pvalues/ulm.rds")
    readouts$ulm <- ulmdf
  }

  if ("ulmB" %in% methods){
    set.seed(rndSeed)
    ulm <- runulm(DARmat = DARmat,
                  matchMtx = mb)

    saveRDS(ulm, "./runATAC_results/raw/ulmB_raw.rds")

    ulmdf <- data.frame(row.names = rownames(ulm[[1]]),
                        score = ulm[[1]]$score,
                        padj = ulm[[1]]$padj,
                        p = ulm[[1]]$p)
    ulmdf <- ulmdf[order(ulmdf$p, -abs(ulmdf$score)),]
    ulmdf$rank <- seq_along(row.names(ulmdf))


    saveRDS(ulmdf, "./runATAC_results/with_pvalues/ulmB.rds")
    readouts$ulmB <- ulmdf
  }
  
  if ("ulmGC" %in% methods){
    set.seed(rndSeed)
    ulmGC <- runulm(DAR = DAR,
                  matchMtx = m)
    
    saveRDS(ulmGC, "./runATAC_results/raw/ulmGC_raw.rds")
    
    ulmGCdf <- data.frame(row.names = rownames(ulmGC[[1]]),
                        score = ulmGC[[1]]$score,
                        padj = ulmGC[[1]]$padj,
                        p = ulmGC[[1]]$p)
    ulmGCdf <- ulmGCdf[order(ulmGCdf$p, -abs(ulmGCdf$score)),]
    ulmGCdf$rank <- seq_along(row.names(ulmGCdf))
    
    
    saveRDS(ulmGCdf, "./runATAC_results/with_pvalues/ulmGC.rds")
    readouts$ulmGC <- ulmGCdf
  }
  
  # Run regreg
  
  if ("regreg" %in% methods){
    set.seed(rndSeed)
    regreg <- runregreg(DARmat = DARmat,
                        matchMtx = mb, alpha=1)
    
    saveRDS(regreg, "./runATAC_results/raw/regreg_raw.rds")
    
    regregdf <- data.frame(row.names = rownames(regreg[[1]]),
                           score = regreg[[1]]$score,
                           padj = regreg[[1]]$FDR,
                           p = regreg[[1]]$p_value)
    regregdf$rank <- seq_along(row.names(regregdf))
    
    saveRDS(regregdf, "./runATAC_results/with_pvalues/regreg.rds")
    readouts$regreg <- regregdf
  }

  if ("regregR" %in% methods){
    set.seed(rndSeed)
    regregR <- runregreg(DARmat = DARmat,
                        matchMtx = mb, alpha=0)
    
    saveRDS(regregR, "./runATAC_results/raw/regregR_raw.rds")
    
    regregRdf <- data.frame(row.names = rownames(regregR[[1]]),
                           score = regregR[[1]]$score,
                           padj = regregR[[1]]$FDR,
                           p = regregR[[1]]$p_value)
    regregRdf$rank <- seq_along(row.names(regregRdf))
    
    saveRDS(regregRdf, "./runATAC_results/with_pvalues/regregR.rds")
    readouts$regregR <- regregRdf
  }

  # run Model-based analysis

  if ("MBA" %in% methods){
    #try({
    set.seed(rndSeed)
    MBA <- runMBA(readlist, pmoi2, dataset$paired)
    saveRDS(MBA, "./runATAC_results/raw/MBA_raw.rds")
    MBA <- MBA[[1]]
    MBAdf <- data.frame(row.names=rownames(MBA),
                        padj=MBA$adj.P.Val,
                        p=MBA$P.Value)
    MBAdf <- MBAdf[order(MBAdf$p),]
    MBAdf$rank =seq_along(row.names(MBAdf))
    saveRDS(MBAdf, "runATAC_results/with_pvalues/MBA.rds")
    readouts$MBA <- MBAdf
    #}, silent=TRUE)
  }


  
  # run BaGFootLike method
  
  if ("BaGFoot" %in% methods){
  #try({
    set.seed(rndSeed)
    BaGFoot <- runBaGFoot(readlist,
                          pmoi2,
                          dataset$paired)
    saveRDS(BaGFoot, "./runATAC_results/raw/BaGFootLike_raw.rds")
    
    BaGFoot <- BaGFoot[[1]]
    BFdf <- data.frame(row.names = rownames(BaGFoot),
                       padj=BaGFoot$FDR,
                       p=BaGFoot$combined.P)
    BFdf <- BFdf[order(BFdf$p),]
    BFdf$rank = seq_along(row.names(BFdf))
    saveRDS(BFdf,"./runATAC_results/with_pvalues/BaGFootLike.rds")
    readouts$BaGFoot <- BFdf
  #}, silent=TRUE)
  }
  
  saveRDS(readouts, "./runATAC_results/others/readouts.rds")
  return(readouts=readouts)
}


runCVariants <- function(datasets){
  library(chromVAR)
  library(limma)
  for(dn in names(datasets)){
    ds <- datasets[[dn]]
    if(!is.null(ds$folder)) dn <- ds$folder
    if(dir.exists(dn)){
      print(dn)
      dev <- readRDS(file.path(dn, "runATAC_results", "raw", "CV_raw.rds"))[[3]]
      assays(dev)$norm <- scale(assays(dev)$z)
      assays(dev)$centered <- scale(assays(dev)$z, scale=FALSE)
      qt <- preprocessCore::normalize.quantiles(assays(dev)$z)
      dimnames(qt) <- dimnames(assays(dev)$z)
      assays(dev)$qt <- qt
      assays(dev)$devnorm <- scale(assays(dev)$deviations)
      assays(dev)$devcentered <- scale(assays(dev)$deviations, scale=FALSE)
      qt <- preprocessCore::normalize.quantiles(assays(dev)$deviations)
      dimnames(qt) <- dimnames(assays(dev)$deviations)
      assays(dev)$devqt <- qt
      
      ass <- c("CV"="z", "CVcentered"="centered", "CVnorm"="norm", "CVqt"="qt",
               "CVdev"="deviations", "CVdevNorm"="devnorm", "CVdevCentered"="devcentered", "CVdevqt"="devqt" )
      design <- c(rep(-1, ncol(dev)/2), rep(1, ncol(dev)/2))
      for(f in names(ass)){
        devMat <- assays(dev)[[ass[[f]]]]
        fit <- eBayes(lmFit(devMat, design))
        CVsel <- topTable(fit, number=Inf)
        CVdf <- data.frame(row.names = rownames(CVsel),
                           logFC = CVsel[, "logFC"], 
                           padj = CVsel[, "adj.P.Val"], 
                           p = CVsel[, "P.Value"],
                           t = CVsel[, "t"])
        CVdf$rank = seq_along(row.names(CVdf))
        saveRDS(CVdf, file.path(dn, "runATAC_results", "with_pvalues", paste0(f, ".rds")))
      }
      dev$condition <- factor(design)
      dd <- chromVAR::differentialDeviations(dev, "condition")
      dd <- dd[order(dd$p_value),]
      colnames(dd) <- c("p","padj")
      dd$rank <- seq_len(nrow(dd))
      saveRDS(dd, file.path(dn, "runATAC_results", "with_pvalues", "CVoriginal.rds"))
      
    }else{
      warning("Could not find dataset ",tf)
    }
  } 
}
