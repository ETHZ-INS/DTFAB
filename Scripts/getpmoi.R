suppressPackageStartupMessages({
  library(MotifDb)
  library(motifmatchr)
  library(memes)
  library(Matrix)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(chromVAR)
  library(data.table)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(universalmotif)
})

getNonRedundantMotifs <- function(format=c("PFMatrix","universal","PWMatrix"),
                                  species=c("Hsapiens","Mmusculus")){
  species <- match.arg(species)
  motifs <- MotifDb::query(MotifDb::MotifDb, c(species,"HOCOMOCO"))
  pat <- paste0("^",species,"-HOCOMOCOv1[0-1]-|_HUMAN.+|_MOUSE.+|core-[A-D]-|secondary-[A-D]-")
  modf <- data.frame(row.names=names(motifs),
                     TF=gsub(pat,"",names(motifs)),
                     grade=gsub(".+\\.","",names(motifs)))
  modf <- modf[order(modf$TF,-as.numeric(grepl("HOCOMOCOv11",row.names(modf))),modf$grade),]
  modf <- modf[!duplicated(modf$TF),]
  motifs <- motifs[row.names(modf)]
  switch(match.arg(format),
         universal=setNames(universalmotif::convert_motifs(motifs), modf$TF),
         PFMatrix=do.call(TFBSTools::PFMatrixList, setNames(
           universalmotif::convert_motifs(motifs, class="TFBSTools-PFMatrix"),
           modf$TF)),
         PWMatrix=do.call(TFBSTools::PWMatrixList, 
                          setNames(universalmotif::convert_motifs(motifs, 
                                                                  class="TFBSTools-PWMatrix"), modf$TF))
  )
}

#' Title
#'
#' @param genome A BSgenome respective to the genome used for the alignment.
#' @param peakpath The path to the merged ATAC-seq peaks from control and treatment conditions.
#' @param spec Species. Either "Hsapiens" or "Mmusculus".
#' @param seqStyle Either "ensembl" or "UCSC" depending on the format of the peak file.
#' @param srcFolder Folder where to find scripts and important objects
#'
#' @return
#' @export
#'
#' @examples
getpmoi <- function(genome,
                    peaks,
                    spec=c("Hsapiens","Mmusculus"),
                    seqStyle=c("ensembl","UCSC"),
                    srcFolder){
  seqStyle <- match.arg(seqStyle)
  # choose the file that contains the correct names from HOCOMOCO v11
  if (spec=="Hsapiens") {
    motifnames <- fread(file.path(srcFolder, "HOCOMOCOv11_core_annotation_HUMAN_mono.tsv"))
  } else if (spec=="Mmusculus") {
    motifnames <- fread(file.path(srcFolder, "HOCOMOCOv11_core_annotation_MOUSE_mono.tsv"))
  }
  
  if(is.character(peaks)){
    peals <- sort(rtracklayer::import(peaks))
  }
  peaks <- keepStandardChromosomes(peaks, 
                                   pruning.mode = "coarse")
  
  seqlevelsStyle(genome) <- seqStyle
  peak_seqs <- get_sequence(peaks, genome)
  
  # Get the motifs in universal format required by memes
  motifs <- getNonRedundantMotifs("universal", species = spec)
  
  BANP_motif <- readRDS(file.path(srcFolder, "BANP.PFMatrix.rds"))
  BANP_universalmotif <- convert_motifs(BANP_motif,
                                        class="universalmotif-universalmotif")
  motifs$BANP <- BANP_universalmotif
  if(!("NR1H3" %in% names(motifs)))
    motifs$NR1H3 <- convert_motifs(
      readRDS(file.path(srcFolder, "NR1H3.PWMatrix.rds")),
      class = "universalmotif-universalmotif")
  
  # Correct inconvential names for motifs
  for (i in seq_along(motifs)){ 
    for (j in seq_along(motifnames$Model)) {
      if (motifs[[i]]@altname == motifnames$Model[[j]])
        names(motifs)[[i]] <- motifnames$`Transcription factor`[[j]]
    }
  }
  
  # Obtain the positions of motif instances which are later required as input for runATAC
  pmoi <- runFimo(peak_seqs, 
                  motifs, 
                  meme_path="/common/meme/bin/", 
                  skip_matched_sequence=TRUE)
  try(saveRDS(pmoi, "./runATAC_results/others/pmoi.rds"), silent=TRUE)
  return(pmoi)
  }