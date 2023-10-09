importFrags <- function(bamPath, which=NULL, annotationStyle="NCBI")
{
  if(is.null(which))
  {
    which <- GRanges(Rle(1:2), 
                     IRanges(start=rep(1,2), width=536870912))
  }
  
  param <- ScanBamParam(what=c('pos', 'qwidth', 'isize'),
                        which=which, 
                        flag=scanBamFlag(isUnmappedQuery=FALSE))
  
  readPairs <- readGAlignmentPairs(bamPath, param=param)
  
  # get fragment coordinates from read pairs
  frags <- GRanges(seqnames(GenomicAlignments::first(readPairs)), 
                   IRanges(start=pmin(start(GenomicAlignments::first(readPairs)), 
                                      start(GenomicAlignments::second(readPairs))), 
                           end=pmax(end(GenomicAlignments::first(readPairs)), 
                                    GenomicAlignments::end(second(readPairs)))))
  
  frags <- granges(frags, use.mcols=TRUE)
  
  # change annotation style
  seqlevelsStyle(frags) <- annotationStyle
  
  return(list(fragments=frags,
              readPairs=readPairs))
}

Probabilities <- function(fragments)
{
  fragDf <- as.data.frame(fragments$fragments)
  nf <- fragDf$width[(fragDf$width<120)]
  mono <- fragDf$width[(fragDf$width>=120 & fragDf$width<300)]
  di <- fragDf$width[(fragDf$width>=300 & fragDf$width<500)]
  tri <- fragDf$width[(fragDf$width>=500)]
  
  prob_nf <- length(nf)/length(fragDf$width)
  prob_mono <- length(mono)/length(fragDf$width)
  prob_di <- length(di)/length(fragDf$width)
  prob_tri <- length(tri)/length(fragDf$width)
  
  probs <- setNames(c(prob_nf, 
                      prob_mono, 
                      prob_di, 
                      prob_tri), 
                    c("nf", 
                      "mono", 
                      "di", 
                      "tri"))
  
  return(probs)
}


# fragments <- importFrags("/mnt/plger/fgerbaldo/BenchmarkTFactivity/Simulation_data/filtered_bam/sorted_bam/SRR16946077.bam")
# 
# fragments
# 
# fragDf <- as.data.frame(fragments$fragments)
# fragDf
# 
# fragDf
# nf <- fragDf$width[(fragDf$width<150)]
# mono <- fragDf$width[(fragDf$width>=150 & fragDf$width<300)]
# di <- fragDf$width[(fragDf$width>=300 & fragDf$width<500)]
# tri <- fragDf$width[(fragDf$width>=500)]
# length(fragDf$width)
# length(nf)+length(mono)+length(di)+length(tri)
# 
# prob_nf <- length(nf)/length(fragDf$width)
# 
# prob_mono <- length(mono)/length(fragDf$width)
# 
# prob_di <- length(di)/length(fragDf$width)
# 
# prob_tri <- length(tri)/length(fragDf$width)
# 
# prob_nf
# prob_mono
# prob_di
# prob_tri
