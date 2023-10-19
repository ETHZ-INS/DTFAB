rankHeatmap <- function(res, rankBreaks=c(1,6,20,50,100,200,400), squeeze=FALSE){
  p <- ggplot(res, aes(reorder(dataset,-sqrt(rank)), reorder(method,-sqrt(rank)),
                  fill=sqrt(rank), label=rank)) + geom_tile() +  geom_text() +
    scale_fill_viridis_c(direction=-1, breaks=sqrt(rankBreaks), 
                         trans="sqrt", labels=c("top",rankBreaks[-1])) +
    labs(fill="Rank of\ntrue TF", y="Methods", x="Datasets") + theme_minimal()
  if(squeeze) p <- p + 
      theme(legend.position="bottom", legend.key.width = unit(1,"cm"),
            axis.text.x=element_text(angle=45, hjust=1, vjust=1))
  p
}

rankHeatmap2 <- function(res){
  rankm <- reshape2::dcast(res, method~dataset, value.var="rank")
  row.names(rankm) <- rankm[,1]
  rankm <- sqrt(as.matrix(rankm[,-1]))
  rankm <- rankm[order(rowMeans(rankm)), order(-colMeans(rankm))]
  type <- sapply(getDatasets(), FUN=function(x) x$type)
  draw(ComplexHeatmap::pheatmap(rankm, cluster_rows=FALSE, cluster_cols=FALSE,
     display_numbers=rankm^2, border=NA, number_color="black",
     color=viridisLite::viridis(100, direction=-1), name="Rank of\ntrue TF",
     column_title="Datasets", column_title_side="bottom",
     annotation_col=as.data.frame(type)[colnames(rankm),,drop=FALSE],
     annotation_colors=list(type=c(deletion="red", dTag="black", ligand="darkblue"))),
     merge=TRUE)
}



sensFDRplot <- function(res){
  res$Sensitivity <- res$trueQ<0.05
  res$FDR[which(is.na(res$FDR))] <- 0
  res2 <- aggregate(res[,c("Sensitivity", "FDR")], by=res[,"method",drop=FALSE], mean)
  ggplot(res2, aes(FDR, Sensitivity, label=method)) + geom_point() + 
    theme_bw() + ggrepel::geom_text_repel(min.segment.length=0)
}

relAUCplot <- function(res, squeeze=FALSE){
  res$auc2 <- gsub("^0\\.",".",round(res$relAUC,2))
  p <- ggplot(res, aes(reorder(dataset,relAUC), reorder(method,relAUC), 
                  fill=relAUC, label=auc2)) + 
    geom_tile() + theme_minimal() + geom_text() +
    scale_fill_viridis_c(option="A") +
    labs(fill="relative\nAUC", y="Methods", x="Datasets") + theme_minimal()
  if(squeeze) p <- p + 
    theme(legend.position="bottom", legend.key.width = unit(1,"cm"),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1))
  p
}

runtimePlot <- function(res){
  res <- res[!is.na(res$elapsed),]
  res$runtime <- pmax(res$elapsed, res$cpu)
  ggplot(res, aes(runtime/60, reorder(method,log10(runtime)))) + 
    geom_boxplot(fill="lightgrey", outlier.size=1) + 
    scale_x_log10() + theme_bw() + labs(x="Runtime (min)",y="")
}

plotAllMetrics <- function(res){
  cowplot::plot_grid(
    rankHeatmap(res, squeeze = TRUE),
    relAUCplot(res, squeeze = TRUE),
    sensFDRplot(res), runtimePlot(res),
    nrow=2, labels="AUTO", scale=0.95, rel_heights = c(5,4)
  )
}

renameMethods <- function(x,renaming=NULL){
  if(is.data.frame(x)){
    stopifnot("method" %in% colnames(x))
    x$method <- renameMethods(x$method, renaming=renaming)
    return(x)
  }
  if(is.null(renaming)) renaming <- c(
    "CVoriginal"="chromVAR", "CV"="chromVAR>limma",
    "CVqt"="chromVAR>Qt>limma", "CVcentered"="chromVAR>center>limma",
    "CVnorm"="chromVAR>scale>limma", ML="MonaLisa.vsOthers",
    MLzero="MonaLisa.vsZero", MLsp="MonaLisa.vsOthers+spearman",
    MLfewerBins="MonaLisa.vsOthers(fewer bins)", MLStabSel="MonaLisa.StabSel",
    MLlower="MonaLisa.vsOthers(smaller zeroBin)", "regreg"="Lasso-lm", "regregR"="Ridge-lm", VIPER="viper(scores)>limma", VIPERb="viper(binary)>limma",
    msVIPER="msViper(scores)", msVIPERb="msViper(binary)", ulmGC="ulm+GC", ulmB="ulm(binary)"
  )
  for(i in names(renaming)) x <- replace(x, x==i, renaming[[i]])
  x <- gsub("decoupleR","decoupleR:",x)
  x <- gsub("::",":",x,fixed=TRUE)
  x
}

filterMainMethods <- function(x, rename=TRUE){
  w <- which(x$method=="CV" | x$method=="chromVAR>scale>limma")
  CVcpu <- setNames(x$cpu[w], x$dataset[w])
  CVelapsed <- setNames(x$elapsed[w], x$dataset[w])
  w <- which(grepl("^CV|^chromVAR",x$method) & is.na(x$cpu))
  x$cpu[w] <- CVcpu[x$dataset[w]]
  x$elapsed[w] <- CVelapsed[x$dataset[w]]
  x$method <- as.character(x$method)
  renames <- c(
    chromVAR = "chromVAR", `chromVAR>Qt>limma` = "chromVAR-adjusted",
    MonaLisa.vsOthers = "MonaLisa", MonaLisa.StabSel = "StabSel", 
    `msViper(scores)` = "msViper", `viper(binary)>limma` = "viper", 
    `decoupleR:consensus` = "decoupleR", GSEA = "GSEA")
  x <- renameMethods(x)
  x <- x[which(x$method %in% names(renames)),]
  renameMethods(x, renames)
}
