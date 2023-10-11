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
  rankm <- rankm[order(-rowMeans(rankm)), order(-colMeans(rankm))]
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
  ggplot(res, aes(runtime/60, method)) + 
    geom_boxplot(fill="lightgrey", outlier.size=1) + 
    scale_x_log10() + theme_bw() + labs(x="Runtime (min)")
}

plotAllMetrics <- function(res){
  cowplot::plot_grid(
    rankHeatmap(res, squeeze = TRUE),
    relAUCplot(res, squeeze = TRUE),
    sensFDRplot(res), runtimePlot(res),
    nrow=2, labels="AUTO", scale=0.95, rel_heights = c(5,4)
  )
}