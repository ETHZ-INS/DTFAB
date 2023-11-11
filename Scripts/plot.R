# rankHeatmap <- function(res, rankBreaks=c(1,6,20,50,100,200,400), squeeze=FALSE){
#   p <- ggplot(res, aes(reorder(dataset,-sqrt(rank)), reorder(method,-sqrt(rank)),
#                   fill=sqrt(rank), label=rank)) + geom_tile() +  geom_text() +
#     scale_fill_viridis_c(direction=-1, breaks=sqrt(rankBreaks), 
#                          trans="sqrt", labels=c("top",rankBreaks[-1])) +
#     labs(fill="Rank of\ntrue TF", y="Methods", x="Datasets") + theme_minimal()
#   if(squeeze) p <- p + 
#       theme(legend.position="bottom", legend.key.width = unit(1,"cm"),
#             axis.text.x=element_text(angle=45, hjust=1, vjust=1))
#   p
# }

rankHeatmap2 <- function(res, rankBreaks=c(1,5,25,75,150,300,600), ..., 
                         column_title="Datasets", doDraw=TRUE,
                         datasetInfo=data.frame(
                          type=sapply(getDatasets(), FUN=function(x) x$type))){
  library(ComplexHeatmap)
  rankm <- reshape2::dcast(res, method~dataset, value.var="rank")
  row.names(rankm) <- rankm[,1]
  rankm <- sqrt(as.matrix(rankm[,-1]))
  ro <- row.names(rankm)[order(rowMeans(rankm, na.rm=TRUE))]
  column_order <- names(sort(-colMeans(rankm, na.rm=TRUE)))
  LFCbased <- grepl("GSEA|ulm|msViper|decoupleR|MonaLisa|diffTF|-lm",
                    row.names(rankm),ignore.case=TRUE)
  ancols <- list(type=c(deletion="darkorange3", dTag="black", ligand="darkslateblue"),
                 LFCbased=c("FALSE"="white", "TRUE"="brown4"))
  colan <- HeatmapAnnotation(df = datasetInfo[colnames(rankm),,drop=FALSE],
                             col=ancols, show_annotation_name=FALSE)
  h <- ComplexHeatmap::Heatmap(rankm, cluster_rows=FALSE, cluster_columns=FALSE,
          cell_fun=function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%1.0f", rankm[i, j]^2), x, y, gp = gpar(fontsize=8))
                      },
     col=viridisLite::viridis(100, direction=-1), name="Rank of\ntrue TF",
     column_title=column_title, row_order=ro, column_order=column_order, ..., 
     left_annotation=rowAnnotation(df=as.data.frame(LFCbased), col=ancols, show_legend=FALSE),
     row_names_side="left", top_annotation=colan,
     heatmap_legend_param=list(at=sqrt(rankBreaks), labels=c("top",rankBreaks[-1])))
  if(!doDraw) return(h)
  draw(h, merge=TRUE)
}


sensFDRplot <- function(res){
  res$Sensitivity <- res$trueQ<0.05
  res$FDR[which(is.na(res$FDR))] <- 0
  if(is.null(res$type)){
    res$type <- ifelse(grepl("GSEA|ulm|msViper|decoupleR|MonaLisa|diffTF|-lm|Lasso",
                             res$method, ignore.case=TRUE), "LFC-based", "sample-wise")
  }
  cols <- setNames(c("#CC6677", "#4477AA"), unique(res$type))
  res2 <- aggregate(res[,c("Sensitivity", "FDR")], by=res[,c("method","type")],
                    na.rm=TRUE, FUN=mean)
  ggplot(res2, aes(FDR, Sensitivity, label=method, colour=type)) + geom_point(show.legend=FALSE) + 
    theme_bw() + ggrepel::geom_text_repel(min.segment.length=0, show.legend=FALSE) +
    scale_color_manual(values=cols) +
    theme(legend.position = "bottom")
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

relAUCplot2 <- function(res, doDraw=TRUE, sqrt=FALSE, column_title="Datasets", 
                        column_order=NULL, row_order=NULL, ..., 
                        datasetInfo=data.frame(
                           type=sapply(getDatasets(), FUN=function(x) x$type))){
  m <- reshape2::dcast(res, method~dataset, value.var="relAUC")
  row.names(m) <- m[,1]
  m <- as.matrix(m[,-1])
  if(is.null(column_order)) column_order <- names(sort(colMeans(m, na.rm=TRUE)))
  if(is.null(row_order)) row_order <- row.names(m)[order(-rowMeans(m, na.rm=TRUE))]

  ancols <- list(type=c(deletion="darkorange3", dTag="black", ligand="darkslateblue"),
                 LFCbased=c("FALSE"="white", "TRUE"="brown4"))
  colan <- HeatmapAnnotation(df = datasetInfo[colnames(m),,drop=FALSE],
                             col=ancols)
  h <- ComplexHeatmap::Heatmap(
      m, cluster_rows=FALSE, cluster_columns=FALSE, ..., column_order=column_order,
      cell_fun=function(j, i, x, y, width, height, fill) {
        v <- gsub("^0","",round(m[i,j],2))
        if(is.na(v)) v <- "."
        if(v=="") v <- "0"
        grid.text(v, x, y, gp=gpar(fontsize=8,
                                   col=ifelse(m[i,j]<0.5,"lightgrey","black")))
      },
      col=viridisLite::magma(100), name="Relative\nnetwork AUC",
      column_title=column_title, row_names_side="left", top_annotation=colan)
  if(!doDraw) return(h)
  draw(h, merge=TRUE)
}

runtimePlot <- function(res, removeNAs=FALSE){
  if(removeNAs) res <- res[!is.na(res$elapsed) & !is.na(res$cpu),]
  res$runtime <- pmax(res$elapsed, res$cpu, na.rm=TRUE)
  ggplot(res, aes(runtime/60, reorder(method,runtime))) + stat_summary() +
    #geom_boxplot(fill="lightgrey", outlier.size=1) + 
    theme_bw() + labs(x="CPU time (min)",y="")
}

plotAllMetrics <- function(res){
  cowplot::plot_grid(
    rankHeatmap(res, squeeze = TRUE),
    relAUCplot(res, squeeze = TRUE),
    sensFDRplot(res), runtimePlot(res),
    nrow=2, labels="AUTO", scale=0.95, rel_heights = c(5,4)
  )
}

renameMethods <- function(x, renaming=NULL){
  if(is.data.frame(x)){
    stopifnot("method" %in% colnames(x))
    x$method <- renameMethods(x$method, renaming=renaming)
    return(x)
  }
  if(is.null(renaming)) renaming <- c(
    "CVoriginal"="chromVAR::differentialDeviations", "CV"="chromVAR(z)>limma",
    "CVqt"="chromVAR(z)>Qt>limma", "CVcentered"="chromVAR(z)>center>limma",
    "CVnorm"="chromVAR(z)>scale>limma", MLzero="monaLisa.vsZero",
    MLsp="monaLisa.vsOthers+spearman", MLfewerBins="monaLisa.vsOthers(fewer bins)",
    MLStabSel="monaLisa.StabSel", MLlower="monaLisa.vsOthers(smaller zeroBin)",
    ML="monaLisa.vsOthers", "regreg"="Lasso-lm", "regregR"="Ridge-lm",
    VIPER="VIPER(scores)>limma", VIPERb="VIPER(binary)>limma",
    msVIPER="msVIPER(scores)", msVIPERb="msVIPER(binary)", ulmGC="ulm+GC", 
    ulmB="ulm(binary)", ulm="ulm(scores)", MBA="insertionModel", GSEA="fGSEA",
    CVdev="chromVAR(deviations)>limma", CVdevCentered="chromVAR(deviations)>center>limma",
    CVdevNorm="chromVAR(deviations)>scale>limma", CVdevqt="chromVAR(deviations)>Qt>limma"
  )
  for(i in names(renaming)) x <- replace(x, x==i, renaming[[i]])
  x <- gsub("decoupleR","decoupleR:",x)
  x <- gsub("::",":",x,fixed=TRUE)
  x
}

filterMethods <- function(x, meth=getMainMethods(), rename=FALSE){
  w <- which(x$method=="CV" | x$method=="chromVAR>scale>limma")
  CVcpu <- setNames(x$cpu[w], x$dataset[w])
  CVelapsed <- setNames(x$elapsed[w], x$dataset[w])
  w <- which(grepl("^CV|^chromVAR",x$method) & is.na(x$cpu))
  x$cpu[w] <- CVcpu[x$dataset[w]]
  x$elapsed[w] <- CVelapsed[x$dataset[w]]
  x$method <- as.character(x$method)
  renames <- meth
  if(rename) x <- renameMethods(x)
  x <- x[which(x$method %in% names(renames)),]
  renameMethods(x, renames)
}

getTopMethods <- function(){
  c(`chromVAR::differentialDeviations` = "chromVAR", 
    `chromVAR(z)>Qt>limma` = "chromVAR-adjusted",
    monaLisa.vsOthers = "monaLisa", monaLisa.StabSel = "monaLisa.StabSel", 
    `msVIPER(scores)` = "msVIPER", `VIPER(binary)>limma` = "VIPER", 
    `decoupleR:consensus` = "decoupleR", fGSEA = "fGSEA")
}

getMainMethods <- function(){
  c(`chromVAR::differentialDeviations` = "chromVAR", 
    `chromVAR(z)>Qt>limma` = "chromVAR-adjusted", diffTF="diffTF",
    monaLisa.vsOthers = "monaLisa", monaLisa.StabSel = "monaLisa.StabSel", 
    `msVIPER(scores)` = "msVIPER", `VIPER(binary)>limma` = "VIPER", 
    `decoupleR:consensus` = "decoupleR", fGSEA = "fGSEA",
    insertionModel="insertionModel", ulmB="ulm","Lasso-lm"="Lasso",
    "BaGFootLike"="BaGFootLike")
}