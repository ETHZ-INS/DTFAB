#' Get benchmark metrics for all dataset
#'
#' @param datasets A list of dataset objects (see getDatasets() )
#' @param path The path in which the datasets' results folders are
#' @param resin The name of the results subfolder
#' @param interactors A named list of interactors per TF
#'
#' @return A data.frame of compiled metrics
compileBenchmark <- function(datasets, rootPath=".", resin="runATAC_results", interactors){
  res <- list()
  for(dn in names(datasets)){
    ds <- datasets[[dn]]
    if(!is.null(ds$folder)) dn <- ds$folder
    print(dn)
    tryCatch(
      res[[dn]] <- getBenchmarkMetrics(ds, file.path(rootPath, dn), interactors=interactors),
      error=function(e) message(e))
  }
  dplyr::bind_rows(res, .id="dataset")
}


#' Get benchmark metrics for one dataset
#'
#' @param dataset A dataset object (see getDatasets() )
#' @param path The path in which the dataset's results are
#' @param resin The name of the results subfolder
#' @param interactors A named list of interactors per TF
#'
#' @return A data.frame of compiled metrics for one dataset
getBenchmarkMetrics <- function(dataset, path=head(dataset$truth,1),
                                resin="runATAC_results", interactors){
  # grab runtimes
  fl <- list.files(file.path(path, resin, "raw"), full=TRUE)
  names(fl) <- gsub("_raw\\.rds$","",basename(fl))
  rt <- sapply(fl, FUN=function(x){
    x <- readRDS(x)
    if(!is.null(names(x)) && "runtime" %in% names(x))
      return(x$runtime)
    return(x[[2]])
  })
  rt <- rbind(elapsed=rt["elapsed",], cpu=colSums(rt[which(row.names(rt)!="elapsed"),]))
  rt <- as.data.frame(t(rt))
  if(any(row.names(rt)=="decoupleR")){
    # grab times for individual methods
    dt <- readRDS(file.path(path, resin, "raw","decoupleR_raw.rds"))$runtime2
    rt <- rbind(rt, data.frame(row.names=paste0("decoupleR",names(dt)),
                               elapsed=dt, cpu=dt))
  }
  # get actual results
  fl <- list.files(file.path(path, resin, "with_pvalues"), full=TRUE)
  names(fl) <- gsub("\\.rds$","",basename(fl))
  res <- lapply(fl, readRDS)
  res <- res[sapply(res, FUN=function(x) isTRUE(nrow(x)>1))]
  res <- dplyr::bind_rows(lapply(res, truth=dataset$truth, interactors=interactors,
                                 FUN=.getBenchmarkMetrics), .id="method")
  res$elapsed <- rt[res$method, "elapsed"]
  res$cpu <- rt[res$method, "cpu"]
  res
}

.getBenchmarkMetrics <- function(x, truth, interactors){
  cofactors <- unique(c(truth, unlist(interactors[truth])))
  allTFs <- row.names(x)
  cofactors <- intersect(cofactors, allTFs)
  cappedNCof <- min(length(cofactors),100)
  optimalAUC <- sum(cumsum(c(rep(1,cappedNCof), rep(0,100-cappedNCof)))/seq_len(100))
  w <- head(which(row.names(x) %in% truth),1)
  x$isCofactor <- row.names(x) %in% cofactors
  auc <- sum(cumsum(head(x$isCofactor,100))/seq_len(100))
  if(is.null(x$padj)) x$padj <- x$adj.P.Val
  xsig <- x[which(x$padj<0.05),,drop=FALSE]
  data.frame(rank=w, trueQ=x[w,"padj"],
             FDR=sum(!xsig$isCofactor)/nrow(xsig),
             AUC=auc, relAUC=auc/optimalAUC)
}


# fetch a list of interactors for all requested TFs
getAllInteractors <- function(datasets, extra=c("CEBPB", "MAZ", "ZNF143", "NR3C1")){
  names(ds2) <- ds2 <- extra
  ds2 <- lapply(ds2, FUN=function(x) list(truth=x, species="h"))
  ints <- lapply(c(datasets,ds2), FUN=function(x){
    tf <- x$truth
    sp <- "9606"
    if(x$species=="m"){
      tf <- tools::toTitleCase(tolower(tf))
      sp <- "10090"
    }
    unique(toupper(unlist(lapply(tf, species=sp, FUN=getInteractors))))
  })
  return(ints)
}
