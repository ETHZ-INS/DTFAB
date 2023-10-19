getInteractors <- function(gene,
                            species=c(9606,10090), # first is human, second is mouse
                            accessKey="ffe973bbe4af1828a826453284dd3d40") # replace the XXX with your access key
{
  stopifnot(length(gene)==1)
  species <- match.arg(species)
  tryCatch({
    d <- read.delim(url(paste0("https://webservice.thebiogrid.org/interactions/?geneList=",gene,"&searchNames=true&taxId=",species,"&includeInteractors=true&accessKey=",accessKey)), header=FALSE)
    d <- read.delim(url(paste0("https://webservice.thebiogrid.org/interactions/?geneList=",gene,"&searchNames=true&taxId=",species,"&includeInteractors=true&accessKey=",accessKey)), header=FALSE)
    d <- d[which(d$V8==gene | d$V9==gene),]
    sort(setdiff(unique(c(d$V8,d$V9)),gene))
  }, error=function(e){
    warning("Could not fetch any interactor for ", gene)
    return(c())
  })
}

getAllInteractors <- function(ds2){
  lapply(setNames(names(ds2),names(ds2)), FUN=function(x){
    if(is.null(ds2[[x]]$truth)){
      truth <- x
    }else{
      truth <- ds2[[x]]$truth
    }
    if(ds2[[x]]$species=="m") truth <- tools::toTitleCase(tolower(truth))
    res <- lapply(truth,
                  species=ifelse(ds2[[x]]$species=="h","9606","10090"),
                  FUN=getInteractors)
    Sys.sleep(1)
    res <- res[which(lengths(res)>0)]
    unique(unlist(res))
  })
}
