getInteractors <- function(gene,
                            species=c(9606,10090), # first is human, second is mouse
                            accessKey="ffe973bbe4af1828a826453284dd3d40") # replace the XXX with your access key
{
  stopifnot(length(gene)==1)
  species <- match.arg(species)
  d <- read.delim(url(paste0("https://webservice.thebiogrid.org/interactions/?geneList=",gene,"&searchNames=true&taxId=",species,"&includeInteractors=true&accessKey=",accessKey)), header=FALSE)
  d <- d[which(d$V8==gene | d$V9==gene),]
  sort(setdiff(unique(c(d$V8,d$V9)),gene))
}
