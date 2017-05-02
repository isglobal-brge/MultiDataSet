#' @describeIn MultiDataSet Get sample names
#' @aliases MultiDataSet-methods
setMethod(
  f = "sampleNames",
  signature = "MultiDataSet",
  definition = function(object) {
    tables <- names(object)
    if (is.null(tables)){
      return(NULL)
    }
    ## CHECK THAT WE RETURN CHARACTER, NOT FACTOR!
    res <- lapply(tables, function(x) as.character(object@phenoData[[x]]$main$id))
    names(res) <- tables
    res
  }
)