#' @describeIn MultiDataSet Retrieve information on experimental phenotypes.
#' @aliases pData
setMethod(
  f = "pData",
  signature = "MultiDataSet", 
  definition = function(object) {
    return(object@phenoData)
  }
)