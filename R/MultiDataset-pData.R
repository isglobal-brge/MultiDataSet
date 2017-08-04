#' @describeIn MultiDataSet Retrieve information on experimental phenotypes
#' @aliases pData
setMethod(
  f = "pData",
  signature = "MultiDataSet", 
  definition = function(object) {
    return(lapply(object@phenoData, function(x) pData(x$main)))
  }
)


#' @describeIn MultiDataSet Retrieve information on experimental phenotypes
#' @aliases phenoData
setMethod(
    f = "phenoData",
    signature = "MultiDataSet", 
    definition = function(object) {
        return(object@phenoData)
    }
)