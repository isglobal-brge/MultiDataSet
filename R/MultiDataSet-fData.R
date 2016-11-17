#' @describeIn MultiDataSet Retrieve information on features.
#' @aliases fData
setMethod(
  f = "fData",
  signature = "MultiDataSet", 
  definition = function(object) {
    return(lapply(object@featureData, pData))
  }
)

#' @describeIn MultiDataSet Retrieve information on features.
#' @aliases featureData
setMethod(
    f = "featureData",
    signature = "MultiDataSet", 
    definition = function(object) {
        return(object@featureData)
    }
)