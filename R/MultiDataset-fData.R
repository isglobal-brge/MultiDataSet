#' @describeIn MultiDataSet Retrieve information on features.
#' @aliases fData
setMethod(
  f = "fData",
  signature = "MultiDataSet", 
  definition = function(object) {
    return(object@featureData)
  }
)