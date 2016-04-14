#' @describeIn MultiDataSet Retrieve all assay data blocks.
#' @aliases assayData
setMethod(
  f = "assayData",
  signature = "MultiDataSet", 
  definition = function(object) {
    return(object@assayData)
  }
)