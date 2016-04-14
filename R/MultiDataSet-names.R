#' @describeIn MultiDataSet Get the names of the slots.
#' @aliases MultiDataSet-methods, names
setMethod(
  f = "names",
  signature = "MultiDataSet",
  definition = function(x) {
    return(names(x@assayData))
  }
)
