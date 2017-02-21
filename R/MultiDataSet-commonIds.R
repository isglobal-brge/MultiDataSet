#' @describeIn MultiDataSet Get the name of the ids common to all datasets
#' @aliases MultiDataSet-methods
setMethod(
  f = "commonIds",
  signature = "MultiDataSet",
  definition = function(object) {
    return(Reduce(intersect, sampleNames(object)))
  }
)
