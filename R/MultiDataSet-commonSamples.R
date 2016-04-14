#' @describeIn MultiDataSet Get a MultiDataSet only with the samples present in all the tables
#' @aliases MultiDataSet-methods
setMethod(
  f = "commonSamples",
  signature = "MultiDataSet",
  definition = function(object) {
    samples <- commonIds(object)
    if (length(samples) == 0){
      stop("There are no samples present in all the datasets")
    }
    return(object[samples, ,])
  }
)
