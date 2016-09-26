#' @describeIn MultiDataSet Retrieve information on feature ranges.
#' @aliases rowRanges
#' @export rowRanges
setMethod(
    f = "rowRanges",
    signature = "MultiDataSet", 
    definition = function(x) {
        return(x@rowRanges)
    }
)