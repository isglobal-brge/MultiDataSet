#' @describeIn MultiDataSet Get the name datasets having rowRanges
#' @aliases MultiDataSet-methods
setMethod(
    f = "rowRangesElements",
    signature = "MultiDataSet",
    definition = function(object) {
        names(object@rowRanges[!is.na(object@rowRanges)])
    }
)