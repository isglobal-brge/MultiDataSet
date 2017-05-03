#' @describeIn MultiDataSet Get the name of the datasets that have rowRanges
#' @aliases MultiDataSet-methods
setMethod(
    f = "rowRangesElements",
    signature = "MultiDataSet",
    definition = function(object) {
        res <- names(object@rowRanges[!is.na(object@rowRanges)])
        if (is.null(res)) {
            res <- character(0)
        }
        return(res)
    }
)