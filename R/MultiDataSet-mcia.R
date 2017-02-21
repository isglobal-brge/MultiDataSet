#' @describeIn MultiDataSet Apply mcia integration method to a MultiDataSet object
#' @aliases MultiDataSet-methods
setMethod(
    f = "w_mcia",
    signature = "MultiDataSet",
    definition = function(object, ...) {
        if(!"omicade4" %in% utils::installed.packages()[,"Package"]) {
            stop("R package 'omicade4' not found. Please, install it before running this function.")
        }
        object <- commonSamples(object)
        omicade4::mcia(as.list(object), ...)
    }
)