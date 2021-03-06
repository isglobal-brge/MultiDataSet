as.list.MultiDataSet <- function(x) {
    ll <- lapply(names(x), function(dtype) {
        elm <- assayDataElementNames(assayData(x)[[dtype]])[1]
        assayDataElement(assayData(x)[[dtype]], elm)
    })
    names(ll) <- names(x)
    return(ll)
}
#' @describeIn MultiDataSet Returns a list with the first matrix of each
#' dataset.
#' @aliases as.list
#' @export
setMethod("as.list", "MultiDataSet", as.list.MultiDataSet)
