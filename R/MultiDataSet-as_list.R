as.list.MultiDataSet <- function(x) {
    ll <- lapply(names(x), function(dtype) {
        elm <- assayDataElementNames(x[[dtype]])[1]
        assayDataElement(x[[dtype]], elm)
    })
    names(ll) <- names(x)
    return(ll)
}
#' @describeIn MultiDataSet Returns a list with the first matrix of each
#' dataset.
#' @aliases as.list
setMethod("as.list", "MultiDataSet", as.list.MultiDataSet)

#' Method to convert \code{MultiDataSet} to a list of matrices.
#'
#' @rdname as_list
#' @aliases as_list
#' @param object \code{MultiDataSet} that will be filled.
#' @return A lsit of matrices, one per dataset
#' @exportMethod as_list
setGeneric("as_list", function(object) standardGeneric("as_list"))
setMethod(
    f = "as_list",
    signature = "MultiDataSet",
    definition = function(object) {
        ll <- lapply(names(object), function(dtype) {
            elm <- assayDataElementNames(object[[dtype]])[1]
            assayDataElement(object[[dtype]], elm)
        })
        names(ll) <- names(object)
        return(ll)
    }
)