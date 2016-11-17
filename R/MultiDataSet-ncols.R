#' @describeIn MultiDataSet Get number of samples of each set
#' @aliases ncols
setMethod(
    f = "ncols",
    signature = c("MultiDataSet"),
    definition = function(object) {
        return(vapply(object@phenoData, nrow, numeric(1)))
    } 
)