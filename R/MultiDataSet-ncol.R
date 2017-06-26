#' @describeIn MultiDataSet Get number of samples of each set
#' @aliases ncol
setMethod(
    f = "ncol",
    signature = c("MultiDataSet"),
    definition = function(x) {
        return(vapply(x@phenoData, nrow, numeric(1)))
    } 
)