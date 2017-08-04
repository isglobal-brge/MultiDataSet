#' @describeIn MultiDataSet Get number of samples of each set
#' @aliases ncol
setMethod(
    f = "ncol",
    signature = c("MultiDataSet"),
    definition = function(x) {
        return(vapply(x@phenoData, function(y) nrow(y[[1]]), numeric(1)))
    } 
)