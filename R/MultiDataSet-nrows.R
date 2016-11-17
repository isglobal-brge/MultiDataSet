#' @describeIn MultiDataSet Get number of features of each set
#' @aliases nrows
setMethod(
    f = "nrows",
    signature = c("MultiDataSet"),
    definition = function(object) {
        return(vapply(object@featureData, nrow, numeric(1)))
    } 
)