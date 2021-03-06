#' @describeIn MultiDataSet Get number of features of each set
#' @aliases nrow
setMethod(
    f = "nrow",
    signature = c("MultiDataSet"),
    definition = function(x) {
        return(vapply(x@featureData, function(y) nrow(y[[1]]), numeric(1)))
    } 
)