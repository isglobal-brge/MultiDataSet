#' @describeIn MultiDataSet Returns the dimensions of the sets
#' @aliases dims
setMethod(
    f = "dims",
    signature = c("MultiDataSet"),
    definition = function(x) {
        dims <- lapply(names(x), function(name) c(Features = nrow(pData(x@featureData[[name]]$main)), 
                                                       Samples =  nrow(pData(x@phenoData[[name]]$main))))
        names(dims) <- names(x)
        return(dims)
    } 
)