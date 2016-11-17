#' @describeIn MultiDataSet Returns the dimensions of the sets
#' @aliases dims
setMethod(
    f = "dims",
    signature = c("MultiDataSet"),
    definition = function(object) {
        dims <- lapply(names(object), function(name) c(Features = nrow(pData(object@featureData[[name]])), 
                                                       Samples =  nrow(pData(object@phenoData[[name]]))))
        names(dims) <- names(object)
        return(dims)
    } 
)