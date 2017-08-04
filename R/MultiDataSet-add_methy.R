#' @describeIn MultiDataSet Method to add a slot of methylation to 
#' \code{MultiDataSet} from a \code{GenomicRatioSet}.
#' @aliases MultiDataSet-methods
setMethod(
    f = "add_methy",
    signature = c("MultiDataSet", "GenomicRatioSet"),
    definition = function(object, methySet, ...) {
        
        object <- add_rse(object, methySet, dataset.type = "methylation")
        
        return(object)
    }
)