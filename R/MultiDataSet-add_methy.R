#' @describeIn MultiDataSet Method to add a slot of methylation to 
#' \code{MultiDataSet} from a \code{GenomicRatioSet}.
#' @aliases MultiDataSet-methods
#' @param methySet \code{GenomicRatioSet} to be used to fill the slot.
#' @param ... Further arguments passed to add_rse or add_se
setMethod(
    f = "add_methy",
    signature = c("MultiDataSet", "GenomicRatioSet"),
    definition = function(object, methySet, ...) {
        
        object <- add_rse(object, methySet, dataset.type = "methylation")
        
        return(object)
    }
)