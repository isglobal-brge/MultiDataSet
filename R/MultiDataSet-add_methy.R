#' @describeIn MultiDataSet Method to add a slot of methylation to \code{MultiDataSet}.
#' @aliases MultiDataSet-methods
#' @param methySet \code{MethylationSet} to be used to fill the slot.
#' @param ... Further arguments passed to \code{add_eset}.
setMethod(
    f = "add_methy",
    signature = c("MultiDataSet", "MethylationSet"),
    definition = function(object, methySet, ...) {
        
        fet <- fData(methySet)
        if (!all(c("position", "chromosome") %in% colnames(fet))){
            stop("fData of methySet must contain columns chromosome and position")
        }
        range <- GenomicRanges::makeGRangesFromDataFrame(fet, seqnames.field = "chromosome", 
                                                         start.field = "position", end.field = "position")

        object <- add_eset(object, methySet, dataset.type = "methylation", GRanges = range, ...)
        
        return(object)
    }
)

#' @describeIn MultiDataSet Method to add a slot of methylation to \code{MultiDataSet}.
#' @aliases MultiDataSet-methods
setMethod(
    f = "add_methy",
    signature = c("MultiDataSet", "RatioSet"),
    definition = function(object, methySet, ...) {
        
        fet <- fData(methySet)
        if (!all(c("position", "chromosome") %in% colnames(fet))){
            stop("fData of methySet must contain columns chromosome and position")
        }

        range <- GenomicRanges::makeGRangesFromDataFrame(fet, seqnames.field = "chromosome", 
                                                         start.field = "position", end.field = "position")
        
        
        object <- add_eset(object, methySet, dataset.type = "methylation", GRanges = range, ...)
        
        return(object)
    }
)
