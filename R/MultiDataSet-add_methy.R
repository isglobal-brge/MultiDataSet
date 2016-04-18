#' @describeIn MultiDataSet Method to add a slot of methylation to \code{MultiDataSet}.
#' @aliases MultiDataSet-methods
#' @param methySet \code{MethylationSet} to be used to fill the slot.
#' @param ... Further arguments passed to \code{add_eset}.
setMethod(
    f = "add_methy",
    signature = c("MultiDataSet", "MethylationSet"),
    definition = function(object, methySet, ...) {
        
        fet <- fData(methySet)
        if (!"position" %in% colnames(fet)){
            stop("methySet must contain a fData with a column called position.")
        }
        colnames(fet)[colnames(fet) == "position"] <- "start"
        
        if (!"end" %in% colnames(fet)){
            fet$end <- fet$start
        }
        range <- GenomicRanges::makeGRangesFromDataFrame(fet)
        
        
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
        if (!"position" %in% colnames(fet)){
            stop("methySet must contain a fData with a column called position.")
        }
        colnames(fet)[colnames(fet) == "position"] <- "start"
        
        if (!"end" %in% colnames(fet)){
            fet$end <- fet$start
        }
        range <- GenomicRanges::makeGRangesFromDataFrame(fet)
        
        
        object <- add_eset(object, methySet, dataset.type = "methylation", GRanges = range, ...)
        
        return(object)
    }
)
