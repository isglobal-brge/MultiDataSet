#' @describeIn MultiDataSet Method to add a slot of expression to \code{MultiDataSet}.
#' @aliases MultiDataSet-methods
#' @param gexpSet \code{ExpressionSet} to be used to fill the slot.
setMethod(
    f = "add_genexp",
    signature = c("MultiDataSet", "ExpressionSet"),
    definition = function(object, gexpSet, ...) {
        
        fet <- fData(gexpSet)
        if (!"end" %in% colnames(fet)){
            fet$end <- fet$start
        }
        range <- GenomicRanges::makeGRangesFromDataFrame(fet)
        names(range) <- featureNames(gexpSet)
        object <- add_eset(object, gexpSet, dataset.type = "expression", GRanges = range, ...)
        
        return(object)
    }
)

#' @describeIn MultiDataSet Method to add a slot of (RNASeq) expression to \code{MultiDataSet}.
#' @aliases MultiDataSet-methods
#' @param rnaSet \code{ExpressionSet} to be used to fill the slot.
setMethod(
    f = "add_rnaseq",
    signature = c("MultiDataSet", "ExpressionSet"),
    definition = function(object, rnaSet, ...) {
        
        fet <- fData(rnaSet)
        if (!"end" %in% colnames(fet)){
            fet$end <- fet$start
        }
        range <- GenomicRanges::makeGRangesFromDataFrame(fet)
        object <- add_eset(object, rnaSet, dataset.type = "rnaseq", GRanges = range, ...)
        
        return(object)
    }
)

