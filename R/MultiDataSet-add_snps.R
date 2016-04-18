#' @describeIn MultiDataSet Method to add a slot of SNPs to \code{MultiDataSet}.
#' @aliases MultiDataSet-methods
#' @param snpSet \code{SnpSet} to be used to fill the slot.
setMethod(
    f = "add_snps",
    signature = c("MultiDataSet", "SnpSet"),
    definition = function(object, snpSet, ...) {
        
        fet <- fData(snpSet)
        if (!"position" %in% colnames(fet)){
            stop("snpSet must contain a fData with a column called position.")
        }
        colnames(fet)[colnames(fet) == "position"] <- "start"
        
        if (!"end" %in% colnames(fet)){
            fet$end <- fet$start
        }
        range <- GenomicRanges::makeGRangesFromDataFrame(fet)
        object <- add_eset(object, snpSet, dataset.type = "snps", GRanges = range, ...)
        return(object)
    }
)
