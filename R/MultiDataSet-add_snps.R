#' @describeIn MultiDataSet Method to add a slot of SNPs to \code{MultiDataSet}.
#' @aliases MultiDataSet-methods
#' @param snpSet \code{SnpSet} to be used to fill the slot.
setMethod(
    f = "add_snps",
    signature = c("MultiDataSet", "SnpSet"),
    definition = function(object, snpSet, ...) {
        
        fet <- fData(snpSet)
        if (!all(c("position", "chromosome") %in% colnames(fet))){
            stop("fData of methySet must contain columns chromosome and position")
        }
        range <- GenomicRanges::makeGRangesFromDataFrame(fet, seqnames.field = "chromosome", 
                                                         start.field = "position", end.field = "position")
        
        object <- add_eset(object, snpSet, dataset.type = "snps", GRanges = range, ...)
        return(object)
    }
)
