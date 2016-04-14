#' @describeIn MultiDataSet Method to add a slot of SNPs to \code{MultiDataSet}.
#' @aliases MultiDataSet-methods
#' @param snpSet \code{SnpSet} to be used to fill the slot.
setMethod(
  f = "add_snps",
  signature = c("MultiDataSet", "SnpSet"),
  definition = function(object, snpSet, ...) {
    object <- add_eset(object, snpSet, dataset.type = "snps", ...)
    return(object)
  }
)
