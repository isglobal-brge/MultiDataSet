#' @describeIn MultiDataSet Method to add an \code{eSet} to \code{MultiDataSet}.
#' @aliases MultiDataSet-methods
#' @param set Object derived from \code{eSet} to be used to fill the slot.
#' @param dataset.name Character with the specific name for this set (NULL by default). It is useful when there 
#' @param warnings Logical to indicate if warnings will be displayed.
#' @param GRanges \code{GenomicRanges} to be included in rowRanges slot. 
setMethod(
    f = "add_eset",
    signature = c("MultiDataSet", "eSet"),
    definition = function(object, set, dataset.type, dataset.name = NULL, warnings = TRUE, overwrite = FALSE, 
                          GRanges) {
        validObject(set)
        dataset.name <- paste(c(dataset.type, dataset.name), collapse = "+")
        
        if(dataset.name %in% names(object)){
            if (!overwrite){
                stop("There is already an object in this slot. Set overwrite = TRUE to overwrite the previous set.")
            }
            if (warnings) {
                warning("Slot '", dataset.name, "' is already set in 'MultiDataSet'. Previous content will be overwritten.")
            }
        }
        
        
        if (missing(GRanges)){
            GRanges <- GenomicRanges::makeGRangesFromDataFrame(fData(set))
            names(GRanges) <- rownames(fData(set))
        } 
        if (!is(GRanges, "GenomicRanges")){
            if (!is.na(GRanges)){
               stop("GRanges should be a GenomicRanges or NA.")
            }
        }
        
        object@assayData[[dataset.name]] <- assayData(set)
        
        pheno <- phenoData(set)
        if (!"id" %in% colnames(pheno)){
            warning("No id column found in pData. The id will be equal to the sampleNames")
            pheno$id <- rownames(pheno)
        }
        object@phenoData[[dataset.name]] <- pheno
        object@featureData[[dataset.name]] <- featureData(set)
        object@rowRanges[[dataset.name]] <- GRanges
        
        returnfunc <- function(env, phe, fet) {
            new(class(set), assayData = env, phenoData = phe, 
                featureData = fet)
        }
        
        object@return_method[[dataset.name]] <- returnfunc
        return(object)
    }
)
