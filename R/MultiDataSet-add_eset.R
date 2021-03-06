#' @describeIn MultiDataSet Method to add an \code{eSet} to \code{MultiDataSet}.
#' @aliases MultiDataSet-methods
#' @param set Object derived from \code{eSet} to be used to fill the slot.
#' @param dataset.name Character with the specific name for this set (NULL by default). It is useful when there 
#' @param sample.tables Character with the names of the slots with sample data besides phenoData.
#' @param feature.tables Character with the names of the slots with feature data besides featureData.
#' @param warnings Logical to indicate if warnings will be displayed.
#' @param GRanges \code{GenomicRanges} to be included in rowRanges slot. 
setMethod(
    f = "add_eset",
    signature = c("MultiDataSet", "eSet"),
    definition = function(object, set, dataset.type, dataset.name = NULL, 
                          sample.tables = "protocolData", feature.tables = NULL,
                          warnings = TRUE, overwrite = FALSE, 
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
        
        

        
        object@assayData[[dataset.name]] <- assayData(set)
        
        pheno <- phenoData(set)
        if (!"id" %in% colnames(pheno)){
            if (warnings) {
                warning("No id column found in pData. The id will be equal to the sampleNames")
            }
            pheno$id <- rownames(pheno)
            
        }
        
        extra <- attributes(set)
        
        object@phenoData[[dataset.name]] <- list(main = pheno)
        if (!is.null(sample.tables)) {
            object@phenoData[[dataset.name]] <- 
                c(object@phenoData[[dataset.name]],  extra[sample.tables])
        }
        
        feat <- featureData(set)
        object@featureData[[dataset.name]] <- list(main = feat)
        if (!is.null(feature.tables)) {
            object@featureData[[dataset.name]] <- 
                c(object@featureData[[dataset.name]], extra[feature.tables])
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
        object@rowRanges[[dataset.name]] <- GRanges
        
        
        extranames <- names(extra)[!names(extra) %in% 
                                       c("assayData", "phenoData", "featureData", "class", sample.tables, feature.tables)]
        extra <- extra[extranames]
        extra <- extra[!sapply(extra, class) == "name"]
        object@extraData[[dataset.name]] <- extra
        
        returnfunc <- function(env, phe, fet, extra) {
            attr <- list(Class = class(set), 
                         assayData = env, 
                         phenoData = phe$main,
                         featureData = fet$main)
            attr <- c(attr, phe[-1], fet[-1], extra)
            do.call("new", attr)
        }
        
        object@return_method[[dataset.name]] <- returnfunc
        return(object)
    }
)
