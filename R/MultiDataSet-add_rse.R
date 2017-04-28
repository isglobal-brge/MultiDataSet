#' @describeIn MultiDataSet Method to add a \code{RangedSummarizedExperiment} to \code{MultiDataSet}.
#' @aliases MultiDataSet-methods
#' @param object \code{MultiDataSet}
#' @param dataset.type Character with the type of data of the omic set (e.g. expression, methylation...)
#' @param overwrite Logical to indicate if the set stored in the slot will be overwritten.
#' @return \code{MultiDataSet}
setMethod(
    f = "add_rse",
    signature = c("MultiDataSet", "RangedSummarizedExperiment"),
    definition = function(object, set, dataset.type, dataset.name, 
                          sample.tables = NULL, feature.tables = "elementMetadata", 
                          warnings = TRUE, overwrite = FALSE) {
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
        
        env <- new("environment")
        lapply(names(SummarizedExperiment::assays(set)), 
               function(x) assign(x, SummarizedExperiment::assays(set)[[x]], env))
        object@assayData[[dataset.name]] <- env
        
        extra <- attributes(set)
        pheno <- Biobase::AnnotatedDataFrame(as.data.frame(SummarizedExperiment::colData(set)))
        if (!"id" %in% colnames(pheno)){
            warning("No id column found in colData. The id will be equal to the sampleNames")
            pheno$id <- rownames(pheno)
        }
        
        object@phenoData[[dataset.name]] <- list(main = pheno)
        if (!is.null(sample.tables)) {
            phetabs <- extra[sample.tables]
            for (tab in names(phetabs)){
                rownames(phetabs[[tab]]) <- rownames(pheno)
            }
            object@phenoData[[dataset.name]] <- c(object@phenoData[[dataset.name]], phetabs)
            
        }
        
        
        feat <- Biobase::AnnotatedDataFrame(as.data.frame(SummarizedExperiment::rowRanges(set)))
        object@featureData[[dataset.name]] <- list(main = feat)
        if (!is.null(feature.tables)) {
            feattabs <- extra[feature.tables]
            for (tab in names(feattabs)){
                rownames(feattabs[[tab]]) <- rownames(feat)
            }
            object@featureData[[dataset.name]] <- c(object@featureData[[dataset.name]], feattabs)
        }
        
        object@rowRanges[[dataset.name]] <- SummarizedExperiment::rowRanges(set)
        
        extra <- attributes(set)
        extranames <- names(extra)[!names(extra) %in% 
                                       c("assays", "colData", "rowRanges", "class", sample.tables, feature.tables)]
        extra <- extra[extranames]
        extra <- extra[!sapply(extra, class) == "name"]
        object@extraData[[dataset.name]] <- extra
        
        returnfunc <- function(env, phe, fet, extra) {
            assays <- SummarizedExperiment::Assays(as.list(env))
            attr <- list(Class = class(set), 
                         assays = assays, 
                         colData = S4Vectors::DataFrame(as(phe$main, "data.frame")),
                         rowRanges = GenomicRanges::makeGRangesFromDataFrame(as(fet$main, "data.frame"), 
                                                                             keep.extra.columns=TRUE))
            attr <- c(attr, phe[-1], fet[-1], extra)
            rownames(attr$elementMetadata) <- NULL
            do.call("new", attr)
        }
        
        object@return_method[[dataset.name]] <- returnfunc
        return(object)
    }
)
