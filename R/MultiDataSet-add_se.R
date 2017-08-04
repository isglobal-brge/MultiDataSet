#' @describeIn MultiDataSet Method to add a \code{SummarizedExperiment} to \code{MultiDataSet}.
#' @aliases MultiDataSet-methods
#' @return \code{MultiDataSet}
setMethod(
    f = "add_se",
    signature = c("MultiDataSet", "SummarizedExperiment"),
    definition = function(object, set, dataset.type, dataset.name, 
                          sample.tables = NULL, feature.tables = "elementMetadata",
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
            object@phenoData[[dataset.name]] <- 
                c(object@phenoData[[dataset.name]], extra[sample.tables])
            
        }
        
        
        feat <- Biobase::AnnotatedDataFrame(as.data.frame(SummarizedExperiment::rowData(set)))
        object@featureData[[dataset.name]] <- list(main = feat)
        if (!is.null(feature.tables)) {
            object@featureData[[dataset.name]] <- 
                c(object@featureData[[dataset.name]], extra[feature.tables])
        }
        
        if (missing(GRanges)){
            GRanges <- GenomicRanges::makeGRangesFromDataFrame(SummarizedExperiment::rowData(set))
            names(GRanges) <- rownames(SummarizedExperiment::rowData(set))
        } 
        if (!is(GRanges, "GenomicRanges")){
            if (!is.na(GRanges)){
                stop("GRanges should be a GenomicRanges or NA.")
            }
        }
        object@rowRanges[[dataset.name]] <- GRanges
        
        
        extranames <- names(extra)[!names(extra) %in% 
                                       c("assays", "colData", "rowRanges", "class", sample.tables, feature.tables)]
        extra <- extra[extranames]
        extra <- extra[!sapply(extra, class) == "name"]
        object@extraData[[dataset.name]] <- extra
        
        returnfunc <- function(env, phe, fet, extra) {
            assays <- SummarizedExperiment::Assays(as.list(env))
            attr <- list(Class = class(set), 
                         assays = assays, 
                         colData = S4Vectors::DataFrame(as(phe$main, "data.frame")))
            attr <- c(attr, phe[-1], fet[-1], extra)
            do.call("new", attr)
        }

        object@return_method[[dataset.name]] <- returnfunc
        return(object)
    }
)
