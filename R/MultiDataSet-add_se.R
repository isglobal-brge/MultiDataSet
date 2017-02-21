#' @describeIn MultiDataSet Method to add a \code{SummarizedExperiment} to \code{MultiDataSet}.
#' @aliases MultiDataSet-methods
#' @return \code{MultiDataSet}
setMethod(
    f = "add_se",
    signature = c("MultiDataSet", "SummarizedExperiment"),
    definition = function(object, set, dataset.type, dataset.name, warnings = TRUE, overwrite = FALSE, 
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
        
        pheno <- Biobase::AnnotatedDataFrame(as.data.frame(SummarizedExperiment::colData(set)))
        if (!"id" %in% colnames(pheno)){
            warning("No id column found in rowRanges. The id will be equal to the sampleNames")
            pheno$id <- rownames(pheno)
        }
        object@phenoData[[dataset.name]] <- pheno
        object@featureData[[dataset.name]] <- 
            Biobase::AnnotatedDataFrame(as.data.frame(SummarizedExperiment::rowData(set)))
        
        if (missing(GRanges)){
            GRanges <- GenomicRanges::makeGRangesFromDataFrame(SummarizedExperiment::colData(set))
            names(GRanges) <- rownames(SummarizedExperiment::rowData(set))
        } 
        if (!is(GRanges, "GenomicRanges")){
            if (!is.na(GRanges)){
                stop("GRanges should be a GenomicRanges or NA.")
            }
        }
        
        object@rowRanges[[dataset.name]] <- GRanges
        
        returnfunc <- function(env, phe, fet) {
            assays <- SummarizedExperiment::Assays(as.list(env))
            new(class(set), assays = assays, colData = S4Vectors::DataFrame(as(phe, "data.frame")), 
                rowRanges = GenomicRanges::makeGRangesFromDataFrame(as(fet, "data.frame"), keep.extra.columns=TRUE), 
                elementMetadata = S4Vectors::DataFrame(matrix(nrow = nrow(assays), ncol = 0 )))
        }
        
        object@return_method[[dataset.name]] <- returnfunc
        return(object)
    }
)
