#' @describeIn MultiDataSet Method to add a \code{matrix} to \code{MultiDataSet}.
#' @aliases MultiDataSet-methods
setMethod(
    f = "add_table",
    signature = c("MultiDataSet", "matrix"),
    definition = function(object, set, dataset.type, dataset.name = NULL, warnings = TRUE, 
                          overwrite = FALSE) {
        dataset.name <- paste(c(dataset.type, dataset.name), collapse = "+")
        
        if(dataset.name %in% names(object)){
            if (!overwrite){
                stop("There is already an object in this slot. Set overwrite = TRUE to overwrite the previous set.")
            }
            if (warnings) {
                warning("Slot '", dataset.name, "' is already set in 'MultiDataSet'. Previous content will be overwritten.")
            }
        }
        
        if (is.null(colnames(set))){
            stop("Set must contain colnames.")
        }
        if (sum(duplicated(colnames(set))) > 0){
            stop("Colnames of set must be unique.")
        }

        if (is.null(rownames(set))){
            stop("Set must contain rownames.")
        }
        if (sum(duplicated(rownames(set))) > 0){
            stop("Rownames of set must be unique.")
        }
        
        env <- new("environment")
        assign("mat", set, env)
        
        object@assayData[[dataset.name]] <- env
        
        pheno <- as(data.frame(id = colnames(set)), "AnnotatedDataFrame")
        rownames(pheno) <- pheno$id

        feats <- as(data.frame(id = rownames(set)), "AnnotatedDataFrame")
        rownames(feats) <- feats$id
        
        object@phenoData[[dataset.name]] <- pheno
        object@featureData[[dataset.name]] <- feats
        object@rowRanges[[dataset.name]] <- NA
        
        returnfunc <- function(env, phe, fet) {
            env$mat
        }
        
        object@return_method[[dataset.name]] <- returnfunc
        return(object)
    }
)
