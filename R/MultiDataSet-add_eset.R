#' @describeIn MultiDataSet Method to add an \code{eSet} to \code{MultiDataSet}.
#' @aliases MultiDataSet-methods
#' @param set Object derived from \code{eSet} to be used to fill the slot.
#' @param dataset.name Character with the specific name for this set (NULL by default). It is useful when there 
#' @param warnings Logical to indicate if warnings will be displayed.
#' @param useGRanges Logical to indicate if the featureData will be transformed to GenomicRanges
setMethod(
  f = "add_eset",
  signature = c("MultiDataSet", "eSet"),
  definition = function(object, set, dataset.type, dataset.name = NULL, warnings = TRUE, overwrite = FALSE, useGRanges = FALSE) {
    #validObject(set)
    dataset.name <- paste(c(dataset.type, dataset.name), collapse = "+")
    
#     if (!"id" %in% varLabels(set)){
#         warning("pData of set should contain a column called id. sampleNames will be used.")
#         set$id <- sampleNames(set)
#     }
    
    if(dataset.name %in% names(object)){
      if (!overwrite){
        stop("There is already an object in this slot. Set overwrite = TRUE to overwrite the previous set.")
      }
      if (warnings) {
        warning("Slot '", dataset.name, "' is already set in 'MultiDataSet'. Previous content will be overwritten.")
      }
    }
    
    
#     if (useGRanges){
#       if (!c("chromosome") %in% fvarLabels(set)){
#         stop("Set must contain a fData with column chromosome.")
#       }
#     
#       if (all(c("position", "start") %in% fvarLabels(set))){
#         stop("Set cannot contain a fData with columns position and start. Only one of them is allowed.")
#       }
#     
#       if (!any(c("position", "start") %in% fvarLabels(set))){
#         stop("Set must contain a fData with columns position or start.")
#       }
#     
#       if ("position" %in% fvarLabels(set)){
#         if (warnings){
#           warning("Column position of set will be renamed to start.")
#         }
#         fvarLabels(set)[fvarLabels(set) == "position"] <- "start"
#       }
#       
#       if (!"end" %in% fvarLabels(set)){
#         fData(set)$end <- fData(set)$start
#       }
#       
#       fdata <- GenomicRanges::makeGRangesFromDataFrame(
#         fData(set), seqnames.field = "chromosome", #start.field = "start", end.field = "end",
#         keep.extra.columns = TRUE)
#       names(fdata) <- rownames(fData(set))
#       
#     } else {
#       warning("The usage of AnnotationDataSet is not recomended and will be removed. Some functionallity can not run properly.")
      fdata <- featureData(set)
    # }
  
    object@assayData[[dataset.name]] <- assayData(set)
    
    pheno <- phenoData(set)
    if (!"id" %in% colnames(pheno)){
      warning("No id column found in pData. The id will be equal to the sampleNames")
      pheno$id <- rownames(pheno)
    }
    object@phenoData[[dataset.name]] <- pheno
    object@featureData[[dataset.name]] <- fdata
    
    returnfunc <- function(env, phe, fet) {
      #fet <- as.data.frame(fet)
      #colnames(fet)[colnames(fet) == "seqnames"] <- "chromosome"
      new(class(set), assayData = env, phenoData = phe, 
          featureData = fet)
    }
    
    object@return_method[[dataset.name]] <- returnfunc
    return(object)
  }
)
