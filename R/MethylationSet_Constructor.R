#' @export
#' @name methylationSet
#' @rdname MethylationSet-class
#' @aliases MethylationSet MethylationSet-methods
#' 
#' @param betas Matrix of beta values
#' @param phenotypes Data.frame or AnnotatedDataFrame with the phenotypes
#' @param annotationDataFrame Data.frame or AnnotatedDataFrame with the annotation of 
#' the methylation sites.
#' @param annoString Character with the name of the annotation used. 
#' @return \code{MethylationSet}
#' @examples
#' showClass("MethylationSet")
methylationSet <- function(betas, phenotypes, annotationDataFrame, annoString = "custom"){
  
  
  if (is(phenotypes, "DataFrame")){
      phenotypes <- data.frame(phenotypes)
  }
  if (is(phenotypes, "data.frame")){
      if (nrow(phenotypes) != ncol(betas)){
          stop("phenotypes must have as rows as columns in betas.")
      }
      rownames(phenotypes) <- colnames(betas)
      phenotypes <- AnnotatedDataFrame(phenotypes)
  }
  if (!is(phenotypes, "AnnotatedDataFrame")){
    stop("phenotypes must be a data.frame or an AnnotatedDataFrame.")
  }

  if (is(annotationDataFrame, "data.frame")){
      if (nrow(annotationDataFrame) != nrow(betas)){
          stop("annotationDataFrame must have the same rows than betas.")
      }
      rownames(annotationDataFrame) <- rownames(betas)
        annotationDataFrame <- AnnotatedDataFrame(annotationDataFrame)
  } 
  if (!is(annotationDataFrame, "AnnotatedDataFrame")){
    stop("annotationDataFrame must be a data.frame or an AnnotatedDataFrame.")
  }
  meth <- assayDataNew(storage.mode = "lockedEnvironment",  meth = betas)
    
  set <- new(Class = "MethylationSet", assayData = meth, phenoData = phenotypes, 
             featureData = annotationDataFrame, annotation = annoString)
  
  return(set)
}