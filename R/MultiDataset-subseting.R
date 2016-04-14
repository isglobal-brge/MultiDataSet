#' @describeIn MultiDataSet Get a set from a slot
#' @aliases MultiDataSet-methods
#' @param x \code{MultiDataSet}
setMethod(
  f = "[[",
  signature = "MultiDataSet",
  definition = function(x, i) {
    if (i %in% names(x)) {
      set <- x@return_method[[i]](x@assayData[[i]], x@phenoData[[i]], x@featureData[[i]])
      validObject(set)
      return(set)
    }
  }
)

#' @describeIn MultiDataSet Subset a MultiDataSet
#' @aliases MultiDataSet-methods [
#' @param i Character corresponding to selected sample names. They should match the id column of phenoData.
#' @param j Character with the name of the selected tables.
#' @param drop Logical indicating if a \code{MultiDataSet} with one table must be converted to the original 
#' set.
#' @param k \code{GenomicRange} used to filter the features.   
setMethod(
    f = "[", 
    signature = "MultiDataSet", 
    definition = function(x, i, j, drop = FALSE, k) {
        if (missing(i) && missing(j) && missing(k)) {
            stop("Specify samples (i), tables (j) or range (k) to subset.")
        }
      
        ## Filtering tables if j is present
        if (!missing(j)) {
            if(is.numeric(j)) {
                j <- names(x)[j]
            }
            if (length(j) == 1 & drop) {
                x <- x@return_method[[j]](x@assayData[[j]], x@phenoData[[j]], x@featureData[[j]])
                validObject(x)
            } else {
                if(sum(j %in% names(x)) != length(j)) {
                    stop("Invalid tables' selection. Given table not present.")
                } else {
                    x@return_method <- x@return_method[j]
                    x@phenoData <- x@phenoData[j]
                    x@featureData <- x@featureData[j]
                    x@assayData <- x@assayData[j]
                }
            }
        }
        ## /
    
        ## If i not mising get common samples alog the full sets
        if(!missing(i)) {
            samp <- sampleNames(x)
            if(length(samp) <= 0)
                stop("Invalid samples' selection. It is an empty object.")
            
            samp$selection <- i
            samp <- Reduce(intersect, samp)
            
            if(length(samp) <= 0)
                stop("Invalid samples' selection. There are no samples.")
            
            i <- samp
            rm(samp)
        }
        ## /
        ## Filter assay data if i or k is present
        if(!missing(i) | !missing(k)) {
            ## assayData; implemented here to avoid function call
            assyD <- list()
            phenD <- list()
            protD <- list()
            featD <- list()
            for(dtype in names(x)) {
#                 if(!missing(k)) {       # k may or may not be present
#                     nfData <- subsetByOverlaps(x@featureData[[data]], k)
#                     sNames <- names(nfData)
#                     
#                     if(!missing(i)) {   # if i and k present, intersection of samples
#                         sNames <- intersect(i, sNames)
#                     }
#                 } else {
#                     sNames <- i
#                 }
#                 
                orig <- assayData(x[[dtype]])
                storage.mode <- Biobase:::assayDataStorageMode(orig)
                assyD[[dtype]] <-
                    switch(storage.mode,
                        environment =,
                        lockedEnvironment = {
                            aData <- new.env(parent=emptyenv())
                        
                            for(nm in ls(orig)){
                                aData[[nm]] <- orig[[nm]][, i, drop = FALSE]
                            }
                            
                            if ("lockedEnvironment" == storage.mode) {
                                Biobase:::assayDataEnvLock(aData)
                            }
                            aData
                        },
                        list = {
                            lapply(orig, function(obj) obj[, i, drop = FALSE])
                        })
                
                phenD[[dtype]] <- x@phenoData[[dtype]][i, , drop = FALSE]
            }
            x@assayData <- assyD
            x@phenoData <- phenD
        }
        ## /
    
        validObject(x)
        return(x)
})