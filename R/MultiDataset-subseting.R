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
#' @param drop ...
#' @param k \code{GenomicRange} used to filter the features.   
setMethod(
    f = "[", 
    signature = "MultiDataSet", # c("MultiDataSet", "character", "ANY", "GRanges"),
    definition = function(x, i, j, k, ..., drop = FALSE) {
        if (missing(i) && missing(j) && missing(k)) {
            stop("Specify samples (i), tables (j) or range (k) to subset.")
        }
        
        ## Filtering tables if j is present
        if (!missing(j)) {
            if(is.numeric(j)) {
                j <- names(x)[j]
            }

                if(sum(j %in% names(x)) != length(j)) {
                    stop("Invalid tables' selection. Given table not present.")
                } else {
                    x@return_method <- x@return_method[j]
                    x@phenoData <- x@phenoData[j]
                    x@featureData <- x@featureData[j]
                    x@assayData <- x@assayData[j]
                    x@rowRanges <- x@rowRanges[j]
                }
            
        }
        ## /
        
        ## If i not mising get common samples alog the full sets
        if(!missing(i)) {
            ## Checking the content of object for possible errors
            samp <- sampleNames(x)
            if(length(samp) <= 0)
                stop("Invalid samples' selection. It is an empty object.")
            
            samp$selection <- i
            samp2 <- Reduce(intersect, samp)
            
            if(length(samp2) <= 0)
                stop("Invalid samples' selection. Given samples are not at any dataset.")
            rm(samp2)
            ## /
            
            ## assayData; implemented here to avoid function call
            assyD <- list()
            phenD <- list()
            for(dtype in names(x)) {
                orig <- assayData(x[[dtype]])
                storage.mode <- Biobase:::assayDataStorageMode(orig)
                ii <- intersect(i, samp[[dtype]])
                assyD[[dtype]] <-
                    switch(storage.mode,
                           environment =,
                           lockedEnvironment = {
                               aData <- new.env(parent=emptyenv())
                               
                               for(nm in ls(orig)){
                                   aData[[nm]] <- orig[[nm]][, ii, drop = FALSE]
                               }
                               
                               if ("lockedEnvironment" == storage.mode) {
                                   Biobase:::assayDataEnvLock(aData)
                               }
                               aData
                           },
                           list = {
                               lapply(orig, function(obj) obj[, ii, drop = FALSE])
                           })
                
                phenD[[dtype]] <- x@phenoData[[dtype]][ii, , drop = FALSE]
            }
            x@assayData <- assyD
            x@phenoData <- phenD
        }
        ## /
        if(!missing(k)) {       # k may or may not be present
            assyD <- list()
            featD <- list()
            rangeD <- list()
            ranges <- rowRanges(x)
            if (sum(is.na(ranges))){
                x <- x[ , names(x)[!is.na(ranges)]]
            }
            for (dtype in names(x)){
                nfData <- GenomicRanges::subsetByOverlaps(x@rowRanges[[dtype]], k)
                fNames <- names(nfData)
                orig <- assayData(x[[dtype]])
                storage.mode <- Biobase:::assayDataStorageMode(orig)
                assyD[[dtype]] <-
                    switch(storage.mode,
                           environment =,
                           lockedEnvironment = {
                               aData <- new.env(parent=emptyenv())
                               
                               for(nm in ls(orig)){
                                   aData[[nm]] <- orig[[nm]][fNames, , drop = FALSE]
                               }
                               
                               if ("lockedEnvironment" == storage.mode) {
                                   Biobase:::assayDataEnvLock(aData)
                               }
                               aData
                           },
                           list = {
                               lapply(orig, function(obj) obj[fNames, , drop = FALSE])
                           })
                featD[[dtype]] <- x@featureData[[dtype]][fNames, , drop = FALSE]
                rangeD[[dtype]] <- x@rowRanges[[dtype]][fNames]
                
            }   
            x@assayData <- assyD
            x@featureData <- featD
            x@rowRanges <- rangeD
        }    
        
        if (length(x) == 1 & drop) {
            x <- x@return_method[[1]](x@assayData[[1]], x@phenoData[[1]], x@featureData[[1]])
        }
        
        validObject(x)
        return(x)
    })