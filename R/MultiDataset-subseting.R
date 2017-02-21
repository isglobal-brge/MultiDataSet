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
#' @param drop If TRUE, sets with no samples or features will be discarded
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
            
            samp2 <- Reduce(union, samp)
            samp2 <- intersect(samp2, i)
            
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
                phen <- x@phenoData[[dtype]]
                sel <- order(match(phen$id, i), na.last = NA)
                assyD[[dtype]] <-
                    switch(storage.mode,
                           environment =,
                           lockedEnvironment = {
                               aData <- new.env(parent=emptyenv())
                               
                               for(nm in ls(orig)){
                                   aData[[nm]] <- orig[[nm]][, sel, drop = FALSE]
                               }
                               
                               if ("lockedEnvironment" == storage.mode) {
                                   Biobase:::assayDataEnvLock(aData)
                               }
                               aData
                           },
                           list = {
                               lapply(orig, function(obj) obj[, sel, drop = FALSE])
                           })
                
                phenD[[dtype]] <- x@phenoData[[dtype]][sel, , drop = FALSE]
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
                nfData <- IRanges::subsetByOverlaps(x@rowRanges[[dtype]], k)
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

#' @describeIn MultiDataSet Filter a subset using feature ids or phenotypes
#' @aliases MultiDataSet-methods
#' @param feat Logical expression indicating features to keep
#' @param phe Logical expression indicating the phenotype of the samples to keep
#' @param keep If FALSE, sets where the expression cannot be evaluated will be discarded.
setMethod(
    f = "subset",
    signature = "MultiDataSet",
    definition = 
        function(x, feat, phe, warnings = TRUE, keep = TRUE) {
            
            if (missing(feat) && missing(phe)) {
                stop("Specify genes id (feat) or phenotypes (phe) to subset.")
            }
            
            if (!missing(feat)){
                assyD <- list()
                featD <- list()
                rangeD <- list()
                ranges <- rowRanges(x)
                
                # Catch the expression
                e <- substitute(feat)
                
                # Get column names that will be used to filter 
                featcols <- all.vars(e)
                noFilteredSets <- character()
                for(dtype in names(x)) {
                    feats <- x@featureData[[dtype]]
                    if (!all(featcols %in% colnames(feats))){
                        noFilteredSets <- c(noFilteredSets, dtype)
                        if (keep) {
                            featD[[dtype]] <- x@featureData[[dtype]]
                            assyD[[dtype]] <- x@assayData[[dtype]]
                            rangeD[[dtype]] <- x@rowRanges[[dtype]]
                        }
                    } else {
                        sel <- eval(e, pData(feats), parent.frame())
                        if (!is.logical(sel))
                            stop("'feat' must be a logical expression")
                        sel <- sel & !is.na(sel)
                        orig <- assayData(x[[dtype]])
                        storage.mode <- Biobase:::assayDataStorageMode(orig)
                        assyD[[dtype]] <-
                            switch(storage.mode,
                                   environment =,
                                   lockedEnvironment = {
                                       aData <- new.env(parent=emptyenv())
                                       
                                       for(nm in ls(orig)){
                                           aData[[nm]] <- orig[[nm]][sel, , drop = FALSE]
                                       }
                                       
                                       if ("lockedEnvironment" == storage.mode) {
                                           Biobase:::assayDataEnvLock(aData)
                                       }
                                       aData
                                   },
                                   list = {
                                       lapply(orig, function(obj) obj[sel, , drop = FALSE])
                                   })
                        featD[[dtype]] <- x@featureData[[dtype]][sel, , drop = FALSE]
                        rangeD[[dtype]] <- x@rowRanges[[dtype]][sel]
                    }
                    
                }
                x@assayData <- assyD
                x@featureData <- featD
                x@rowRanges <- rangeD
                
                if (length(noFilteredSets) && noFilteredSets == names(x)){
                    stop("feat expression could not be applied to any of the sets.")
                }
                
                if (warnings & length(noFilteredSets)) {
                    warn <- "The following sets could not be filtered by feature id"
                    if (keep){
                        warning(paste0(warn, ": ", paste(noFilteredSets, collapse = ", ")))
                    } else{
                        warning(paste(warn, "and have been discarded:", paste(noFilteredSets, collapse = ", ")))
                    }
                }
            }
            
            if (!missing(phe)){
                assyD <- list()
                phenD <- list()
                
                # Catch the expression
                e <- substitute(phe)
                
                # Get column names that will be used to filter 
                phenocols <- all.vars(e)
                noFilteredSets <- character()
                for(dtype in names(x)) {
                    phen <- x@phenoData[[dtype]]
                    if (!all(phenocols %in% colnames(phen))){
                        noFilteredSets <- c(noFilteredSets, dtype)
                        if (keep) {
                            phenD[[dtype]] <- x@phenoData[[dtype]]
                            assyD[[dtype]] <- x@assayData[[dtype]]
                        }
                    } else {
                        sel <- eval(e, pData(phen), parent.frame())
                        if (!is.logical(sel))
                            stop("'phe' must be a logical expression")
                        sel <- sel & !is.na(sel)
                        orig <- assayData(x[[dtype]])
                        storage.mode <- Biobase:::assayDataStorageMode(orig)
                        assyD[[dtype]] <-
                            switch(storage.mode,
                                   environment =,
                                   lockedEnvironment = {
                                       aData <- new.env(parent=emptyenv())
                                       
                                       for(nm in ls(orig)){
                                           aData[[nm]] <- orig[[nm]][, sel, drop = FALSE]
                                       }
                                       
                                       if ("lockedEnvironment" == storage.mode) {
                                           Biobase:::assayDataEnvLock(aData)
                                       }
                                       aData
                                   },
                                   list = {
                                       lapply(orig, function(obj) obj[, sel, drop = FALSE])
                                   })
                        
                        phenD[[dtype]] <- x@phenoData[[dtype]][sel, , drop = FALSE]
                    }
                    
                }
                x@assayData <- assyD
                x@phenoData <- phenD
                
                if (length(noFilteredSets) && noFilteredSets == names(x)){
                    stop("phe expression could not be applied to any of the sets.")
                }
                
                if (warnings & length(noFilteredSets)) {
                    warn <- "The following sets could not be filtered by phenotype"
                    if (keep){
                        warning(paste(warn, ":", paste(noFilteredSets, collapse = ", ")))
                    } else{
                        warning(paste(warn, "and have been discarded:", paste(noFilteredSets, collapse = ", ")))
                    }
                }
            }
            
            
            validObject(x)
            return(x)
        })
