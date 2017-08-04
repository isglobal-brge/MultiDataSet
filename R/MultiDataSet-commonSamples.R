#' @describeIn MultiDataSet Get a MultiDataSet only with the samples present in all the tables
#' @param unify.names Logical indicating if sample names of the sets should be unified.
#' @aliases MultiDataSet-methods
setMethod(
  f = "commonSamples",
  signature = "MultiDataSet",
  definition = function(object, unify.names = FALSE) {
    samples <- commonIds(object)
    if (length(samples) == 0){
      stop("There are no samples present in all the datasets")
    }
    object <- object[samples, ]
    if (unify.names){
        assyD <- list()
        phenD <- list()
        
        dups <- vapply(sampleNames(object), function(x) sum(duplicated(x)), numeric(1))
        if (any(dups != 0)){
            stop("Sample names cannot be unified. ", paste(names(dups)[dups != 0], collapse = ", "), " contain(s) duplicated sample names.")
        }
        
        for(dtype in names(object)) {
            orig <- assayData(object)[[dtype]]
            storage.mode <- Biobase:::assayDataStorageMode(orig)
            ids <- object@phenoData[[dtype]]$main$id
            
            assyD[[dtype]] <-
                switch(storage.mode,
                       environment =,
                       lockedEnvironment = {
                           aData <- new.env(parent=emptyenv())
                           
                           for(nm in ls(orig)){
                               aData[[nm]] <- orig[[nm]]
                               colnames(aData[[nm]]) <- ids
                               
                           }
                           
                           if ("lockedEnvironment" == storage.mode) {
                               Biobase:::assayDataEnvLock(aData)
                           }
                           aData
                       },
                       list = {
                           lapply(orig, function(obj) {
                               aData[[nm]] <- orig[[nm]]
                               colnames(aData[[nm]]) <- ids
                               aData[[nm]]
                           })
                       })
            
            phenD[[dtype]] <- object@phenoData[[dtype]]
            for (tab in names(phenD[[dtype]])){
                rownames(phenD[[dtype]][[tab]]) <- ids
            }
        }
        object@assayData <- assyD
        object@phenoData <- phenD
    }
        return(object)
  }
)
