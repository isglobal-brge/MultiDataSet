#' @describeIn ResultSet Getter to obtain the raw \code{data.frame} from
#' association and integration analysis.
#' @param contrast Numeric matrix with the contrasts used to perform the analyses
#' @param fNames Character vector with the names of the fData columns that will 
#' be added to the results data.frame. 
#' @param ... Further arguments passed to \link[limma]{topTable}
setMethod(
    f = "getAssociation",
    signature = "ResultSet",
    definition = function(object, rid = 1, coef = 2, contrast = NULL, 
                          fNames = NULL, ...) {
        
        res <- object@results[[rid]]$result
        
        if(is(res, "MArrayLM")) {
            fit <- res; rm(res)
            if (!is.null(contrast)){
                fit <- limma::contrasts.fit(fit, contrast)
            }
            fit <- limma::eBayes(fit)
            res <- limma::topTable(fit, coef = coef, number = Inf, confint = TRUE, ...)
            res$SE <- sqrt(fit$s2.post) * fit$stdev.unscaled
            
            ## Add fData to results
            if (!is.null(fNames)){
                fData <- object@fData[[1]]
                
                if (!all(fNames %in% colnames(fData))){
                    stop("All fNames must be present in ResultSet fData.")
                }
                
                res <- cbind(res, fData[rownames(res), fNames])  
            }
        }
        
        return(res)
})
