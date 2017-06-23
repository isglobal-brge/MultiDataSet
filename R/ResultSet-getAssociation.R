#' @describeIn ResultSet Getter to obtain the raw \code{data.frame} from
#' association and integration analysis.
#' @param sort Indicates of the result should be ordered by P-Value
setMethod(
    f = "getAssociation",
    signature = "ResultSet",
    definition = function(object, rid, coef = 2, contrasts = NULL, 
                          fNames = c("chromosome", "start", "end", "genesymbol"), ...) {
        
        res <- object@results[[rid]]$result
        
        if(class(res) == "MArrayLM") {
            if (!is.null(contrasts)){
                res <- limma::contrasts.fit(res, contrasts)
            }
            res <- limma::eBayes(res)
            res <- limma::topTable(res, coef = coef, number = Inf, ...)
            
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
