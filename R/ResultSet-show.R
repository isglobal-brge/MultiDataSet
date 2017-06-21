setMethod(
    f = "show",
    signature="ResultSet",
    definition = function(object) {
        cat("Object of class 'ResultSet'\n", sep="")
        cat(" . created with:", object@fun_origin, "\n")
        if(object@fun_origin == "association") {
            cat("    . via:", paste0(object@options$names, collapse=" and "), "\n")
        }
        cat(" . sva: ", ifelse(object@options$sva, "yes", "no"), "\n")
        if(object@fun_origin == "crossomics") {
            cat("    . method:", object@options$method, " (", object@options$package, ")\n")
        }
        cat(" . #results:", length(object@results), "( error:", sum(!is.na(sapply(object@results, "[[", "error"))), ")\n")
        if (length(object@fData) == 1){
            nr <- nrow(object@fData[[1]])
            nc <- ncol(object@fData[[1]])
            cat(" . featureData: ", nr, " probes x ", nc, " variables\n", sep="")
        } else {
            cat(" . featureData: ", length(object@fData), "\n")
            for(nm in names(object@fData)) {
                nr <- nrow(object@fData[[nm]])
                nc <- ncol(object@fData[[nm]])
                cat("    . ", nm, ": ", nr, "x", nc, "\n", sep="")
            }
        }
    }
)
