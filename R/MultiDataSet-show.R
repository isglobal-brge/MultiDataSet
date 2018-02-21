setMethod(
    f = "show",
    signature="MultiDataSet",
    definition = function(object) {
        cat("Object of class 'MultiDataSet'\n")
        cat(" . assayData:", length(object@assayData), "elements\n")
        for(name in names(object)) {
            cat("    . ", name, ": ", nrow(object@featureData[[name]]$main), " features, ",
                nrow(object@phenoData[[name]]$main), " samples \n", sep = "")
        }
        
        cat( " . featureData:\n")
        for(name in names(object)) {
            cat("    . ", name, ": ", nrow(object@featureData[[name]]$main), " rows, ",
                length(varLabels(object@featureData[[name]]$main)), " cols ", sep = "")
            
            nms <- varLabels(object@featureData[[name]]$main)
            cat(ifelse(length(nms) == 0, "", paste0("(", .wrpvec(nms), ")")), "\n", sep = "")
        }
        cat( " . rowRanges:\n")
        for(name in names(object)) {
            cat("    . ", name, ": ", 
                ifelse(is(object@rowRanges[[name]], "GenomicRanges"), "YES\n", "NO\n"), sep="")
        }
        
        
        cat( " . phenoData:\n")
        for(name in names(object)) {
            cat("    . ", name, ": ", nrow(object@phenoData[[name]]$main), " samples, ",
                ncol(object@phenoData[[name]]$main), " cols", sep = "")
            
            nms <- colnames(object@phenoData[[name]]$main)
            cat(" (", ifelse(is.null(nms), "", .wrpvec(nms)), ")\n", sep="")
        }
    }
    
)

.wrpvec <- function(vector_names) {
    if (length(vector_names) < 2) {
        return(vector_names)
    } else if (length(vector_names) == 2) {
        return(paste0(vector_names[1], ", ", vector_names[2]))
    } else {
        nn <- length(vector_names) - 1
        return(paste0(vector_names[1], ", ..., ", vector_names[nn]))
    }
}
