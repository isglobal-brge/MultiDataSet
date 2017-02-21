#' @describeIn MultiDataSet Apply iClusterPlus clustering method to a MultiDataSet object
#' @aliases MultiDataSet-methods
#' @param commonSamples Logical to indicate if common samples are selected
setMethod(
    f = "w_iclusterplus",
    signature = "MultiDataSet",
    definition = function(object, commonSamples = TRUE, ...) {
        if(!"iClusterPlus" %in% utils::installed.packages()[,"Package"]) {
            stop("R package 'iClusterPlus' not found. Please, install it al run again this function.")
        }

        if(length(object) > 4) {
            stop("'iClusterPlus' only allows four datasets.")
        }
        
        # Obtain the common samples between datsets
        if(commonSamples) {
            object <- commonSamples(object)
        }
        
        # Put the names iClusterPlus requires on the datasets
        datasets <- lapply(as.list(object), t)
        names(datasets) <- paste("dt", 1:length(datasets), sep="")
        
        # Generate the "type" argument for iClusterPlus
        datasets[["type"]] <- sapply(names(object), function(nm) {
            ifelse(startsWith(nm, "expression"), "gaussian",
                   ifelse(startsWith(nm, "methylation"), "gaussian", "multinomial"))
        })
        names(datasets[["type"]]) <- paste("dt", 1:(length(datasets) -1), sep="")
        
        # Call iClusterPlus with the generated arguments and the user's arguments
        do.call(iClusterPlus::iClusterPlus, c(datasets, list(...)))
    }
)