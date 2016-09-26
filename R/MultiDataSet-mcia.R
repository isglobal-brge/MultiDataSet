setMethod(
    f = "w_mcia",
    signature = "MultiDataSet",
    definition = function(object, ...) {
        if(!"omicade4" %in% installed.packages()[,"Package"]) {
            stop("R package 'omicade4' not found. Please, install it al run again this function.")
        }
        require(omicade4)
        
        object <- commonSamples(object)
        omicade4::mcia(as_list(object), ...)
    }
)