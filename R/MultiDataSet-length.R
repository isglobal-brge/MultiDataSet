#' @describeIn MultiDataSet Returns the number of sets into the object.
#' @aliases length
setMethod(
    f = "length",
    signature = c("MultiDataSet"),
    definition = function(x) {
        return(length(names(x)))
    } 
)