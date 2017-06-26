#' @describeIn ResultSet Returns the names of the variables of the models used in
#' a \code{ResultSet}.
setMethod(
    f = "varLabels",
    signature="ResultSet",
    definition = function(object) {
        return(lapply(object@results, function(x) {
            if(class(x$result) == "MArrayLM") {
                colnames(x$result$design)
            } else {
                NULL
            }
        }))
    })