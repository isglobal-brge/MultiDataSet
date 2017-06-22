#' @describeIn ResultSet Returns the names of the analyses
#' stored in the \code{ResultSet}.
setMethod(
    f = "rid",
    signature = "ResultSet",
    definition = function(object) {
        return(names(object@results))
    }
)
