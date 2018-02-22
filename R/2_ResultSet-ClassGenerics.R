#' Class ResultSet
#'
#' Class \code{ResultSet} used to encapsulate results from \code{MEAL} and 
#' \code{omicrexposome}.
#'
#' @name ResultSet
#' @aliases ResultSet-class
#' @rdname ResultSet-class
#' @exportClass ResultSet
#' @slot fun_origin Character containing the function that creates the object.
#' @slot results List containing the results of the association/integration.
#' @slot fData List containing the feature-data of the original objects.
#' @slot options list of options used to create the \code{ResultSet}.
#' @return An object of class \code{ResultSet}
setClass(
    Class = "ResultSet",
    representation = representation(
        fun_origin = "character",
        results = "list",
        fData = "list",
        options = "list"
    )
)

## ------------------------------------------------------------------------- ##

#' Method to extrat feature result from a ResultSet
#'
#' Homologous methods from \code{limma}, \code{getAssociation} resturns a
#' \code{data.frame} with the \code{logFC} and  \code{PValue} per
#' featrue for the selcted \code{coef} and for given result (\code{rid}).
#'
#' @name getAssociation
#' @rdname getAssociation-methods
#' @aliases getAssociation
#' @param object A \code{\link{ResultSet}} object.
#' @param rid The name or index of the result to be extracted.
#' @param coef (default \code{2}) Index of the coefficient to be extracted.
#' @param contrast (default \code{1}) When \code{code} corresponds to a
#' multicategorical variable, contasr selects the comparison.
#' @param fNames (default \code{c("chromosome", "start", "end", "genesymbol")})
#' Corresponds to the columns selected from \code{fData} that will be
#' incorporated into the resulting \code{data.frame}.
#' @param ... Further arguments passed to \link[limma]{topTable}
#' @return A \code{data.frame} with the result of the association study,
#' including P-Value and Fold Change.
#' @examples
#' data(rset)
#' getAssociation(rset, rid=1, fNames = c("chromosome", "position"))
#' @export getAssociation
setGeneric("getAssociation", function(object, rid = 1, coef = 2, contrast = NULL, 
       fNames = NULL, ...)
    standardGeneric("getAssociation")
)

#' Method to get the options sued to create the ResultSet
#'
#' Method that returns a list with the options used to create the
#' \code{ResultSet}.
#'
#' @name opt
#' @rdname opt-methods
#' @aliases opt
#' @param object A \code{\link{ResultSet}} object.
#' @return A list with the options used to create the \code{ResultSet}.
#' @examples
#' data(rset)
#' opt(rset)
#' @export opt
setGeneric("opt", function(object)
    standardGeneric("opt")
)
