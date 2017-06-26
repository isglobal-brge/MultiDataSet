#' @export 
#' @name create_resultset
#' @rdname ResultSet-class
#' @aliases ResultSet-methods
#' 
#' @param fOrigin Chracter with the function used to run the analysis.
#' @param lResults List with the results
#' @param fData List with the feature data.
#' @param lOptions List with additional options
#' @examples
#' create_resultset("hello", list(), list(), list())
create_resultset <- function(fOrigin, lResults, fData, lOptions = list()){
    new("ResultSet",
        fun_origin = fOrigin,
        results = lResults,
        fData = fData,
        options = lOptions
    )
}
