#' @export 
#' @name create_resultset
#' @rdname ResultSet-class
#' @aliases ResultSet-methods
#' 
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
