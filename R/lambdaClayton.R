#' Lambda Calculation for a vector of P-Values
#'
#' Implementation of Clayton's lambda score for a vector of P-Values
#'
#' @author Juan R. Gonzalez
#' @param x Vector of P-Value
#' @param trim (default \code{0.5})
#' @return A lambda value, indicating the inflation/deflation of the
#' analysis.
#' @examples 
#' lambdaClayton(runif(30))
#' @export
lambdaClayton <- function(x, trim = 0.5) {
    xx <- qnorm(1-x)^2
    N <- length(xx)
    obsvd <- sort(xx, na.last=NA)
    expctd <- qchisq(p = (1:N)/(N + 1), 1)
    Nu <- floor(trim * N)
    lambda <- mean(obsvd[1:Nu])/mean(expctd[1:Nu])
    lambda
}
