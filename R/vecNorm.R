#----------------------------------------------------------------
# Roxy package build comments:
#' vecNorm(x)
#'
#' This function calculates the norm (magnitude, Euclidian length) of a vector x.
#'
#' @param x A vector.  Defaults to NULL.
#'
#' @keywords
#'
#' @return A scalar, the norm of vector x.  
#'
#' @examples
#' x <- c(1,3,5,7,9)
#' vecNorm(x)
#'
#' @author Jennifer Starling
#'
#' @export
#----------------------------------------------------------------


vecNorm <- function(x) sqrt(sum(x^2))