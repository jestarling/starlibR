#----------------------------------------------------------------
# Roxy package build comments:
#' vecAngle(x,y)
#'
#' This function calculates the angle in degrees between two vectors.
#'
#' @param x A vector.  Defaults to NULL.
#' @param y A vector.  Defaults to NULL.
#'
#' @keywords
#'
#' @return A scalar, the angle between the two vectors.
#'
#' @examples
#' x <- c(1,3,5,7,9)
#' y <- c(2,4,6,8,10)
#' vecAngle(x,y)
#'
#' @author Jennifer Starling
#'
#' @export
#----------------------------------------------------------------

vecAngle <- function(x,y) acos(sum(x %*% y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)) ))
