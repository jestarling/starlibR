# Compute the vector cross product between x and y, and return the components
# indexed by i.  X and Y must be vectors of length 3.

#----------------------------------------------------------------
# Roxy package build comments:
#' CrossProduct3D(a,b,i)
#'
#' This function calculates the cross product of two vectors, each of 
#' length 3.  Specifying i = 1, 2 or 3 returns a single element of the cross-product
#' vector, based on the i index.  Leaving i blank returns the entire cross-product vector.
#'
#' @param a A vector of length 3.
#' @param b A vector of length 3.
#' @param i An index value of 1, 2, or 3.
#'
#' @keywords
#'
#' @return A vector (length 3) containing the cross-product of the two input vectors.
#'
#' @examples
#' a <- c(2,3,4)
#' b <- c(5,6,7)
#' x <- CrossProduct3D(a,b)
#'
#' @author Jennifer Starling
#'
#' @export
#----------------------------------------------------------------

CrossProduct3D <- function(x, y, i=1:3) {
  # Project inputs into 3D, since the cross product only makes sense in 3D.
  To3D <- function(x) head(c(x, rep(0, 3)), 3)
  x <- To3D(x)
  y <- To3D(y)

  # Indices should be treated cyclically (i.e., index 4 is "really" index 1, and
  # so on).  Index3D() lets us do that using R's convention of 1-based (rather
  # than 0-based) arrays.
  Index3D <- function(i) (i - 1) %% 3 + 1

  # The i'th component of the cross product is:
  # (x[i + 1] * y[i + 2]) - (x[i + 2] * y[i + 1])
  # as long as we treat the indices cyclically.
  return (x[Index3D(i + 1)] * y[Index3D(i + 2)] -
          x[Index3D(i + 2)] * y[Index3D(i + 1)])
}
