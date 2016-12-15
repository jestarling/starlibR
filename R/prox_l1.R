#' Proximal operator of the scaled L1 norm.

#'

#' Computes the proximal operator of the L1 norm: \deqn{h(x) = \lambda ||x||_1

#' ,} where \eqn{\lambda} is a scaling factor.

#'

#' @param x The input vector

#' @param t The step size

#' @param opts List of parameters, which must include:

#'   \itemize{\code{opts$lambda}, the scaling factor of the L1 norm.}

#'

#' @return The proximal of \eqn{h} at {x} with step size \eqn{t}, given by

#'   \deqn{prox_h(x,t) = argmin_u [ t h(u) + 1/2 || x - u ||^2 ]}.

#'

#' @export

#' @examples prox.l1(c(1,3,-2), 1.5, list(lambda=1))



soft.thresholding <- function(x, lambda){

  # Computes the soft thresholding estimator

  # ----------------------------------------

  # Args: 

  #   - x: vector of the observations

  #   - lambda: penalization parameter (threshold)

  # Returns: 

  #   - theta: the soft thresholding estimator

  # ------------------------------------------

  theta <- sign(x) * pmax(rep(0, length(x)), abs(x) - lambda)

  return (theta)

}



