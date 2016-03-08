#----------------------------------------------------------------
# Roxy package build comments:
#' myBoxCox(x) A Box-Cox transformation function.
#'
#' This function calculates the theta value for the y^theta box-cox transformation
#' to obtain normally distributed transformed data.
#'
#' @param x A vector of data to be transformed to normality.
#'
#' @keywords transformation box-cox normal
#'
#' @return x A list containing two objects.
#' @return x$theta The theta value to use for the box-cox transformation.
#' @return x$directions A text string reminder for how to perform the transformation using theta.
#'
#' @examples
#' x_lognorm <- rlnorm(100, meanlog = 0, sdlog = 1)
#' x_norm <- myBoxCox(x)
#'
#' @author Jennifer Starling
#'
#' @export

#----------------------------------------------------------------

myBoxCox <- function(x){

n = length(x)
y = abs(x)
yt0 = log(x)
s = sum(yt0)
varyt0 = var(yt0)
Lt0 = -1*s - .5*n*(log(2*pi*varyt0)-1)
th = 0
Lt = 0
t = -3.01
i = 0

while(t < 3){
	t = t+.001
	i = i+1
	th[i] = t
	yt = (x^t -1)/t
	varyt = var(yt,na.rm=TRUE)
	Lt[i] = (t-1)*s - .5*n*(log(2*pi*varyt)-1)
	if(abs(th[i])<1.0e-10)Lt[i]<-Lt0
	if(abs(th[i])<1.0e-10)th[i]<-0
}

# The following outputs the values of the likelihood and theta and yields
# the value of theta where likelihood is a maximum.

out = cbind(th,Lt)
Ltmax= max(Lt)
imax= which(Lt==max(Lt))
theta= th[imax]

x <- list()
x$theta <- theta
x$directions <- "Tranformation: Use data^theta to obtain normally transformed data."
return(x)	#Return theta value for power transformation.

#Take y^thmax to obtain normally transformed data.
}