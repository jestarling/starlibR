#Dirichlet Log-likelihood function:
dirichlet.loglik <- function(data,alpha){	
	N <- nrow(data)						#Number of (j) observations.
	logxi <- colMeans(log(data))		#logxi = (1/N)*sum(log(x_ij)
	logl <- N* (lgamma(sum(alpha)) - sum(lgamma(alpha)) + sum((alpha-1)*logxi))

	return(logl)	
}