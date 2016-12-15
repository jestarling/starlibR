linear.regression <- function(X,y){

	n = length(y)
	X = as.matrix(cbind(rep(1,n),X))
	beta_hat <- solve(t(X) %*% X) %*% t(X) %*% y

	#Calculate y_hat and RSS.
	y_hat <- X %*% beta_hat
	RSS <- t(y - X %*% beta_hat) %*% (y - X %*% beta_hat)

	n <- nrow(X)	#Number of obs.
	p <- 4			#Number of predictors.
	sigma2_hat <- RSS / (n-p-1)		#Estimate of sigma^2.
	sigma_hat <- sqrt(sigma2_hat)	#Estimate of sigma.

	#Estimated standard error and CI for each coefficient:
	beta_hat_se <- sqrt( sigma2_hat * diag(solve(t(X) %*% X)) ) #Standard errors.

	#95% CI for beta_hat's.
	lb <- beta_hat -  1.96 * beta_hat_se
	ub <- beta_hat + 1.96 * beta_hat_se
	sig <- ifelse( (lb < 0) & (ub > 0), 0,1)

	#Display CI, including star for betas not equal to zero based on the CI.
	ci = data.frame(lower.bound=lb,upper.bound=ub,sig.coefs=sig)
	
	return(list(beta_hat = beta_hat, 
				RSS = RSS, 
				sigma_hat = sigma_hat, 
				beta_hat_se = beta_hat_se, 
				ci = ci,
				n=n,
				p = ncol(X)-1))
} #End myLinearReg function.
