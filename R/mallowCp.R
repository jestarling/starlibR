#---FUNCTION-------------------------
#Custom function to calculate mallow's cp.  (Using normalized version.)
mallowCp <- function(X,y){
	model = myLinearReg(X,y)
	n = model$n		#Number of obs.
	p = model$p		#Number of predictors.
	
	RSS = model$RSS	#RSS for fitted model.
	sig2 = model$sigma_hat^2	#Sigma^2_hat for fitted model.
	cp = (1/n) * (RSS + 2*p*sig2)	#Mallow's cp for fitted model.
	
	return(cp)
}