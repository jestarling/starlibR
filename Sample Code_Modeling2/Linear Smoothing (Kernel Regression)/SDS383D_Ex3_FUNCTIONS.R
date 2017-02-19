#SDS 383D - Exercise 3
#Functions
#Jennifer Starling
#Feb 15, 2017

#-------------------------------------------------------------
### Kernel-Related Functions (Data simulation, linear smoother)

sim_noisy_data = function(x,f,sig2){
	#FUNCTION: 	Simulates noisy data from some nonlinear function f(x).
	#INPUTS:	x = independent observations.
	#			f = function f(x) to simulate from
	#			sig2 = variance of the e ~ N(0,sig2) noise.
	#OUTPUTS:	x = generated x values.
	#			y = generated y = f(x) + e values.

	fx = f(x)
	e = rnorm(length(x),0,sqrt(sig2))
	return(y = fx+e)
}

linear_smoother = function(x,y,x_star,h=1,K){
	#FUNCTION: 	Linear smoothing function for kernel regression.
	#INPUTS:	x = a scalar or vector of regression covariates.
	#			x_star = scalar or vector of new x values for prediction.
	#			h = a positive bandwidth.
	#			K = a kernel function.  Default is set to Gaussian kernel.
	#OUTPUT:	yhat = a scalar or vector of smoothed x values.
	
	yhat=0	#Initialize yhat.
	
	for (i in 1:length(x_star)){
		w = (1/h) * K((x-x_star[i])/h) #Calculates weights.
		w = w / sum(w)					#Normalize weights.
		yhat[i] = crossprod(w,y)
	}
	return(yhat)	
}	

#-------------------------------------------------------------
### Kernel Functions

K_uniform = function(x){
	#FUNCTION: 	Uniform kernel.
	#INPUTS:	x = a scalar or vector of values.
	#OUTPUT:	k = a scalar or vector of smoothed x values.
	
	k = .5 * ifelse(abs(x)<=1,rep(1,length(x)),rep(0,length(x)))
	return(k)
}

K_gaussian = function(x){
	#FUNCTION: Gaussian kernel.
	#INPUTS:	x = a scalar or vector of values.
	#OUTPUT:	k = a scalar or vector of smoothed x values.
	
	k = (1/sqrt(2*pi)) * exp(-x^2/2)
	return(k)
}

#-------------------------------------------------------------
### Cross-Validation Functions (tuning bandwidth)

tune_h = function(test,train,K,h){
	#FUNCTION: 	Function to tune bandwidth h for 
	#				specified test/train data sets and 
	#				specified kernel K.
	#INPUTS:	test = a test data set. Must have two cols, x and y.
	#			train = a training data set. Must have two cols, x and y.
	#			K = the kernel function
	#			h = a scalar or vector of bandwidths to try.
	#OUTPUTS:	pred_err_test = prediction error for testing data for each h.
	#			fhat_test = predicted values for testing data's x vals.
	
	#Extract training and test x and y vectors.
	x = train[,1]
	y = train[,2]
	x_star = test[,1]
	y_star = test[,2]
	
	#Calculate predicted points for test data set, and prediction error.
	fhat_star = linear_smoother(x,y,x_star,h,K)
	
	#Calculate predicted points for test data set, and prediction error.
	pred_err_test = sum(fhat_star-y_star)^2 / (length(x))
	
	#Return function outputs:
	return(list(fhat_star=fhat_star, pred_err_test = pred_err_test))
}