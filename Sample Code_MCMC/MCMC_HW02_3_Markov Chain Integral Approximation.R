#MCMC Homework 2
#Problem 3
#Feb 9, 2017
#Jennifer Starling

#Approximate the integral 
#I = \int_{-\infty}^{\infty} (1 + x^4) * exp(-(1/2)*x^2) 
#using a Monte Carlo chain sample (Xn), 
#given by X_{n+1} ~ N(rho*X_n, 1-rho^2).
#First Xn value is X1 ~ N(0,1).

#------------------------------------------------
mcmc_integral_approx = function(n,rho){
	#FUNCTION: MCMC Integral Approximation Function
	#INPUTS:	n = length of markov chain
	#			rho = correlation for X_{n+1} ~ N(rho*X_n, 1-rho^2)	
	#OUTPUTS:	Ihat_N = integral approximation
	#			x = markov chain
	
	#Initialize chain.
	x = rnorm(1,0,1)
	
	#Iterate through markov chain.
	for (i in 2:n){
		#Previous x value, x_n.
		x0 = x[i-1]	
		
		#Set new x value, x_{n+1}.
		x[i] = rnorm(1,rho*x0,sqrt(1-rho^2))
	}
	
	#Calculate gx values and Ihat_N integral approximation.
	gx = sqrt(2*pi) / (1 + x^4)
	Ihat_N = (1/n) * sum(gx)
	
	#Return function output.
	return(list(x=x,Ihat_N=Ihat_N))
}
#------------------------------------------------

#Run for a one-off case to obtain a single Ihat_N approximation.
Ihat_N = mcmc_integral_approx(n=1000,rho=0)
Ihat_N$Ihat_N  #Shoudl be approx 1.696

#Try various values of rho.  Run B=1000 estimates of Ihat_N at each rho,
#and estimate var(Ihat_N) to see which rho gives best variance.

B=10000								#Number of Ihat_N estimates for ea rho.
rho = seq(-1,1,by=.25)				#Vector of rho values to try.
Ihat_N_seq = rep(0,B)				#Temp vector to hold all B Ihat_N.
var_Ihat_N = rep(0,length(rho))		#Vector to hold Ihat_N variances.

#Iterate through rho values.
for (r in 1:length(rho)){
	#Iterate through b values for each rho.
	for (b in 1:B){
		#Calculate Ihat_N_seq estimates.
		Ihat_N_seq[b] =  mcmc_integral_approx(n=1000,rho=rho[r])$Ihat_N
	}
	
	#Ihat_N variance estimate for each rho.
	var_Ihat_N[r] = var(Ihat_N_seq)
}			

cbind(rho,var_Ihat_N)	#Display variance results for each rho.		

#Plot rho vs variance.
pdf(file='/Users/jennstarling/UTAustin/2017S_MCMC/Homework/Homework 02/Figures/Rho.pdf')
plot(rho,var_Ihat_N,type='l',col='blue',lwd=2,xlab='rho',ylab='Var(Ihat_N)',main='MC Integral Approximation Variance',
	ylim=c(0,.002))
dev.off()
