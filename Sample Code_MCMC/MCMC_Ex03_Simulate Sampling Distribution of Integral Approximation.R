#MCMC
#Exercise 3
#Jan 28, 2017

# 1. Sampling Distribution of IHAT_N integral approximation.

#   Suppose I = integral of gx*fx dx from 0 to 1.
#	where fx = uniform density, and gx = x^2.
#	Var Ihat_N = (1/N) * ( integral of gx^x - (integral of gx)^2) = sigma^2/n.
#	Verify by simulation that sqrt(N)(Ihat-I) ~ N(0,sigma^2).

N = 1000			#Set sample size.
fx = runif(N,0,1)	#Random sample of uniform(0,1) realizations.
gx = fx^2

I = 1/3				#True value of integral of x^2 from 0 to 1.

#Simulate a single sigma2_hat value.
Ihat_N = (1/N) * sum(gx)
Var_Ihat_N = (1/N) * ((1/N)*sum(gx^2) - ((1/N)*sum(gx))^2)
sigma2_hat = Var_Ihat_N * N

#Simulate B=10000 Ihat_N values.  Estimate variance, calculate sigma2_hat.
B = 10000
Ihat_vec = rep(0,B)

for (b in 1:B){
	fx = runif(N,0,1)	#Random sample of uniform(0,1) realizations.
	gx = fx^2
	Ihat_vec[b] = (1/N) * sum(gx)	
}

var_boot = var(Ihat_vec)	#Estimate bootstrapped var_Ihat_N.
sigma2_boot = var_boot*N	#Estimate bootstrapped sigma2_hat.

#Standardize the Ihat_N results and check out the sampling distribution.
Standardized_Ihat_N = sqrt(N) * (Ihat_vec - I)
mean(Standardized_Ihat_N)
var(Standardized_Ihat_N)

#Var we expected to see based on Var_Ihat_N formula provided:
sigma2_hat

hist(Standardized_Ihat_N,main='Sampling distribution of standardized Ihat_N estimates')

