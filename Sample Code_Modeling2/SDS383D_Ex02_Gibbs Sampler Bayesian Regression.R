#Stats Modeling 2
#Exercise 02 - Gibbs Sampler for Bayesian Linear Model

#Housekeeping
rm(list=ls())

#Full hierarchical model:
#	(y | B,w,L) ~ N(XB,(wL)^{-1})
#	L = diag(l1,...,ln) #Lambda_i values.
#	li ~ iid Gamma(h/2,h/2)
#	(B|w) ~ N(m,(wK)^{-1})
#	w ~ Gamma(d/2,eta/2)

#Conditional posteriors:
#	p(B|y,w,L) ~ N(m_star,(wK_star)^{-1})
#	p(w|y,L) ~ Gamma(d_star/2,eta_star/2)
#	p(li | y,B,w) ~ Gamma( (h+1)/2, (w(y_i - x_i'B)^2 + h)/2 )

#Posterior parameters:
#	d_star = d + n
#	K_star = (X^T %*% L %*% X + K)
#	m_star = (K_star)^{-1} (X^T %*% L %*% y + K %*% m) 
		   = (X^T %*% L %*% X + K)^{-1} %*% (X^T %*% L %*% y + K %*% m)
# 	eta_star = eta + y^T %*% y + m^T %*% K %*% M - m_star^T %*% K_star %*% m_star

#----------------------------------------
gibbs = function(iter=11000,burn=1000,thin=2){
	#PURPOSE: Gibbs Sampler for model described above.
	#INPUTS:	iter = number of posterior samples to generate.
	#			burn = number of burn-in values to discard.
	#			thin = thinning to reduce autocorr of chain.
	#OUTPUT:	gibbs = matrix of values sampled from posterior.  Rows = samples.
	
	#Require mvtnorm for sampling from multivariate normal.
	require(mvtnorm)
	
	# Set up objects to hold sampled posterior results.
	mat = matrix(0,nrow=iter, ncol = p + 1 + n)
		# Since lambda is a diagonal matrix, storing just diagonal 
		# elements diag(lambda_1...lambda_n), as 
		# part of mat instead of storing list of matrices.
		
	#Initialize each element of chain.
	lambda = rep(.1,n)
	beta = rep(0,p)
	omega = 1
	
	#Iterate through sampler.
	for (i in 1:iter){
		
		#1. Update beta & relevant posterior parameters.
		XtLX = t(X) %*% diag(lambda) %*% X #Precache to save time.
		K_star = XtLX + K
		m_star = solve(K_star) %*% (t(X) %*% diag(lambda) %*% y + K %*% m)
		
		beta = rmvnorm(1,m_star,(1/omega) * solve(K_star))
		
		#2. Update omega & relevant posterior parameters.
		XtLX = t(X) %*% diag(lambda) %*% X #Precache to save time.
		d_star = d + n
		eta_star = t(y) %*% diag(lambda) %*% y + eta + t(m) %*% K %*% m - t(m_star) %*% K_star %*% m_star 
		#Note: m_star, K_star do not change from beta update; no need to update here.
		
		omega = rgamma(1,d_star/2,rate=eta_star/2)
		
		#3. Update lambda.
		lambda = rgamma(n,(h+1)/2, rate = (omega * (y-X %*% t(beta))^2 + h)/2)
		
		#Add values to results matrix.
		mat[i,] = c(beta,omega,lambda)
	}
	
	#Set up variable names for matrix of results.
	colnames(mat) = c(unlist(strsplit(paste("beta",1:p,sep="")," ")),
		"omega",
		unlist(strsplit(paste("lambda",1:n,sep="")," ")))
	
	#Keep only values after burn-in period.
	mat = mat[(burn+1):iter,]
	
	#Thin values.
	mat = mat[seq(1,nrow(mat),by=2),]
	
	#Return results.
	return(mat)
}
#----------------------------------------

### Read data file & set up model variables (n,p,etc).
data = read.csv('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/statsmod/Course-Data/gdpgrowth.csv',header=T)

n = nrow(data)
X = cbind(rep(1,n),data$DEF60)
y = data$GR6096
p = ncol(X)

#Initialize prior parameters.
m <- rep(0, p)
d = .01
eta = .01
K = diag(c(.01,.01))

h = .01

### Run sampler to obtain samples from joint posterior p(beta,omega,lambda|y).
output = gibbs(iter=11000,burn=1000,thin=2)

### Analysis:

#Plot beta_hat values to show chain stabilized.
pdf(file='/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-02/Figures/Gibbs_beta_stabilized.pdf')
par(mfrow=c(2,1))
plot(1:nrow(output),output[,1],type='l',main='Beta1 Chain',xlab='iter',ylab=expression(beta[1]))
plot(1:nrow(output),output[,2],type='l',main='Beta2 Chain',xlab='iter',ylab=expression(beta[2]))
dev.off()

#Posterior mean of beta.
beta_post_mean = colMeans(output[,1:2])
beta_post_mean

#Omega posterior mean.
mean(output[,3])

#Posterior mean of sigma^2
sigma2_post_mean = 1/mean(output[,3])
sigma2_post_mean

### Lambda_i precision plot.

#Calculate posterior mean of 1/lambda_i for each country (each row X_i is a country).
#We notice that five of the 1/lambda_i values are larger than the others. 
#These correspond to observations with lower precision, ie higher variance.
lam_post_mean_inv = 1/colMeans(output[,4:ncol(output)])
idx_outliers = which(lam_post_mean_inv >= sort(lam_post_mean_inv,decreasing=T)[5])
country_outliers = as.character(data$CODE[idx_outliers])

pdf(file='/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-02/Figures/Lambda_i_Inverse.pdf')
plot(1:n,lam_post_mean_inv,main= 'Posterior Mean of lambda for h=1',xlab='Countries',ylab=expression(1/lambda[i]))
abline(h=sort(lam_post_mean_inv,decreasing=T)[5] - .05,col='blue') #Plots h-line below top 5 points.
dev.off()

#Heavy-tailed Bayesian regression line w outliers labeled.
pdf(file='/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-02/Figures/Gibbs_HeavyTail_Bayes_LM.pdf')
plot(X[,2],y,pch=19,col='black',xlab='Defense Spending',ylab='GDP Growth Rate',main='Heavy-Tailed Bayesian Regression Line (h=1)')
abline(a=beta_post_mean[1], b=beta_post_mean[2],lwd=2,col='blue')
points(X[idx_outliers,2],y[idx_outliers],pch=19,col='red')
text(X[idx_outliers,2],y[idx_outliers],labels=country_outliers,pos=4,offset=.5)
dev.off()

#This plot confirms that these observations are those furthest from the regression line.
#The regression line is more robust to these outliers than the previous lines (Frequentiest and Bayesian).

### Analysis of choice of h under vague priors:

#How does choice h affect things?  Used h=1 previously.  Try a few different h values.

#Look at 1/lambda_i precision plots for a variety of h values.
pdf(file='/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-02/Figures/H_Compare.pdf')
par(mfrow=c(2,3))
H = c(.01,.1,1,2,5,10)
for (g in H){
	h = g
	output = gibbs(iter=11000,burn=1000,thin=2)
	lam_post_mean_inv = 1/colMeans(output[,4:ncol(output)])
	plot(1:n,lam_post_mean_inv,main= paste('Lambda Post. Mean, h=',h),xlab='Countries',ylab=expression(1/lambda[i]))
	abline(h=sort(lam_post_mean_inv,decreasing=T)[5] - .01,col='blue') #Plots h-line below top 5 points.
}
dev.off()