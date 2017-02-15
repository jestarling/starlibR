#Stats Modeling 2 - Exercise 2
#Jennifer Starling
#Feb 6, 2017

#Housekeeping
rm(list=ls())

#Read data file.
data = read.csv('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/statsmod/Course-Data/gdpgrowth.csv',header=T)

#Load libraries.
library(mvtnorm)

#Consider a linear model with intercept for GDP growth rate (GR6096) versus
#its level of defense spending as a fraction of GDP (DEF60).

#Fit the Bayesian linear model to this data set, choosing Lambda = I, and
#something diagonal and vague for the prior precision matrix K = diag(k1,k2).

#Inspect the fitted line (graphically).  Are you happy with the fit?  Why/not?

#---------------------------------------
#Set up model variables.
n = nrow(data)
X = cbind(rep(1,n),data$DEF60)
y = data$GR6096
p = ncol(X)

#####################################
###   Frequentist linear model.   ###
#####################################

freq_lm = lm(y~X-1)
summary(freq_lm)
beta_hat = freq_lm$coefficients

#Visually inspect fit.
#pdf(file='/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-02/LaTeX Files/Figures/Frequentist_Linear_Model.pdf')
plot(X[,2],y,pch=19,col='black',xlab='X',ylab='y',main='Frequentist Linear Model')
abline(a=beta_hat[1], b=beta_hat[2],lwd=2,col='blue')
#dev.off()

################################################
###   Bayesian linear model (vague prior).   ###
################################################

# 	(y | B, sigma2) ~ N(XB, (wL)^{-1})
#	(B|w) ~ N(m,(wK)^{-1})
#	w ~ Gamma(d/2,eta/2)

#Set up L matrix and some vague prior parameters.
L = diag(n)	#Begin with Lambda = I.

m <- rep(0, p)
d = .01
eta = .01
K = diag(c(.01,.01))

#Precache X^T %*% L %*% X for use in hyperparameter calculation.
XtLX = t(X) %*% L %*% X

#Update posterior parameters.
d_star = d+n
K_star = XtLX + K
m_star = solve(K_star) %*% (t(X) %*% L %*% y + K %*% m)
eta_star = eta + t(y) %*% y + t(m) %*% K %*% m - t(m_star) %*% K_star %*% m_star

# Calculate posterior beta_hat estimate.  This is posterior mean of beta.
# From Part C, p(beta|y) ~ t(m*,d*,K)
beta_hat_post = m_star

#Plot frequenstist and Bayesian lm results over data for comparison.
pdf(file='/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-02/Figures/Freq_Bayes_LM_Compare.pdf')
plot(X[,2],y,pch=19,col='black',
	main='Frequentist & Bayesian Linear Models',
	xlab='Defense Spending',ylab='GDP Growth Rate')
abline(a=beta_hat[1], b=beta_hat[2],lwd=2,col='blue')
abline(a=beta_hat_post[1],b=beta_hat_post[2],lwd=2,col='firebrick3')
legend('topright',legend=c("Freq LM","Bayes LM"),lwd=2,lty=1,col=c('blue','firebrick3'))
dev.off()