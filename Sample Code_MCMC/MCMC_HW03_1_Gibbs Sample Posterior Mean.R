#MCMC Homework 3
#Problem 3
#Feb 9, 2017
#Jennifer Starling

# Joint posterior is f(a,b) \propto a^4 * b^6 * exp(-a-b-3*a*b)
# for a>0, b>0.
# Find full conditionals f(a|b) and f(b|a).
# Implement chain and evaluate I = int a*b*f(a,b) da db (posterior means).

#Set chain length and thinning/burn-in period parameters.
B=11000
burnin = 1000
thin = 2

#Initialize joint posterior samples a and b to begin chain.  
a=c(1,rep(0,B-1))
b=c(1,rep(0,B-1))

#Run chain.
for (i in 2:B){
	a[i] = rgamma(1, 5, 1+3*b[i-1])
	b[i] = rgamma(1, 7, 1+3*a[i])
}

#Burn and thin.
a = a[-c(1:burnin)]
b = b[-c(1:burnin)]

a = a[seq(1,length(a),by=thin)]
b = b[seq(1,length(b),by=thin)]

#Plot resulting Markov Chain to show stationary, along with histogram.
pdf(file='/Users/jennstarling/UTAustin/2017S_MCMC/Homework/Homework 03/Figures/1_Markov_Chain.pdf')
par(mfrow=c(2,1))
plot(a,type='l',main='Chain Trace: a')
plot(b,,type='l',main='Chain Trace: b')
dev.off()

#Calculate posterior means. (Evaluate integral.)
mean(a*b)