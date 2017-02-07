#MCMC Exercise 4
#Problem 3
#Feb 7, 2017
#Jennifer Starling

#Let f(x) be the density: f(x) \propto exp(-2 * abs(x)^4), -inf < x < inf.
#Using importance sampling, estimate I = integral_{-inf}^{inf} cosh(x)*f(x) dx.

#Target density.
f = function(x) exp(-2 * abs(x)^4)
g = function(x) ((1+exp(2*x)) * exp(-2 * abs(x)^4))/2

n = 1000	#1000 samples

#Sampling density is exp(1).
x = rexp(n,1)
gx = g(x)

#Approximate integral value.
Ihat_N = 2 * (1/n) * sum(gx)
Ihat_N #Dispay result: was 1.7433

#Visualize target density. (Confirming symmetric and can get away w splitting integral.)
pdf(file='/Users/jennstarling/UTAustin/2017S_MCMC/Exercises/Exercise-04/Tex/Figures/Prob3_fx.pdf')
x = seq(-10,10,by=.01)
plot(x,f(x),type='l',main='Visualizing target density')
dev.off()