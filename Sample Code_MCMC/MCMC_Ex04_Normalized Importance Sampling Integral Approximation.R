#MCMC Exercise 4
#Problem 3
#Feb 7, 2017
#Jennifer Starling

#Let f(x) be the density: f(x) \propto exp(-2 * abs(x)^4), -inf < x < inf.
#Using importance sampling, estimate I = integral_{-inf}^{inf} cosh(x)*f(x) dx.

#Set up appropriate functions.
f_star = function(x) exp(-2 * abs(x)^4)	#Target density up to a constant of proportionality.
g = function(x) cosh(x)					#g(x) density.
h = function(x) dnorm(x,0,1)			#h(x) density.  (Will match h(x) used to generate iid x_i's.)
h_gen = function(n) rnorm(n,0,1)	

#----------------------------------------------------
norm_import_samp = function(n=1000,f_star,g,h,h_gen){
	#FUNCTION: Normalized Importance Sampling.
	#PURPOSE: 	Estimates integral I = E_f(g(x)) = int g(x)*f(x)dx
	#INPUTS:	f_star: function, the target density up to a constant of proportionality.
	#			g: function, the function multiplied by f(x) in the integral.
	#			h: function, the pdf of the function used to generate h(x_i).  	Example: dnorm(x,0,1)
	#			h_gen: function, generates realizations from the h density.		Example: rnorm(n,0,1)
	#OUTPUTS: 	Integral estimate Ihat_N.
	
	#Generate n idd observations from h(x).
	x = h_gen(n)
	
	#Approximate denominator; integral of f_star(x) dx.
	Ihat_denom = (1/n) * sum(f_star(h(x)) / h(x))

	#Approximate numerator; (1/N) sum_{i=1}^{n} g(x)f*(x) / h(x)
	Ihat_num = (1/n) * sum(g(hx)*f_star(h(x)) / h(x))

	Ihat_N = Ihat_num / Ihat_denom
	Ihat_N
	
	return(Ihat_N)
}
#----------------------------------------------------

#Run function and display results for single Ihat_N estimate.
Ihat_N = norm_import_samp(n=1000,f_star,g,h,h_gen)
Ihat_N

#Seems to have high variability.  Run B=10000 times for mean.
B=1000
Ihat_N_vec = rep(0,B)

for (b in 1:B){
	Ihat_N_vec[b] = norm_import_samp(n=1000,f_star,g,h,h_gen)
}

mean(Ihat_N_vec)
#Dispay result: was 1.7786, which is close to true integral value for this example.

#NOTE: How to pick hx? Want something similar in shape to fx.
#This one looks fairly normal.
x = seq(-10,10,by=.01)
y = f_star(x)
plot(x,y,type='l',col='blue')