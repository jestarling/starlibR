#MCMC Exercise 01
#R Code
#Jennifer Starling
#Jan 2017

Ihat_fct = function(N,reps){
	
	Ihat = rep(0,reps)
	N = 1000
	
	for (i in 1:reps){
		x = rcauchy(N,0,1)
		Ihat[i] = pi * (1/N) * sum(x)
	}
	
	var = var(Ihat)
	return(list(var=var,Ihat=Ihat))
}

results = Ihat_fct(N=1000,reps=1000)

results$var
hist(results$Ihat)

par(mfrow=c(3,4))
vars = rep(0,12)
for (j in 1:12){
	results = Ihat_fct(N=1000,reps=1000)
	hist(results$Ihat)
	vars[j] = results$var
}

vars

#----------------------------------

#Problem 6 ctd.

N = 1000
a = 2 #shape
b = 2 #scale

Ihat_fctn = function(N,reps){
	
	Ihats = rep(0,reps)
	
	save_g = matrix(0,nrow=N,ncol=reps)
	
	for (i in 1:reps){
		x = rgamma(N, shape=a, scale=b)

		gx = x / ( (1+x^2) * (b^a)/(gamma(a)) * x^(a-1) * exp(-b*x))

		Ihat[i] = (1/N) * (sum(gx) - sum(gx))
		save_g[,i] = gx
	}
	
	var = var(Ihat)
	return(list(Ihat=Ihat,var=var,gx = save_g))
}

results = Ihat_fctn(N=1000,reps=1000)

results$var
hist(results$Ihat)

par(mfrow=c(3,4))
vars = rep(0,12)
for (j in 1:12){
	results = Ihat_fctn(N=1000,reps=1000)
	hist(results$Ihat,freq=T)
	vars[j] = results$var
}

round(range(apply(results$gx,2,range)),2)

vars

