#ADAPTIVE REJECTION SAMPLING FUNCTION

#Same code as Exercise 3, just with different functions.
#More general, as no 'a' parameter for h function.

#Set up original function fx, hx=logfx, and hpx=d/dx (hx).
f  = function(x) x^3 * exp(-2*x) * (1-exp(-2*x))^3
h  = function(x) 3*log(x) - 2*x + 3*log(1-exp(-2*x))
hp = function(x) 3/x - 2 + 6/(exp(2*x)-1)

#Inspect log-concave h for good starting points.
pdf(file='/Users/jennstarling/UTAustin/2017S_MCMC/Exercises/Exercise-04/Problem_1_ARS_logfx.pdf')
x = seq(0,10,by=.01)
hx = h(x)
plot(x,hx,type='l',col='blue',main='logfx to eyeball starting points')
dev.off()

#Call function.
fx_samp = adaptive_rej_sampling(f,h,hp,n=10000,c(.1,7),xlb=-Inf,xub=Inf)

#Plot results.
pdf(file='/Users/jennstarling/UTAustin/2017S_MCMC/Exercises/Exercise-04/Problem_1_ARS.pdf')
fx_samp$p_accept
hist(fx_samp$fx_sample,freq=F,main='Adaptive Rejection Sampling',xlab='Samples of f(x)')
x=seq(0,10,by=.001)
lines(x,3*f(x),col='blue')
dev.off()

##############################
###   MAIN ARS FUNCTION:   ###
##############################

adaptive_rej_sampling = function(f,h,hp,n=10,x,xlb=0,xub=1){
	#FUNCTION: Performs adaptive rejection sampling.
	#
	#INPUTS:
	#	f = target function f(x).  Must be log-concave.
	#	h = log of target function, logf(x).
	#	hp = first derivative of log of target function, d/dx log(f(x))
	#	M = max of f(x).  (Must be a bounded function.)
	#	n = sample size.  (Returned sample will be lower size, bc some rejections.)
	#	x = vector of starting points.  Must be length 2 or longer.
	#	xlb = lower bound of domain for f(x).  Can be -Inf.
	#	xub = upper bound of domain for f(x).  Can be Inf.
	#
	#OUTPUTS:
	#	fx_sample = vector of realizations sampled from f(x).
	#	p_accept = prob(acceptance); accepts/n.
	
	fx_sample = numeric() 	#Empty vector to store sampled values from fx.
	accepts = 0		#For tracking P(Acceptance).

	for(ns in 1:n){
		
		#1. Update envelope: piecewise linear envelope in log-scale.
		new_envelope = update_envelope(x,h,hp,xlb,xub)
		
		#-----------------------
		
		#2. Calculate areas of piecewise exponential functions based on
		#exponentiating each tangent line over its range z_p to z_p+1.
		Areas = calculate_exp_areas(new_envelope)
		
		#-----------------------
		
		#3. Randomly select an area, based on uniform sampling, weighted by areas.
		A_idx = sample(1:length(Areas),size=1,replace=F,prob=Areas/sum(Areas))
		#A_idx = sample(1:length(Areas),size=1,replace=F)
		
		#-----------------------
		
		#4. Sample an x* from the truncated exponential distibution from selected area.
		hx = new_envelope$hx[A_idx]
		hpx = new_envelope$hpx[A_idx]		#lambda (rate) for the selected exponential.
		z1 = new_envelope$z[A_idx]			#Lower bound for truncated exponential.  Shifts to zero for sampling.
		z2 = new_envelope$z[A_idx+1]		#Upper bound for truncated exponential.

			#Note: we only know how to sample from exponential at zero.  We must do two things:
			#1. Check if slope hpx_star is +/-.
			#2. Adjust x* realization by shifting it.
		
		#Shift and rotate the exponential sample as needed.
		if(hpx<0){
			#If slope negative, generate using neg slope, and shift result.
			x_star = rtexp(n=1,m=-hpx,t=z2-z1)
			x_star = z1 + x_star
		} else{
			#If slope positive, generate using slope, and shift/flip result.
			x_star = rtexp(n=1,m=hpx,t=-(z1-z2))	
			x_star = z2 - x_star
		}
		
		#-----------------------
		
		#5. Acceptance test:
		u = runif(1,0,1)	#Generate one unif(0,1) realization.
		
		#Set up the proposal function.  This is the exponential function selected
		#for sampling, above, for the selected area.
		x1 = new_envelope$x[A_idx]		#For calculation of value of tangent line at x_star.
		b =  new_envelope$b[A_idx]		#The slope for the calculation of the tangent line value at x_star.
		
		log_g = hp(x1) * x_star + b	#Value of piecewise linear envelope at x_star.
		g = exp(log_g)					#Value of piecewise exponential at x_star.
	
		#Accept x_star as from f(x) if fx/gx <= a uniform draw.
		#	If accept, then add x_star to fx_sample vector.
		#	If reject, then add x_star to x vector of points to sample, and start over.
		
		if (u <= f(x_star) / g){
			fx_sample = c(fx_sample,x_star)
			accepts = accepts + 1
		} else{
			x = c(x,x_star)
		}
		
	} #Start over for another sample until reach n tries.
	
	return(list(fx_sample=fx_sample,p_accept = accepts/n))
}

#############################
###   HELPER FUNCTIONS:   ###
#############################

update_envelope = function(x,h,hp,xlb,xub){
	#FUNCTION: Update piece-wise linear envelope of logf(x).
	#
	#INPUTS:
	#	x = vector of points (unsorted)
	#	h = function for log(f(x)).
	#	hp = function for d/dx log(f(x)).
	#	xdomain = lower and upper bounds for x.
	#
	#OUTPUTS:
	#	x = vector of points, sorted
	#	hx = h(x) = log(f(x)) values for each x.
	#	z = intersection points for tangents (in between each x value)
	#	m = slope for tangent lines.
	#	b = y-intercept for tangent lines
	
	x = sort(x)
	hx = h(x)	#Values of logf(x) at each x.
	hpx = hp(x)	#Slopes

	m = hpx			#slopes
	b = hx - hpx*x	#intercepts
	
	#Calculate points where tangents meet.
	z = rep(0,length(x)-1)	
	
	for (i in 1:length(z)){
		z[i] = (hx[i+1] - hx[i] - x[i+1]*hpx[i+1] + x[i]*hpx[i]) / (hpx[i] - hpx[i+1])
	}
	
	z = c(xlb,z,xub)	#Add in bounds
	
	return(list(x=x,hx=hx,hpx=hpx,m=m,b=b,z=z))
}

calculate_exp_areas = function(envelope){
	#FUNCTION: Calculate area for each part of piecewise
	#			exponential functions created by exponentiating 
	#			piecewise linear envelope.
	#INPUTS:
	#	envelope: Takes in a piecewise linear envelope function, as created by
	#				update_envelope().
	#OUTPUTS:	
	#	areas: A vector of areas, A1 to Ak, where k = length(z).
	
	x = envelope$x	#Length of x is number of areas to calculate.
	z = envelope$z	#Extract z elements from envelope.
	hx = envelope$hx
	hpx = envelope$hpx				
	
	Areas = rep(0,length(x)) #Empty vector to hold areas.
	
	#Loop through areas.
	for (i in 1:length(Areas)){
		Areas[i] = (exp(hx[i] - hpx[i]*x[i]) / hpx[i]) * (exp(hpx[i]*z[i+1]) - exp(hpx[i]*z[i]))	
	}	
	
	return(Areas)
}

rtexp = function(n,m,t){
	#PURPOSE: Draws n random samples from inverse cdf of truncated exponential.
	#n = number of samples.
	#m = rate
	#t = level of truncation; x value at which to truncate.
	u = runif(n)
	itex = -log(1-u*(1-exp(-t*m)))/m
	return(itex)
}

#NOT USED.  Embedded in update_envelope function.
tangent_intersect = function(x1,x2,hx1,hx2,hpx1,hpx2){
	#PURPOSE: Returns x value where two tangent lines intersect.
	#x1, x2 = two points
	#hx1, hx2 = function values evaluated at x1 and x2 for some h(x).
	#hpx1,hpx2 = derivative values evaluated at x1 and x2 for h'(x).
	
	xtan = (hx2 - hx1 - hx2*hpx2 + x1*hpx1) / (hpx1-hpx2)
	return(xtan)
}

#NOT USED.
itexp = function(u,m,t){
	#PURPOSE: Computes the inverse cdf for the truncated exponential at a specified quantile u.
	#u = quantile
	#m = rate
	#t = level of truncation; x value at which to truncate.
	return(-log(1-u*(1-exp(-t*m)))/m)
}