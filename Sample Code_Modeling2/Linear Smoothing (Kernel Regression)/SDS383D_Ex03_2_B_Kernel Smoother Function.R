#Stats Modeling 2
#Exercise 3
#Curve Fitting by Linear Smoothing - Part B 
#Jennifer Starling
#Feb 15, 2017

#------------------------------------------------------------
###   	Simulate noisy data from a non-linear function,
###		subtract the sample means from the simulated x and y,
###		and use the smoother function to fit the kernel smoother
###		for some choice of h.  Plot estimated functions for a 
###		range of bandwidths large enough to yield noticeable changes
###		in the qualitative behavior of the prediction functions.
#------------------------------------------------------------

### Environment setup.

#Housekeeping.
rm(list=ls())

#Load functions.
source("/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-03/R Files/SDS383D_Ex3_FUNCTIONS.R")

#------------------------------------------------------------
### Data generation.

# Simulate x and y from a non-linear function x.  
x = runif(100,0,1)	#Generate a sample of x values.
f = function(x) sin(2*pi*x)	#Set the non-linear function.

#Generate noisy data, and extract x and y.
y = sim_noisy_data(x,f,sig2=.75)

#Subtract means from simulated x and y.
x = x - mean(x)
y = y - mean(y)

#Plot noisy data, with means subtracted.
pdf(file='/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-03/Figures/noisydata.pdf')
plot(x,y,main='Noisy data',xlab='x',ylab='y')
dev.off()

#Set up a vector of x_star values, and estimate function value at these points.
x_star = seq(min(x),max(x),by=.01)
fhat_star = linear_smoother(x,y,x_star,h=.01,K=K_gaussian)

plot(x,y,main='Gaussian Kernel Smoother')
lines(x_star,fhat_star,col='red')

#------------------------------------------------------------
### Plot linear smoothing output for various bandwidths.

#Set up bandwidths to try, and corresponding colors.
H = c(.01,.1,.5,1,5)
col = rainbow(length(H))

#Open pdf file for two plots.
pdf(file='/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-03/Figures/Smoothers.pdf')
par(mfrow=c(2,1))

#Plot for a variety of bandwidths using Gaussian kernel.
plot(x,y,main='Gaussian Kernel Smoother')
for (i in 1:length(H)){
	fhat_star = linear_smoother(x,y,x_star,h=H[i],K=K_gaussian)
	lines(x_star,fhat_star,col=col[i],lwd=2)
}
legend('topleft',legend=paste("h=",H,sep=''),lwd=2,lty=1,col=col,bg='white')

#Plot for a variety of bandwidths using uniform kernel.
plot(x,y,main='Uniform Kernel Smoother')
for (i in 1:length(H)){
	fhat_star = linear_smoother(x,y,x_star,h=H[i],K=K_uniform)
	lines(x_star,fhat_star,col=col[i],lwd=2)
}
legend('topleft',legend=paste("h=",H,sep=''),lwd=2,lty=1,col=col,bg='white')

dev.off() #Close pdf file for two plots.


