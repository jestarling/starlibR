#Stats Modeling 1
#HW 04

################################################
###  SDS 383C - Stats Mod 1 - Problem 3b-ii  ###
################################################

#EM Algorithm for:

#	y = (y1, y2, y3)^T
#	observed y = (38, 34, 125)^T

#	y are observed counts from Multinomial(.5-t/2, t/4, .5+t/4)

myEM = function(y,theta0,gamma0,iter){
	
	#Note: gamma represents Z4.
	
	#Save each y value as its own variable.
	y1 = y[1]
	y2 = y[2]
	y3 = y[3]
	
	#Initialize vectors to hold iterations of parameters.
	theta = rep(0,iter)
	theta[1] = theta0
	
	gamma = rep(0,iter)
	gamma[1] = gamma0
	
	for(i in 1:iter){
		#Expectation Step
		gamma[i+1] = y3 * (theta[i] / (2 + theta[i]))
		
		#Maximization Step
		theta[i+1] = (y2 + gamma[i+1]) / (y1 + y2 + gamma[i+1]) 
	}
	
	return(list(thetahat=theta[iter], theta=theta, gamma=gamma, iter=iter, t=1:(iter+1)))
}

#Call function with our y values and some sensible starting guesses for the parameter estimates.
y = c(38,34,125)

gamma0 = 50
theta0 = .5

results = myEM(y,theta0,gamma0,iter=10)

#Output theta.hat value and plot theta_t versus iterations t.
results

jpeg(file='/Users/jennstarling/UTAustin/2016_Fall_SDS 383C_Statistical Modeling 1/Homework/HW 04/3-b-ii.jpg')
plot(results$t,results$theta,type='l',col='blue', xlab='iterations (t)', ylab='theta',
	main=paste('EM Convergence to theta = ',round(results$thetahat,4)))
dev.off()