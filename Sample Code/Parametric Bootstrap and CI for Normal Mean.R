###########################################################################
##       SDS 383C - Stats Mod 1 - Fall 2016 - PROBLEM 2-c-ii:            ##
###########################################################################

#Parametric bootstrap for var(theta).
B <- 1000
mu_boot <- rep(0,B)
theta_boot <- rep(0,B)	#Initialize vector to hold bootstrap theta estimates.

mu = 5
sigma = 1
n = 100

th_hat <- exp(mean(rnorm(n,mu,sigma)))

for (i in 1:B){
	#Draw a sample from the distribution, since parametric.
	y_boot <- rnorm(n,mu,sigma)
	
	#Calculate statistics from the bootstrap draw.
	mu_boot[i] <- mean(y_boot)
	theta_boot[i] <- exp(mu_boot[i])
}

#Display bootstrapped variance.
var_theta_boot <- var(theta_boot)
var_theta_boot

#Display 95% CI: (Using percentile interval)
quantile(theta_boot,c(.025,.975))

#Display 95% CI: (Using normal interval, for fun.  This interval is wider.)
th_hat + c(-1,1) * 1.96 * sqrt(var_theta_boot)
