#Inversion of digamma function using 
#"Estimating a Dirichlet Distribution" by Thomas P. Minka, 2000, Appendix 3 method.
inv.digamma <- function(y){
	
  #Initialize psi inverse iterations.	
  psiinv_y <- array(NA, dim=c(6,length(y)))

  #Set gamma to -Psi(1)
  gamma <- -digamma(1)
  
  #Initialize psiinv_y values.
  for (i in 1:length(y)){
    if (y[i] >= -2.22)
      psiinv_y[1,i] <- exp(y[i]) + 0.5
    else
      psiinv_y[1,i] <- -1/(y[i]+gamma)
  }
  
  #Iterate through Newton's algorithm to update psiinv_y values.
  #Only 6 iterations needed, per author.
  for (i in 2:6){
    psiinv_y[i,] <- ( psiinv_y[i-1,] - 
		(digamma(psiinv_y[i-1,]) - y)/(trigamma(psiinv_y[i-1,])) )
  }
  
  return (psiinv_y[6,])
}