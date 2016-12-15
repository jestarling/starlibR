
###################################
###        PROBLEM :3           ###
###################################

#Read in data.
bcdata = read.table(file='/Users/jennstarling/UTAustin/2016_Fall_SDS 383C_Statistical Modeling 1/Homework/HW 03/breast-cancer.data')
colnames(bcdata) = c( 'CT','CS','CSh','MA','SECS','BN','BC','NN','Mit','class')

#Recode class as 1's and 0's.  (Malignant = 4, Benign = 2)
bcdata$class = ifelse(bcdata$class==4,1,0)

#Save unscaled data set for Gaussian
bcdata_unscaled = bcdata

#Scale just the x columns.
#bcdata=as.matrix(cbind(scale(bcdata[,1:9]),bcdata[,10,drop=F]))

#Set up percentages of training set to sample.
train.pct = c(.01,.02,.03,.125,.625,1)

#Set up matrix to hold NB and LR prediction errors for each train size for each of 5 iterations.
pred.err.nb = matrix(0,nrow=5,ncol=length(train.pct))
pred.err.nbg = matrix(0,nrow=5,ncol=length(train.pct))
pred.err.lr = matrix(0,nrow=5,ncol=length(train.pct))

#Loop through 5 random train/test splits, using train=2/3 of obs.
for (i in 1:5){
	
	#Randomly split training/test data, using 2/3 obs for training.
	split = split.data(bcdata,train.pct=.66)
	
	#Save unscaled data for Naive Bayes.
	test_unscaled = split$test
	train_unscaled=split$train
	
	#Save scaled data for LR.  (Scale X cols only, not Y.)
	test = as.matrix(cbind(scale(split$test[,1:9]),split$test[,10,drop=F]))
	train = as.matrix(cbind(scale(split$train[,1:9]),split$train[,10,drop=F]))
	
	#Loop through each training set percent size.
	for (j in 1:length(train.pct)){
		
		#Sample a train.pct-sized subset of the training obs to use as training set.
		idx = sample(1:nrow(train),size=train.pct[j]*nrow(train),replace=F)
		tr = train[idx,]
		tr_unscaled = train_unscaled[idx,]
		
		#-------------------------------
		#LOGISTIC REGRESSION/IRLS STEPS:
		#Run IRLS for training data subset.
		output.irls = irls(X=tr[,1:9],y=tr[,10])
		
		#Predict IRLS yhat values on full test set.
		logodds = test[,1:9] %*% output.irls$beta_hat
		odds = exp(logodds)
		yhat.prob = odds/(1+odds)
		yhat=ifelse(yhat.prob>.5,1,0)
		
		#Calculate prediction error.
		pred.err.lr[i,j] = sum(yhat!=test[,10])/length(yhat)
		#-------------------------------
		#NAIVE BAYES STEPS:
		
		#Multinomial NB
		output.nb.multinom = naiveBayesMultinomial(data=tr_unscaled,testdata=test_unscaled)
		pred.err.nb[i,j] = output.nb.multinom$pred.err
		
		#Gaussian NB
		output.nb.gauss = naiveBayesGaussian(data=tr_unscaled,testdata=test_unscaled)
		pred.err.nbg[i,j] = output.nb.gauss$pred.err
		#-------------------------------
		
	} #end j loop through train.pct vector
} #end i loop through 5 random test/train splits

#Calculate average prediction error for each method for each training size pct.
avg.pred.err.lr = colMeans(pred.err.lr)
avg.pred.err.nb = colMeans(pred.err.nb)
avg.pred.err.nbg = colMeans(pred.err.nbg)

jpeg(file='/Users/jennstarling/UTAustin/2016_Fall_SDS 383C_Statistical Modeling 1/Homework/HW 03/nb_lr_prederr.jpg')
plot(train.pct,avg.pred.err.nb,type='l',xlab='train.pct',ylab='pred.error',
	main='Naive Bayes & Logistic Regression (IRLS)',xlim=c(.05,1))
lines(train.pct,avg.pred.err.lr,type='l',col='blue',lty=2)
legend(.7,.25,c('lr','nb'),col=c('black','blue'),lty=c(1,2))
dev.off()

#Note: I calculated Gaussian Naive Bayes also just for fun.  Not shown on plot, since paper used multinomial.
#Function is at bottom of R file.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Data splitting function:
#	Inputs: A data set and a percent of the data set to use as training data.
#	Output: A train and a test data set.
split.data = function(data,train.pct){
	
	#Split test and training data.  Use 2/3 of each class for train.
	class0.idx = which(data[,"class"]==0)
	class1.idx = which(data[,"class"]==1)

	#Number of obs in each class.
	nclass0 = table(data[,"class"])[1]
	nclass1 = table(data[,"class"])[2]

	#Sample 2/3 of the samples from each class.
	train0.idx = sample(class0.idx,nclass0*train.pct)
	train1.idx = sample(class1.idx,nclass1*train.pct)

	#Construct test and training sets.
	train = bcdata[c(train0.idx,train1.idx),]
	test = bcdata[-c(train0.idx,train1.idx),]
	
	return(list(train=train,test=test))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Iteratively Reweighted Least Squares Function
#(Binomial Distribution)
#	Inputs:
#		X = design matrix
#		y = vector of 1/0 response values
#		maxiter = maximum number of iterations allowed
#		conv = tolerance level for convergence check
#	Outputs:
#		loglik = loglikelihood function (should be maximized)
#		beta_hat = estimated coefficients
#		iter = iterations completed
#		converged = convergence indicator

irls = function(X,y,maxiter=500,conv=1E-14){
	
	#------------------------------------------------------------------	
	#Binomial Loglikelihood function. 
		#Inputs: 	X = design matrix 
		#			Y = vector of 1, 0 values. 
		#   		m = sample size vector m.
		#			w = probabilities vector.  (Diagonal of W matrix.)
		#Output: Binomial loglikelihood value.
	logl <- function(X,y,beta,m,w){
		logl <- sum(y*log(w+1E-4) + (m-y)*log(1-w+1E-4)) #Calculate log-likelihood.
		return(logl)	
	}
	#------------------------------------------------------------------
	
	#Indicator for whether convergence met.
	converged <- 0		
	
	#Set m vector (sample sizes) as a vector of 1's.
	m = rep(1,nrow(X))	
	
	#1. Initialize matrix to hold beta vector for each iteration.
	betas <- matrix(0,nrow=maxiter+1,ncol=ncol(X)) 
	betas[1,] <- rep(0,ncol(X))	#Initialize beta vector to 0 to start.
	
	#1a. Initialize weights.
	w = 1 / (1 + exp(-X %*% betas[1,]))
	
	#2. Initialize values for log-likelihood.
	loglik <- rep(0,maxiter) 	#Initialize vector to hold loglikelihood fctn.
	loglik[1] <- logl(X,y,betas[1,],m,w)

	#3. Perform Newton's method iterations.	
	for (i in 2:maxiter){
		
		#Calculate probabilities vector w_i.
		w = 1 / (1 + exp(-X %*% betas[i-1,]))
		
		#Set up W matrix, using w as diagonal elements.	
		W = diag(as.vector(w))
		
		#Adding a small regularization penalty to avoid singularity.
		#Using lambda=1E-5.
		lambdaI = diag(1E-5,nrow=ncol(X))
		
		#Update betas.
		betas[i,]  = solve(t(X) %*% W %*% X + lambdaI) %*% t(X) %*% W %*% y
		
		#Update loglikelihood.
		loglik[i] = logl(X,y,betas[i,],m,w)

		#Check if convergence met: If yes, exit loop.
		if (i>2 && abs(loglik[i]-loglik[i-1])/abs(loglik[i-1]+1E-3) < conv ){
			converged=1;
			break;
		}
	} #end IRLS iterations.
	
	return(list(beta_hat=betas[i,],iter=i,logl = loglik[1:i],converged=converged))
} #End IRLS function.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

naiveBayesGaussian = function(data,testdata){
	#Data must be unscaled.
	
	n = nrow(data)						#Number of obs.
	p = ncol(data)-1					#Number of predictors.
	k = length(unique(data[,"class"]))	#Number of response class buckets.
	xlevs = 10							#Number of levels of each x variable.
	
	#Calculate P(Y=y)
	PY = (table(data[,"class"])/n)[2]
	
	#Calculate means and variances for each class for each predictor.

	#Separate data into Y=0 and Y=1 to make conditional probs easier to calculate. 
	data.y1 = data[which(data[,"class"]==1),]
	data.y0 = data[which(data[,"class"]==0),]

	#Calculate Gaussian means and variances for each class, for each predictor Xi.
	x.means.y1 = colMeans(data.y1)
	x.means.y0 = colMeans(data.y0)
	 
	x.vars.y1 = colVars(as.matrix(data.y1))
	x.vars.y0 = colVars(as.matrix(data.y0))
	
	#Prob matrices will hold P(Xi | Y) for each predictor (col), with 10 buckets (rows) for categories of predictor.
	#X = {1,2,3,4,5,6,7,8,9,10}
	PX.given.Y0 = matrix(0,nrow=10,ncol=9)
	PX.given.Y1 = matrix(0,nrow=10,ncol=9)
	
	#Fill in probabilities:  (Adding 1 before normalizing to avoid 0 probabilities.)
	for (i in 1:p){
		PX.given.Y1[,i] = dnorm(1:10,x.means.y1[i],sqrt(x.vars.y1[i]))
		PX.given.Y0[,i] = dnorm(1:10,x.means.y0[i],sqrt(x.vars.y0[i]))
		
	}
	
	#Calculate probabilities for new data:
	Y.test = rep(0,nrow(testdata))
	
	for (i in 1:nrow(testdata)){
		xlevels = unlist(testdata[i,-10])	#Extract X values for observation i.
		
		#Extract P(X|Y=1) and P(X|Y=0) values for x permutation and calculate product of probabilities.
		#(Using product due to the Naive Bayes conditional independence assumption.)
		PX.g.Y1 = prod(PX.given.Y1[cbind(xlevels,1:9)])		
		PX.g.Y0 = prod(PX.given.Y0[cbind(xlevels,1:9)])
		
		ratio = (PY*PX.g.Y1) / ((1-PY)*PX.g.Y0)	#Ratio of P(Y=1|X) to P(Y=0|X), using bayes rule.
		
		#Calculate predicted yhat value.
		Y.test[i] = ifelse(ratio > 1, 1, 0)	
	}
	
	#Calculate prediction error.
	pred.err = sum(Y.test!=test[,10])/length(Y.test)
	
	#Return function outputs.
	return(list(yhat=Y.test,pred.err=pred.err))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Naive Bayes Function:
naiveBayesMultinomial = function(data,testdata){
	#Data must be unscaled.
	
	n = nrow(data)						#Number of obs.
	p = ncol(data)-1					#Number of predictors.
	k = length(unique(data[,"class"]))	#Number of response class buckets.
		
	#Calculate P(Y=y)
	PY = (table(data[,"class"])/n)[2]
	
	#Calculate P(X=x|Y) for both Y=1 and Y=0
	
	#a. Separate data into Y=0 and Y=1 to make conditional probs easier to calculate. 
	data.y0 = data[which(data[,"class"]==0),]
	data.y1 = data[which(data[,"class"]==1),]
	
	#b. Prob matrices will hold P(Xi | Y) for each predictor (col), with 10 buckets (rows) for categories of predictor.
	#X = {1,2,3,4,5,6,7,8,9,10}
	PX.given.Y0 = matrix(0,nrow=10,ncol=9)
	PX.given.Y1 = matrix(0,nrow=10,ncol=9)
	
	#Fill in probabilities:  (Adding 1 before normalizing to avoid 0 probabilities.)
	for (i in 1:p){
		PX.given.Y0[,i] = (table(factor(data.y0[,i],lev=1:10))+1)/n	
		PX.given.Y1[,i] = (table(factor(data.y1[,i],lev=1:10))+1)/n
	}
	
	#Calculate probabilities for new data:
	Y.test = rep(0,nrow(testdata))
	
	for (i in 1:nrow(testdata)){
		xlevels = unlist(testdata[i,-10])	#Extract X values for observation i.
		
		#Extract P(X|Y=1) and P(X|Y=0) values for x permutation and calculate product of probabilities.
		#(Using product due to the Naive Bayes conditional independence assumption.)
		PX.g.Y1 = prod(PX.given.Y1[cbind(xlevels,1:9)])		
		PX.g.Y0 = prod(PX.given.Y0[cbind(xlevels,1:9)])
		
		ratio = (PY*PX.g.Y1) / ((1-PY)*PX.g.Y0)	#Ratio of P(Y=1|X) to P(Y=0|X), using bayes rule.
		
		#Calculate predicted yhat value.
		Y.test[i] = ifelse(ratio > 1, 1, 0)	
	}
	
	#Calculate prediction error.
	pred.err = sum(Y.test!=test[,10])/length(Y.test)
	
	#Return function outputs.
	return(list(yhat=Y.test,pred.err=pred.err))
}

naiveBayesMultinomial(data=train,testdata=test)


naiveBayesGaussian(data=train,testdata=test)
