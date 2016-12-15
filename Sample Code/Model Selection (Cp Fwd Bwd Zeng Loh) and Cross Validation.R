#Stats Modeling 1, HW 3
#Jennifer Starling
#25 Sept 2016

rm(list=ls()) #Begin with clean workspace.
setwd('/Users/jennstarling/UTAustin/starlib/data/')

library(leaps) #For Forward and Backward stepwise selection.
library("Hmisc") #For error bars in plot.

###################################
###        PROBLEM 1:           ###
###################################

#Read data:
data <- read.csv(file='mpg_data.csv',header=T)
attach(data)

#-------------------------------------------------------------------
#Part A: Fit a multiple linear regression model to predict MPG.

X <- data[,c("VOL","HP","SP","WT")]
y <- data$MPG

mpgLM <- myLinearReg(X,y)
mpgLM

#For fun, can check my results against the built-in lm function.
mylm <- lm(MPG ~ VOL + HP + SP + WT,data=data,model=T)
summary(mylm)

#---FUNCTION-------------------------
#My custom linear regression function.
#	Inputs: X = matrix of predictors (no intercept col).
#			y = vector of responses.
#	Output: Linear model details.
myLinearReg <- function(X,y){

	n = length(y)
	X = as.matrix(cbind(rep(1,n),X))
	beta_hat <- solve(t(X) %*% X) %*% t(X) %*% y

	#Calculate y_hat and RSS.
	y_hat <- X %*% beta_hat
	RSS <- t(y - X %*% beta_hat) %*% (y - X %*% beta_hat)

	n <- nrow(X)	#Number of obs.
	p <- 4			#Number of predictors.
	sigma2_hat <- RSS / (n-p-1)		#Estimate of sigma^2.
	sigma_hat <- sqrt(sigma2_hat)	#Estimate of sigma.

	#Estimated standard error and CI for each coefficient:
	beta_hat_se <- sqrt( sigma2_hat * diag(solve(t(X) %*% X)) ) #Standard errors.

	#95% CI for beta_hat's.
	lb <- beta_hat -  1.96 * beta_hat_se
	ub <- beta_hat + 1.96 * beta_hat_se
	sig <- ifelse( (lb < 0) & (ub > 0), 0,1)

	#Display CI, including star for betas not equal to zero based on the CI.
	ci = data.frame(lower.bound=lb,upper.bound=ub,sig.coefs=sig)
	
	return(list(beta_hat = beta_hat, 
				RSS = RSS, 
				sigma_hat = sigma_hat, 
				beta_hat_se = beta_hat_se, 
				ci = ci,
				n=n,
				p = ncol(X)-1))
} #End myLinearReg function.
#---END FUNCTION---------------------

#-------------------------------------------------------------------
#Part B: Use Mallow Cp to select a best sub-model.  To search through models, use 
# (i) Forward stepwize, and (ii) Backward stepwise.  Summarize findings.

#---FUNCTION-------------------------
#Custom function to calculate mallow's cp.  (Using normalized version.)
myCp <- function(X,y){
	model = myLinearReg(X,y)
	n = model$n		#Number of obs.
	p = model$p		#Number of predictors.
	
	RSS = model$RSS	#RSS for fitted model.
	sig2 = model$sigma_hat^2	#Sigma^2_hat for fitted model.
	cp = (1/n) * (RSS + 2*p*sig2)	#Mallow's cp for fitted model.
	
	return(cp)
}
#---END FUNCTION---------------------

#---FUNCTION-------------------------
#Custom function for forward stepwise.
myFwdStepwise <- function(X,y){
	
	Xwhole = X				#X matrix with all values.
	Xnew = NULL				#Start with null model.  (myLinearReg adds intercept)
	keeps = character()		#Empty char vector to hold names of kept vars.
	
	#Loop through possible numbers of predictors to add.
	for (i in 1:ncol(X)){ 
		
		cp_xnew = myCp(Xnew,y)		#Benchmark cp for xnew model.
		cp = rep(0,ncol(Xwhole))	#Vector to hold CP values for remaining predictors.
		vars = names(Xwhole)		#Vector of variable names still remaining to add.
		
		
		#Loop through each column of remaining X variables.
		for (j in 1:ncol(Xwhole)){
			#Calculate cp for each potential new model.
			cp[j] = myCp(cbind(Xnew,Xwhole[,j]),y)
		}
		
		#If none of the new models have lower cp than the old model,
		#break and return the old model.
		if (min(cp) >= cp_xnew){
			break;
		}
			
		#Pick the Xwhole column with the loweset cp.
		#Remove it from Xwhole and add it to Xnew.
		moveVar = vars[which(cp==min(cp))]	#Pick the variable to move.
		Xnew = cbind(Xnew,Xwhole[,moveVar])
		Xwhole = Xwhole[,-which(cp==min(cp)),drop=F]	
		keeps = c(keeps,moveVar)	#Keep track of vars added to model.
	}
	
	colnames(Xnew) = keeps
	final_cp = myCp(Xnew,y)
	return(list(X=Xnew,vars=keeps,cp = final_cp))	
} #End fwd stepwise function.
#---END FUNCTION---------------------

#---FUNCTION-------------------------
#Custom function for backward stepwise.
myBackwardStepwise <- function(X,y){
	
	Xwhole = X				#X matrix with all values.
	
	#Loop through possible numbers of predictors to add.
	for (i in 1:ncol(X)){ 
		
		cp_xnew = myCp(Xwhole,y)	#Benchmark cp for xnew model.
		cp = rep(0,ncol(Xwhole))	#Vector to hold CP values for remaining predictors.
		vars = names(Xwhole)		#Vector of variable names still remaining to add.
		
		#Loop through each column of remaining X variables.
		for (j in 1:ncol(Xwhole)){
			#Calculate cp for each potential new model.
			cp[j] = myCp(cbind(Xwhole[,-j]),y)
		}
		
		#If none of the new models have lower cp than the old model,
		#break and return the old model.
		if (min(cp) >= cp_xnew){
			break;
		}
			
		#Pick the Xwhole column with the loweset cp to remove.
		#Remove it from Xwhole and add it to Xnew.
		removeVar = vars[which(cp==min(cp))]	#Pick the variable to move.

		Xwhole = Xwhole[,-which(cp==min(cp)),drop=F]	
	}
	
	final_cp = myCp(Xwhole,y)
	return(list(X=Xwhole,vars=colnames(Xwhole),cp = final_cp))
		
} #End backward stepwise function.
#---END FUNCTION---------------------

#(i) Forward stepwise:
myFwd  = myFwdStepwise(X,y)
myFwd$vars
myFwd$cp

#Can compare results with the leaps package.  (Note - they use different Cp version.)
regfit.fwd = regsubsets(MPG ~ VOL + HP + SP + WT,data=data,method="forward")
#Leaps package comparison:
cbind(summary(regfit.fwd)$which, Cp=summary(regfit.fwd)$cp)

#(ii) Backward stepwise:
myBck  = myBackwardStepwise(X,y)
myBck$vars
myBck$cp

#Can compare results with the leaps package.  (Note - they use different Cp version.)
regfit.bwd = regsubsets(MPG ~ VOL + HP + SP + WT,data=data,method="backward")
#Leaps package comparison:
cbind(summary(regfit.bwd)$which,Cp=summary(regfit.bwd)$cp)

#-------------------------------------------------------------------
#Part C: Zheng-Loh Model Selection.

#Fit full model and obtain z-values (t-values, in this case) for each covariate.
mylm <- lm(MPG ~ VOL + HP + SP + WT,data=data,model=T)

#Save t-values in decreasing abs value order for the coefficients, excluding the intercept:
tvals <- summary(mylm)$coefficients[-1,3] 
tvals_ordered <- sort(abs(tvals),decreasing=T)

var_added <- names(tvals_ordered)
nmodels <- length(var_added)

#Set up data frame for output:
output <- data.frame(var_added, j=1:nmodels,rss=rep(0,nmodels),sigma_hat=rep(0,nmodels),jstar= rep(0,nmodels))

#Loop through model sizes, calculating the j's, adding one variable at a time.
for (j in 1:length(tvals_ordered)){
	
	#Set up the formula using only the first j variables.
	vars_kept <- noquote(paste(var_added[1:j],collapse="+"))
	formula = as.formula(paste("MPG ~ ",vars_kept,sep=""))
	
	#Fit a linear model using these j variables.
	temp_lm <- lm(formula,data=data)
	
	#Extract RSS and sigma_hat from the linear model.
	rss = output$rss[j] = sum(temp_lm$residuals^2)
	sigma = output$sigma_hat[j] <- summary(temp_lm)$sigma
	
	#Calculate jstar value.  Model with lowest jstar value is best.
	output$jstar[j] = rss + j * sigma^2 * log(nrow(data))
}

output #Display output
plot(1:4,output$jstar,type='l',col='blue',main='Zheng-Loh Model Selection',
	xlab='Model complexity (# predictors)', ylab='j.star')
	

###################################
###        PROBLEM 2:           ###
###################################

#Load data set.
pdata <- read.table(file='prostate_data.txt')

#Split data set into test and train.  
#Only train data set will be used for choosing the best model size and best model.
pdata_train <- pdata[pdata$train==T,-10]
pdata_test <- pdata[pdata$train==F,-10]

#-------------------------------------------------------------------
#Part A: Carry out a best-subset linear regression analysis.  Compute AIC, BIC, five and
# ten-fold cross-validation estimates of prediction error.  (Write your own code for cross-val.)

#---------------------	
#Choosing subset size based on AIC and BIC.
#Shows best model of each predictor size, including which predictors in that model.
regfit.bestsubsets <- regsubsets(lpsa ~ lcavol + lweight + age + 
	lbph + svi + lcp + gleason + pgg45, 
	data=pdata_train, method='exhaustive')

#---------------------	
#Cross-validation to pick best subset size using 10-fold CV.

nfolds <- 10

#Initialize output to hold prediction error and its SE for each k # of predictors.
output1 <- data.frame(k=1:8,pred_err=rep(0,8),se=rep(0,8))

#Loop through predictor sizes and call cross-val function.
for(k in 1:8){
	cv = myCrossVal(pdata_train,k,cv_folds=nfolds)
	output1$pred_err[k] = cv$mean_pred_test_err
	output1$se[k] = sqrt(cv$var_pred_test_err)/sqrt(nrow(pdata_train)/nfolds)
}

#Display cross-val results to choose best model size.
output1

#Save 10-fold error for output:
cvPredErr10 = output1$pred_err

#---------------------		
#Cross-validation to pick best subset size using 5-fold CV.

nfolds <- 5

#Initialize output to hold prediction error and its SE for each k # of predictors.
output2 <- data.frame(k=1:8,pred_err=rep(0,8),se=rep(0,8))

#Loop through predictor sizes and call cross-val function.
for(k in 1:8){
	cv = myCrossVal(pdata_train,k,cv_folds=nfolds)
	output2$pred_err[k] = cv$mean_pred_test_err
	output2$se[k] = sqrt(cv$var_pred_test_err)/sqrt(nrow(pdata_train)/nfolds)
}

#Display cross-val results to choose best model size.
output2

#Save 5-fold error for output:
cvPredErr5 = output2$pred_err

#---------------------
#Output all results from above.

cbind( 	summary(regfit.bestsubsets)$which,
		rss=summary(regfit.bestsubsets)$rss,
		aic=summary(regfit.bestsubsets)$cp,
		bic=summary(regfit.bestsubsets)$bic,
		cvPredErr10,
		cvPredErr5)

#---------------------
#Output cross-validation plots.	
#Plot predicted test error for each k size, with error bars.	

#10-fold CV plot:

plot(output1$k, output1$pred_err, type="n",ylim=c(.15,1.2),main='CV for 10-fold Cross-Validation',
	xlab='Number of predictors(k)',ylab='Mean Prediction Error (red line = 1 StDev Rule)')
with (
  data = output1
  , expr = errbar(k, pred_err, pred_err+se, pred_err-se, add=T, pch=1, cap=.1)
)

#Add dotted blue line to show 1-standard-deviation rule to pick most parsimonious model.
abline(h=output[which(output$pred_err==min(output$pred_err)),]$pred_err+output[which(output$pred_err==min(output$pred_err)),]$se,col='red')

	

#5-fold CV plot:
plot(output2$k, output2$pred_err, type="n",ylim=c(.15,1.2),main='CV for 10-fold Cross-Validation',
	xlab='Number of predictors(k)',ylab='Mean Prediction Error (red line = 1 StDev Rule)')
with (
  data = output2
  , expr = errbar(k, pred_err, pred_err+se, pred_err-se, add=T, pch=1, cap=.1)
)

#Add dotted blue line to show 1-standard-deviation rule to pick most parsimonious model.
abline(h=output[which(output$pred_err==min(output$pred_err)),]$pred_err+output[which(output$pred_err==min(output$pred_err)),]$se,col='red')


#---------------------
#After using CV to select most parsimonious model size, do best-subset to fit the appropriate model(s).
#The most parsimonious model has k=2 predictors.
#Now do best subsets on 2 predictors and display the best model:

best.cv.model <- regsubsets(lpsa ~ lcavol + lweight + age + lbph + svi + lcp + gleason + pgg45, 
	data=pdata_train, method='exhaustive',nvmax=2)
coefi <- coef(best.cv.model,id=2)	
coefi

#---FUNCTION-------------------------
#Cross-validation Function: Calculates predicted test error for a user-input number of predictors k.
myCrossVal <- function(data,k,cv_folds=10){
	#data = holds predictors and response.
	#k = number of parameters to include in model.
	#folds = number of cross-val folds.
	
	#Randomly shuffle data.
	data = data[sample(nrow(data)),]
	
	#Create 'folds' number of equally sized folds.
	folds <- cut(seq(1,nrow(data)),breaks=cv_folds,labels=F)
	
	#Initialize vector to hold prediction error for each cv fold iteration.
	pred_test_error <- rep(0,cv_folds)
	
	#Perform cross-validation.
	for (i in 1:cv_folds){
		
		#Split up data using folds.
		testIndices <- which(folds==i,arr.ind=T)
		testData <- data[testIndices, ]
		trainData <- data[-testIndices, ]
		
		#Do best subset selection for that number of predictors.  New model for each fold.
		best.fit.k <- regsubsets(lpsa~.,data=trainData,method="exhaustive",nvmax=k)
		
		#Compute test set error.
		test.mat <- model.matrix(lpsa~.,data=testData) #Set up model matrix.
		coefi <- coef(best.fit.k,id=k) #Extract coefficients from fitted model with k predictors.
		pred <- test.mat[,names(coefi)] %*% coefi #Calculate predicted values on test data. yhat=X*betahat.
		pred_test_error[i] <- mean((testData$lpsa-pred)^2)		
	}
	
	#Return average predicted test error for each k value.
	return(list(mean_pred_test_err=mean(pred_test_error),var_pred_test_err=var(pred_test_error)))
}
#---END FUNCTION---------------------