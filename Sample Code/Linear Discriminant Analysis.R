
setwd('/Users/jennstarling/UTAustin/starlib/data')
library(MASS)		#For use of lda() function.

#Set up colors to use for scatter plots.
cl = c('black','dark red','red','orange','yellow','green','light blue','blue','pink','purple','grey')

#Read in data.
train = read.table(file='vowel.train.txt', header=T,sep=',')
test = read.table(file='vowel.test.txt', header=T,sep=',')

#Get rid of superfluous first 'row name' columns.
train = train[,-1]
test = test[,-1]

#Scale and center data.
train = as.data.frame(cbind(train[,1], scale(train[,-1],scale=F))) #Center data only
test = as.data.frame(cbind(test[,1], scale(test[,-1],scale=F)))	#Center data only.

#Fix column headers (getting erased during scaling).
colnames(train)[1] = colnames(test)[1] = 'y'	

#Extract X matrices for calculating projections.
X.train = as.matrix(train[,-1])
X.test = as.matrix(test[,-1])

#Perform LDA.
mylda = lda(y~., data=train)

#Save projection matrix V.
V = mylda$scaling

#Calculate projected data.
projected.data = X.train %*% V

#Recreate plot 4.8 from ELS.  
#(Scatter plots of various pairs of projections.)
par(mfrow=c(2,2), oma=c(0,0,2,0))
plot(projected.data[,1], projected.data[,3], col=cl[train[,1]], xlab='Coordinate 1', ylab='Coordinate 3',cex=.75)
plot(projected.data[,2], projected.data[,3], col=cl[train[,1]], xlab='Coordinate 2', ylab='Coordinate 3',cex=.75)
plot(projected.data[,1], projected.data[,7], col=cl[train[,1]], xlab='Coordinate 1', ylab='Coordinate 7',cex=.75)
plot(projected.data[,9], projected.data[,10], col=cl[train[,1]], xlab='Coordinate 9', ylab='Coordinate 10',cex=.75)
title("Linear Discriminant Analysis", outer=TRUE)

#------------------------------------------------------------------------
#Recreate plot from 4.11 from ELS.
#(Scatter plot of projection with coords 1 and 2, with decision boundaries plotted.)
par(mfrow=c(1,1))
plot(projected.data[,1], projected.data[,2], col=cl[train[,1]], xlab='Canonical Coordinate 1', ylab='Canonical Coordinate 2',cex=.75)

#Per Purna, do not worry about adding decision boundaries.  

###########################################
###  Extra Credit: Reproduce Plot 4.10  ###
###########################################

dims = 1:10	#Vector to hold dimensions.
test.err = rep(0,10)	#Vector to hold training error for each dim.
train.err = rep(0,10)	#Vector to hold test error for each dim.

for (d in 1:10){	#Loop through dimensions.
	
	#-------------------------
	#Training Error:
	mlda = lda(y~.,data=train)
	lda.pred = predict(mylda,newdata=train,dimen=d)
	train.err[d] = mean(lda.pred$class!=train$y)
		
	#-------------------------
	#Calc test error.
	mylda = lda(y~., data=train)	#Fit lda model on train set.
	lda.pred = predict(mylda,newdata=test,dimen=d)
	test.err[d] = mean(lda.pred$class!=test$y)	
}

#Plot results:
plot(dims,test.err,type='l',col='red',xlab='Dimension',ylab='Misclassification Error',
	main='Figure 4.10 ESL',ylim=c(.2,1))
points(dims,train.err,type='l',col='blue')
legend(6,.75,c('Train Error','Test Error'),lty=c(1,1), col=c('blue','red'))

