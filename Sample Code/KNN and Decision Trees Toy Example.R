
###################################
###   	KNN				###
###################################

library(class)	#Contains knn() function.

library(dplyr)
library(lubridate)
set.seed(100)

setwd('/Users/jennstarling/Desktop/Integra')
stocks <- read.csv('stocks.csv')


#----------------------
#DATA PREP:

#Put all data before the year 2014 into the training set, and the rest into the test set.
stocks$Date <- ymd(stocks$Date)
stocksTrain <- year(stocks$Date) < 2014

#Now, we need to build the training set. It will consist of the prices of stocks of Apple, Google, and Microsoft on the previous day. For this, we can use the lag function in dplyr.

predictors <- cbind(lag(stocks$Apple, default = 210.73), lag(stocks$Google, default = 619.98), lag(stocks$MSFT, default = 30.48))

#Since for the very first value (corresponding to January 4, 2010), the lag function has nothing to compare it to, it will default to NA. To avoid this, I set the default values for each stock to be its value on the previous business day (December 31, 2009).

#----------------------
#FIT THE KNN MODEL AND APPLY TO TRAINING DATA:

prediction <- knn(predictors[stocksTrain, ], predictors[!stocksTrain, ], stocks	
	$Increase[stocksTrain], k = 1)

#We can see it’s accuracy using table.
table(prediction, stocks$Increase[!stocksTrain])

summary(prediction) #Shows number of yes/no classifications.

#Measure accuracy of classification.
mean(prediction == stocks$Increase[!stocksTrain])

#----------------------
#VARY K VALUES TO ASSESS ACCURACY

#We can use a for loop to see how the algorithm performs for different values of k.
accuracy <- rep(0, 10)
k <- 1:10
for(x in k){
  prediction <- knn(predictors[stocksTrain, ], predictors[!stocksTrain, ],
                    stocks$Increase[stocksTrain], k = x)
  accuracy[x] <- mean(prediction == stocks$Increase[!stocksTrain])
}

plot(k, accuracy, type = 'b')

###################################
###   	RANDOM FORESTS		###
###################################


## randomForest
library(randomForest)
fit.rf = randomForest(frmla, data=raw)
print(fit.rf)
importance(fit.rf)
plot(fit.rf)
plot( importance(fit.rf), lty=2, pch=16)
lines(importance(fit.rf))

imp = importance(fit.rf)
impvar = rownames(imp)[order(imp[, 1], decreasing=TRUE)]
op = par(mfrow=c(1, 3))

for (i in seq_along(impvar)) {
partialPlot(fit.rf, raw, impvar[i], xlab=impvar[i],
main=paste("Partial Dependence on", impvar[i]),
ylim=c(0, 1))
}

###################################
###   	CLASSIFICATION TREE	###
###################################

# Classification Tree with rpart library(rpart)  # grow tree  fit <- rpart(Kyphosis ~ Age + Number + Start,    method="class", data=kyphosis)  printcp(fit) # display the results  plotcp(fit) # visualize cross-validation results  summary(fit) # detailed summary of splits  # plot tree  plot(fit, uniform=TRUE,     main="Classification Tree for Kyphosis") text(fit, use.n=TRUE, all=TRUE, cex=.8)  # create attractive postscript plot of tree  post(fit, file = "c:/tree.ps",     title = "Classification Tree for Kyphosis")
   click to view
# prune the tree  pfit<- prune(fit, cp=   fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"])  # plot the pruned tree  plot(pfit, uniform=TRUE,     main="Pruned Classification Tree for Kyphosis") text(pfit, use.n=TRUE, all=TRUE, cex=.8) post(pfit, file = "c:/ptree.ps",     title = "Pruned Classification Tree for Kyphosis")

###################################
###   	REGRESSION TREE		###
###################################

# Regression Tree Example library(rpart)  # grow tree  fit <- rpart(Mileage~Price + Country + Reliability + Type,     method="anova", data=cu.summary)  printcp(fit) # display the results  plotcp(fit) # visualize cross-validation results  summary(fit) # detailed summary of splits  # create additional plots  par(mfrow=c(1,2)) # two plots on one page  rsq.rpart(fit) # visualize cross-validation results     # plot tree  plot(fit, uniform=TRUE,     main="Regression Tree for Mileage ") text(fit, use.n=TRUE, all=TRUE, cex=.8)  # create attractive postcript plot of tree  post(fit, file = "c:/tree2.ps",     title = "Regression Tree for Mileage ")
   click to view
# prune the tree  pfit<- prune(fit, cp=0.01160389) # from cptable     # plot the pruned tree  plot(pfit, uniform=TRUE,     main="Pruned Regression Tree for Mileage") text(pfit, use.n=TRUE, all=TRUE, cex=.8) post(pfit, file = "c:/ptree2.ps",     title = "Pruned Regression Tree for Mileage")