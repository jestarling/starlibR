library(MASS)
library(neuralnet)

#######################
###  THE DATA SET   ###
#######################

#We are going to use the Boston dataset in the MASS package.
#The Boston dataset is a collection of data about housing values in the suburbs of Boston. 
#Our goal is to predict the median value of owner-#occupied homes (medv) 
#using all the other continuous variables available.

set.seed(500)
data <- Boston

#First we need to check that no datapoint is missing, otherwise we need to fix the dataset.
#There is no missing data.
apply(data,2,function(x) sum(is.na(x)))

#--------------------------------
#Split data randomly into a test and training set.
#First a linear regression model, and test it on the test set.  (using gml instead of lm, as is useful for cross-val.)

index <- sample(1:nrow(data),round(0.75*nrow(data)))
train <- data[index,]
test <- data[-index,]
lm.fit <- glm(medv~., data=train)
summary(lm.fit)
pr.lm <- predict(lm.fit,test)
MSE.lm <- sum((pr.lm - test$medv)^2)/nrow(test)

#Display Mean Squared Error for model.
MSE.lm  

#--------------------------------
#NEURAL NETWORK PREP STEPS (scaling/centering data)
#Before fitting a neural network, some preparation need to be done. Neural networks are not that easy to train and tune.
#As a first step, we are going to address data preprocessing.

#It is good practice to normalize your data before training a neural network. I cannot emphasize enough how important this step #is: depending on your dataset, avoiding normalization may lead to useless results or to a very difficult training process (most #of the times the algorithm will not converge before the number of maximum iterations allowed). You can choose different methods #to scale the data (z-normalization, min-max scale, etc…). I chose to use the min-max method and scale the data in the interval #[0,1]. Usually scaling in the intervals [0,1] or [-1,1] tends to give better results.
#We therefore scale and split the data before moving on:

maxs <- apply(data, 2, max) 
mins <- apply(data, 2, min)

scaled <- as.data.frame(scale(data, center = mins, scale = maxs - mins))

train_ <- scaled[index,]
test_ <- scaled[-index,]

#--------------------------------
#FITTING THE NEURAL NETWORK
#As far as I know there is no fixed rule as to how many layers and neurons to use although there are several more or less accepted rules of #thumb. Usually, if at all necessary, one hidden layer is enough for a vast numbers of applications. As far as the number of neurons is #concerned, it should be between the input layer size and the output layer size, usually 2/3 of the input size. At least in my brief experience #testing again and again is the best solution since there is no guarantee that any of these rules will fit your model best.

#Since this is a toy example, we are going to use 2 hidden layers with this configuration: 13:5:3:1. The input layer has 13 inputs, the two #hidden layers have 5 and 3 neurons and the output layer has, of course, a single output since we are doing regression.

#Let’s fit the net:
library(neuralnet)

#Set up formula statement
n <- names(train_)
f <- as.formula(paste("medv ~", paste(n[!n %in% "medv"], collapse = " + ")))

#Fit the net.
nn <- neuralnet(f,data=train_,hidden=c(5,3),linear.output=T)

#Some notes on fitting the network.
#1. For some reason the formula y~. is not accepted in the neuralnet() function. You need to first write the formula and then pass it as an argument in the fitting function.

#2. The hidden argument accepts a vector with the number of neurons for each hidden layer, while the argument linear.output is used to specify whether we want to do regression linear.output=TRUE or classification linear.output=FALSE

#--------------------------------
#PLOTTING THE RESULTS:

plot(nn)

#This is the graphical representation of the model with the weights on each connection:

#The black lines show the connections between each layer and the weights on 
#each connection while the blue lines show the bias term added in each step. 
#The bias can be thought as the intercept of a linear model.
#The net is essentially a black box so we cannot say that much about the fitting, 
#the weights and the model. Suffice to say that the training algorithm has converged and 
#therefore the model is ready to be used.

#--------------------------------
#PREDICTING MEDV USING THE NEURAL NETWORK:

#Now we can try to predict the values for the test set and calculate the MSE. 
#Remember that the net will output a normalized prediction, 
#so we need to scale it back in order to make a meaningful comparison (or just a simple prediction).

pr.nn <- compute(nn,test_[,1:13]) 

pr.nn_ <- pr.nn$net.result*(max(data$medv)-min(data$medv))+min(data$medv)

test.r <- (test_$medv)*(max(data$medv)-min(data$medv))+min(data$medv)

#Calculate MSE for neural network.
MSE.nn <- sum((test.r - pr.nn_)^2)/nrow(test_)

#Compare the two MSEs
print(paste(MSE.lm,MSE.nn))

#Apparently the net is doing a better work than the linear #model at predicting medv. Once again, be careful because this #result depends on the train-test split performed above. Below, #after the visual plot, we are going to perform a fast cross #validation in order to be more confident about the results.

#A first visual approach to the performance of the network and the linear model on the test set is plotted below.

par(mfrow=c(1,3))

plot(test$medv,pr.nn_,col='red',main='Real vs predicted NN',pch=18,cex=0.7)
abline(0,1,lwd=2)
legend('bottomright',legend='NN',pch=18,col='red', bty='n')

plot(test$medv,pr.lm,col='blue',main='Real vs predicted lm',pch=18, cex=0.7)
abline(0,1,lwd=2)
legend('bottomright',legend='LM',pch=18,col='blue', bty='n', cex=.95)

plot(test$medv,pr.nn_,col='red',main='Real vs predicted NN',pch=18,cex=0.7)
points(test$medv,pr.lm,col='blue',pch=18,cex=0.7)
abline(0,1,lwd=2)
legend('bottomright',legend=c('NN','LM'),pch=18,col=c('red','blue'))

#By visually inspecting the plot we can see that the predictions made by the 
#neural network are (in general) more concetrated around the line 
#(a perfect alignment with the line would indicate a MSE of 0 and thus an 
#ideal perfect prediction) than those made by the linear model.

#--------------------------------
# A (FAST) CROSS-VALIDATION:

#While there are different kind of cross validation methods, the basic idea 
#is repeating the following process a number of time:
#train-test split
#		Do the train-test split
#		Fit the model to the train set
#		Test the model on the test set
#		Calculate the prediction error
#		Repeat the process K times
#Then by calculating the average error we can get a grasp of how the model is doing.
#We are going to implement a fast cross validation using a for loop for the neural network and the cv.glm() function in the boot #package for the linear model.
#As far as I know, there is no built-in function in R to perform cross validation on this kind of neural network, if you do know #such a function, please let me know in the comments. 

#Here is the 10 fold cross validated MSE for the linear model:

library(boot)
set.seed(200)
lm.fit <- glm(medv~.,data=data)
cv.glm(data,lm.fit,K=10)$delta[1]

#Now 10-fold CV for the net. 
#Note that I am splitting the data in this way: 90% train set and 10% test set in a random way for 10 times. 
#I am also initializing a progress bar using the plyr library because 
#I want to keep an eye on the status of the process since the fitting of the neural network may take a while.

set.seed(450)
cv.error <- NULL
k <- 10

library(plyr) 
pbar <- create_progress_bar('text')
pbar$init(k)

for(i in 1:k){
    index <- sample(1:nrow(data),round(0.9*nrow(data)))
    train.cv <- scaled[index,]
    test.cv <- scaled[-index,]
    
    nn <- neuralnet(f,data=train.cv,hidden=c(5,2),linear.output=T)
    
    pr.nn <- compute(nn,test.cv[,1:13])
    pr.nn <- pr.nn$net.result*(max(data$medv)-min(data$medv))+min(data$medv)
    
    test.cv.r <- (test.cv$medv)*(max(data$medv)-min(data$medv))+min(data$medv)
    
    cv.error[i] <- sum((test.cv.r - pr.nn)^2)/nrow(test.cv)
    
    pbar$step()
}

#After a while, the process is done, we calculate the average MSE and plot the results as a boxplot
mean(cv.error)
cv.error

#The code for the box plot of the CV error:
boxplot(cv.error,xlab='MSE CV',col='cyan',
        border='blue',names='CV error (MSE)',
        main='CV error (MSE) for NN',horizontal=TRUE)

#-------------------------------------
#FINAL NOTE ON INTERPRETABILITY:

#Neural networks resemble black boxes a lot: explaining their outcome is 
#much more difficult than explaining the outcome of simpler model such as a linear model. 
#Therefore, depending on the kind of application you need, you might want to take 
#into account this factor too. Furthermore, as you have seen above, extra care is 
#needed to fit a neural network and small changes can lead to different results.

