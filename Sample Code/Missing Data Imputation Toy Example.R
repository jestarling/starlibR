#Missing data imputation toy example.

library(mice)

#-------------------------------------------------

#DATA SETUP & EXPLORATION:

data <- airquality
data[4:10,3] <- rep(NA,7)
data[1:5,4] <- NA

#Show the number of NA data points for each predictor variable.
summary(data)

#There are two types of missing data:
#MCAR: missing completely at random. This is the desirable scenario in case of missing data.
#MNAR: missing not at random. Missing not at random data is a more serious issue and in this case it might be wise to check the data gathering process further and try to understand why the information is missing. For instance, if most of the people in a survey did not answer a certain question, why did they do that? Was the question unclear?
#Assuming data is MCAR, too much missing data can be a problem too. Usually a safe maximum threshold is 5% of the total for large datasets. If missing #data for a certain feature or sample is more than 5% then you probably should leave that feature or sample out. We therefore check for features (columns) #and samples (rows) where more than 5% of the data is missing using a simple function

pMiss <- function(x){sum(is.na(x))/length(x)}
apply(data,2,pMiss)
apply(data,1,pMiss)

#Can see if there is any pattern in the missing data, or if it is missing at random.
md.pattern(data)

#Better visual representation of missing data using VIM package.
library(VIM)
aggr_plot <- aggr(data, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(data), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))

#The plot helps us understanding that almost 70% of the samples are not missing any information, 22% are missing the Ozone value, and the remaining ones #show other missing patterns. Through this approach the situation looks a bit clearer in my opinion.
#Another (hopefully) helpful visual approach is a special box plot:

marginplot(data[c(1,2)])

#Obviously here we are constrained at plotting 2 variables at a time only, but nevertheless we can gather some interesting insights.
#The red box plot on the left shows the distribution of Solar.R with Ozone missing while the blue box plot shows the distribution of #the remaining datapoints. Likewhise for the Ozone box plots at the bottom of the graph.
#If our assumption of MCAR data is correct, then we expect the red and blue box plots to be very similar.

#-------------------------------------------------

#PERFORMING THE IMPUTATION
tempData <- mice(data,m=5,maxit=50,meth='pmm',seed=500)
summary(tempData)

	?mice #shows you the available methods of imputation.
	# pmm = predictive mean matching.
	# norm = bayesian linear regression
	# norm.nob = linear regression ignoring model error
	# norm.boot = linear regression using bootstrap
	# mean = linear regression, predicted values.
	# l2.norm = unconditional mean imputation.
	# logreg = logistic regression (factor, 2-levels)
	# lda = linear discriminant analysis (factor, categories)
	# cart  = classification/regression trees.
	# rf = random forest
	# ... and others.

#NOTE: m=5 refers to the number of imputed datasets. Five is the default value.

#RETURN COMPLETED DATA SET WITH IMPUTED VALUES.
completedData <- complete(tempData,1)

#The missing values have been replaced with the imputed values in the first of the five datasets. 
#If you wish to use another one, just change the second parameter in the complete() function.

#-------------------------------------------------

#SANITY CHECK IMPUTED DATA:
#If you would like to check the imputed data, for instance for the variable Ozone, you need to enter the following line of code
tempData$imp$Ozone

#The output shows the imputed data for each observation (first column left) within each imputed dataset (first row at the top).

#If you need to check the imputation method used for each variable, mice makes it very easy to do
tempData$meth

#-------------------------------------------------

#EXPLORING RESULTS OF IMPUTATION FOR VALIDITY

#Let’s compare the distributions of original and imputed data using a some useful plots.
#First of all we can use a scatterplot and plot Ozone against all the other variables
xyplot(tempData,Ozone ~ Wind+Temp+Solar.R,pch=18,cex=1)

	#What we would like to see is that the shape of the magenta points (imputed) matches the 
	#shape of the blue ones (observed). The matching shape tells us that the imputed 
	#values are indeed “plausible values”.

#Another helpful plot is the density plot:

densityplot(tempData)

	#The density of the imputed data for each imputed dataset is showed in magenta 
	#while the density of the observed data is showed in blue. Again, under our 
	#previous assumptions we expect the distributions to be similar.

#Another useful visual take on the distributions can be obtained using the stripplot() function that shows the distributions of the variables as individual points

stripplot(tempData, pch = 20, cex = 1.2)

#-------------------------------------------------

#POOLING MODELS FROM EACH IMPUTED DATA SET

#Suppose that the next step in our analysis is to fit a linear model to the data. 
#You may ask what imputed dataset to choose. 
#The mice package makes it again very easy to fit a a model to each of 
#the imputed dataset and then pool the results together

modelFit1 <- with(tempData,lm(Temp~ Ozone+Solar.R+Wind))
summary(pool(modelFit1))
