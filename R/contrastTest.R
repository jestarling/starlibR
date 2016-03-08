#CONTRASTS TESTING

#Housekeeping

contrastTest <- function(data=NULL){
	
	library(car)
	library(MASS)
	set.seed(1)

	#If a data set is not provided, load and use the contrast_data set in starlib.
	if(is.null(data)){
		data(contrast_data,package='starlib')
		data <- contrast_data
	}

	x <- factor(data$x)
	
	contr <- list() #Value to hold output models.

#NOTE: Original data source is:
#data <- read.csv("/Users/jennstarling/TAMU/STAT 645 (Biostats)/Data Sets/contrast_test.csv",header=T)

	#Investigate means of y:
	meany <- mean(y) #Mean of overall y values
	meanA <- mean(y[which(x=='A')]) #Group mean for x=A
	meanB <- mean(y[which(x=='B')]) #Group mean for x=B
	meanC <- mean(y[which(x=='C')]) #Group mean for x=C
	grand_mean <- mean(c(meanA,meanB,meanC))
	
	contr$means <- data.frame(meany=meany,meanA=meanA,meanB=meanB,meanC=meanC,grandMean=grand_mean)

#Insert some blank space for fomatting:
writeLines('
	')
	
print(paste('Mean of vector of y values: ',meany,sep=""))
print(paste('Mean of y values where x=A ',meanA,sep=""))
print(paste('Mean of y values where x=B ',meanB,sep=""))
print(paste('Mean of y values where x=C ',meanC,sep=""))
print(paste('Mean of grand mean (mean of meanA, meanB, meanC): ',grand_mean,sep=""))

#-----------------------------------------------------------------
# DEFAULT CONTRAST:

#Display default contrast for X:
print(paste('Default contasts for x:'), contrasts(x))
writeLines('
	')
contrasts(x)

	contr$default_contrast <- contrasts(x)

#Linear model, with and w/o intercept, with default contrast:
my.lm1 <- lm(y~x,data=data)
my.lm2 <- lm(y~x-1,data=data)

	contr$default_lm_withIntercept <- my.lm1
	contr$default_lm_noIntercept <- my.lm2

writeLines('
	Linear models, with (1) and without (2) an intercept:
	')
	
summary(my.lm1)
summary(my.lm2)

writeLines('
	DEFAULT CONTRAST LINEAR MODEL RESULTS:
	For lm1 (model with intercept), treats A as baseline category.
	Intercept = mean(x=A values) 
	xB Coeff = mean(X=B values) - mean(X=A values) 
	xC Coeff = mean(X=C values) - mean(X=A values) 
	
	For lm2, ie model without intercept (ie intercept = 0): 
	xA Coeff = mean(X=A values) 
	xB Coeff = mean(X=B values) 
	xC Coeff = mean(X=C values)')

#For lm1, ie model with intercept (treats A as baseline category):
#	Intercept = mean(X=A values)
#	xB Coeff = mean(X=B values) - mean(X=A values)
#	xC Coeff = mean(X=C values) - mean(X=A values)

#For lm2, ie model without intercept:
#	No intercept, ie Intercept = 0
#	xA Coeff = mean(X=A values)
#	xB Coeff = mean(X=B values)
#	xC Coeff = mean(X=C values)

#-----------------------------------------------------------------
# CONTR.SUM() CONTRAST:

#Set contrast for x to contr.sum()
contrasts(x) <- contr.sum(nlevels(x))
print(paste('Contrasts for x using contr.sum() :'),contrasts(x))
contrasts(x)

	contr$contrsum_contrast <- contrasts(x)

#Linear model, with and w/o intercept, with updated contrasts:
my.lm1 <- lm(y~x,data=data,contrasts=list(x=contr.sum(nlevels(x)))) #Syntax 1
my.lm1 <- lm(y~x,data=data,contrasts=list(x="contr.sum")) #Syntax 2

my.lm2 <- lm(y~x-1,data=data,contrasts=list(x='contr.sum"'))

	contr$contrsum_lm_withIntercept <- my.lm1
	contr$contrsum_lm_noIntercept <- my.lm2

 #Intercept = grand mean
 #meanA = Intercept - x1 coeff
 #meanB = Intercept - x2 coeff
 #meanC = Intercept - (x1 coeff + x2 coeff)
 
#NOTE: The grand mean = global mean ie mean(y) here, bc data is balanced.
#Balanced means that all groups have same # obs.
#If not balanced, then the intercept is not mean(y) anymore - it is the 
#mean of the group means, like an average of averages.
#This interpretation is a bit weird, but the coeffs still have similar meanings,
#and can still do LR and F-Tests and significance tests for coeffs.

print('
Linear models, with (1) and without (2) an intercept:
')
summary(my.lm1)
summary(my.lm2)

writeLines('
	SUM.CONTR() CONTRAST LINEAR MODEL RESULTS:
	Intercept = grand mean
	meanA = Intercept - x1 coeff
	meanB = Intercept - x2 coeff
	Intercept - (x1 coeff + x2 coeff)
	
	NOTE REGARDING BALANCED DATA AND THE GRAND MEAN: 
	The grand mean = global mean ie mean(y) for simulated data, since data is balanced.
	Balanced means that all groups have same # obs.
	If not balanced, then the intercept is not mean(y) anymore - it is the 
		mean of the group means, like an average of averages.
	This interpretation is a bit weird, but the coeffs still have similar meanings,
		and can still do LR and F-Tests and significance tests for coeffs.
	')

	return(contr)
}