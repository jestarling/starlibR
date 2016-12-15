#Summary statistic/exploratory data stuff in R.

#Working with iris data.
data = iris

#---------------------------------------
#SUMMARY OF DATA:
class(iris)

dim(data)	#Dimensions

names(data)	#Column names.

str(data) 	#Displays structure.

attributes(data)	#Displays col names, row names, and class.

head(data)
tail(data,10)

summary(data) #5-num spread for continuous, frequency of factor levels

head(data)
tail(data,10)

summary(data)

library(Hmisc)
describe(data[,1:5])

#---------------------------------------
# EXPLORING A CONTINUOUS VARIABLE

variable = data$Sepal.Length

range(variable)
quantile(variable)
quantile(variable,c(0,.3,.6,.9,1)) #Specific quantiles
var(variable)	#Variance
hist(variable)	#Histogram

#Normal QQ plot.
qqnorm(variable)
qqline(variable)

#Density plot.
plot(density(variable))

#Histogram with density estimator.
hist(variable,freq=F)
lines(density(variable),col='blue')

#Scatter plot.
library(calibrate)
names = rownames(data)
idx = 1:nrow(data)

plot(variable)
textxy(idx,variable, labs=names, cex=.7,col='blue')


#---------------------------------------
# EXPLORING A CATEGORICAL VARIABLE
variable = iris$Species

#Frequency table:
table(variable)

#Pie chart
pie(table(variable))

#Bar chart
barplot(table(variable))

#---------------------------------------
# EXPLORING PREDICTORS TOGETHER (CORRELATION STRUCTURE, ETC)

#Cov & corr matrices for first 4 (continuous) predictors.
cov(data[,1:4]) 
cor(data[,1:4])

#Aggregate stats of one variable by another.
aggregate(Sepal.Length ~ Species, summary, data=data)

#Side by side BoxPlots
