#This code gives a demonstration of how to plot an estimated density funciton
#and shdae an area under the curve.

densityShadeDemo = function()
{
x <- rnorm(5000,0,1) #Generate standard normal observations
t <- 2	#Threshold value - shade all values >= than this value.

hist(x,prob=T)

#-------------------------------------------------------------------
#Freq histogram is default, where each bar is a count of number of times
#value is seen.
#prob=T scales the hist so that each bar corresponds to a probability, ie
#scale so total area under bars = 1.

#A density is basically a histogram with smoothed curve over it.
#Areas under the curve correspond to probabilities.

#density() gives the KDE (kernel density estimator) for x.
#This prints out info about the density.
dens <- density(x) #Fits a smooth curve to a distribution.
dens$x[1:10] #Shows first ten x coordinates for how to draw the density estimate curve.
dens$y[1:10] #Shows first ten y coordinates for how to draw the density estimate curve.

#Can add the density curve to the histogram plot:
hist(x,prob=T)
lines(dens,col='blue',lwd=2)  #This is the same as doing lines(dens$x,dens$y,col='blue',lwd=2)

#-------------------------------------------------------------------
#Now, plot the density and shade part of the area under it.
plot(density(x)) #or plot(d)

#Say we want to shade the areas where value >=t.
#There is a function called polygon we can use.  It makes lots of little rectangles and shades them.
#If you give it enough points, it will look smooth.  Have to specify width and height of each bin.

polygon( c(2, dens$x[dens$x>=t], max(dens$x)), #All x points that define region
	c(0, dens$y[dens$x>=2],0), #All y points that define region.
	col='red')

#-------------------------------------------------------------------
#This is another way to do this:
plot(dens,main = 'Density plot',
	xlab='density',ylab='scores') #Density plot of the B scores.

#Adens the grey shading for scores >= t.
x1 <- min(which(dens$x >= t))
x2<- max(which(dens$x <= max(dens$x)))
with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col="gray"))

}
