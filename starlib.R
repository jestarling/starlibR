#This is the code to load up the starlib package.

library(devtools)

setwd("/Users/jennstarling/TAMU/starlib")
install("/Users/jennstarling/TAMU/starlib") 

# Github varstar package install.
#library(devtools)
#install_github(‘jstarling1/starlib’,’jstarling1’)

library(starlib)		#Load custom starlib package for all custom functions.
data(package='starlib')	#View available starlib data sets.
ls("package:starlib")	#View all functions in starlib package.
