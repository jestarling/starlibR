#This is the code to load up the starlib package from my local machine.

library(devtools)

setwd("/Users/jennstarling/UT Austin/starlib")
install("/Users/jennstarling/UT Austin/starlib") 

# Github varstar package install.
#library(devtools)
#install_github(‘jstarling1/starlib’,’jstarling1’)

library(starlib)		#Load custom starlib package for all custom functions.
data(package='starlib')	#View available starlib data sets.
ls("package:starlib")	#View all functions in starlib package.

#This is the code to install the starlib package from github.
library(devtools)
install_github(repo="jstarling1/starlib", username = NULL, ref = "master", 
    subdir = NULL, auth_token = github_pat(quiet), host = "api.github.com", 
    force = FALSE, quiet = FALSE)

#Update package documentation in R (roxygen):
document()