# Demo: Set up a new project using starUtils.

# First, create the project directory and create an R Project within that directory.  
#Then, open the project, and execute the following code.

# Install package.
library(devtools)
install_github('jestarling/starUtils')

# Load package.
library(starUtils)

# Verify root is correct.
find_root(is_rstudio_project)

# Create project directories.
setUpProject("JenniferStarling")

# Verify it worked.
dir()
