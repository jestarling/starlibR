% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/myBoxCox.R
\name{myBoxCox}
\alias{myBoxCox}
\title{myBoxCox(x) A Box-Cox transformation function.}
\usage{
myBoxCox(x)
}
\arguments{
\item{x}{A vector of data to be transformed to normality.}
}
\value{
x A list containing two objects.

x$theta The theta value to use for the box-cox transformation.

x$directions A text string reminder for how to perform the transformation using theta.
}
\description{
This function calculates the theta value for the y^theta box-cox transformation
to obtain normally distributed transformed data.
}
\examples{
x_lognorm <- rlnorm(100, meanlog = 0, sdlog = 1)
x_norm <- myBoxCox(x)

}
\author{
Jennifer Starling
}
\keyword{box-cox}
\keyword{normal}
\keyword{transformation}

