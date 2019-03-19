\name{fit_gamma}
\alias{fit_gamma}
\title{Fit a gamma distribution to observed data}
\description{
	Fits a gamma distribution to observed data.
}
\usage{
	fit_gamma(x)
}
\arguments{
  \item{x}{a vector of observed values.}
}
\details{
	Uses the sample mean (m) and variance (v) to calculate shape and scale parameters of a gamma distribution as follows:
	scale = v/m
	shape = m^2/v
	
	Useful for deciding appropriate parameters for simulations.
}
\value{
	Named list of parameters.
}
\examples{
data <- rgamma(1000, shape=1, scale=10);
params <- fit_gamma(data)
round(unlist(params)) == c(1,10) # TRUE TRUE
}
\keyword{gamma, simulation}
