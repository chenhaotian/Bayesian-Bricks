% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Gaussian_Inference.r
\name{dPosteriorPredictive.GaussianInvWishart}
\alias{dPosteriorPredictive.GaussianInvWishart}
\title{Posterior predictive density function of a "GaussianInvWishart" object}
\usage{
\method{dPosteriorPredictive}{GaussianInvWishart}(obj, x, LOG = TRUE, ...)
}
\arguments{
\item{obj}{A "GaussianInvWishart" object.}

\item{x}{matrix, or the ones that can be converted to matrix, each row of x is an observation.}

\item{LOG}{Return the log density if set to "TRUE".}

\item{...}{Additional arguments to be passed to other inherited types.}
}
\value{
A numeric vector of the same length as nrow(x), the posterior predictive density.
}
\description{
Generate the the density value of the posterior predictive distribution of the following structure:
    \deqn{x \sim Gaussian(mu,Sigma)}
    \deqn{Sigma \sim InvWishart(v,S)}
mu is known. Gaussian() is the Gaussian distribution. See \code{?dGaussian} and \code{?dInvWishart} for the definition of the distributions.\cr
The model structure and prior parameters are stored in a "GaussianInvWishart" object. \cr
Posterior predictive density is p(x|v,S,mu).
}
\examples{
obj <- GaussianInvWishart(gamma=list(mu=c(-1.5,1.5),v=3,S=diag(2)))
x <- rGaussian(100,mu = c(-1.5,1.5),Sigma = matrix(c(0.1,0.03,0.03,0.1),2,2))
xNew <- rGaussian(100,mu = c(-1.5,1.5),Sigma = matrix(c(0.1,0.03,0.03,0.1),2,2))
ss <- sufficientStatistics(obj=obj,x=x,foreach = FALSE)
## update piror with x
posterior(obj=obj,ss = ss)
## use the posterior to calculate the probability of observing each xNew
dPosteriorPredictive(obj = obj,x = xNew,LOG = TRUE)
}
\references{
Gelman, Andrew, et al. Bayesian data analysis. CRC press, 2013.

MARolA, K. V., JT KBNT, and J. M. Bibly. Multivariate analysis. AcadeInic Press, Londres, 1979.
}
\seealso{
\code{\link{GaussianInvWishart}}, \code{\link{dPosteriorPredictive.GaussianInvWishart}}, \code{\link{marginalLikelihood.GaussianInvWishart}}
}
