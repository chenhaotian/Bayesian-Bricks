% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Gaussian_Inference.r
\name{GaussianGaussian}
\alias{GaussianGaussian}
\title{Create objects of type "GaussianGaussian".}
\usage{
GaussianGaussian(
  objCopy = NULL,
  ENV = parent.frame(),
  gamma = list(Sigma = 1, m = 0, S = 1)
)
}
\arguments{
\item{objCopy}{an object of type "GaussianGaussian". If "objCopy" is not NULL, the function create a new "GaussianGaussian" object by copying the content from objCopy, otherwise this new object will be created by using "ENV" and "gamma". Default NULL.}

\item{ENV}{environment, specify where the object will be created.}

\item{gamma}{list, a named list of parameters, gamma=list(Sigma,m,S). Where gamma$Sigma is the known covariance matrix of x, gamma$m and gamma$S are the prior mean and covariance matrix of mu.}
}
\value{
An object of class "GaussianGaussian".
}
\description{
Create an object of type "GaussianGaussian", which represents the Gaussian and Gaussian conjugate structure:
     \deqn{x \sim Gaussian(mu,Sigma)}
     \deqn{mu \sim Gaussian(m,S)}
Where Sigma is known. Gaussian() is the Gaussian distribution. See \code{?dGaussian} for the definition of Gaussian distribution.\cr
The created object will be used as a place for recording and accumulating information in the related inference/sampling functions such as posterior(), posteriorDiscard(), MAP(), marginalLikelihood(), dPosteriorPredictive(), rPosteriorPredictive() and so on.
}
\examples{
obj <- GaussianGaussian(gamma=list(Sigma=matrix(c(2,1,1,2),2,2),m=c(0.2,0.5),S=diag(2)))
obj #print the content
}
\references{
Gelman, Andrew, et al. Bayesian data analysis. CRC press, 2013.
}
\seealso{
\code{\link{posterior.GaussianGaussian}},\code{\link{posteriorDiscard.GaussianGaussian}},\code{\link{MAP.GaussianGaussian}},\code{\link{MPE.GaussianGaussian}},\code{\link{marginalLikelihood.GaussianGaussian}},\code{\link{rPosteriorPredictive.GaussianGaussian}},\code{\link{dPosteriorPredictive.GaussianGaussian}}.
}
