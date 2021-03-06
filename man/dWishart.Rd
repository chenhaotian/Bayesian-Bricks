% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Gamma_Inference.r
\name{dWishart}
\alias{dWishart}
\title{Density function of Wishart distribution}
\usage{
dWishart(x, df, rate, LOG = TRUE)
}
\arguments{
\item{x}{matrix, a symmetric positive-definite matrix.}

\item{df}{numeric, the degree of freedom.}

\item{rate}{matrix, a symmetric positive-definite matrix, the 'rate', or 'inverse-scale' parameter. The 'rate' parameter in Wishart is the 'scale' parameter in InvWishart}

\item{LOG}{logical, return log density of LOG=TRUE, default TRUE.}
}
\value{
A numeric vector, the density values.
}
\description{
For a random matrix x, the density function of Wishart distribution is defined as:
    \deqn{(2^{(df p)/2} Gamma_p(df/2) |rate|^{-df/2})^{-1} |x|^{(df-p-1)/2} exp(-1/2 tr(x rate))}
Where x is a pxp symmetric positive definite matrix, Gamma_p() is the multivariate Gamma function of dimension p.
}
\examples{
##generate a symmetric positive-definite matrix
x <- crossprod(matrix(rnorm(15),5,3))
rate <- crossprod(matrix(rnorm(15),5,3)) #the prior inverse-scale of x
dWishart(x,df = 5,rate = rate,LOG = TRUE)
dWishart(x,df = 5,rate = rate,LOG = FALSE)
}
\references{
Wishart, John. "The generalized product moment distribution in samples from a normal multivariate population." Biometrika (1928): 32-52.

MARolA, K. V., JT KBNT, and J. M. Bibly. Multivariate analysis. AcadeInic Press, Londres, 1979.
}
