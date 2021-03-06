% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Gaussian_Inference.r
\name{sufficientStatistics_Weighted.LinearGaussianGaussian}
\alias{sufficientStatistics_Weighted.LinearGaussianGaussian}
\title{Weighted sufficient statistics of a "LinearGaussianGaussian" object}
\usage{
\method{sufficientStatistics_Weighted}{LinearGaussianGaussian}(obj, x, w, A, b = NULL, foreach = FALSE, ...)
}
\arguments{
\item{obj}{A "LinearGaussianGaussian" object.}

\item{x}{matrix, Gaussian samples, when x is a matrix, each row is a sample of dimension ncol(x). when x is a vector, x is length(x) samples of dimension 1.}

\item{w}{numeric, sample weights.}

\item{A}{matrix or list. when x is a \eqn{N x 1} matrix, A must be a matrix of \eqn{N x dimz}, dimz is the dimension of z; When x is a \eqn{N x dimx} matrix, where \eqn{dimx > 1}, A can be either a list or a matrix. When A is a list, \eqn{A = {A_1,A_2,...A_N}} is a list of \eqn{dimx x dimz} matrices. If A is a single \eqn{dimx x dimz} matrix, it will be replicated N times into a length N list}

\item{b}{matrix, when x is a \eqn{N x 1} matrix, b must also be a \eqn{N x 1} matrix or length N vector; When x is a \eqn{N x dimx} matrix, where \eqn{dimx > 1}, b can be either a matrix or a vector. When b is a matrix, \eqn{b={b_1^T,...,b_N^T}} is a \eqn{N x dimx} matrix, each row is a transposed vector. When b is a length \eqn{dimx} vector, it will be transposed into a row vector and replicated N times into a \eqn{N x dimx} matrix. When b = NULL, it will be treated as a vector of zeros. Default NULL.}

\item{foreach}{logical, specifying whether to return the sufficient statistics for each observation. Default FALSE.}

\item{...}{Additional arguments to be passed to other inherited types.}
}
\value{
If foreach=TRUE, will return a list of sufficient statistics for each row of x, otherwise will return the sufficient statistics of x as a whole.
}
\description{
For following model structure:
     \deqn{x \sim Gaussian(A z + b, Sigma)}
     \deqn{z \sim Gaussian(m,S)}
Where Sigma is known. A is a \eqn{dimx x dimz} matrix, x is a \eqn{dimx x 1} random vector, z is a \eqn{dimz x 1} random vector, b is a \eqn{dimm x 1} vector. Gaussian() is the Gaussian distribution. See \code{?dGaussian} for the definition of Gaussian distribution.\cr
For weight vector w and one dimensional observations: x is a vector of length N, or a \eqn{N x 1} matrix, each row is an observation, must satisfy \eqn{nrow(x)=length(w)}; A is a \eqn{N x dimz} matrix; b is a length N vector. The sufficient statistics are:
\itemize{
 \item  SA = \eqn{A^T (A w) / Sigma}
 \item SAx = \eqn{A^T ((x-b) w) / Sigma}
}
For weight vector w and \eqn{dimx} dimensional observations: x must be a \eqn{N x m} matrix, each row is an observation, must satisfy \eqn{nrow(x)=length(w)}; A can be either a list or a matrix. When A is a list, \eqn{A = {A_1,A_2,...A_N}} is a list of \eqn{dimx x dimz} matrices. If A is a single \eqn{dimx x dimz} matrix, it will be replicated N times into a length N list; b can be either a matrix or a vector. When b is a matrix, \eqn{b={b_1^T,...,b_N^T}} is a \eqn{N x dimx} matrix, each row is a transposed vector. When b is a length \eqn{dimx} vector, it will be transposed into a row vector and replicated N times into a \eqn{N x dimx} matrix. The sufficient statistics are:
\itemize{
 \item   SA = \eqn{sum_{i=1:N} w_i A_i^T Sigma^{-1} A_i}
 \item  SAx = \eqn{sum_{i=1:N} w_i A_i^T Sigma^{-1} (x_i-b_i)}
}
}
\examples{
obj <- LinearGaussianGaussian(gamma=list(Sigma=matrix(c(2,1,1,2),2,2),m=c(0.2,0.5,0.6),S=diag(3)))
x <- rGaussian(100,mu = runif(2),Sigma = diag(2))
w <- runif(100)
A <- matrix(runif(6),2,3)
b <- runif(2)
sufficientStatistics_Weighted(obj,x=x,w=w,A=A,b=b)
sufficientStatistics_Weighted(obj,x=x,w=w,A=A,b=b,foreach = TRUE)
}
\references{
Murphy, Kevin P. Machine learning: a probabilistic perspective. MIT press, 2012.
}
\seealso{
\code{\link{LinearGaussianGaussian}}, \code{\link{sufficientStatistics.LinearGaussianGaussian}}
}
