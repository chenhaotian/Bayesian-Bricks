


#' @title Create objects of type '"BasicBayesian"'.
#' @description
#' A Basic Bayesian Object is with following conditional dependency structure: 
#'      \deqn{theta|gamma \sim H(gamma)}
#'      \deqn{X|theta \sim F(theta)}
#' Where H(gamma) is usually called "the prior distribution", F(theta) is called "the observation distribution". Objects of type "LinearGaussianGaussian", "GaussianGaussian", "GaussianInvWishart", "GaussianNIW", "GaussianNIG", "CatDirichlet" and "CatDP" are all "BasicBayesian"s.
#' 
#' @seealso \code{\link{LinearGaussianGaussian}} for Linear Gaussian and Gaussian conjugate structure, \code{\link{GaussianGaussian}} for Gaussian-Gaussian conjugate structure, \code{\link{GaussianInvWishart}} for Gaussian-Inverse-Wishart conjugate structure, \code{\link{GaussianNIW}} for Gaussian-NIW conjugate structure, \code{\link{GaussianNIG}} for Gaussian-NIG conjugate structure, \code{\link{CatDirichlet}} for Categorical-Dirichlet conjugate structure, \code{\link{CatDP}} for Categorical-DP conjugate structure ...
#' @param ENV The environment where you want to create the BasicBayesian object
#' @return An object of class "BasicBayesian".
#' @export
BasicBayesian <- function(ENV=parent.frame()){
    object <- new.env(parent = ENV)
    
    object$gamma <- list()
    object$H <- character(0)
    object$F <- character(0)

    class(object) <- c("BasicBayesian","BayesianBrick")
    return(object)

}

#' @title Get sample sufficient statistics
#' @description
#' This is a generic function that will generate the sufficient statistics of a given Bayesian bricks object. i.e. \cr
#' for the model structure:
#'      \deqn{theta|gamma \sim H(gamma)}
#'      \deqn{x|theta \sim F(theta)}
#' get the sufficient statistics T(x). \cr
#' For a given sample set x, each row of x is an observation, and a Bayesian bricks object obj. \code{sufficientStatistics()} return the sufficient statistics for different model structures:
#' \subsection{class(obj)="LinearGaussianGaussian"}{
#'      \deqn{x \sim Gaussian(A z + b, Sigma)}
#'      \deqn{z \sim Gaussian(m,S)}
#' The sufficient statistics are:
#' \itemize{
#'  \item   SA = \eqn{sum_{i=1:N} A_i^T Sigma^{-1} A_i}
#'  \item  SAx = \eqn{sum_{i=1:N} A_i^T Sigma^{-1} (x_i-b_i)}
#' }
#' See \code{?sufficientStatistics.LinearGaussianGaussian} for details.
#' }
#' \subsection{class(obj)="GaussianGaussian"}{
#' Where
#' \deqn{x \sim Gaussian(mu,Sigma)}
#' \deqn{mu \sim Gaussian(m,S)}
#' Sigma is known.
#' The sufficient statistics are:
#' \itemize{
#' \item N: the effective number of samples.
#' \item xsum: the row sums of the samples.
#' }
#' See \code{?sufficientStatistics.GaussianGaussian} for details.
#' }
#' \subsection{class(obj)="GaussianInvWishart"}{
#' Where
#' \deqn{x \sim Gaussian(mu,Sigma)}
#' \deqn{Sigma \sim InvWishart(v,S)}
#' mu is known. \cr
#' The sufficient statistics are:
#' \itemize{
#' \item N: the effective number of samples.
#' \item xsum: the sample scatter matrix centered on the mean vector.
#' }
#' See \code{?sufficientStatistics.GaussianInvWishart} for details.
#' }
#' \subsection{class(obj)="GaussianNIW"}{
#' Where
#' \deqn{x \sim Gaussian(mu,Sigma)}
#' \deqn{Sigma \sim InvWishart(v,S)}
#' \deqn{mu \sim Gaussian(m,Sigma/k)}
#' The sufficient statistics are:
#' \itemize{
#' \item N: the effective number of samples.
#' \item xsum: the row sums of the samples.
#' \item S: the uncentered sample scatter matrix.
#' }
#' See \code{?sufficientStatistics.GaussianNIW} for details.
#' }
#' \subsection{class(obj)="GaussianNIG"}{
#' Where
#' \deqn{x \sim Gaussian(X beta,sigma^2)}
#' \deqn{sigma^2 \sim InvGamma(a,b)}
#' \deqn{beta \sim Gaussian(m,sigma^2 V)}
#' X is a row vector, or a design matrix where each row is an obervation.
#' The sufficient statistics are:
#' \itemize{
#' \item N: the effective number of samples.
#' \item SXx: covariance of X and x
#' \item SX: the uncentered sample scatter matrix.
#' \item Sx: the variance of x
#' }
#' See \code{?sufficientStatistics.GaussianNIG} for details.
#' }
#' \subsection{class(obj)="CatDirichlet"}{
#' Where
#' \deqn{x \sim Categorical(pi)}
#' \deqn{pi \sim Dirichlet(alpha)}
#' The sufficient statistics of CatDirichlet object can either be x itself, or the counts of the unique labels in x.\cr
#' See \code{?sufficientStatistics.CatDirichlet} for details.
#' }
#' \subsection{class(obj)="CatDP"}{
#' Where
#' \deqn{x \sim Categorical(pi)}
#' \deqn{pi \sim DirichletProcess(alpha)}
#' The sufficient statistics of CatDP object can either be x itself, or the counts of the unique labels in x.\cr
#' See \code{?sufficientStatistics.CatDP} for details.
#' }
#' \subsection{class(obj)="DP"}{
#' Where
#'      \deqn{pi|alpha \sim DP(alpha,U)}
#'     \deqn{z|pi \sim Categorical(pi)}
#'      \deqn{theta_z|psi \sim H0(psi)}
#'      \deqn{x|theta_z,z \sim F(theta_z)}
#' The sufficient statistics of "DP" object is the same sufficient statistics of the "BasicBayesian" inside the "DP".
#' See \code{?sufficientStatistics.DP} for details.
#' }
#' \subsection{class(obj)="HDP"}{
#' Where
#'      \deqn{G_j|gamma \sim DP(gamma,U), j = 1:J}
#'      \deqn{pi_j|G_j,alpha \sim DP(alpha,G_j)}
#'      \deqn{z|pi_j \sim Categorical(pi_j)}
#'      \deqn{k|z,G_j \sim Categorical(G_j), if z is a sample from the base measure G}
#'      \deqn{theta_k|psi \sim H0(psi)}
#' The sufficient statistics of "HDP" object is the same sufficient statistics of the "BasicBayesian" inside the "HDP".
#' See \code{?sufficientStatistics.HDP} for details.
#' }
#' \subsection{class(obj)="HDP2"}{
#' Where
#'      \deqn{G_m |eta \sim DP(eta,U), m = 1:M}
#'      \deqn{G_mj|gamma,G_m \sim DP(gamma,G_m), j = 1:J_m}
#'      \deqn{pi_mj|G_mj,alpha \sim DP(alpha,G_mj)}
#'      \deqn{z|pi_mj \sim Categorical(pi_mj)}
#'      \deqn{k|z,G_mj \sim Categorical(G_mj), if z is a sample from the base measure G_mj}
#'      \deqn{u|k,G_m \sim Categorical(G_m), if k is a sample from the base measure G_m}
#'      \deqn{theta_u|psi \sim H0(psi)}
#'      \deqn{x|theta_u,u \sim F(theta_u)}
#' The sufficient statistics of "HDP2" object is the same sufficient statistics of the "BasicBayesian" inside the "HDP2".
#' See \code{?sufficientStatistics.HDP2} for details.
#' }
#'
#' @seealso \code{\link{sufficientStatistics.LinearGaussianGaussian}} for Linear Gaussian and Gaussian conjugate structure, \code{\link{sufficientStatistics.GaussianGaussian}} for Gaussian-Gaussian conjugate structure, \code{\link{sufficientStatistics.GaussianInvWishart}} for Gaussian-Inverse-Wishart conjugate structure, \code{\link{sufficientStatistics.GaussianNIW}} for Gaussian-NIW conjugate structure, \code{\link{sufficientStatistics.GaussianNIG}} for Gaussian-NIG conjugate structure, \code{\link{sufficientStatistics.CatDirichlet}} for Categorical-Dirichlet conjugate structure, \code{\link{sufficientStatistics.CatDP}} for Categorical-DP conjugate structure ...
#' @param obj a "BayesianBrick" object used to select a method.
#' @param x a set of samples.
#' @param ... further arguments passed to or from other methods.
#' @return An object of corresponding sufficient statistics class, such as "ssGaussian"
#' @export
#' @examples
#' x <- rGaussian(10,mu = 1,Sigma = 1)
#' obj <- GaussianNIW()                    #an GaussianNIW object
#' sufficientStatistics(obj=obj,x=x)
sufficientStatistics <- function(obj,x,...) UseMethod("sufficientStatistics")

#' @title Get weighted sample sufficient statistics
#' @description
#' This is a generic function that will generate the weighted sufficient statistics of a given "BayesianBrick" object. That is, for the model structure:
#'      \deqn{theta|gamma \sim H(gamma)}
#'      \deqn{x|theta \sim F(theta)}
#' get the weighted sufficient statistics T(x).
#' For a given sample set x, each row of x is an observation, the sample weights w, and a Bayesian bricks object obj. \code{sufficientStatistics_Weighted()} return the weighted sufficient statistics for different model structures:
#' \subsection{class(obj)="LinearGaussianGaussian"}{
#'      \deqn{x \sim Gaussian(A z + b, Sigma)}
#'      \deqn{z \sim Gaussian(m,S)}
#' The sufficient statistics are:
#' \itemize{
#'  \item   SA = \eqn{sum_{i=1:N} w_i A_i^T Sigma^{-1} A_i}
#'  \item  SAx = \eqn{sum_{i=1:N} w_i A_i^T Sigma^{-1} (x_i-b_i)}
#' }
#' See \code{?sufficientStatistics.LinearGaussianGaussian} for details.
#' }
#' \subsection{class(obj)="GaussianGaussian"}{
#' Where
#' \deqn{x \sim Gaussian(mu,Sigma)}
#' \deqn{mu \sim Gaussian(m,S)}
#' Sigma is known.
#' The sufficient statistics are:
#' \itemize{
#' \item N: the effective number of samples.
#' \item xsum: the row sums of the samples.
#' }
#' See \code{?sufficientStatistics_Weighted.GaussianGaussian} for details.
#' }
#' \subsection{class(obj)="GaussianInvWishart"}{
#' Where
#' \deqn{x \sim Gaussian(mu,Sigma)}
#' \deqn{Sigma \sim InvWishart(v,S)}
#' mu is known. \cr
#' The sufficient statistics are:
#' \itemize{
#' \item N: the effective number of samples.
#' \item xsum: the sample scatter matrix centered on the mean vector.
#' }
#' See \code{?sufficientStatistics_Weighted.GaussianInvWishart} for details.
#' }
#' \subsection{class(obj)="GaussianNIW"}{
#' Where
#' \deqn{x \sim Gaussian(mu,Sigma)}
#' \deqn{Sigma \sim InvWishart(v,S)}
#' \deqn{mu \sim Gaussian(m,Sigma/k)}
#' The sufficient statistics are:
#' \itemize{
#' \item N: the effective number of samples.
#' \item xsum: the row sums of the samples.
#' \item S: the uncentered sample scatter matrix.
#' }
#' See \code{?sufficientStatistics_Weighted.GaussianNIW} for details.
#' }
#' \subsection{class(obj)="GaussianNIG"}{
#' Where
#' \deqn{x \sim Gaussian(X beta,sigma^2)}
#' \deqn{sigma^2 \sim InvGamma(a,b)}
#' \deqn{beta \sim Gaussian(m,sigma^2 V)}
#' X is a row vector, or a design matrix where each row is an obervation.
#' The sufficient statistics are:
#' \itemize{
#' \item N: the effective number of samples.
#' \item SXx: covariance of X and x
#' \item SX: the uncentered sample scatter matrix.
#' \item Sx: the variance of x
#' }
#' See \code{?sufficientStatistics_Weighted.GaussianNIG} for details.
#' }
#' \subsection{class(obj)="CatDirichlet"}{
#' Where
#' \deqn{x \sim Categorical(pi)}
#' \deqn{pi \sim Dirichlet(alpha)}
#' The sufficient statistics of CatDirichlet object can either be x itself, or the counts of the unique labels in x.\cr
#' See \code{?sufficientStatistics_Weighted.CatDirichlet} for details.
#' }
#' \subsection{class(obj)="CatDP"}{
#' Where
#' \deqn{x \sim Categorical(pi)}
#' \deqn{pi \sim DirichletProcess(alpha)}
#' The sufficient statistics of CatDP object can either be x itself, or the counts of the unique labels in x.\cr
#' See \code{?sufficientStatistics_Weighted.CatDP} for details.
#' }
#' \subsection{class(obj)="DP"}{
#' Where
#'      \deqn{pi|alpha \sim DP(alpha,U)}
#'     \deqn{z|pi \sim Categorical(pi)}
#'      \deqn{theta_z|psi \sim H0(psi)}
#'      \deqn{x|theta_z,z \sim F(theta_z)}
#' The sufficient statistics of "DP" object is the same sufficient statistics of the "BasicBayesian" inside the "DP".
#' See \code{?sufficientStatistics_Weighted.DP} for details.
#' }
#' \subsection{class(obj)="HDP"}{
#' Where
#'      \deqn{G_j|gamma \sim DP(gamma,U), j = 1:J}
#'      \deqn{pi_j|G_j,alpha \sim DP(alpha,G_j)}
#'      \deqn{z|pi_j \sim Categorical(pi_j)}
#'      \deqn{k|z,G_j \sim Categorical(G_j), if z is a sample from the base measure G}
#'      \deqn{theta_k|psi \sim H0(psi)}
#' The sufficient statistics of "HDP" object is the same sufficient statistics of the "BasicBayesian" inside the "HDP".
#' See \code{?sufficientStatistics_Weighted.HDP} for details.
#' }
#' \subsection{class(obj)="HDP2"}{
#' Where
#'      \deqn{G_m |eta \sim DP(eta,U), m = 1:M}
#'      \deqn{G_mj|gamma,G_m \sim DP(gamma,G_m), j = 1:J_m}
#'      \deqn{pi_mj|G_mj,alpha \sim DP(alpha,G_mj)}
#'      \deqn{z|pi_mj \sim Categorical(pi_mj)}
#'      \deqn{k|z,G_mj \sim Categorical(G_mj), if z is a sample from the base measure G_mj}
#'      \deqn{u|k,G_m \sim Categorical(G_m), if k is a sample from the base measure G_m}
#'      \deqn{theta_u|psi \sim H0(psi)}
#'      \deqn{x|theta_u,u \sim F(theta_u)}
#' The sufficient statistics of "HDP2" object is the same sufficient statistics of the "BasicBayesian" inside the "HDP2".
#' See \code{?sufficientStatistics_Weighted.HDP2} for details.
#' }
#' @seealso \code{\link{sufficientStatistics_Weighted.LinearGaussianGaussian}} for Linear Gaussian and Gaussian conjugate structure, \code{\link{sufficientStatistics_Weighted.GaussianGaussian}} for Gaussian-Gaussian conjugate structure, \code{\link{sufficientStatistics_Weighted.GaussianInvWishart}} for Gaussian-Inverse-Wishart conjugate structure, \code{\link{sufficientStatistics_Weighted.GaussianNIW}} for Gaussian-NIW conjugate structure, \code{\link{sufficientStatistics_Weighted.GaussianNIG}} for Gaussian-NIG conjugate structure, \code{\link{sufficientStatistics_Weighted.CatDirichlet}} for Categorical-Dirichlet conjugate structure, \code{\link{sufficientStatistics_Weighted.CatDP}} for Categorical-DP conjugate structure ...
#' @param obj a "BayesianBrick" object used to select a method.
#' @param x a set of samples.
#' @param w numeric, sample weights.
#' @param ... further arguments passed to or from other methods.
#' @return An object of corresponding sufficient statistics class, such as "ssGaussian"
#' @export
#' @examples
#' \donttest{
#' x <- rGaussian(10,mu = 1,Sigma = 1)
#' w <- runif(10)
#' obj <- GaussianNIW()                    #an GaussianNIW object
#' sufficientStatistics_Weighted(obj=obj,x=x,w=w)
#' }
sufficientStatistics_Weighted <- function(obj,x,w,...) UseMethod("sufficientStatistics_Weighted")

#' @title update the prior distribution with sufficient statistics
#' @description update the prior distribution with sufficient statistics
#' @param obj A "BayesianBrick" object used to select a method.
#' @param ... further arguments passed to or from other methods.
#' @return None, or an error message if the update fails.
posterior_bySufficientStatistics <- function(obj,...) UseMethod("posterior_bySufficientStatistics")

#' @title update the prior distribution with sufficient statistics
#' @description
#' This is a generic function that will update the prior distribution of a "BayesianBrick" object by adding information of the observation's sufficient statistics. i.e. for the model structure:
#'      \deqn{theta|gamma \sim H(gamma)}
#'      \deqn{x|theta \sim F(theta)}
#' update gamma to gamma_posterior by adding the information of x to gamma.\cr
#' For a given sample set x or it's sufficient statistics ss, and a Bayesian bricks object obj, \code{posterior()} will update the posterior parameters in obj for different model structures:
#' \subsection{class(obj)="LinearGaussianGaussian"}{
#'      \deqn{x \sim Gaussian(A z + b, Sigma)}
#'      \deqn{z \sim Gaussian(m,S)}
#' \code{posterior()} will update m and S in obj.
#' See \code{?posterior.LinearGaussianGaussian} for details.
#' }
#' \subsection{class(obj)="GaussianGaussian"}{
#' Where
#' \deqn{x \sim Gaussian(mu,Sigma)}
#' \deqn{mu \sim Gaussian(m,S)}
#' Sigma is known.
#' \code{posterior()} will update m and S in obj.
#' See \code{?posterior.GaussianGaussian} for details.
#' }
#' \subsection{class(obj)="GaussianInvWishart"}{
#' Where
#' \deqn{x \sim Gaussian(mu,Sigma)}
#' \deqn{Sigma \sim InvWishart(v,S)}
#' mu is known.
#' \code{posterior()} will update v and S in obj.
#' See \code{?posterior.GaussianInvWishart} for details.
#' }
#' \subsection{class(obj)="GaussianNIW"}{
#' Where
#' \deqn{x \sim Gaussian(mu,Sigma)}
#' \deqn{Sigma \sim InvWishart(v,S)}
#' \deqn{mu \sim Gaussian(m,Sigma/k)}
#' \code{posterior()} will update m, k, v and S in obj.
#' See \code{?posterior.GaussianNIW} for details.
#' }
#' \subsection{class(obj)="GaussianNIG"}{
#' Where
#' \deqn{x \sim Gaussian(X beta,sigma^2)}
#' \deqn{sigma^2 \sim InvGamma(a,b)}
#' \deqn{beta \sim Gaussian(m,sigma^2 V)}
#' X is a row vector, or a design matrix where each row is an obervation.
#' \code{posterior()} will update m, V, a and b in obj.
#' See \code{?posterior.GaussianNIG} for details.
#' }
#' \subsection{class(obj)="CatDirichlet"}{
#' Where
#' \deqn{x \sim Categorical(pi)}
#' \deqn{pi \sim Dirichlet(alpha)}
#' \code{posterior()} will update alpha in obj.
#' See \code{?posterior.CatDirichlet} for details.
#' }
#' \subsection{class(obj)="CatDP"}{
#' Where
#' \deqn{x \sim Categorical(pi)}
#' \deqn{pi \sim DirichletProcess(alpha)}
#' \code{posterior()} will update alpha in obj.
#' See \code{?posterior.CatDP} for details.
#' }
#' \subsection{class(obj)="DP"}{
#' Where
#'      \deqn{pi|alpha \sim DP(alpha,U)}
#'     \deqn{z|pi \sim Categorical(pi)}
#'      \deqn{theta_z|psi \sim H0(psi)}
#'      \deqn{x|theta_z,z \sim F(theta_z)}
#' \code{posterior()} will update alpha and psi in obj.
#' See \code{?posterior.DP} for details.
#' }
#' \subsection{class(obj)="HDP"}{
#' Where
#'      \deqn{G_j|gamma \sim DP(gamma,U), j = 1:J}
#'      \deqn{pi_j|G_j,alpha \sim DP(alpha,G_j)}
#'      \deqn{z|pi_j \sim Categorical(pi_j)}
#'      \deqn{k|z,G_j \sim Categorical(G_j), if z is a sample from the base measure G}
#'      \deqn{theta_k|psi \sim H0(psi)}
#' \code{posterior()} will update gamma, alpha and psi in obj.
#' See \code{?posterior.HDP} for details.
#' }
#' \subsection{class(obj)="HDP2"}{
#' Where
#'      \deqn{G_m |eta \sim DP(eta,U), m = 1:M}
#'      \deqn{G_mj|gamma,G_m \sim DP(gamma,G_m), j = 1:J_m}
#'      \deqn{pi_mj|G_mj,alpha \sim DP(alpha,G_mj)}
#'      \deqn{z|pi_mj \sim Categorical(pi_mj)}
#'      \deqn{k|z,G_mj \sim Categorical(G_mj), if z is a sample from the base measure G_mj}
#'      \deqn{u|k,G_m \sim Categorical(G_m), if k is a sample from the base measure G_m}
#'      \deqn{theta_u|psi \sim H0(psi)}
#'      \deqn{x|theta_u,u \sim F(theta_u)}
#' \code{posterior()} will update eta, gamma, alpha and psi in obj.
#' See \code{?posterior.HDP2} for details.
#' }
#'
#' @seealso \code{\link{posterior.LinearGaussianGaussian}} for Linear Gaussian and Gaussian conjugate structure, \code{\link{posterior.GaussianGaussian}} for Gaussian-Gaussian conjugate structure, \code{\link{posterior.GaussianInvWishart}} for Gaussian-Inverse-Wishart conjugate structure, \code{\link{posterior.GaussianNIW}} for Gaussian-NIW conjugate structure, \code{\link{posterior.GaussianNIG}} for Gaussian-NIG conjugate structure, \code{\link{posterior.CatDirichlet}} for Categorical-Dirichlet conjugate structure, \code{\link{posterior.CatDP}} for Categorical-DP conjugate structure ...
#' @param obj A "BayesianBrick" object used to select a method.
#' @param ... further arguments passed to or from other methods.
#' @return None, or an error message if the update fails.
#' @export
posterior <- function(obj,...) UseMethod("posterior")

#' @title update the prior distribution with sufficient statistics
#' @description update the prior distribution with sufficient statistics
#' @param obj A "BayesianBrick" object used to select a method.
#' @param ... other parameters.
#' @return None, or an error message if the update fails.
posteriorDiscard_bySufficientStatistics <- function(obj,...) UseMethod("posteriorDiscard_bySufficientStatistics")

#' @title update the prior distribution with sufficient statistics
#' @description
#' Contrary to posterior(), posteriorDiscard() a generic function that will update the prior distribution of a "BayesianBrick" object by removing the information provided by the observation's sufficient statistics. i.e. for the model structure:
#'      \deqn{theta|gamma \sim H(gamma)}
#'      \deqn{x|theta \sim F(theta)}
#' update gamma to gamma_posterior by removing the information of x from gamma.
#' For a given sample set x or it's sufficient statistics ss, and a Bayesian bricks object obj, \code{posteriorDiscard()} will update the posterior parameters in obj for different model structures:
#' \subsection{class(obj)="LinearGaussianGaussian"}{
#'      \deqn{x \sim Gaussian(A z + b, Sigma)}
#'      \deqn{z \sim Gaussian(m,S)}
#' \code{posteriorDiscard()} will update m and S in obj.
#' See \code{?posteriorDiscard.LinearGaussianGaussian} for details.
#' }
#' \subsection{class(obj)="GaussianGaussian"}{
#' Where
#' \deqn{x \sim Gaussian(mu,Sigma)}
#' \deqn{mu \sim Gaussian(m,S)}
#' Sigma is known.
#' \code{posteriorDiscard()} will update m and S in obj.
#' See \code{?posteriorDiscard.GaussianGaussian} for details.
#' }
#' \subsection{class(obj)="GaussianInvWishart"}{
#' Where
#' \deqn{x \sim Gaussian(mu,Sigma)}
#' \deqn{Sigma \sim InvWishart(v,S)}
#' mu is known.
#' \code{posteriorDiscard()} will update v and S in obj.
#' See \code{?posteriorDiscard.GaussianInvWishart} for details.
#' }
#' \subsection{class(obj)="GaussianNIW"}{
#' Where
#' \deqn{x \sim Gaussian(mu,Sigma)}
#' \deqn{Sigma \sim InvWishart(v,S)}
#' \deqn{mu \sim Gaussian(m,Sigma/k)}
#' \code{posteriorDiscard()} will update m, k, v and S in obj.
#' See \code{?posteriorDiscard.GaussianNIW} for details.
#' }
#' \subsection{class(obj)="GaussianNIG"}{
#' Where
#' \deqn{x \sim Gaussian(X beta,sigma^2)}
#' \deqn{sigma^2 \sim InvGamma(a,b)}
#' \deqn{beta \sim Gaussian(m,sigma^2 V)}
#  X is a row vector, or a design matrix where each row is an obervation.
#' \code{posteriorDiscard()} will update m, V, a and b in obj.
#' See \code{?posteriorDiscard.GaussianNIG} for details.
#' }
#' \subsection{class(obj)="CatDirichlet"}{
#' Where
#' \deqn{x \sim Categorical(pi)}
#' \deqn{pi \sim Dirichlet(alpha)}
#' \code{posteriorDiscard()} will update alpha in obj.
#' See \code{?posteriorDiscard.CatDirichlet} for details.
#' }
#' \subsection{class(obj)="CatDP"}{
#' Where
#' \deqn{x \sim Categorical(pi)}
#' \deqn{pi \sim DirichletProcess(alpha)}
#' \code{posteriorDiscard()} will update alpha in obj.
#' See \code{?posteriorDiscard.CatDP} for details.
#' }
#' \subsection{class(obj)="DP"}{
#' Where
#'      \deqn{pi|alpha \sim DP(alpha,U)}
#'     \deqn{z|pi \sim Categorical(pi)}
#'      \deqn{theta_z|psi \sim H0(psi)}
#'      \deqn{x|theta_z,z \sim F(theta_z)}
#' \code{posteriorDiscard()} will update alpha and psi in obj.
#' See \code{?posteriorDiscard.DP} for details.
#' }
#' \subsection{class(obj)="HDP"}{
#' Where
#'      \deqn{G_j|gamma \sim DP(gamma,U), j = 1:J}
#'      \deqn{pi_j|G_j,alpha \sim DP(alpha,G_j)}
#'      \deqn{z|pi_j \sim Categorical(pi_j)}
#'      \deqn{k|z,G_j \sim Categorical(G_j), if z is a sample from the base measure G}
#'      \deqn{theta_k|psi \sim H0(psi)}
#' \code{posteriorDiscard()} will update gamma, alpha and psi in obj.
#' See \code{?posteriorDiscard.HDP} for details.
#' }
#' \subsection{class(obj)="HDP2"}{
#' Where
#'      \deqn{G_m |eta \sim DP(eta,U), m = 1:M}
#'      \deqn{G_mj|gamma,G_m \sim DP(gamma,G_m), j = 1:J_m}
#'      \deqn{pi_mj|G_mj,alpha \sim DP(alpha,G_mj)}
#'      \deqn{z|pi_mj \sim Categorical(pi_mj)}
#'      \deqn{k|z,G_mj \sim Categorical(G_mj), if z is a sample from the base measure G_mj}
#'      \deqn{u|k,G_m \sim Categorical(G_m), if k is a sample from the base measure G_m}
#'      \deqn{theta_u|psi \sim H0(psi)}
#'      \deqn{x|theta_u,u \sim F(theta_u)}
#' \code{posteriorDiscard()} will update eta, gamma, alpha and psi in obj.
#' See \code{?posteriorDiscard.HDP2} for details.
#' }
#'
#' @seealso \code{\link{posteriorDiscard.LinearGaussianGaussian}} for Linear Gaussian and Gaussian conjugate structure, \code{\link{posteriorDiscard.GaussianGaussian}} for Gaussian-Gaussian conjugate structure, \code{\link{posteriorDiscard.GaussianInvWishart}} for Gaussian-Inverse-Wishart conjugate structure, \code{\link{posteriorDiscard.GaussianNIW}} for Gaussian-NIW conjugate structure, \code{\link{posteriorDiscard.GaussianNIG}} for Gaussian-NIG conjugate structure, \code{\link{posteriorDiscard.CatDirichlet}} for Categorical-Dirichlet conjugate structure, \code{\link{posteriorDiscard.CatDP}} for Categorical-DP conjugate structure ...
#' @param obj A "BayesianBrick" object used to select a method.
#' @param ... further arguments passed to or from other methods.
#' @return None, or an error message if the update fails.
#' @export
posteriorDiscard <- function(obj,...) UseMethod("posteriorDiscard")

#' @title Get the Maximum A Posteriori(MAP) estimate of a "BayesianBrick" object
#' @description
#' This is a generic function that will generate the MAP estimate of a given "BayesianBrick" object. For the model structure:  \cr
#'      \deqn{theta|gamma \sim H(gamma)}
#'      \deqn{x|theta \sim F(theta)}
#' MAP estimate of theta is theta_MAP = argmax_theta p(theta|gamma,x).
#' For a given Bayesian bricks object obj, the MAP estimate will be:
#' \subsection{class(obj)="LinearGaussianGaussian"}{
#' Where
#'      \deqn{x \sim Gaussian(A z + b, Sigma)}
#'      \deqn{z \sim Gaussian(m,S)}
#' \code{MAP()} will return the MAP estimate of z.
#' See \code{?MAP.LinearGaussianGaussian} for details.
#' }
#' \subsection{class(obj)="GaussianGaussian"}{
#' Where
#' \deqn{x \sim Gaussian(mu,Sigma)}
#' \deqn{mu \sim Gaussian(m,S)}
#' Sigma is known.
#' \code{MAP()} will return the MAP estimate of mu.
#' See \code{?MAP.GaussianGaussian} for details.
#' }
#' \subsection{class(obj)="GaussianInvWishart"}{
#' Where
#' \deqn{x \sim Gaussian(mu,Sigma)}
#' \deqn{Sigma \sim InvWishart(v,S)}
#' mu is known.
#' \code{MAP()} will return the MAP estimate of Sigma.
#' See \code{?MAP.GaussianInvWishart} for details.
#' }
#' \subsection{class(obj)="GaussianNIW"}{
#' Where
#' \deqn{x \sim Gaussian(mu,Sigma)}
#' \deqn{Sigma \sim InvWishart(v,S)}
#' \deqn{mu \sim Gaussian(m,Sigma/k)}
#' \code{MAP()} will return the MAP estimate of mu and Sigma.
#' See \code{?MAP.GaussianNIW} for details.
#' }
#' \subsection{class(obj)="GaussianNIG"}{
#' Where
#' \deqn{x \sim Gaussian(X beta,sigma^2)}
#' \deqn{sigma^2 \sim InvGamma(a,b)}
#' \deqn{beta \sim Gaussian(m,sigma^2 V)}
#' X is a row vector, or a design matrix where each row is an obervation.
#' \code{MAP()} will return the MAP estimate of beta and sigma^2.
#' See \code{?MAP.GaussianNIG} for details.
#' }
#' \subsection{class(obj)="CatDirichlet"}{
#' Where
#' \deqn{x \sim Categorical(pi)}
#' \deqn{pi \sim Dirichlet(alpha)}
#' \code{MAP()} will return the MAP estimate of pi.
#' See \code{?MAP.CatDirichlet} for details.
#' }
#' \subsection{class(obj)="CatDP"}{
#' Where
#' \deqn{x \sim Categorical(pi)}
#' \deqn{pi \sim DirichletProcess(alpha)}
#' \code{MAP()} will return the MAP estimate of pi.
#' See \code{?MAP.CatDP} for details.
#' }
#' @seealso \code{\link{MAP.LinearGaussianGaussian}} for Linear Gaussian and Gaussian conjugate structure, \code{\link{MAP.GaussianGaussian}} for Gaussian-Gaussian conjugate structure, \code{\link{MAP.GaussianInvWishart}} for Gaussian-Inverse-Wishart conjugate structure, \code{\link{MAP.GaussianNIW}} for Gaussian-NIW conjugate structure, \code{\link{MAP.GaussianNIG}} for Gaussian-NIG conjugate structure, \code{\link{MAP.CatDirichlet}} for Categorical-Dirichlet conjugate structure, \code{\link{MAP.CatDP}} for Categorical-DP conjugate structure ...
#' @param obj A "BayesianBrick" object used to select a method.
#' @param ... further arguments passed to or from other methods.
#' @return A list of the MAP estimates
#' @export
MAP <- function(obj,...) UseMethod("MAP")

#' @title Get the Mean Posterior Estimate(MPE) of a "BayesianBrick" object
#' @description
#' This is a generic function that will generate the MPE estimate of a given "BayesianBrick" object. i.e. for the model structure: \cr
#'      \deqn{theta|gamma \sim H(gamma)}
#'      \deqn{x|theta \sim F(theta)}
#' MPE estimate of theta is theta_MPE = E(theta|gamma,x), E() is the expectation function.
#' For a given Bayesian bricks object obj, the MPE estimate will be:
#' \subsection{class(obj)="LinearGaussianGaussian"}{
#' Where
#'      \deqn{x \sim Gaussian(A z + b, Sigma)}
#'      \deqn{z \sim Gaussian(m,S)}
#' \code{MPE()} will return the MPE estimate of z.
#' See \code{?MPE.LinearGaussianGaussian} for details.
#' }
#' \subsection{class(obj)="GaussianGaussian"}{
#' Where
#' \deqn{x \sim Gaussian(mu,Sigma)}
#' \deqn{mu \sim Gaussian(m,S)}
#' Sigma is known.
#' \code{MPE()} will return the MPE estimate of mu.
#' See \code{?MPE.GaussianGaussian} for details.
#' }
#' \subsection{class(obj)="GaussianInvWishart"}{
#' Where
#' \deqn{x \sim Gaussian(mu,Sigma)}
#' \deqn{Sigma \sim InvWishart(v,S)}
#' mu is known.
#' \code{MPE()} will return the MPE estimate of Sigma.
#' See \code{?MPE.GaussianInvWishart} for details.
#' }
#' \subsection{class(obj)="GaussianNIW"}{
#' Where
#' \deqn{x \sim Gaussian(mu,Sigma)}
#' \deqn{Sigma \sim InvWishart(v,S)}
#' \deqn{mu \sim Gaussian(m,Sigma/k)}
#' \code{MPE()} will return the MPE estimate of mu and Sigma.
#' See \code{?MPE.GaussianNIW} for details.
#' }
#' \subsection{class(obj)="GaussianNIG"}{
#' Where
#' \deqn{x \sim Gaussian(X beta,sigma^2)}
#' \deqn{sigma^2 \sim InvGamma(a,b)}
#' \deqn{beta \sim Gaussian(m,sigma^2 V)}
#' X is a row vector, or a design matrix where each row is an obervation.
#' \code{MPE()} will return the MPE estimate of beta and sigma^2.
#' See \code{?MPE.GaussianNIG} for details.
#' }
#' \subsection{class(obj)="CatDirichlet"}{
#' Where
#' \deqn{x \sim Categorical(pi)}
#' \deqn{pi \sim Dirichlet(alpha)}
#' \code{MPE()} will return the MPE estimate of pi.
#' See \code{?MPE.CatDirichlet} for details.
#' }
#' \subsection{class(obj)="CatDP"}{
#' Where
#' \deqn{x \sim Categorical(pi)}
#' \deqn{pi \sim DirichletProcess(alpha)}
#' \code{MPE()} will return the MPE estimate of pi.
#' See \code{?MPE.CatDP} for details.
#' }
#' @seealso \code{\link{MPE.LinearGaussianGaussian}} for Linear Gaussian and Gaussian conjugate structure, \code{\link{MPE.GaussianGaussian}} for Gaussian-Gaussian conjugate structure, \code{\link{MPE.GaussianInvWishart}} for Gaussian-Inverse-Wishart conjugate structure, \code{\link{MPE.GaussianNIW}} for Gaussian-NIW conjugate structure, \code{\link{MPE.GaussianNIG}} for Gaussian-NIG conjugate structure, \code{\link{MPE.CatDirichlet}} for Categorical-Dirichlet conjugate structure, \code{\link{MPE.CatDP}} for Categorical-DP conjugate structure ...
#' @param obj A "BayesianBrick" object used to select a method.
#' @param ... further arguments passed to or from other methods.
#' @return A list of MPE estimates
#' @export
MPE <- function(obj,...) UseMethod("MPE")

#' @title Get the density from the posterior distribution.
#' @description
#' This is a generic function that will generate the the density value of the posterior distribution. i.e. for the model structure: \cr
#'      \deqn{theta|gamma \sim H(gamma)}
#'      \deqn{x|theta \sim F(theta)}
#' get the probability density/mass from the distribution \eqn{theta \sim H(gamma)}.
#' For a given Bayesian bricks object obj and an observation of theta, \code{dPosterior()} will calculate the density value for different model structures:
#' \subsection{class(obj)="LinearGaussianGaussian"}{
#' Where
#'      \deqn{x \sim Gaussian(A z + b, Sigma)}
#'      \deqn{z \sim Gaussian(m,S)}
#' \code{dPosterior()} will return p(theta|m,S)
#' See \code{?dPosterior.LinearGaussianGaussian} for details.
#' }
#' \subsection{class(obj)="GaussianGaussian"}{
#' Where
#' \deqn{x \sim Gaussian(mu,Sigma)}
#' \deqn{mu \sim Gaussian(m,S)}
#' Sigma is known.
#' \code{dPosterior()} will return p(mu|m,S)
#' See \code{?dPosterior.GaussianGaussian} for details.
#' }
#' \subsection{class(obj)="GaussianInvWishart"}{
#' Where
#' \deqn{x \sim Gaussian(mu,Sigma)}
#' \deqn{Sigma \sim InvWishart(v,S)}
#' mu is known.
#' \code{dPosterior()} will return p(Sigma|v,S)
#' See \code{?dPosterior.GaussianInvWishart} for details.
#' }
#' \subsection{class(obj)="GaussianNIW"}{
#' Where
#' \deqn{x \sim Gaussian(mu,Sigma)}
#' \deqn{Sigma \sim InvWishart(v,S)}
#' \deqn{mu \sim Gaussian(m,Sigma/k)}
#' \code{dPosterior()} will return p(mu,Sigma|m,k,v,S)
#' See \code{?dPosterior.GaussianNIW} for details.
#' }
#' \subsection{class(obj)="GaussianNIG"}{
#' Where
#' \deqn{x \sim Gaussian(X beta,sigma^2)}
#' \deqn{sigma^2 \sim InvGamma(a,b)}
#' \deqn{beta \sim Gaussian(m,sigma^2 V)}
#' X is a row vector, or a design matrix where each row is an obervation.
#' \code{dPosterior()} will return p(beta,sigma^2|m,V,a,b)
#' See \code{?dPosterior.GaussianNIG} for details.
#' }
#' \subsection{class(obj)="CatDirichlet"}{
#' Where
#' \deqn{x \sim Categorical(pi)}
#' \deqn{pi \sim Dirichlet(alpha)}
#' \code{dPosterior()} will return p(pi|alpha)
#' See \code{?dPosterior.CatDirichlet} for details.
#' }
#' @seealso \code{\link{dPosterior.LinearGaussianGaussian}} for Linear Gaussian and Gaussian conjugate structure, \code{\link{dPosterior.GaussianGaussian}} for Gaussian-Gaussian conjugate structure, \code{\link{dPosterior.GaussianInvWishart}} for Gaussian-Inverse-Wishart conjugate structure, \code{\link{dPosterior.GaussianNIW}} for Gaussian-NIW conjugate structure, \code{\link{dPosterior.GaussianNIG}} for Gaussian-NIG conjugate structure, \code{\link{dPosterior.CatDirichlet}} for Categorical-Dirichlet conjugate structure ...
#' @param obj A "BayesianBrick" object used to select a method.
#' @param ... further arguments passed to or from other methods.
#' @return numeric, the density value
#' @export
dPosterior <- function(obj,...) UseMethod("dPosterior")

#' @title Generate random samples from the posterior distribution
#' @description
#' This is a generic function that will generate random samples from the posterior distribution. i.e. for the model structure: \cr
#'      \deqn{theta|gamma \sim H(gamma)}
#'      \deqn{x|theta \sim F(theta)}
#' generate random sampels of theta from the distribution \eqn{theta \sim H(gamma)}.
#' For a given Bayesian bricks object obj, \code{rPosterior()} will generate random samples for different model structures:
#' \subsection{class(obj)="LinearGaussianGaussian"}{
#' Where
#'      \deqn{x \sim Gaussian(A z + b, Sigma)}
#'      \deqn{z \sim Gaussian(m,S)}
#' \code{rPosterior()} will generate random samples from Gaussian(m,S)
#' See \code{?rPosterior.LinearGaussianGaussian} for details.
#' }
#' \subsection{class(obj)="GaussianGaussian"}{
#' Where
#' \deqn{x \sim Gaussian(mu,Sigma)}
#' \deqn{mu \sim Gaussian(m,S)}
#' Sigma is known.
#' \code{rPosterior()} will generate random samples from Gaussian(m,S)
#' See \code{?rPosterior.GaussianGaussian} for details.
#' }
#' \subsection{class(obj)="GaussianInvWishart"}{
#' Where
#' \deqn{x \sim Gaussian(mu,Sigma)}
#' \deqn{Sigma \sim InvWishart(v,S)}
#' mu is known.
#' \code{rPosterior()} will generate random samples from InvWishart(v,S)
#' See \code{?rPosterior.GaussianInvWishart} for details.
#' }
#' \subsection{class(obj)="GaussianNIW"}{
#' Where
#' \deqn{x \sim Gaussian(mu,Sigma)}
#' \deqn{Sigma \sim InvWishart(v,S)}
#' \deqn{mu \sim Gaussian(m,Sigma/k)}
#' \code{rPosterior()} will generate random samples from NIW(m,k,v,S)
#' See \code{?rPosterior.GaussianNIW} for details.
#' }
#' \subsection{class(obj)="GaussianNIG"}{
#' Where
#' \deqn{x \sim Gaussian(X beta,sigma^2)}
#' \deqn{sigma^2 \sim InvGamma(a,b)}
#' \deqn{beta \sim Gaussian(m,sigma^2 V)}
#' X is a row vector, or a design matrix where each row is an obervation.
#' \code{rPosterior()} will generate random samples from NIG(m,V,a,b)
#' See \code{?rPosterior.GaussianNIG} for details.
#' }
#' \subsection{class(obj)="CatDirichlet"}{
#' Where
#' \deqn{x \sim Categorical(pi)}
#' \deqn{pi \sim Dirichlet(alpha)}
#' \code{rPosterior()} will generate random samples from Dirichlet(alpha)
#' See \code{?rPosterior.CatDirichlet} for details.
#' }
#' @seealso \code{\link{rPosterior.LinearGaussianGaussian}} for Linear Gaussian and Gaussian conjugate structure, \code{\link{rPosterior.GaussianGaussian}} for Gaussian-Gaussian conjugate structure, \code{\link{rPosterior.GaussianInvWishart}} for Gaussian-Inverse-Wishart conjugate structure, \code{\link{rPosterior.GaussianNIW}} for Gaussian-NIW conjugate structure, \code{\link{rPosterior.GaussianNIG}} for Gaussian-NIG conjugate structure, \code{\link{rPosterior.CatDirichlet}} for Categorical-Dirichlet conjugate structure ...
#' @param obj A "BayesianBrick" object used to select a method.
#' @param ... further arguments passed to or from other methods.
#' @return numeric, the density value
#' @export
rPosterior <- function(obj,...) UseMethod("rPosterior")

#' @title Get the marginal likelihood of a "BayesianBrick" object
#' @description
#' This is a generic function that will generate the marginal likelihood of a set of observations conditioned on a given "BayesianBrick" object. i.e. for the model structure: \cr
#'      \deqn{theta|gamma \sim H(gamma)}
#'      \deqn{x|theta \sim F(theta)}
#' Marginal likelihood is p(x|gamma), p() is the probability density/mass function for continuous/discrete x.
#' For a given Bayesian bricks object obj and a sample set x, \code{marginalLikelihood_bySufficientStatistics()} will calculate the marginal likelihood for different model structures:
#' \subsection{class(obj)="LinearGaussianGaussian"}{
#' Where
#'      \deqn{x \sim Gaussian(A z + b, Sigma)}
#'      \deqn{z \sim Gaussian(m,S)}
#' \code{marginalLikelihood_bySufficientStatistics()} will return p(x|m,S,A,b,Sigma)
#' See \code{?marginalLikelihood_bySufficientStatistics.LinearGaussianGaussian} for details.
#' }
#' \subsection{class(obj)="GaussianGaussian"}{
#' Where
#' \deqn{x \sim Gaussian(mu,Sigma)}
#' \deqn{mu \sim Gaussian(m,S)}
#' Sigma is known.
#' \code{marginalLikelihood_bySufficientStatistics()} will return p(x|m,S,Sigma)
#' See \code{?marginalLikelihood_bySufficientStatistics.GaussianGaussian} for details.
#' }
#' \subsection{class(obj)="GaussianInvWishart"}{
#' Where
#' \deqn{x \sim Gaussian(mu,Sigma)}
#' \deqn{Sigma \sim InvWishart(v,S)}
#' mu is known.
#' \code{marginalLikelihood_bySufficientStatistics()} will return p(x|mu,v,S)
#' See \code{?marginalLikelihood_bySufficientStatistics.GaussianInvWishart} for details.
#' }
#' \subsection{class(obj)="GaussianNIW"}{
#' Where
#' \deqn{x \sim Gaussian(mu,Sigma)}
#' \deqn{Sigma \sim InvWishart(v,S)}
#' \deqn{mu \sim Gaussian(m,Sigma/k)}
#' \code{marginalLikelihood_bySufficientStatistics()} will return p(x|m,k,v,S)
#' See \code{?marginalLikelihood_bySufficientStatistics.GaussianNIW} for details.
#' }
#' \subsection{class(obj)="GaussianNIG"}{
#' Where
#' \deqn{x \sim Gaussian(X beta,sigma^2)}
#' \deqn{sigma^2 \sim InvGamma(a,b)}
#' \deqn{beta \sim Gaussian(m,sigma^2 V)}
#' X is a row vector, or a design matrix where each row is an obervation.
#' \code{marginalLikelihood_bySufficientStatistics()} will return p(x,X|m,V,a,b)
#' See \code{?marginalLikelihood_bySufficientStatistics.GaussianNIG} for details.
#' }
#' \subsection{class(obj)="CatDirichlet"}{
#' Where
#' \deqn{x \sim Categorical(pi)}
#' \deqn{pi \sim Dirichlet(alpha)}
#' \code{marginalLikelihood_bySufficientStatistics()} will return p(x|alpha)
#' See \code{?marginalLikelihood_bySufficientStatistics.CatDirichlet} for details.
#' }
#' \subsection{class(obj)="CatDP"}{
#' Where
#' \deqn{x \sim Categorical(pi)}
#' \deqn{pi \sim DirichletProcess(alpha)}
#' \code{marginalLikelihood_bySufficientStatistics()} will return p(x|alpha)
#' See \code{?marginalLikelihood_bySufficientStatistics.CatDP} for details.
#' }
#' @seealso \code{\link{marginalLikelihood_bySufficientStatistics.LinearGaussianGaussian}} for Linear Gaussian and Gaussian conjugate structure, \code{\link{marginalLikelihood_bySufficientStatistics.GaussianGaussian}} for Gaussian-Gaussian conjugate structure, \code{\link{marginalLikelihood_bySufficientStatistics.GaussianInvWishart}} for Gaussian-Inverse-Wishart conjugate structure, \code{\link{marginalLikelihood_bySufficientStatistics.GaussianNIW}} for Gaussian-NIW conjugate structure, \code{\link{marginalLikelihood_bySufficientStatistics.GaussianNIG}} for Gaussian-NIG conjugate structure, \code{\link{marginalLikelihood_bySufficientStatistics.CatDirichlet}} for Categorical-Dirichlet conjugate structure, \code{\link{marginalLikelihood_bySufficientStatistics.CatDP}} for Categorical-DP conjugate structure ...
#' @param obj A "BayesianBrick" object used to select a method.
#' @param ... further arguments passed to or from other methods.
#' @return numeric, the marginal likelihood
#' @export
marginalLikelihood_bySufficientStatistics <- function(obj,...) UseMethod("marginalLikelihood_bySufficientStatistics")

#' @title Get the marginal likelihood of a "BayesianBrick" object
#' @description
#' This is a generic function that will generate the marginal likelihood of a set of observations conditioned on a given "BayesianBrick" object. i.e. for the model structure: \cr
#'      \deqn{theta|gamma \sim H(gamma)}
#'      \deqn{x|theta \sim F(theta)}
#' Marginal likelihood is p(x|gamma), p() is the probability density/mass function for continuous/discrete x.
#' For a given Bayesian bricks object obj and a sample set x, \code{marginalLikelihood()} will calculate the marginal likelihood for different model structures:
#' \subsection{class(obj)="LinearGaussianGaussian"}{
#' Where
#'      \deqn{x \sim Gaussian(A z + b, Sigma)}
#'      \deqn{z \sim Gaussian(m,S)}
#' \code{marginalLikelihood()} will return p(x|m,S,A,b,Sigma)
#' See \code{?marginalLikelihood.LinearGaussianGaussian} for details.
#' }
#' \subsection{class(obj)="GaussianGaussian"}{
#' Where
#' \deqn{x \sim Gaussian(mu,Sigma)}
#' \deqn{mu \sim Gaussian(m,S)}
#' Sigma is known.
#' \code{marginalLikelihood()} will return p(x|m,S,Sigma)
#' See \code{?marginalLikelihood.GaussianGaussian} for details.
#' }
#' \subsection{class(obj)="GaussianInvWishart"}{
#' Where
#' \deqn{x \sim Gaussian(mu,Sigma)}
#' \deqn{Sigma \sim InvWishart(v,S)}
#' mu is known.
#' \code{marginalLikelihood()} will return p(x|mu,v,S)
#' See \code{?marginalLikelihood.GaussianInvWishart} for details.
#' }
#' \subsection{class(obj)="GaussianNIW"}{
#' Where
#' \deqn{x \sim Gaussian(mu,Sigma)}
#' \deqn{Sigma \sim InvWishart(v,S)}
#' \deqn{mu \sim Gaussian(m,Sigma/k)}
#' \code{marginalLikelihood()} will return p(x|m,k,v,S)
#' See \code{?marginalLikelihood.GaussianNIW} for details.
#' }
#' \subsection{class(obj)="GaussianNIG"}{
#' Where
#' \deqn{x \sim Gaussian(X beta,sigma^2)}
#' \deqn{sigma^2 \sim InvGamma(a,b)}
#' \deqn{beta \sim Gaussian(m,sigma^2 V)}
#' X is a row vector, or a design matrix where each row is an obervation.
#' \code{marginalLikelihood()} will return p(x,X|m,V,a,b)
#' See \code{?marginalLikelihood.GaussianNIG} for details.
#' }
#' \subsection{class(obj)="CatDirichlet"}{
#' Where
#' \deqn{x \sim Categorical(pi)}
#' \deqn{pi \sim Dirichlet(alpha)}
#' \code{marginalLikelihood()} will return p(x|alpha)
#' See \code{?marginalLikelihood.CatDirichlet} for details.
#' }
#' \subsection{class(obj)="CatDP"}{
#' Where
#' \deqn{x \sim Categorical(pi)}
#' \deqn{pi \sim DirichletProcess(alpha)}
#' \code{marginalLikelihood()} will return p(x|alpha)
#' See \code{?marginalLikelihood.CatDP} for details.
#' }
#' @seealso \code{\link{marginalLikelihood.LinearGaussianGaussian}} for Linear Gaussian and Gaussian conjugate structure, \code{\link{marginalLikelihood.GaussianGaussian}} for Gaussian-Gaussian conjugate structure, \code{\link{marginalLikelihood.GaussianInvWishart}} for Gaussian-Inverse-Wishart conjugate structure, \code{\link{marginalLikelihood.GaussianNIW}} for Gaussian-NIW conjugate structure, \code{\link{marginalLikelihood.GaussianNIG}} for Gaussian-NIG conjugate structure, \code{\link{marginalLikelihood.CatDirichlet}} for Categorical-Dirichlet conjugate structure, \code{\link{marginalLikelihood.CatDP}} for Categorical-DP conjugate structure ...
#' @param obj A "BayesianBrick" object used to select a method.
#' @param ... further arguments passed to or from other methods.
#' @return numeric, the marginal likelihood
#' @export
marginalLikelihood <- function(obj,...) UseMethod("marginalLikelihood")

#' @title Get the density value of the posterior predictive distribution
#' @description
#' This is a generic function that will generate the the density value of the posterior predictive distribution. i.e. for the model structure: \cr
#'      \deqn{theta|gamma \sim H(gamma)}
#'      \deqn{x|theta \sim F(theta)}
#' get the probability density/mass of the posterior predictive distribution of a new sample x_new: p(x_new|gamma).
#' For a given Bayesian bricks object obj and a new sample x, \code{dPosteriorPredictive()} will calculate the marginal likelihood for different model structures:
#' \subsection{class(obj)="LinearGaussianGaussian"}{
#' Where
#'      \deqn{x \sim Gaussian(A z + b, Sigma)}
#'      \deqn{z \sim Gaussian(m,S)}
#' \code{dPosteriorPredictive()} will return p(x|m,S,A,b,Sigma)
#' See \code{?dPosteriorPredictive.LinearGaussianGaussian} for details.
#' }
#' \subsection{class(obj)="GaussianGaussian"}{
#' Where
#' \deqn{x \sim Gaussian(mu,Sigma)}
#' \deqn{mu \sim Gaussian(m,S)}
#' Sigma is known.
#' \code{dPosteriorPredictive()} will return p(x|m,S,Sigma)
#' See \code{?dPosteriorPredictive.GaussianGaussian} for details.
#' }
#' \subsection{class(obj)="GaussianInvWishart"}{
#' Where
#' \deqn{x \sim Gaussian(mu,Sigma)}
#' \deqn{Sigma \sim InvWishart(v,S)}
#' mu is known.
#' \code{dPosteriorPredictive()} will return p(x|mu,v,S)
#' See \code{?dPosteriorPredictive.GaussianInvWishart} for details.
#' }
#' \subsection{class(obj)="GaussianNIW"}{
#' Where
#' \deqn{x \sim Gaussian(mu,Sigma)}
#' \deqn{Sigma \sim InvWishart(v,S)}
#' \deqn{mu \sim Gaussian(m,Sigma/k)}
#' \code{dPosteriorPredictive()} will return p(x|m,k,v,S)
#' See \code{?dPosteriorPredictive.GaussianNIW} for details.
#' }
#' \subsection{class(obj)="GaussianNIG"}{
#' Where
#' \deqn{x \sim Gaussian(X beta,sigma^2)}
#' \deqn{sigma^2 \sim InvGamma(a,b)}
#' \deqn{beta \sim Gaussian(m,sigma^2 V)}
#' X is a row vector, or a design matrix where each row is an obervation.
#' \code{dPosteriorPredictive()} will return p(x,X|m,V,a,b)
#' See \code{?dPosteriorPredictive.GaussianNIG} for details.
#' }
#' \subsection{class(obj)="CatDirichlet"}{
#' Where
#' \deqn{x \sim Categorical(pi)}
#' \deqn{pi \sim Dirichlet(alpha)}
#' \code{dPosteriorPredictive()} will return p(x|alpha)
#' See \code{?dPosteriorPredictive.CatDirichlet} for details.
#' }
#' \subsection{class(obj)="CatDP"}{
#' Where
#' \deqn{x \sim Categorical(pi)}
#' \deqn{pi \sim DirichletProcess(alpha)}
#' \code{dPosteriorPredictive()} will return p(x|alpha)
#' See \code{?dPosteriorPredictive.CatDP} for details.
#' }
#' @seealso \code{\link{dPosteriorPredictive.LinearGaussianGaussian}} for Linear Gaussian and Gaussian conjugate structure, \code{\link{dPosteriorPredictive.GaussianGaussian}} for Gaussian-Gaussian conjugate structure, \code{\link{dPosteriorPredictive.GaussianInvWishart}} for Gaussian-Inverse-Wishart conjugate structure, \code{\link{dPosteriorPredictive.GaussianNIW}} for Gaussian-NIW conjugate structure, \code{\link{dPosteriorPredictive.GaussianNIG}} for Gaussian-NIG conjugate structure, \code{\link{dPosteriorPredictive.CatDirichlet}} for Categorical-Dirichlet conjugate structure, \code{\link{dPosteriorPredictive.CatDP}} for Categorical-DP conjugate structure ...
#' @param obj A "BayesianBrick" object used to select a method.
#' @param ... further arguments passed to or from other methods.
#' @return numeric, the density value
#' @export
dPosteriorPredictive <- function(obj,...) UseMethod("dPosteriorPredictive")

#' @title Generate random samples from the posterior predictive distribution
#' @description
#' This is a generic function that will generate random samples from the posterior predictive distribution. i.e. for the model structure: \cr
#'      \deqn{theta|gamma \sim H(gamma)}
#'      \deqn{x|theta \sim F(theta)}
#' generate x_new from the posterior predictive distribution of x|gamma.
#' For a given Bayesian bricks object obj, \code{rPosteriorPredictive()} will generate random samples from different model structures:
#' \subsection{class(obj)="LinearGaussianGaussian"}{
#' Where
#'      \deqn{x \sim Gaussian(A z + b, Sigma)}
#'      \deqn{z \sim Gaussian(m,S)}
#' \code{rPosteriorPredictive()} will generate samples from the distribution of x|m,S,A,b,Sigma
#' See \code{?rPosteriorPredictive.LinearGaussianGaussian} for details.
#' }
#' \subsection{class(obj)="GaussianGaussian"}{
#' Where
#' \deqn{x \sim Gaussian(mu,Sigma)}
#' \deqn{mu \sim Gaussian(m,S)}
#' Sigma is known.
#' \code{rPosteriorPredictive()} will generate samples from the distribution of x|m,S,Sigma
#' See \code{?rPosteriorPredictive.GaussianGaussian} for details.
#' }
#' \subsection{class(obj)="GaussianInvWishart"}{
#' Where
#' \deqn{x \sim Gaussian(mu,Sigma)}
#' \deqn{Sigma \sim InvWishart(v,S)}
#' mu is known.
#' \code{rPosteriorPredictive()} will generate samples from the distribution of x|mu,v,S
#' See \code{?rPosteriorPredictive.GaussianInvWishart} for details.
#' }
#' \subsection{class(obj)="GaussianNIW"}{
#' Where
#' \deqn{x \sim Gaussian(mu,Sigma)}
#' \deqn{Sigma \sim InvWishart(v,S)}
#' \deqn{mu \sim Gaussian(m,Sigma/k)}
#' \code{rPosteriorPredictive()} will generate samples from the distribution of x|m,k,v,S
#' See \code{?rPosteriorPredictive.GaussianNIW} for details.
#' }
#' \subsection{class(obj)="GaussianNIG"}{
#' Where
#' \deqn{x \sim Gaussian(X beta,sigma^2)}
#' \deqn{sigma^2 \sim InvGamma(a,b)}
#' \deqn{beta \sim Gaussian(m,sigma^2 V)}
#' X is a row vector, or a design matrix where each row is an obervation.
#' \code{rPosteriorPredictive()} will generate samples from the distribution of x,X|m,V,a,b
#' See \code{?rPosteriorPredictive.GaussianNIG} for details.
#' }
#' \subsection{class(obj)="CatDirichlet"}{
#' Where
#' \deqn{x \sim Categorical(pi)}
#' \deqn{pi \sim Dirichlet(alpha)}
#' \code{rPosteriorPredictive()} will generate samples from the distribution of x|alpha
#' See \code{?rPosteriorPredictive.CatDirichlet} for details.
#' }
#' \subsection{class(obj)="CatDP"}{
#' Where
#' \deqn{x \sim Categorical(pi)}
#' \deqn{pi \sim DirichletProcess(alpha)}
#' \code{rPosteriorPredictive()} will generate samples from the distribution of x|alpha
#' See \code{?rPosteriorPredictive.CatDP} for details.
#' }
#' @seealso \code{\link{rPosteriorPredictive.GaussianNIW}} for Gaussian-NIW conjugate structure, \code{\link{rPosteriorPredictive.GaussianNIG}} for Gaussian-NIG conjugate structure, \code{\link{rPosteriorPredictive.CatDirichlet}} for Categorical-Dirichlet conjugate structure, \code{\link{rPosteriorPredictive.CatDP}} for Categorical-DP conjugate structure ...
#' @param obj A "BayesianBrick" object used to select a method.
#' @param n integer, specify the number of samples to be generated.
#' @param ... further arguments passed to or from other methods.
#' @return a matrix or vector or list of random samples, depends on the type of 'obj'.
#' @export
rPosteriorPredictive <- function(obj,n,...) UseMethod("rPosteriorPredictive")



#' Print the content of an BasicBayesian object
#' @param x An object of the target type.
#' @param ... other parameters passed to print.
#' @return None.
#' @export
print.BasicBayesian <- function(x,...){
    cat("Prior distribution: ",x$H,"\n")
    cat("Observation distribution: ",x$F,"\n")        
    cat("Parameters for prior distribution:\n")
    print(x$gamma)
}

#' a internal version of "is", only for this package
#' @param object an object to be test with
#' @param class2 the class you want to test
#' @return logical, return TRUE if 'class2' is in 'class(object)'.
#' @noRd
.is <- function(object,class2){
    isTRUE(class2 %in% class(object))
}

#' log multivariate gamma function
#' @param a numeric, a numeric vector to be evaluated.
#' @param p numeric, the dimension.
#' @return numeric, the log value of the multivariate gamma function.
#' @noRd
.lmvgamma <- function(a,p){
    sapply(a,function(ai){
        p*(p-1)/4*log(pi)+ sum(lgamma(ai+(1-(1:p))/2))
    },simplify = TRUE)
}

#' determinant of a matrix
#' @param m numeric or matrix. When m is a numeric vector, it's treated as a diagonal matrix.
#' @return logical, return TRUE if 'class2' is in 'class(object)'.
#' @noRd
.det <- function(m){
    if(is.matrix(m)){
        base::det(m)
    }else if(is.vector(m)){
        base::prod(m)
    }
}

#' @title log sum exp
#' @description For each row l of a matrix x, calculate log(sum(exp(l))).
#' @param x matrix, the values in x a all logged. If x is a numeric vector, it will be converted to a matrix with 1 row.
#' @return numeric, the logsumexp of each row of x.
#' @export
#' @examples
#' \donttest{
#' ## Normalize the rows of x to make them sum up to 1
#' x <- matrix(runif(6,-1000,-20),nrow=2)
#' x <- x-logsumexp(x)
#' x <- exp(x)
#' rowSums(x)
#' }
logsumexp <- function(x){
    if(is.vector(x)) x <- matrix(x,nrow =1L)
    else if(!is.matrix(x)) x <- as.matrix(x)
    x <- as.matrix(x)
    a <- apply(x,1,max)
    log(rowSums(exp(x-a)))+a
}

#' @title Inverse of a positive definite symmetric matrix
#' @description Use Cholesky decomposition to calculate the inverse of S,  where \eqn{S = A'A}, A is a upper diagonal matrix. inv(S) = inv(A)inv(A)'.
#' @param S a symmetric positive definitive matrix.
#' @param returnUpper logical, return inv(A) if returnUpper=TRUE,return inv(S) if returnUpper=FALSE, default FALSE.
#' @return A matrix, the inverse of "S".
#' @export
#' @examples
#' Sigma = matrix(c(2,1,1,2),2,2)      # some positive definite symmetric matrix
#' pdsInverse(Sigma)                   # get inv(Sigma)
#' pdsInverse(Sigma,returnUpper=TRUE)  # get inv(A), where Sigma=A'A, A is upper triangle
pdsInverse <- function(S,returnUpper=FALSE){
    if(missing(S)) stop("'S' not specified!")
    if(length(S)==1L){
        if(returnUpper==FALSE)
            return(1/S)
        else
            return(sqrt(S))
    }else{
        A <- chol(S)                              #Cholesky decomposition S=A'A
        Ai <- backsolve(A,diag(nrow(S)))          #A_i = A^(-1)
        if(returnUpper==FALSE)
            return(crossprod(t(Ai)))
        else
            return(Ai)
    }
}

#' @title determinant of a positive definite symmetric matrix
#' @description Use Cholesky decomposition to calculate the determinant of S, where \eqn{S = A'A}, A is a upper diagonal matrix. det(S) = det(A)*det(A).
#' @param S a symmetric positive definitive matrix.
#' @param LOG logical, return log(det(S)) if TRUE, return det(S) if FALSE, default FALSE.
#' @return A matrix, the determinant of "S".
#' @export
#' @examples
#' Sigma = matrix(c(1.3,1,1,2),2,2)      # some positive definite symmetric matrix
#' pdsDeterminant(Sigma)                 # get inv(Sigma)
#' pdsDeterminant(Sigma,LOG=TRUE)        # get inv(Sigma)
#' pdsInverse(Sigma,returnUpper=TRUE)  # get inv(A), where Sigma=A'A, A is upper triangle
pdsDeterminant <- function(S,LOG=FALSE){
    if(missing(S)) stop("'S' not specified!")
    if(length(S)==1L){
        return(S)
    }else{
        A <- chol(S)                              #Cholesky decomposition S=A'A
        if(LOG){
            return(sum(log(diag(A))*2))
        }else{
            return(prod(diag(A))^2)
        }
    }
}

#' @title additional release questions
#' @description This is a list of additional questions you want devtools::release() to ask when releasing your package.
release_questions <- function(){
  c("Have you set the correct version number?",
    "Have you removed the irrelevant code blocks, such as 'if(FALSE){...}', from your R files?",
    "Have you add all the files in ./R folder to DESCRIPTION?",
    "Have you removed the unfinished examples from the vignette?",
    "Have you add all the references to README and function documentation?",
    "Have you replaced all of the dontrun{} with donttest{}?",
    "Have you rebuilt README.md from README.raw.md, using publish_markdown.r? Note that you have to rebuild README.md whenever you modified README.raw.md.",
    "Have you removed all the par() settings in the examples? Check your documentations by run \"grep 'par(' ./man/*\" from the terminal."
    )
}
