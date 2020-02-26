


#' Create objects of type '"BasicBayesian"'.
#'
#' A Basic Bayesian Object is with following conditional dependency structure: \cr
#'      theta|gamma ~ H(gamma) \cr
#'      X|theta ~ F(theta) \cr
#' Where H(gamma) is usually called "the prior distribution", F(theta) is called "the observation distribution". Objects of type "GaussianNIW", "GaussianNIG", "CatDirichlet", and "CatDP" are all "BasicBayesian"s.
#' 
#' @seealso \code{\link{GaussianNIW}} for Gaussian-NIW conjugate structure, \code{\link{GaussianNIG}} for Gaussian-NIG conjugate structure, \code{\link{CatDirichlet}} for Categorica-Dirichlet conjugate structure, \code{\link{CatDP}} for Categorica-DP conjugate structure ...
#' @param ENV The environment where you want to create the BasicBayesian object
#' @return An object of class "BasicBayesian"
#' @export
BasicBayesian <- function(ENV=parent.frame()){
    object <- new.env(parent = ENV)
    
    object$gamma <- list()
    object$H <- character(0)
    object$F <- character(0)

    class(object) <- c("BasicBayesian","BayesianBrick")
    return(object)

}

#' Get sample sufficient statistics
#'
#' This is a generatic function that will generate the sufficient statistics of a given bayesian bricks object. i.e.
#' for the model structure: \cr
#'      theta|gamma ~ H(gamma) \cr
#'      X|theta ~ F(theta) \cr
#' get the sufficient statistics T(X)
#'
#' @seealso \code{\link{sufficientStatistics.GaussianNIW}} for Gaussian-NIW conjugate structure, \code{\link{sufficientStatistics.GaussianNIG}} for Gaussian-NIG conjugate structure, \code{\link{sufficientStatistics.CatDirichlet}} for Categorica-Dirichlet conjugate structure, \code{\link{sufficientStatistics.CatDP}} for Categorica-DP conjugate structure ...
#' @param obj A "BayesianBrick" object used to select a method.
#' @param ... further arguments passed to or from other methods.
#' @return An object of corresponding sufficient statistics class, such as "ssGaussian"
#' @export
#' @examples
#' x <- rGaussian(10,mu = 1,Sigma = 1)
#' obj <- GaussianNIW()                    #an GaussianNIW object
#' sufficientStatistics(obj=obj,x=x)
sufficientStatistics <- function(obj,...) UseMethod("sufficientStatistics")

#' Get weighted sample sufficient statistics
#'
#' This is a generatic function that will generate the weighted sufficient statistics of a given "BayesianBrick" object. That is, for the model structure: \cr
#'      theta|gamma ~ H(gamma) \cr
#'      X|theta ~ F(theta) \cr
#' get the weigeted sufficient statistics T(X)
#'
#' @seealso \code{\link{sufficientStatistics_Weighted.GaussianNIW}} for Gaussian-NIW conjugate structure, \code{\link{sufficientStatistics_Weighted.GaussianNIG}} for Gaussian-NIG conjugate structure, \code{\link{sufficientStatistics_Weighted.CatDirichlet}} for Categorica-Dirichlet conjugate structure, \code{\link{sufficientStatistics_Weighted.CatDP}} for Categorica-DP conjugate structure ...
#' @param obj A "BayesianBrick" object used to select a method.
#' @param ... further arguments passed to or from other methods.
#' @return An object of corresponding sufficient statistics class, such as "ssGaussian"
#' @export
#' @examples
#' \dontrun{
#' x <- rGaussian(10,mu = 1,Sigma = 1)
#' w <- runif(10)
#' obj <- GaussianNIW()                    #an GaussianNIW object
#' sufficientStatistics_Weighted(obj=obj,x=x,w=w)
#' }
sufficientStatistics_Weighted <- function(obj,...) UseMethod("sufficientStatistics_Weighted")

#' update the prior distribution with sufficient statistics
#' @param obj A "BayesianBrick" object used to select a method.
#' @param ... further arguments passed to or from other methods.
posterior_bySufficientStatistics <- function(obj,...) UseMethod("posterior_bySufficientStatistics")

#' update the prior distribution with sufficient statistics
#'
#' This is a generatic function that will update the prior distribution of a "BayesianBrick" object by adding information of the observation's sufficient statistics. i.e. for the model structure: \cr
#'      theta|gamma ~ H(gamma) \cr
#'      X|theta ~ F(theta) \cr
#' update gamma to gamma_posterior by adding the information of X to gamma
#'
#' @seealso \code{\link{posterior.GaussianNIW}} for Gaussian-NIW conjugate structure, \code{\link{posterior.GaussianNIG}} for Gaussian-NIG conjugate structure, \code{\link{posterior.CatDirichlet}} for Categorica-Dirichlet conjugate structure, \code{\link{posterior.CatDP}} for Categorica-DP conjugate structure ...
#' @param obj A "BayesianBrick" object used to select a method.
#' @param ... further arguments passed to or from other methods.
#' @return None, or an error message if the update fails.
#' @export
posterior <- function(obj,...) UseMethod("posterior")

#' update the prior distribution with sufficient statistics
#' @param obj A "BayesianBrick" object used to select a method.
#' @param ... other parameters.
posteriorDiscard_bySufficientStatistics <- function(obj,...) UseMethod("posteriorDiscard_bySufficientStatistics")

#' update the prior distribution with sufficient statistics
#'
#' Contrary to posterior(), posteriorDiscard() a generatic function that will update the prior distribution of a "BayesianBrick" object by removing the information provided by the observation's sufficient statistics. i.e. for the model structure: \cr
#'      theta|gamma ~ H(gamma) \cr
#'      X|theta ~ F(theta) \cr
#' update gamma to gamma_posterior by removing the information of X from gamma
#'
#' @seealso \code{\link{posterior.GaussianNIW}} for Gaussian-NIW conjugate structure, \code{\link{posterior.GaussianNIG}} for Gaussian-NIG conjugate structure, \code{\link{posterior.CatDirichlet}} for Categorica-Dirichlet conjugate structure, \code{\link{posterior.CatDP}} for Categorica-DP conjugate structure ...
#' @param obj A "BayesianBrick" object used to select a method.
#' @param ... further arguments passed to or from other methods.
#' @return None, or an error message if the update fails.
#' @export
posteriorDiscard <- function(obj,...) UseMethod("posteriorDiscard")

#' Get the MAP estimate of a "BayesianBrick" object
#'
#' This is a generatic function that will generate the MAP estimate of a given "BayesianBrick" object. For the model structure:  \cr
#'      theta|gamma ~ H(gamma) \cr
#'      X|theta ~ F(theta) \cr
#' get theta_MAP
#'
#' @seealso \code{\link{MAP.GaussianNIW}} for Gaussian-NIW conjugate structure, \code{\link{MAP.GaussianNIG}} for Gaussian-NIG conjugate structure, \code{\link{MAP.CatDirichlet}} for Categorica-Dirichlet conjugate structure, \code{\link{MAP.CatDP}} for Categorica-DP conjugate structure ...
#' @param obj A "BayesianBrick" object used to select a method.
#' @param ... further arguments passed to or from other methods.
#' @return A list of the MAP estimates
#' @export
MAP <- function(obj,...) UseMethod("MAP")

#' Get the MPE estimate of a "BayesianBrick" object
#'
#' This is a generatic function that will generate the MPE estimate of a given "BayesianBrick" object. i.e. for the model structure: \cr
#'      theta|gamma ~ H(gamma) \cr
#'      X|theta ~ F(theta) \cr
#' get theta_MPE
#'
#' @seealso \code{\link{MPE.GaussianNIW}} for Gaussian-NIW conjugate structure, \code{\link{MPE.GaussianNIG}} for Gaussian-NIG conjugate structure, \code{\link{MPE.CatDirichlet}} for Categorica-Dirichlet conjugate structure, \code{\link{MPE.CatDP}} for Categorica-DP conjugate structure ...
#' @param obj A "BayesianBrick" object used to select a method.
#' @param ... further arguments passed to or from other methods.
#' @return A list of MPE estimates
#' @export
MPE <- function(obj,...) UseMethod("MPE")

#' Get the marginal likelihood of a "BayesianBrick" object
#'
#' This is a generatic function that will generate the marginal likelihood of a set of observations conditioned on a given "BayesianBrick" object. i.e. for the model structure: \cr
#'      theta|gamma ~ H(gamma) \cr
#'      X|theta ~ F(theta) \cr
#' get p(X|gamma)
#'
#' @seealso \code{\link{marginalLikelihood_bySufficientStatistics.GaussianNIW}} for Gaussian-NIW conjugate structure, \code{\link{marginalLikelihood_bySufficientStatistics.GaussianNIG}} for Gaussian-NIG conjugate structure, \code{\link{marginalLikelihood_bySufficientStatistics.CatDirichlet}} for Categorica-Dirichlet conjugate structure, \code{\link{marginalLikelihood_bySufficientStatistics.CatDP}} for Categorica-DP conjugate structure ...
#' @param obj A "BayesianBrick" object used to select a method.
#' @param ... further arguments passed to or from other methods.
#' @return numeric, the marginal likelihood
#' @export
marginalLikelihood_bySufficientStatistics <- function(obj,...) UseMethod("marginalLikelihood_bySufficientStatistics")

#' Get the marginal likelihood of a "BayesianBrick" object
#'
#' This is a generatic function that will generate the marginal likelihood of a set of observations conditioned on a given "BayesianBrick" object. i.e. for the model structure: \cr
#'      theta|gamma ~ H(gamma) \cr
#'      X|theta ~ F(theta) \cr
#' get p(X|gamma)
#'
#' @seealso \code{\link{marginalLikelihood.GaussianNIW}} for Gaussian-NIW conjugate structure, \code{\link{marginalLikelihood.GaussianNIG}} for Gaussian-NIG conjugate structure, \code{\link{marginalLikelihood.CatDirichlet}} for Categorica-Dirichlet conjugate structure, \code{\link{marginalLikelihood.CatDP}} for Categorica-DP conjugate structure ...
#' @param obj A "BayesianBrick" object used to select a method.
#' @param ... further arguments passed to or from other methods.
#' @return numeric, the marginal likelihood
#' @export
marginalLikelihood <- function(obj,...) UseMethod("marginalLikelihood")

#' Get the density value of the posterior predictive distribution
#'
#' This is a generatic function that will generate the the density value of the posterior predictive distribution. i.e. for the model structure: \cr
#'      theta|gamma ~ H(gamma) \cr
#'      X|theta ~ F(theta) \cr
#' get p(x_new|gamma)
#'
#' @seealso \code{\link{dPosteriorPredictive.GaussianNIW}} for Gaussian-NIW conjugate structure, \code{\link{dPosteriorPredictive.GaussianNIG}} for Gaussian-NIG conjugate structure, \code{\link{dPosteriorPredictive.CatDirichlet}} for Categorica-Dirichlet conjugate structure, \code{\link{dPosteriorPredictive.CatDP}} for Categorica-DP conjugate structure ...
#' @param obj A "BayesianBrick" object used to select a method.
#' @param ... further arguments passed to or from other methods.
#' @return numeric, the density value
#' @export
dPosteriorPredictive <- function(obj,...) UseMethod("dPosteriorPredictive")

#' Generate random samples from the posterior predictive distribution
#'
#' This is a generatic function that will generate random samples from the posterior predictive distribution. i.e. for the model structure: \cr
#'      theta|gamma ~ H(gamma) \cr
#'      X|theta ~ F(theta) \cr
#' generate x_new from G(x|gamma), G() is the posterior predictive distribution
#'
#' @seealso \code{\link{rPosteriorPredictive.GaussianNIW}} for Gaussian-NIW conjugate structure, \code{\link{rPosteriorPredictive.GaussianNIG}} for Gaussian-NIG conjugate structure, \code{\link{rPosteriorPredictive.CatDirichlet}} for Categorica-Dirichlet conjugate structure, \code{\link{rPosteriorPredictive.CatDP}} for Categorica-DP conjugate structure ...
#' @param obj A "BayesianBrick" object used to select a method.
#' @param ... further arguments passed to or from other methods.
#' @return a matrix or vector of random samples
#' @export
rPosteriorPredictive <- function(obj,...) UseMethod("rPosteriorPredictive")

#' the density of observation distribution
#' @param obj An object of the target type.
#' @param ... other parameters.
dObservationDistribution <- function(obj,...) UseMethod("dObservationDistribution")

#' random number generation of the observation distribution
#' @param obj An object of the target type.
#' @param ... other parameters.
rObservationDistribution <- function(obj,...) UseMethod("rObservationDistribution")

#' Print the content of an BayesicBasyeisn object
#' @param x An object of the target type.
#' @param ... other parameters passed to print.
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
.is <- function(object,class2){
    isTRUE(class2 %in% class(object))
}

#' log sum exp
#' 
#' For each row l of a matrix x, calculate log(sum(exp(l))).
#' 
#' @param x matrix, the values in x a all logged. If x is a numeric vector, it will be converted to a matrix with 1 row.
#' @return numeric, the logsumexp of each row of x.
#' @export
#' @examples
#' \dontrun{
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

