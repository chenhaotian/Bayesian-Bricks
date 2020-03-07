


#' Create objects of type '"BasicBayesian"'.
#'
#' A Basic Bayesian Object is with following conditional dependency structure: \cr
#'      theta|gamma ~ H(gamma) \cr
#'      X|theta ~ F(theta) \cr
#' Where H(gamma) is usually called "the prior distribution", F(theta) is called "the observation distribution". Objects of type "GaussianNIW", "GaussianNIG", "CatDirichlet", and "CatDP" are all "BasicBayesian"s.
#' 
#' @seealso \code{\link{GaussianNIW}} for Gaussian-NIW conjugate structure, \code{\link{GaussianNIG}} for Gaussian-NIG conjugate structure, \code{\link{CatDirichlet}} for Categorical-Dirichlet conjugate structure, \code{\link{CatDP}} for Categorical-DP conjugate structure ...
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
#' This is a generic function that will generate the sufficient statistics of a given Bayesian bricks object. i.e.
#' for the model structure: \cr
#'      theta|gamma ~ H(gamma) \cr
#'      X|theta ~ F(theta) \cr
#' get the sufficient statistics T(X)
#'
#' @seealso \code{\link{sufficientStatistics.GaussianNIW}} for Gaussian-NIW conjugate structure, \code{\link{sufficientStatistics.GaussianNIG}} for Gaussian-NIG conjugate structure, \code{\link{sufficientStatistics.CatDirichlet}} for Categorical-Dirichlet conjugate structure, \code{\link{sufficientStatistics.CatDP}} for Categorical-DP conjugate structure ...
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
#' This is a generic function that will generate the weighted sufficient statistics of a given "BayesianBrick" object. That is, for the model structure: \cr
#'      theta|gamma ~ H(gamma) \cr
#'      X|theta ~ F(theta) \cr
#' get the weighted sufficient statistics T(X)
#'
#' @seealso \code{\link{sufficientStatistics_Weighted.GaussianNIW}} for Gaussian-NIW conjugate structure, \code{\link{sufficientStatistics_Weighted.GaussianNIG}} for Gaussian-NIG conjugate structure, \code{\link{sufficientStatistics_Weighted.CatDirichlet}} for Categorical-Dirichlet conjugate structure, \code{\link{sufficientStatistics_Weighted.CatDP}} for Categorical-DP conjugate structure ...
#' @param obj A "BayesianBrick" object used to select a method.
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
sufficientStatistics_Weighted <- function(obj,...) UseMethod("sufficientStatistics_Weighted")

#' update the prior distribution with sufficient statistics
#' @param obj A "BayesianBrick" object used to select a method.
#' @param ... further arguments passed to or from other methods.
#' @return None, or an error message if the update fails.
posterior_bySufficientStatistics <- function(obj,...) UseMethod("posterior_bySufficientStatistics")

#' update the prior distribution with sufficient statistics
#'
#' This is a generic function that will update the prior distribution of a "BayesianBrick" object by adding information of the observation's sufficient statistics. i.e. for the model structure: \cr
#'      theta|gamma ~ H(gamma) \cr
#'      X|theta ~ F(theta) \cr
#' update gamma to gamma_posterior by adding the information of X to gamma
#'
#' @seealso \code{\link{posterior.GaussianNIW}} for Gaussian-NIW conjugate structure, \code{\link{posterior.GaussianNIG}} for Gaussian-NIG conjugate structure, \code{\link{posterior.CatDirichlet}} for Categorical-Dirichlet conjugate structure, \code{\link{posterior.CatDP}} for Categorical-DP conjugate structure ...
#' @param obj A "BayesianBrick" object used to select a method.
#' @param ... further arguments passed to or from other methods.
#' @return None, or an error message if the update fails.
#' @export
posterior <- function(obj,...) UseMethod("posterior")

#' update the prior distribution with sufficient statistics
#' @param obj A "BayesianBrick" object used to select a method.
#' @param ... other parameters.
#' @return None, or an error message if the update fails.
posteriorDiscard_bySufficientStatistics <- function(obj,...) UseMethod("posteriorDiscard_bySufficientStatistics")

#' update the prior distribution with sufficient statistics
#'
#' Contrary to posterior(), posteriorDiscard() a generic function that will update the prior distribution of a "BayesianBrick" object by removing the information provided by the observation's sufficient statistics. i.e. for the model structure: \cr
#'      theta|gamma ~ H(gamma) \cr
#'      X|theta ~ F(theta) \cr
#' update gamma to gamma_posterior by removing the information of X from gamma
#'
#' @seealso \code{\link{posterior.GaussianNIW}} for Gaussian-NIW conjugate structure, \code{\link{posterior.GaussianNIG}} for Gaussian-NIG conjugate structure, \code{\link{posterior.CatDirichlet}} for Categorical-Dirichlet conjugate structure, \code{\link{posterior.CatDP}} for Categorical-DP conjugate structure ...
#' @param obj A "BayesianBrick" object used to select a method.
#' @param ... further arguments passed to or from other methods.
#' @return None, or an error message if the update fails.
#' @export
posteriorDiscard <- function(obj,...) UseMethod("posteriorDiscard")

#' Get the Maximum A Posteriori(MAP) estimate of a "BayesianBrick" object
#'
#' This is a generic function that will generate the MAP estimate of a given "BayesianBrick" object. For the model structure:  \cr
#'      theta|gamma ~ H(gamma) \cr
#'      X|theta ~ F(theta) \cr
#' MAP esitmate of theta is theta_MAP = argmax_theta p(theta|gamma,x).
#'
#' @seealso \code{\link{MAP.GaussianNIW}} for Gaussian-NIW conjugate structure, \code{\link{MAP.GaussianNIG}} for Gaussian-NIG conjugate structure, \code{\link{MAP.CatDirichlet}} for Categorical-Dirichlet conjugate structure, \code{\link{MAP.CatDP}} for Categorical-DP conjugate structure ...
#' @param obj A "BayesianBrick" object used to select a method.
#' @param ... further arguments passed to or from other methods.
#' @return A list of the MAP estimates
#' @export
MAP <- function(obj,...) UseMethod("MAP")

#' Get the Mean Posterior Estimate(MPE) of a "BayesianBrick" object
#'
#' This is a generic function that will generate the MPE estimate of a given "BayesianBrick" object. i.e. for the model structure: \cr
#'      theta|gamma ~ H(gamma) \cr
#'      X|theta ~ F(theta) \cr
#' MPE esitmate of theta is theta_MAP = E(theta|gamma,x), E() is the expectation function.
#'
#' @seealso \code{\link{MPE.GaussianNIW}} for Gaussian-NIW conjugate structure, \code{\link{MPE.GaussianNIG}} for Gaussian-NIG conjugate structure, \code{\link{MPE.CatDirichlet}} for Categorical-Dirichlet conjugate structure, \code{\link{MPE.CatDP}} for Categorical-DP conjugate structure ...
#' @param obj A "BayesianBrick" object used to select a method.
#' @param ... further arguments passed to or from other methods.
#' @return A list of MPE estimates
#' @export
MPE <- function(obj,...) UseMethod("MPE")

#' Get the marginal likelihood of a "BayesianBrick" object
#'
#' This is a generic function that will generate the marginal likelihood of a set of observations conditioned on a given "BayesianBrick" object. i.e. for the model structure: \cr
#'      theta|gamma ~ H(gamma) \cr
#'      X|theta ~ F(theta) \cr
#' Marginal likelihood is p(X|gamma), p() is the probability density/mass function for contunious/discrete X.
#'
#' @seealso \code{\link{marginalLikelihood_bySufficientStatistics.GaussianNIW}} for Gaussian-NIW conjugate structure, \code{\link{marginalLikelihood_bySufficientStatistics.GaussianNIG}} for Gaussian-NIG conjugate structure, \code{\link{marginalLikelihood_bySufficientStatistics.CatDirichlet}} for Categorical-Dirichlet conjugate structure, \code{\link{marginalLikelihood_bySufficientStatistics.CatDP}} for Categorical-DP conjugate structure ...
#' @param obj A "BayesianBrick" object used to select a method.
#' @param ... further arguments passed to or from other methods.
#' @return numeric, the marginal likelihood
#' @export
marginalLikelihood_bySufficientStatistics <- function(obj,...) UseMethod("marginalLikelihood_bySufficientStatistics")

#' Get the marginal likelihood of a "BayesianBrick" object
#'
#' This is a generic function that will generate the marginal likelihood of a set of observations conditioned on a given "BayesianBrick" object. i.e. for the model structure: \cr
#'      theta|gamma ~ H(gamma) \cr
#'      X|theta ~ F(theta) \cr
#' Marginal likelihood is p(X|gamma), p() is the probability density/mass function for contunious/discrete X.
#'
#' @seealso \code{\link{marginalLikelihood.GaussianNIW}} for Gaussian-NIW conjugate structure, \code{\link{marginalLikelihood.GaussianNIG}} for Gaussian-NIG conjugate structure, \code{\link{marginalLikelihood.CatDirichlet}} for Categorical-Dirichlet conjugate structure, \code{\link{marginalLikelihood.CatDP}} for Categorical-DP conjugate structure ...
#' @param obj A "BayesianBrick" object used to select a method.
#' @param ... further arguments passed to or from other methods.
#' @return numeric, the marginal likelihood
#' @export
marginalLikelihood <- function(obj,...) UseMethod("marginalLikelihood")

#' Get the density value of the posterior predictive distribution
#'
#' This is a generic function that will generate the the density value of the posterior predictive distribution. i.e. for the model structure: \cr
#'      theta|gamma ~ H(gamma) \cr
#'      X|theta ~ F(theta) \cr
#' get the probability density/mass of the posterior predictive distribution of a new sample x_new: p(x_new|gamma).
#'
#' @seealso \code{\link{dPosteriorPredictive.GaussianNIW}} for Gaussian-NIW conjugate structure, \code{\link{dPosteriorPredictive.GaussianNIG}} for Gaussian-NIG conjugate structure, \code{\link{dPosteriorPredictive.CatDirichlet}} for Categorical-Dirichlet conjugate structure, \code{\link{dPosteriorPredictive.CatDP}} for Categorical-DP conjugate structure ...
#' @param obj A "BayesianBrick" object used to select a method.
#' @param ... further arguments passed to or from other methods.
#' @return numeric, the density value
#' @export
dPosteriorPredictive <- function(obj,...) UseMethod("dPosteriorPredictive")

#' Generate random samples from the posterior predictive distribution
#'
#' This is a generic function that will generate random samples from the posterior predictive distribution. i.e. for the model structure: \cr
#'      theta|gamma ~ H(gamma) \cr
#'      X|theta ~ F(theta) \cr
#' generate x_new from the posterior predictive distribution of x|gamma.
#'
#' @seealso \code{\link{rPosteriorPredictive.GaussianNIW}} for Gaussian-NIW conjugate structure, \code{\link{rPosteriorPredictive.GaussianNIG}} for Gaussian-NIG conjugate structure, \code{\link{rPosteriorPredictive.CatDirichlet}} for Categorical-Dirichlet conjugate structure, \code{\link{rPosteriorPredictive.CatDP}} for Categorical-DP conjugate structure ...
#' @param obj A "BayesianBrick" object used to select a method.
#' @param ... further arguments passed to or from other methods.
#' @return a matrix or vector of random samples
#' @export
rPosteriorPredictive <- function(obj,...) UseMethod("rPosteriorPredictive")



#' Print the content of an BasicBasyesian object
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

#' @title additional release questions
#' @description This is a list of additional questions you want devtools::release() to ask when releasing your package.
release_questions <- function() {
  c("Have you set the correct version number?",
    "Have you removed the irrelevant code blocks?",
    "Have you add all R files to DESCRIPTION?",
    "Have you removed the unfinished lines from vignette?",
    "Have you add all the references to README and function documentation?",
    "Have you add all the references?",
    "Have you built README.md from README.raw.md?")
}
