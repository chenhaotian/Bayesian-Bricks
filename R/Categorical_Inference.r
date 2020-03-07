
#' @include Bayesian_Bricks.r



countFreq <- function(x,uniqx=NULL){
    if(!is.factor(x)){
        if(!is.null(uniqx)) x <- factor(x,levels=uniqx)
        else x <- factor(x)
    }
    as.numeric(tabulate(x,nbins = length(levels(x))))
}
countFreq_Weighted <- function(x,uniqx=NULL,w){
    if(length(x)!=length(w)) stop("length of x and w don't match")
    if(!is.factor(x)){
        if(!is.null(uniqx)) x <- factor(x,levels=uniqx)
        else x <- factor(x)
    }
    vapply(split(x=w,f=x),sum,FUN.VALUE = numeric(1),USE.NAMES = FALSE)
}

#' Random generation for Categorical distribution
#'
#' Generate random integer samples from a Categorical distribution.
#'
#' @seealso \code{\link{dCategorical}}
#' @param n integer, number of samples.
#' @param p numeric, probabilities.
#' @return An integer vector of length n.
#' @export
#' @examples
#' rCategorical(n=20,p=c(1,2))
rCategorical <- function(n,p){
    sample.int(n=length(p), size = n, replace = TRUE, prob = p)
}

#' Probability mass function for Categorical distribution
#'
#' Calculate probability masses for integer valued Categorical random samples.
#'
#' @seealso \code{\link{rCategorical}}
#' @param x integer, categorical samples.
#' @param p numeric, probabilities.
#' @return A numeric vector of the same length of 'x'.
#' @export
#' @examples
#' \donttest{
#' dCategorical(x=c(1L,2L,1L),p=c(1,2))
#' }
dCategorical <- function(x,p){
    p[x]
}

#' Random generation for Dirichelt distribution
#'
#' Generate random samples from Dirichlet distribution.
#'
#' @seealso \code{\link{dDir}}
#' @param n integer, number of samples.
#' @param alpha numeric, Dirichlet parameter.
#' @return A numeric matrix of n rows and length(alpha) columns.
#' @export
#' @examples
#' rDir(5,c(1,2,3)) #generate 5 samples with parameters c(1,2,3)
#' @import stats
rDir <- function(n, alpha){
    l <- length(alpha)
    x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
    sm <- x %*% rep(1, l)
    x/as.vector(sm)
}

#' Density function for Dirichelt distribution
#'
#' Calculate the densities of a given set of Dirichlet samples.
#'
#' @seealso \code{\link{rDir}}
#' @param x matrix or numeric vector, if matrix every row of x is an observation, if numeric vector, it's the same as a matrix with only one row.
#' @param alpha numeric, Dirichlet parameter.
#' @param LOG logical, return the log density if set to "TRUE".
#' @return A numeric vector of density values.
#' @export
#' @examples
#' x <- rDir(5,c(1,2,3)) #generate 5 samples with parameters c(1,2,3)
#' dDir(x,c(1,2,3))
#' dDir(x,c(1,2,3),LOG=TRUE)
dDir <- function(x, alpha,LOG=FALSE){
    if(!is.matrix(x))
        if(is.numeric(x))
            x <- t(x)
        else
            stop("Error in ddir():x must be a matrix or an integer/numeric vector")
    if(length(alpha)==1) alpha <- rep(alpha,ncol(x))
    if(ncol(x)!=length(alpha))
        stop("Error in ddir(): columns of x and length of alpha don't match")
    out <- apply(x,1,function(l){
        s <- (alpha - 1) * log(l)
        s <- ifelse(alpha == 1 & l == 0, -Inf, s)
        sum(s)
    })
    logD <- sum(lgamma(alpha)) - lgamma(sum(alpha))
    out <- out- logD
    if(!LOG) out <- exp(out)
    return(out)
}

#' Create objects of type "CatDirichlet".
#'
#' Create an object of type "CatDirichlet", which represents the Categorical-Dirichlet(Multinomial-Dirichlet) conjugate structure: \cr
#'      pi|alpha ~ Dir(alpha) \cr
#'      x|pi ~ Categorical(pi) \cr
#' This object will be used as a place for recording and accumulating information in the related inference/sampling functions such as posterior(), posteriorDiscard(), MAP() and so on. A categorical distribution is defined on a set of unique labels, usually these labels are integers, they can also be characters and factors.
#'
#' @seealso \code{\link{posterior.CatDirichlet}},\code{\link{posteriorDiscard.CatDirichlet}},\code{\link{MAP.CatDirichlet}},\code{\link{marginalLikelihood.CatDirichlet}} ...
#' @param objCopy an object of type "CatDirichlet". If "objCopy" is not NULL, the function create a new "CatDirichlet" object by copying the content from objCopy, otherwise this new object will be created by using "ENV" and "gamma". Default NULL.
#' @param ENV environment, specify where the object will be created.
#' @param gamma list, a named list of parameters, gamma=list(alpha,uniqueLabels). Where gamma$alpha is a numeric vector specifying the parameters of the Dirichlet distribution, gamma$uniqueLables is a integer/character vector specifying the unique category labels of the Categorical distribution.
#' @return An object of class "CatDirichlet".
#' @export
#' @examples
#' obj <- CatDirichlet(gamma=list(alpha=c(1,2,1),uniqueLabels = letters[1:3]))
#' obj #print the content
#' @references Murphy, Kevin P. Machine learning: a probabilistic perspective. MIT press, 2012.
CatDirichlet <- function(objCopy=NULL,ENV=parent.frame(),gamma=list(alpha=1,uniqueLabels=1L)){
    object <- BasicBayesian(ENV = ENV)

    if(!is.null(objCopy)){
        if(!.is(objCopy,"CatDirichlet")) stop("'objCopy' must be of class 'CatDirichlet'")
        object$gamma <- objCopy$gamma
        object$H <- objCopy$H
        object$F <- objCopy$F
    }else{
        if(!missing(gamma))
            if((!is.list(gamma)) |
               (!all(names(gamma) %in% c("alpha","uniqueLabels"))))
                stop("gamma must be list(alpha,uniqueLabels)")
        if(length(gamma$alpha)==1L) gamma$alpha <- rep(gamma$alpha,length(gamma$uniqueLabels))
        if(length(gamma$alpha)!=length(gamma$uniqueLabels)) stop("length of 'alpha' and 'uniqueLabels' in 'gamma' don't match!")
        object$gamma <- gamma
        object$H <- "Dirichlet"
        object$F <- "Categorical"
    }
    
    class(object) <- c("CatDirichlet",class(object))
    return(object)
}

#' Sufficient statistics of a "CatDirichlet" object
#'
#' For following Categorical-Dirichlet model structure: \cr
#'      pi|alpha ~ Dir(alpha) \cr
#'      x|pi ~ Categorical(pi) \cr
#' the sufficient statistics of a set of samples x is the effective counts of each uniqueLabel in x. i.e. T(x)[i] = sum(uniqueLabels[i]%in%x). \cr
#' Unique values of x must be in obj$gamma$uniqueLabels, where "obj" is a "CatDirichlet" object, see examples below.
#'
#' @seealso \code{\link{CatDirichlet}}, \code{\link{sufficientStatistics_Weighted.CatDirichlet}} 
#' @param obj A "CatDirichlet" object.
#' @param x numeric,integer or character, samples of the Categorical distribution.
#' @param foreach logical, specifying whether to return the sufficient staistics for each observation. Default FALSE.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return An object of class "ssCat", the sufficient statistics of a set of categorical samples. Or an object of the same class as x if foreach=TRUE.
#' @export
#' @examples
#' obj <- CatDirichlet(gamma=list(alpha=runif(26,1,2),uniqueLabels = letters))
#' x <- sample(letters,size = 20,replace = TRUE)
#' w <- runif(20)
#' sufficientStatistics(obj=obj,x=x)       #return the counts of each unique label
#' sufficientStatistics_Weighted(obj=obj,x=x,w=w) #return the weighted counts of each unique lable
#' sufficientStatistics(obj=obj,x=x,foreach = TRUE) #return the sample itself
#' @references Murphy, Kevin P. Machine learning: a probabilistic perspective. MIT press, 2012.
sufficientStatistics.CatDirichlet <- function(obj,x,foreach=FALSE,...){
    if(missing(x)) stop("'x' must be specified")
    if(!is.vector(x)) x <- as.vector(x)
    if(foreach){
        x
    }else{
        ss <- countFreq(x=x,uniqx = obj$gamma$uniqueLabels)
        class(ss) <- "ssCat"
        ss
    }
}

#' Weighted sufficient statistics of a "CatDirichlet" object
#'
#' For following Categorical-Dirichlet model structure: \cr
#'      pi|alpha ~ Dir(alpha) \cr
#'      x|pi ~ Categorical(pi) \cr
#' the sufficient statistics of a set of samples x and weights w are: \cr
#'   the effective counts (in this case the sum of the weight w) of each uniqueLabel in x \cr
#' Unique values of x must be in obj$gamma$uniqueLabels, where "obj" is a "CatDirichlet" object, see examples below.
#' 
#' @seealso @seealso \code{\link{sufficientStatistics.CatDirichlet}} \code{\link{CatDirichlet}}
#' @param obj A "CatDirichlet" object.
#' @param x numeric,integer or character, samples of the Categorical distribution.
#' @param w numeric, sample weights.
#' @param foreach logical, specifying whether to return the sufficient staistics for each observation. Default FALSE.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return An object of class "ssCat", the sufficient statistics of a set of categorical samples. Or an object of the same class as x if foreach=TRUE.
#' @export
#' @examples
#' obj <- CatDirichlet(gamma=list(alpha=runif(26,1,2),uniqueLabels = letters))
#' x <- sample(letters,size = 20,replace = TRUE)
#' w <- runif(20)
#' sufficientStatistics(obj=obj,x=x)       #return the counts of each unique label
#' sufficientStatistics_Weighted(obj=obj,x=x,w=w) #return the weighted counts of each unique lable
#' @references Murphy, Kevin P. Machine learning: a probabilistic perspective. MIT press, 2012.
sufficientStatistics_Weighted.CatDirichlet <- function(obj,x,w,foreach=FALSE,...){
    if(missing(x)|missing(w)) stop("'x' and 'w' must be both specified")
    if(!is.vector(x)) x <- as.vector(x)
    if(!is.vector(w)) w <- as.vector(w)
    if(length(x)!=length(w)) stop("length of 'x' and 'w' don't match!")
    if(foreach){
        x
    }else{
        ss <- countFreq_Weighted(x=x,uniqx = obj$gamma$uniqueLabels,w=w)
        class(ss) <- "ssCat"
        ss
    }
}

#' Update a "CatDirichlet" object with sample sufficient statistics
#'
#' For the model structure: \cr
#'      pi|alpha ~ Dirichlet(alpha) \cr
#'      x|pi ~ Categorical(pi) \cr
#' update alpha by adding the information of newly observed samples x. \cr
#' The model structure and prior parameters are stored in a "CatDirichlet" object, the prior parameters in this object will be updated after running this function.
#'
#' @seealso \code{\link{CatDirichlet}},\code{\link{posteriorDiscard.CatDirichlet}}
#' @param obj A "CatDirichlet" object.
#' @param ss Sufficient statistics of x. In Categorical-Dirichlet case the sufficient statistic of sample x can be either x itself, of an "ssCat" object generated by the function sufficientStatistics.CatDirichlet().
#' @param w Sample weights, default NULL.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return None. the gamma stored in "obj" will be updated based on "ss".
#' @export
#' @examples
#' obj <- CatDirichlet(gamma=list(alpha=rep(1,26),uniqueLabels = letters))
#' x <- sample(letters,size = 20,replace = TRUE)
#' w <- runif(20)
#' posterior(obj=obj,ss=x)
#' obj
#' posteriorDiscard(obj=obj,ss=x)
#' obj
#' ## weighted sample
#' posterior(obj=obj,ss=x,w=w)
#' obj
#' posteriorDiscard(obj=obj,ss=x,w=w)
#' obj
#' @references Murphy, Kevin P. Machine learning: a probabilistic perspective. MIT press, 2012.
posterior.CatDirichlet <- function(obj,ss,w=NULL,...){
    if(missing(ss)) stop("'ss' must be specified")
    if(.is(ss,"ssCat")){
        posterior_bySufficientStatistics.CatDirichlet(obj=obj,ss=ss)
        invisible(return())
    }
    if(!is.vector(ss)) ss <- as.vector(ss)
    idx <- match(ss,obj$gamma$uniqueLabels,nomatch = NA)
    if(anyNA(idx)) stop("Detect un-recorded label!")
    if(!is.null(w)){
        if(!is.vector(w)) w <- as.vector(w)
        if(length(ss)!=length(w)) stop("length of 'ss' and 'w' don't match!")
        if(length(ss)==1L){
            obj$gamma$alpha[idx] <- obj$gamma$alpha[idx]+w
        }else{
            for(i in 1L:length(ss)) obj$gamma$alpha[idx[i]] <- obj$gamma$alpha[idx[i]]+w[i]
        }
    }else{
        if(length(ss)==1L){
            obj$gamma$alpha[idx] <- obj$gamma$alpha[idx]+1
        }else{
            for(id in idx) obj$gamma$alpha[id] <- obj$gamma$alpha[id]+1
        }
    }
}

#' Update a "CatDirichlet" object with sample sufficient statistics
#' @param obj A "CatDirichlet" object.
#' @param ss Sufficient statistics of x. In Categorical-Dirichlet case the sufficient statistic of sample x can be either x itself, of an "ssCat" object generated by the function sufficientStatistics.CatDirichlet().
#' @param ... Additional arguments to be passed to other inherited types.
#' @return None. the gamma stored in "obj" will be updated based on "ss".
#' @references Murphy, Kevin P. Machine learning: a probabilistic perspective. MIT press, 2012.
posterior_bySufficientStatistics.CatDirichlet <- function(obj,ss,...){
    if(missing(ss)) stop("'ss' must be specified")
    if(!.is(ss,"ssCat")) stop("'ss' must be of class 'ssCat', you need to use sufficientStatistics() to generate 'ssCat' objects")
    
    if(length(ss)!=length(obj$gamma$alpha)) stop("length 'ss' and dirichlet parameters don't match")
    obj$gamma$alpha <- obj$gamma$alpha + as.numeric(ss)
}

#' Update a "CatDirichlet" object with sample sufficient statistics
#'
#' Contrary to posterior(), this function will update alpha by removing the information of observed samples x for the model structure: \cr
#'      pi|alpha ~ Dirichlet(alpha) \cr
#'      x|pi ~ Categorical(pi) \cr
#' the model structure and prior parameters are stored in a "CatDirichlet" object, the prior parameters in this object will be updated after running this function.
#'
#' @seealso \code{\link{CatDirichlet}},\code{\link{posterior.CatDirichlet}}
#' @param obj A "CatDirichlet" object.
#' @param ss Sufficient statistics of x. In Categorical-Dirichlet case the sufficient statistic of sample x can be either x itself, of an "ssCat" object generated by the function sufficientStatistics.CatDirichlet().
#' @param w Sample weights,default NULL.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return None. the prior parameters stored in "obj" will be updated with the information in "ss".
#' @export
#' @examples
#' obj <- CatDirichlet(gamma=list(alpha=rep(1,26),uniqueLabels = letters))
#' x <- sample(letters,size = 20,replace = TRUE)
#' w <- runif(20)
#' posterior(obj=obj,ss=x)
#' obj
#' posteriorDiscard(obj=obj,ss=x)
#' obj
#' ## weighted sample
#' posterior(obj=obj,ss=x,w=w)
#' obj
#' posteriorDiscard(obj=obj,ss=x,w=w)
#' obj
#' @references Murphy, Kevin P. Machine learning: a probabilistic perspective. MIT press, 2012.
posteriorDiscard.CatDirichlet <- function(obj,ss,w=NULL,...){
    if(missing(ss)) stop("'ss' must be specified")
    if(.is(ss,"ssCat")){
        posteriorDiscard_bySufficientStatistics.CatDirichlet(obj = obj,ss = ss)
        invisible(return())
    }
    if(!is.vector(ss)) ss <- as.vector(ss)
    idx <- match(ss,obj$gamma$uniqueLabels,nomatch = NA)
    if(anyNA(idx)) stop("Detect un-recorded label!")
    if(!is.null(w)){
        if(!is.vector(w)) w <- as.vector(w)
        if(length(ss)!=length(w)) stop("length of 'ss' and 'w' don't match!")
        if(length(ss)==1L){
            obj$gamma$alpha[idx] <- obj$gamma$alpha[idx]-w
        }else{
            for(i in 1L:length(ss)) obj$gamma$alpha[idx[i]] <- obj$gamma$alpha[idx[i]]-w[i]
        }
    }else{
        if(length(ss)==1L){
            obj$gamma$alpha[idx] <- obj$gamma$alpha[idx]-1
        }else{
            for(id in idx) obj$gamma$alpha[id] <- obj$gamma$alpha[id]-1
        }
    }

}

#' Update the prior Dirichlet distribution with sample sufficient statistics
#' @seealso \code{\link{posteriorDiscard}}
#' @param obj A "CatDirichlet" object.
#' @param ss Sufficient statistics of x. In Categorical-Dirichlet case the sufficient statistic of sample x can be either x itself, of an "ssCat" object generated by the function sufficientStatistics.CatDirichlet().
#' @param ... Additional arguments to be passed to other inherited types.
#' @return None. the gamma stored in "obj" will be updated based on "ss".
#' @references Murphy, Kevin P. Machine learning: a probabilistic perspective. MIT press, 2012.
posteriorDiscard_bySufficientStatistics.CatDirichlet <- function(obj,ss,...){
    if(missing(ss)) stop("'ss' must be specified")
    if(!.is(ss,"ssCat")) stop("'ss' must be of class 'ssCat', you need to use sufficientStatistics() to generate 'ssCat' objects")
    if(length(ss)!=length(obj$gamma$alpha)) stop("length 'ss' and dirichlet parameters don't match")
    obj$gamma$alpha <- obj$gamma$alpha - as.numeric(ss)
}

#' MAP estimate of a "CatDirichlet" object
#'
#' Generate the MAP estimate of "pi" in following Categorical-Dirichlet structure: \cr
#'      pi|alpha ~ Dirichlet(alpha) \cr
#'      x|pi ~ Categorical(pi) \cr
#' the model structure and prior parameters are stored in a "CatDirichlet" object. \cr
#' MAP is pi_MAP = argmax p(pi|alpha,x).
#'
#' @seealso \code{\link{CatDirichlet}}
#' @param obj A "CatDirichlet" object.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return A numeric vector, the MAP estimate of "pi".
#' @export
#' @examples
#' obj <- CatDirichlet(gamma=list(alpha=rep(1,26),uniqueLabels = letters))
#' x <- sample(letters,size = 20,replace = TRUE)
#' w <- runif(20)
#' posterior(obj=obj,ss=x,w=w)
#' MAP(obj)
#' @references Murphy, Kevin P. Machine learning: a probabilistic perspective. MIT press, 2012.
MAP.CatDirichlet <- function(obj,...){
    tmp <- obj$gamma$alpha-1
    tmp/sum(tmp)
}

#' MPE of a "CatDirichlet" object
#'
#' Generate the MPE of "pi" in following Categorical-Dirichlet structure: \cr
#'      pi|alpha ~ Dirichlet(alpha) \cr
#'      x|pi ~ Categorical(pi) \cr
#' the model structure and prior parameters are stored in a "CatDirichlet" object. \cr
#' MPE is pi_MPE = E(pi|alpha,x), E() is the expectation function.
#'
#' @seealso \code{\link{CatDirichlet}}
#' @param obj A "CatDirichlet" object.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return A numeric vector, the MPE of "pi".
#' @export
#' @examples
#' obj <- CatDirichlet(gamma=list(alpha=rep(1,26),uniqueLabels = letters))
#' x <- sample(letters,size = 20,replace = TRUE)
#' w <- runif(20)
#' posterior(obj=obj,ss=x,w=w)
#' MPE(obj)
#' @references Murphy, Kevin P. Machine learning: a probabilistic perspective. MIT press, 2012.
MPE.CatDirichlet <- function(obj,...){
    obj$gamma$alpha/sum(obj$gamma$alpha)
}

#' Marginal likelihood of a "CatDirichlet" object
#'
#' Generate the marginal likelihood of the following model structure: \cr
#'      pi|alpha ~ Dirichlet(alpha) \cr
#'      x|pi ~ Categorical(pi) \cr
#' the model structure and prior parameters are stored in a "CatDirichlet" object. \cr
#' Marginal likelihood is the likelihood of x|alpha.
#'
#' @seealso @seealso \code{\link{CatDirichlet}}, \code{\link{marginalLikelihood_bySufficientStatistics.CatDirichlet}}
#' @param obj A "CatDirichlet" object.
#' @param x numeric/integer/character vector, observed Categorical samples.
#' @param LOG Return the log density if set to "TRUE".
#' @param ... Additional arguments to be passed to other inherited types.
#' @return numeric, the marginal likelihood.
#' @export
#' @examples
#' obj <- CatDirichlet(gamma=list(alpha=runif(26,1,2),uniqueLabels = letters))
#' x <- sample(letters,size = 20,replace = TRUE)
#' marginalLikelihood(obj=obj,x=x,LOG = TRUE) #marginal likelihood
#' ss <- sufficientStatistics(obj = obj,x=x)
#' marginalLikelihood_bySufficientStatistics(obj=obj,ss = ss,LOG = TRUE)
#' @references Murphy, Kevin P. Machine learning: a probabilistic perspective. MIT press, 2012.
marginalLikelihood.CatDirichlet <- function(obj,x,LOG=TRUE,...){
    if(missing(x)) stop("'x' must be specified")
    if(!is.vector(x)) x <- as.vector(x)
    ss <- countFreq(x=x,uniqx = obj$gamma$uniqueLabels)
    class(ss) <- "ssCat"
    marginalLikelihood_bySufficientStatistics.CatDirichlet(obj=obj,ss=ss,LOG = LOG)
}

#' Marginal likelihood of a "CatDirichlet" object, usnig sufficient statistics
#'
#' Generate the marginal likelihood of a set of observations of the following model structure: \cr
#'      pi|alpha ~ Dirichlet(alpha) \cr
#'      x|pi ~ Categorical(pi) \cr
#' the model structure and prior parameters are stored in a "CatDirichlet" object. \cr
#' Marginal likelihood is the likelihood of x|alpha
#'
#' @seealso @seealso \code{\link{CatDirichlet}}, \code{\link{marginalLikelihood.CatDirichlet}}
#' @param obj A "CatDirichlet" object.
#' @param ss Sufficient statistics of x. In Categorical-Dirichlet case the sufficient statistic of sample x can be either x itself, of an "ssCat" object generated by the function sufficientStatistics.CatDirichlet().
#' @param LOG Return the log density if set to "TRUE".
#' @param ... Additional arguments to be passed to other inherited types.
#' @return numeric, the marginal likelihood.
#' @export
#' @examples
#' obj <- CatDirichlet(gamma=list(alpha=runif(26,1,2),uniqueLabels = letters))
#' x <- sample(letters,size = 20,replace = TRUE)
#' marginalLikelihood(obj=obj,x=x,LOG = TRUE) #marginal likelihood
#' ss <- sufficientStatistics(obj = obj,x=x)
#' marginalLikelihood_bySufficientStatistics(obj=obj,ss = ss,LOG = TRUE)
#' @references Murphy, Kevin P. Machine learning: a probabilistic perspective. MIT press, 2012.
marginalLikelihood_bySufficientStatistics.CatDirichlet <- function(obj,ss,LOG=TRUE,...){
    if(missing(ss)) stop("'ss' must be specified")
    if(!.is(ss,"ssCat")){
        ss <- sufficientStatistics.CatDirichlet(obj = obj,x=ss,foreach = FALSE)
    }
    if(length(ss)!=length(obj$gamma$alpha)) stop("length 'ss' and dirichlet parameters don't match")
    a <- obj$gamma$alpha+as.numeric(ss)
    aa <- sum(a)
    aa0 <- sum(obj$gamma$alpha)
    logp <- lgamma(aa0) - lgamma(aa) + sum(lgamma(a)-lgamma(obj$gamma$alpha))
    if(!LOG) logp <- exp(logp)
    logp
}

#' Posterior predictive density function of a "CatDirichlet" object
#'
#' Generate the the density value of the posterior predictive distribution of the following structure: \cr
#'      pi|alpha ~ Dirichlet(alpha) \cr
#'      x|pi ~ Categorical(pi) \cr
#' the model structure and prior parameters are stored in a "CatDirichlet" object. \cr
#' Posterior predictive is a distribution of x|alpha.
#'
#' @seealso @seealso \code{\link{CatDirichlet}}, \code{\link{dPosteriorPredictive.CatDirichlet}}, \code{\link{marginalLikelihood.CatDirichlet}}
#' @param obj A "CatDirichlet" object.
#' @param x numeric/integer/character vector, observed Categorical samples.
#' @param LOG Return the log density if set to "TRUE".
#' @param ... Additional arguments to be passed to other inherited types.
#' @return A numeric vector, the posterior predictive density.
#' @export
#' @examples
#' obj <- CatDirichlet(gamma=list(alpha=runif(26,1,2),uniqueLabels = letters))
#' x <- sample(letters,size = 20,replace = TRUE)
#' ## res1 and res2 should provide the same result
#' res1 <- dPosteriorPredictive(obj = obj,x=x,LOG = TRUE)
#' res2 <- numeric(length(x))
#' for(i in seq_along(x)) res2[i] <- marginalLikelihood(obj=obj,x=x[i],LOG = TRUE)
#' @references Murphy, Kevin P. Machine learning: a probabilistic perspective. MIT press, 2012.
dPosteriorPredictive.CatDirichlet <- function(obj,x,LOG=TRUE,...){
    if(missing(x)) stop("'x' must be specified")
    if(!is.vector(x)) x <- as.vector(x)
    probs <- obj$gamma$alpha/sum(obj$gamma$alpha)
    out <- probs[match(x,obj$gamma$uniqueLabels)]
    if(LOG) out <- log(out)
    out
}

#' Posterior predictive random generation of a "CatDirichlet" object
#'
#' Generate random samples from the posterior predictive distribution of the following structure: \cr
#'      pi|alpha ~ Dirichlet(alpha) \cr
#'      x|pi ~ Categorical(pi) \cr
#' the model structure and prior parameters are stored in a "CatDirichlet" object \cr
#' posterior predictive is a distribution of x|alpha
#'
#' @seealso @seealso \code{\link{CatDirichlet}}, \code{\link{dPosteriorPredictive.CatDirichlet}}
#' @param obj A "CatDirichlet" object.
#' @param n integer, number of samples.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return A vector of the same type as obj$gamma$uniqueLabels.
#' @export
#' @examples
#' obj <- CatDirichlet(gamma=list(alpha=runif(26,1,2),uniqueLabels = letters))
#' rPosteriorPredictive(obj=obj,n=200)
#' @references Murphy, Kevin P. Machine learning: a probabilistic perspective. MIT press, 2012.
rPosteriorPredictive.CatDirichlet <- function(obj,n,...){
    if(missing(n)) stop("'n' must be specified")
    n <- as.integer(n)
    obj$gamma$uniqueLabels[sample.int(n=length(obj$gamma$uniqueLabels), size = n, replace = TRUE,prob = obj$gamma$alpha/sum(obj$gamma$alpha))]
}
