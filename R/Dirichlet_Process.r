





#' @include Bayesian_Bricks.r Categorical_Inference.r Gaussian_Inference.r Gamma_Inference.r


#' @title a plus b with NA values
#' @description make NA+NA=0 and NA+(not NA)=0 instead of NA+NA=NA and NA+(not NA)=NA
#' @param e1 numeric
#' @param e2 numeric
#' @return numeric, the sum of "e1" and "e2". The length of the returned vector is the same length as the longest of "e1" or "e2".
`%plus%` <- function(e1,e2){
    e1[is.na(e1)] <- 0
    e2[is.na(e2)] <- 0
    e1+e2
}



#' @title Create objects of type "CatDP".
#' @description
#' Create an object of type "CatDP", which represents the Categorical-Dirichlet-Process(Multinomial-Dirichlet-Process) conjugate structure on positive integers:
#'      \deqn{pi|alpha \sim DP(alpha,U)}
#'      \deqn{x|pi \sim Categorical(pi)}
#' where DP(alpha,U) is a Dirichlet Process on positive integers, alpha is the "concentration parameter" of the Dirichlet Process, U is the "base measure" of this Dirichlet process, it is an uniform distribution on all positive integers.Categorical() is the Categorical distribution. See \code{dCategorical} for the definition of the Categorical distribution. \cr
#' In the case of CatDP, x can only be positive integers. \cr
#' This object will be used as a place for recording and accumulating information in the related inference/sampling functions such as posterior(), posteriorDiscard(), MAP(), marginalLikelihood(), dPosteriorPredictive(), rPosteriorPredictive() and so on.
#'
#' @seealso \code{\link{posterior.CatDP}},\code{\link{posteriorDiscard.CatDP}},\code{\link{MAP.CatDP}},\code{\link{marginalLikelihood.CatDP}}, \code{\link{dPosteriorPredictive.CatDP}}, \code{\link{rPosteriorPredictive.CatDP}}
#' @param objCopy an object of type "CatDP". If "objCopy" is not NULL, the function create a new "CatDP" object by copying the content from objCopy, otherwise this new object will be created by using "ENV" and "gamma". Default NULL.
#' @param ENV environment, specify where the object will be created.
#' @param gamma list, a named list of parameters, gamma=list(alpha). Where gamma$alpha is a numeric value specifying the concentration parameter of the Dirichlet Process.
#' @return An object of class "CatDP".
#' @export
#' @examples
#' obj <- CatDP(gamma=list(alpha=2))
#' obj #print the content
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
CatDP <- function(objCopy=NULL,ENV=parent.frame(),gamma=list(alpha=1)){
    object <- BasicBayesian(ENV = ENV)
    if(!is.null(objCopy)){
        if(!.is(objCopy,"CatDP")) stop("'objCopy' must be of class 'CatDP'")
        object$gamma <- objCopy$gamma
        object$H <- objCopy$H
        object$F <- objCopy$F
    }else{
        if(!missing(gamma))
            if((!is.list(gamma)) |
               (!all(names(gamma) %in% c("alpha"))))
                stop("gamma must be of list(alpha)")
        if(length(gamma$alpha)>1) stop("gamma$alpha must be a numeric value of length 1!")
        object$gamma <- gamma
        object$gamma$N <- 0             # number of observed samples
        object$gamma$maxLabel <- 0L     # the length of nk (the largest integer label)
        object$gamma$newLabel <- 1L     # the new label (the label of if draw from H0)
        object$gamma$emptyLabel <- integer(0)
        object$gamma$nk <- numeric(0)   # sample induced partition frequencies
        object$gamma$prop <- numeric(0) # probs of the sample induced partitions
        object$gamma$pH0 <- 1           # probability of drawing from H0
        object$gamma$pFreq <- 0         # probability of drawing from the frequency table        
        object$H <- "DP"
        object$F <- "Cat"
    }
    class(object) <- c("CatDP",class(object))
    return(object)
}

#' @title Sufficient statistics of a "CatDP" object
#' @description
#' For following model structure:
#'      \deqn{pi|alpha \sim DP(alpha,U)}
#'      \deqn{x|pi \sim Categorical(pi)}
#' where DP(alpha,U) is a Dirichlet Process on positive integers, alpha is the "concentration parameter" of the Dirichlet Process, U is the "base measure" of this Dirichlet process, it is an uniform distribution on all positive integers.Categorical() is the Categorical distribution. See \code{dCategorical} for the definition of the Categorical distribution. \cr
#' In the case of CatDP, x can only be positive integers. \cr
#' The model structure and prior parameters are stored in a "CatDP" object. \cr
#' The sufficient statistics of a set of samples x is:
#' \itemize{
#'    \item unique positive integer values
#'    \item effective counts of the unique positive integers
#' }
#' @seealso \code{\link{CatDP}}, \code{\link{sufficientStatistics_Weighted.CatDP}} 
#' @param obj A "CatDP" object.
#' @param x integer, the elements of the vector must all greater than 0, the samples of a Categorical distribution.
#' @param foreach logical, specifying whether to return the sufficient statistics for each observation. Default FALSE.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return An object of class "ssCatDP", the sufficient statistics of a set of categorical samples. Or an integer vector same as x if foreach=TRUE.
#' @export
#' @examples
#' x <- sample(1L:10L,size = 4,replace = TRUE)
#' obj <- CatDP()
#' ## an object of class "ssCatDP", which contains the counts of each unique integer
#' sufficientStatistics(obj=obj,x=x)
#' ## will return x itself
#' sufficientStatistics(obj=obj,x=x,foreach = TRUE)
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
sufficientStatistics.CatDP <- function(obj,x,foreach=FALSE,...){
    if(missing(x)) stop("'x' must be specified")
    if(!is.vector(x)) x <- as.vector(x)
    if(!is.integer(x)) stop("'x' must be a integer vector!")
    if(foreach){
        x
    }else{
        partitionsLabel <- unique(x)
        ss <- list(partitionsLabel=partitionsLabel,
                   freq=countFreq(x=x,uniqx = partitionsLabel))
        class(ss) <- "ssCatDP"
        ss
    }
}

#' @title Weighted sufficient statistics of a "CatDP" object
#' @description
#' For following model structure:
#'      \deqn{pi|alpha \sim DP(alpha,U)}
#'      \deqn{x|pi \sim Categorical(pi)}
#' where DP(alpha,U) is a Dirichlet Process on positive integers, alpha is the "concentration parameter" of the Dirichlet Process, U is the "base measure" of this Dirichlet process, it is an uniform distribution on all positive integers.Categorical() is the Categorical distribution. See \code{dCategorical} for the definition of the Categorical distribution. \cr
#' In the case of CatDP, x can only be positive integers. \cr
#' The model structure and prior parameters are stored in a "CatDP" object. \cr
#' The sufficient statistics of a set of samples x is:
#' \itemize{
#'    \item unique positive integer values
#'    \item effective counts of the unique positive integers
#' }
#' @seealso \code{\link{CatDP}}, \code{\link{sufficientStatistics.CatDP}} 
#' @param obj A "CatDP" object.
#' @param x integer, the elements of the vector must all greater than 0, the samples of a Categorical distribution.
#' @param w numeric, sample weights
#' @param foreach logical, specifying whether to return the sufficient statistics for each observation. Default FALSE.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return An object of class "ssCatDP", the sufficient statistics of a set of categorical samples. Or an integer vector same as x if foreach=TRUE.
#' @export
#' @examples
#' x <- sample(1L:10L,size = 4,replace = TRUE)
#' obj <- CatDP()
#' w <- runif(4)
#' ## return an object of class "ssCatDP" contains the weighted counts of each unique integer
#' sufficientStatistics_Weighted(obj=obj,x=x,w=w)
#' ## return x itself, no matter what w is
#' sufficientStatistics_Weighted(obj=obj,x=x,w=w,foreach = TRUE)
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
sufficientStatistics_Weighted.CatDP <- function(obj,x,w,foreach=FALSE,...){
    if(missing(x)|missing(w)) stop("'x' and 'w' must be both specified")
    if(!is.vector(x)) x <- as.vector(x)
    if(!is.integer(x)) stop("'x' must be a integer vector!")
    if(!is.vector(w)) w <- as.vector(w)
    if(length(x)!=length(w)) stop("length of 'x' and 'w' don't match!")
    if(foreach){
        x
    }else{
        partitionsLabels <- unique(x)
        ss <- list(partitionsLabels=partitionsLabels,
                   freq=countFreq_Weighted(x=x,uniqx = partitionsLabels,w=w))
        class(ss) <- "ssCatDP"
        ss
    }
}

#' @title Update a "CatDP" object with sample sufficient statistics
#' @description
#' For the model structure:
#'      \deqn{pi|alpha \sim DP(alpha,U)}
#'      \deqn{x|pi \sim Categorical(pi)}
#' where DP(alpha,U) is a Dirichlet Process on positive integers, alpha is the "concentration parameter" of the Dirichlet Process, U is the "base measure" of this Dirichlet process, it is an uniform distribution on all positive integers.Categorical() is the Categorical distribution. See \code{dCategorical} for the definition of the Categorical distribution. \cr
#' In the case of CatDP, x can only be positive integers. \cr
#' Update prior knowledge by adding the information of newly observed samples x. The model structure and prior parameters are stored in a "CatDP" object, the prior parameters in this object will be updated after running this function.
#'
#' @seealso \code{\link{CatDP}},\code{\link{posteriorDiscard.CatDP}},\code{\link{sufficientStatistics.CatDP}}
#' @param obj A "CatDP" object.
#' @param ss Sufficient statistics of x. In Categorical-DP case the sufficient statistic of sample x can either be an object of type "ssCatDP" generated by sufficientStatistics(), or x itself(if x is a integer vector with all positive values).
#' @param w Sample weights, default NULL.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return None. the model stored in "obj" will be updated based on "ss".
#' @export
#' @examples
#' ## generate some integer samples
#' x <- sample(1L:10L,size = 40,replace = TRUE)
#' obj <- CatDP()
#' obj2 <- CatDP()
#' obj3 <- CatDP()
#' ## update CatDP object with sufficient statistics
#' ss <- sufficientStatistics(obj=obj,x=x)
#' posterior(obj = obj,ss = ss)
#' ## or, update with x itself
#' posterior(obj = obj2,ss = x)
#' ## or, update with x itself, one by one
#' for(xx in x) posterior(obj = obj3,ss = xx)
#' ## obj, obj2, obj3 should be the same:
#' obj
#' obj2
#' obj3
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
posterior.CatDP <- function(obj,ss,w=NULL,...){
    if(missing(ss)) stop("'ss' must be specified")
    if(.is(ss,"ssCatDP")){
        posterior_bySufficientStatistics.CatDP(obj=obj,ss=ss)
    }else{
        if(!is.vector(ss)) ss <- as.vector(ss)
        if(!is.integer(ss)) stop("'ss' must be a 'ssCatDP' object or an integer vector!")
        if(!is.null(w)){
            if(!is.vector(w)) w <- as.vector(w)
            if(length(ss)!=length(w)) stop("length of 'ss' and 'w' don't match!")
            if(length(ss)==1L){
                if(is.na(obj$gamma$nk[ss])) obj$gamma$nk[ss] <- w
                else obj$gamma$nk[ss] <- obj$gamma$nk[ss]+w
                Diff <- ss-obj$gamma$maxLabel
                if(Diff>0L){
                    if(Diff>1L){
                        obj$gamma$emptyLabel <- c(obj$gamma$emptyLabel,(obj$gamma$maxLabel+1L):(ss-1L))
                        obj$gamma$nk[obj$gamma$emptyLabel] <- 0
                    }
                    obj$gamma$maxLabel <- ss
                }
                if(length(obj$gamma$emptyLabel)>0) obj$gamma$emptyLabel <- setdiff(obj$gamma$emptyLabel,ss)
                if(obj$gamma$newLabel==ss){
                    if(length(obj$gamma$emptyLabel)>0) obj$gamma$newLabel <- min(obj$gamma$emptyLabel)
                    else obj$gamma$newLabel <- obj$gamma$maxLabel+1L
                    
                }
                obj$gamma$N <- obj$gamma$N+w
            }else{
                tmpLabel <- unique(ss)
                tmpFreq <- countFreq_Weighted(x=ss,uniqx = tmpLabel,w=w)
                obj$gamma$nk[tmpLabel] <- obj$gamma$nk[tmpLabel] %plus% tmpFreq
                obj$gamma$nk[is.na(obj$gamma$nk)] <- 0
                obj$gamma$emptyLabel <- which(obj$gamma$nk==0)
                obj$gamma$maxLabel <- length(obj$gamma$nk)
                obj$gamma$newLabel <- min(obj$gamma$emptyLabel,obj$gamma$maxLabel+1L)
                obj$gamma$N <- obj$gamma$N+sum(tmpFreq)
            }
        }else{
            if(length(ss)==1L){
                if(is.na(obj$gamma$nk[ss])) obj$gamma$nk[ss] <- 1
                else obj$gamma$nk[ss] <- obj$gamma$nk[ss]+1
                Diff <- ss-obj$gamma$maxLabel
                if(Diff>0L){
                    if(Diff>1L){
                        obj$gamma$emptyLabel <- c(obj$gamma$emptyLabel,(obj$gamma$maxLabel+1L):(ss-1L))
                        obj$gamma$nk[obj$gamma$emptyLabel] <- 0
                    }
                    obj$gamma$maxLabel <- ss
                }
                if(length(obj$gamma$emptyLabel)>0) obj$gamma$emptyLabel <- setdiff(obj$gamma$emptyLabel,ss)
                if(obj$gamma$newLabel==ss){
                    if(length(obj$gamma$emptyLabel)>0) obj$gamma$newLabel <- min(obj$gamma$emptyLabel)
                    else obj$gamma$newLabel <- obj$gamma$maxLabel+1L
                }
                obj$gamma$N <- obj$gamma$N+1
            }else{
                tmpLabel <- unique(ss)
                tmpFreq <- countFreq(x=ss,uniqx = tmpLabel)
                obj$gamma$nk[tmpLabel] <- obj$gamma$nk[tmpLabel] %plus% tmpFreq
                obj$gamma$nk[is.na(obj$gamma$nk)] <- 0
                obj$gamma$emptyLabel <- which(obj$gamma$nk==0)
                obj$gamma$maxLabel <- length(obj$gamma$nk)
                obj$gamma$newLabel <- min(obj$gamma$emptyLabel,obj$gamma$maxLabel+1L)
                obj$gamma$N <- obj$gamma$N+sum(tmpFreq)
            }
        }
        obj$gamma$prop <- obj$gamma$nk / obj$gamma$N
        obj$gamma$pH0 <- obj$gamma$alpha/(obj$gamma$alpha+obj$gamma$N)
        obj$gamma$pFreq <- 1-obj$gamma$pH0
    }
}

#' @title Update a "CatDP" object with sample sufficient statistics
#' @description Update a "CatDP" object with sample sufficient statistics
#' @param obj A "CatDP" object.
#' @param ss Sufficient statistics of x. In Categorical-DP case the sufficient statistic of sample x can either be an object of type "ssCatDP" generated by sufficientStatistics(), or x itself(if x is a integer vector with all positive values).
#' @param ... Additional arguments to be passed to other inherited types.
#' @return None. the model stored in "obj" will be updated based on "ss".
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
posterior_bySufficientStatistics.CatDP <- function(obj,ss,...){
    if(missing(ss)) stop("'ss' must be specified")
    if(!.is(ss,"ssCatDP")) stop("'ss' must be of class 'ssCatDP', you need to use sufficientStatistics() to generate 'ssCatDP' objects")
    obj$gamma$nk[ss$partitionsLabel] <- obj$gamma$nk[ss$partitionsLabel] %plus% ss$freq
    obj$gamma$nk[is.na(obj$gamma$nk)] <- 0
    obj$gamma$maxLabel <- length(obj$gamma$nk)
    obj$gamma$emptyLabel <- which(obj$gamma$nk==0)
    obj$gamma$newLabel <- min(obj$gamma$emptyLabel,obj$gamma$maxLabel+1L)
    obj$gamma$N <- obj$gamma$N+sum(ss$freq)
    obj$gamma$prop <- obj$gamma$nk / obj$gamma$N
    obj$gamma$pH0 <- obj$gamma$alpha/(obj$gamma$alpha+obj$gamma$N)
    obj$gamma$pFreq <- 1-obj$gamma$pH0
}

#' @title Update a "CatDP" object with sample sufficient statistics
#' @description
#' For the model structure:
#'      \deqn{pi|alpha \sim DP(alpha,U)}
#'      \deqn{x|pi \sim Categorical(pi)}
#' where DP(alpha,U) is a Dirichlet Process on positive integers, alpha is the "concentration parameter" of the Dirichlet Process, U is the "base measure" of this Dirichlet process, it is an uniform distribution on all positive integers.Categorical() is the Categorical distribution. See \code{dCategorical} for the definition of the Categorical distribution. \cr
#' In the case of CatDP, x can only be positive integers. \cr
#' Contrary to posterior(), this function will update the prior knowledge by removing the information of observed samples x. The model structure and prior parameters are stored in a "CatDP" object, the prior parameters in this object will be updated after running this function.
#' @seealso \code{\link{CatDP}},\code{\link{posterior.CatDP}}
#' @param obj A "CatDP" object.
#' @param ss Sufficient statistics of x. In Categorical-DP case the sufficient statistic of sample x can either be an object of type "ssCatDP" generated by sufficientStatistics(), or x itself(if x is a integer vector with all positive values).
#' @param w Sample weights,default NULL.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return None. the model stored in "obj" will be updated with the information in "ss".
#' @export
#' @examples
#' ## generate some integer samples
#' x <- sample(1L:10L,size = 40,replace = TRUE)
#' obj <- CatDP()
#' ss <- sufficientStatistics(obj=obj,x=x)
#' posterior(obj = obj,ss = ss)
#' obj2 <- CatDP(objCopy = obj)            #create obj2 contains the sames info as obj
#' obj3 <- CatDP(objCopy = obj)            #create obj3 contains the sames info as obj
#' ## discard by samples
#' posteriorDiscard(obj = obj,ss = x)
#' ## or discard by samples, one by one
#' for(xx in x) posteriorDiscard(obj = obj2,ss = xx)
#' ## or discard by sufficient statistics
#' posteriorDiscard(obj = obj3,ss = ss)
#' ## obj, obj2 and obj3 should be the same:
#' obj
#' obj2
#' obj3
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
posteriorDiscard.CatDP <- function(obj,ss,w=NULL,...){
    if(missing(ss)) stop("'ss' must be specified")
    if(.is(ss,"ssCatDP")){
        posteriorDiscard_bySufficientStatistics.CatDP(obj=obj,ss=ss)
    }else{
        if(!is.vector(ss)) ss <- as.vector(ss)
        if(!is.integer(ss)) stop("'ss' must be a 'ssCatDP' object or an integer vector!")
        if(!is.null(w)){
            if(!is.vector(w)) w <- as.vector(w)
            if(length(ss)!=length(w)) stop("length of 'ss' and 'w' don't match!")
            if(length(ss)==1L){
                if(is.na(obj$gamma$nk[ss])) stop("No way to discard label that has never been observed before")
                else if(obj$gamma$nk[ss]<=0) stop("No way to discard label that is already zero")
                else obj$gamma$nk[ss] <- obj$gamma$nk[ss]-w
                if(obj$gamma$nk[ss]<=0){
                    obj$gamma$emptyLabel <- c(obj$gamma$emptyLabel,ss)
                    obj$gamma$newLabel <- min(obj$gamma$newLabel,ss)
                }
                obj$gamma$N <- obj$gamma$N-w
            }else{
                tmpLabel <- unique(ss)
                tmpFreq <- countFreq_Weighted(x=ss,uniqx = tmpLabel,w=w)
                obj$gamma$nk[tmpLabel] <- obj$gamma$nk[tmpLabel] - tmpFreq
                if(anyNA(obj$gamma$nk)) stop("No way to discard label that has never been observed before")
                obj$gamma$emptyLabel <- which(obj$gamma$nk==0)
                obj$gamma$newLabel <- min(obj$gamma$emptyLabel,obj$gamma$maxLabel+1L)
                obj$gamma$N <- obj$gamma$N-sum(tmpFreq)
            }
        }else{
            if(length(ss)==1L){
                if(is.na(obj$gamma$nk[ss])) stop("No way to discard label that has never been observed before")
                else if(obj$gamma$nk[ss]<=0) stop("No way to discard label that is already zero")
                else obj$gamma$nk[ss] <- obj$gamma$nk[ss]-1
                if(obj$gamma$nk[ss]<=0){
                    obj$gamma$emptyLabel <- c(obj$gamma$emptyLabel,ss)
                    obj$gamma$newLabel <- min(obj$gamma$newLabel,ss)
                }
                obj$gamma$N <- obj$gamma$N-1
            }else{
                tmpLabel <- unique(ss)
                tmpFreq <- countFreq(x=ss,uniqx = tmpLabel)
                obj$gamma$nk[tmpLabel] <- obj$gamma$nk[tmpLabel] - tmpFreq
                if(anyNA(obj$gamma$nk)) stop("No way to discard label that has never been observed before")
                obj$gamma$emptyLabel <- which(obj$gamma$nk==0)
                obj$gamma$newLabel <- min(obj$gamma$emptyLabel,obj$gamma$maxLabel+1L)
                obj$gamma$N <- obj$gamma$N-sum(tmpFreq)
            }
        }
        if(obj$gamma$N>0)
            obj$gamma$prop <- obj$gamma$nk / obj$gamma$N
        else
            obj$gamma$prop <- rep(0,obj$gamma$maxLabel)
        obj$gamma$pH0 <- obj$gamma$alpha/(obj$gamma$alpha+obj$gamma$N)
        obj$gamma$pFreq <- 1-obj$gamma$pH0
    }
}

#' @title Update a "CatDP" object with sample sufficient statistics
#' @param obj A "CatDP" object.
#' @param ss Sufficient statistics of x. In Categorical-DP case the sufficient statistic of sample x can either be an object of type "ssCatDP" generated by sufficientStatistics(), or x itself(if x is a integer vector with all positive values).
#' @param ... Additional arguments to be passed to other inherited types.
#' @return None. the model stored in "obj" will be updated based on "ss".
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
posteriorDiscard_bySufficientStatistics.CatDP <- function(obj,ss,...){
    if(missing(ss)) stop("'ss' must be specified")
    if(!.is(ss,"ssCatDP")) stop("'ss' must be of class 'ssCatDP', you need to use sufficientStatistics() to generate 'ssCatDP' objects")

    obj$gamma$nk[ss$partitionsLabel] <- obj$gamma$nk[ss$partitionsLabel] - ss$freq
    if(anyNA(obj$gamma$nk)) stop("No way to discard label that has never been observed before")
    obj$gamma$maxLabel <- length(obj$gamma$nk)
    obj$gamma$emptyLabel <- which(obj$gamma$nk==0)
    obj$gamma$newLabel <- min(obj$gamma$emptyLabel,obj$gamma$maxLabel+1L)
    
    obj$gamma$N <- obj$gamma$N-sum(ss$freq)
    if(obj$gamma$N>0)
        obj$gamma$prop <- obj$gamma$nk / obj$gamma$N
    else
        obj$gamma$prop <- rep(0,obj$gamma$maxLabel)
    obj$gamma$pH0 <- obj$gamma$alpha/(obj$gamma$alpha+obj$gamma$N)
    obj$gamma$pFreq <- 1-obj$gamma$pH0
}

#' @title Maximum A Posteriori(MAP) estimate of a "CatDP" object
#' @description
#' Generate the MAP estimate of "pi" in following model structure:
#'      \deqn{pi|alpha \sim DP(alpha,U)}
#'      \deqn{x|pi \sim Categorical(pi)}
#' where DP(alpha,U) is a Dirichlet Process on positive integers, alpha is the "concentration parameter" of the Dirichlet Process, U is the "base measure" of this Dirichlet process, it is an uniform distribution on all positive integers.Categorical() is the Categorical distribution. See \code{dCategorical} for the definition of the Categorical distribution. \cr
#' In the case of CatDP, x can only be positive integers. \cr
#' The model structure and prior parameters are stored in a "CatDP" object. \cr
#' The MAP estimate of pi is pi_MAP = argmax_pi p(pi|alpha,x).
#' @seealso \code{\link{CatDP}}
#' @param obj A "CatDP" object.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return numeric.
#' @export
#' @examples
#' x <- sample(1L:10L,size = 40,replace = TRUE)
#' obj <- CatDP()
#' posterior(obj = obj,ss = x)
#' MAP(obj)
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
MAP.CatDP <- function(obj,...){
    out <- obj$gamma$nk
    out[obj$gamma$newLabel] <- 0L
    out <- out+obj$gamma$alpha-1
    out/sum(out)
}

#' @title Mean Posterior Estimate(MPE) of a "CatDP" object
#' @description
#' Generate the MPE estimate of "pi" in following model structure:
#'      \deqn{pi|alpha \sim DP(alpha,U)}
#'      \deqn{x|pi \sim Categorical(pi)}
#' where DP(alpha,U) is a Dirichlet Process on positive integers, alpha is the "concentration parameter" of the Dirichlet Process, U is the "base measure" of this Dirichlet process, it is an uniform distribution on all positive integers.Categorical() is the Categorical distribution. See \code{dCategorical} for the definition of the Categorical distribution. \cr
#' In the case of CatDP, x can only be positive integers. \cr
#' The model structure and prior parameters are stored in a "CatDP" object. \cr
#' The MPE of pi is pi_MPE = E(pi|alpha,x), E() is the expectation function.
#' @seealso \code{\link{CatDP}}
#' @param obj A "CatDP" object.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return numeric.
#' @export
#' @examples
#' x <- sample(1L:10L,size = 40,replace = TRUE)
#' obj <- CatDP()
#' posterior(obj = obj,ss = x)
#' MPE(obj)
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
MPE.CatDP <- function(obj,...){
    out <- obj$gamma$nk
    out[obj$gamma$newLabel] <- 0L
    out <- out+obj$gamma$alpha
    out/sum(out)
}

#' @title Marginal likelihood of a "CatDP" object
#' @description
#' Generate the marginal likelihood of the following model structure:
#'      \deqn{pi|alpha \sim DP(alpha,U)}
#'      \deqn{x|pi \sim Categorical(pi)}
#' where DP(alpha,U) is a Dirichlet Process on positive integers, alpha is the "concentration parameter" of the Dirichlet Process, U is the "base measure" of this Dirichlet process, it is an uniform distribution on all positive integers.Categorical() is the Categorical distribution. See \code{dCategorical} for the definition of the Categorical distribution. \cr
#' In the case of CatDP, x can only be positive integers. \cr
#' The model structure and prior parameters are stored in a "CatDP" object. \cr
#' Marginal likelihood = p(x|alpha).
#' @seealso \code{\link{CatDP}}, \code{\link{marginalLikelihood_bySufficientStatistics.CatDP}}
#' @param obj A "CatDP" object.
#' @param x integer, the elements of the vector must all greater than 0, the samples of a Categorical distribution.
#' @param LOG Return the log density if set to "TRUE".
#' @param ... Additional arguments to be passed to other inherited types.
#' @return numeric, the marginal likelihood.
#' @export
marginalLikelihood.CatDP <- function(obj,x,LOG=TRUE,...){
    if(missing(x)) stop("'x' must be specified")
    if(!is.vector(x)) x <- as.vector(x)
    if(!is.integer(x)) stop("'x' must be a integer vector!")
    stop("marginalLikelihood method for class 'CatDP' is not implemented yet")
}

#' @title Marginal likelihood of a "CatDP" object, using sufficient statistics
#' @description
#' Generate the marginal likelihood of the following model structure:
#'      \deqn{pi|alpha \sim DP(alpha,U)}
#'      \deqn{x|pi \sim Categorical(pi)}
#' where DP(alpha,U) is a Dirichlet Process on positive integers, alpha is the "concentration parameter" of the Dirichlet Process, U is the "base measure" of this Dirichlet process, it is an uniform distribution on all positive integers.Categorical() is the Categorical distribution. See \code{dCategorical} for the definition of the Categorical distribution. \cr
#' In the case of CatDP, x can only be positive integers. \cr
#' The model structure and prior parameters are stored in a "CatDP" object. \cr
#' Marginal likelihood = p(x|alpha).
#' @seealso \code{\link{CatDP}}, \code{\link{marginalLikelihood_bySufficientStatistics.CatDP}}
#' @param obj A "CatDP" object.
#' @param ss Sufficient statistics of x. In Categorical-DP case the sufficient statistic of sample x can either be an object of type "ssCatDP" generated by sufficientStatistics(), or x itself(if x is a integer vector with all positive values).
#' @param LOG Return the log density if set to "TRUE".
#' @param ... Additional arguments to be passed to other inherited types.
#' @return numeric, the marginal likelihood.
#' @export
marginalLikelihood_bySufficientStatistics.CatDP <- function(obj,ss,LOG=TRUE,...){
    if(missing(ss)) stop("'ss' must be specified")
    if(!.is(ss,"ssCatDP")) stop("'ss' must be of class 'ssCatDP', you need to use sufficientStatistics() to generate 'ssCatDP' objects")
    stop("marginalLikelihood_bySufficientStatistics method for class 'CatDP' is not implemented yet")
}

#' @title Posterior predictive density function of a "CatDP" object
#' @description
#' Generate the the density value of the posterior predictive distribution of the following structure:
#'      \deqn{pi|alpha \sim DP(alpha,U)}
#'      \deqn{x|pi \sim Categorical(pi)}
#' where DP(alpha,U) is a Dirichlet Process on positive integers, alpha is the "concentration parameter" of the Dirichlet Process, U is the "base measure" of this Dirichlet process, it is an uniform distribution on all positive integers.Categorical() is the Categorical distribution. See \code{dCategorical} for the definition of the Categorical distribution. \cr
#' In the case of CatDP, x can only be positive integers. \cr
#' The model structure and prior parameters are stored in a "CatDP" object. \cr
#' Posterior predictive density is p(x|alpha).
#' @seealso \code{\link{CatDP}}, \code{\link{dPosteriorPredictive.CatDP}}, \code{\link{marginalLikelihood.CatDP}}
#' @param obj A "CatDP" object.
#' @param x integer, the elements of the vector must all greater than 0, the samples of a Categorical distribution.
#' @param LOG Return the log density if set to "TRUE".
#' @param ... Additional arguments to be passed to other inherited types.
#' @return A numeric vector, the posterior predictive density.
#' @export
#' @examples
#' x <- sample(1L:10L,size = 40,replace = TRUE)
#' obj <- CatDP()
#' ss <- sufficientStatistics(obj=obj,x=x)
#' posterior(obj = obj,ss = ss)
#' dPosteriorPredictive(obj = obj,x=1L:11L,LOG = FALSE)
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
dPosteriorPredictive.CatDP <- function(obj,x,LOG=TRUE,...){
    if(missing(x)) stop("'x' must be specified")
    if(!is.vector(x)) x <- as.vector(x)
    if(!is.integer(x)) stop("'x' must be a integer vector!")
    out <- obj$gamma$prop[x]
    out <- vapply(out,function(p){
        if(is.na(p)) return(obj$gamma$pH0)
        else if(p==0) return(obj$gamma$pH0)
        else return(obj$gamma$pFreq*p)
    },FUN.VALUE = numeric(1),USE.NAMES = FALSE)
    if(LOG) out <- log(out)
    out
}

#' @title Generate random samples from the posterior predictive distribution of a "CatDP" object
#' @description
#' Generate random samples from the posterior predictive distribution of the following structure:
#'      \deqn{pi|alpha \sim DP(alpha,U)}
#'      \deqn{x|pi \sim Categorical(pi)}
#' where DP(alpha,U) is a Dirichlet Process on positive integers, alpha is the "concentration parameter" of the Dirichlet Process, U is the "base measure" of this Dirichlet process, it is an uniform distribution on all positive integers.Categorical() is the Categorical distribution. See \code{dCategorical} for the definition of the Categorical distribution. \cr
#' In the case of CatDP, x can only be positive integers. \cr
#' The model structure and prior parameters are stored in a "CatDP" object. \cr
#' Posterior predictive distribution is the distribution of x|alpha.
#' @seealso \code{\link{CatDP}}, \code{\link{dPosteriorPredictive.CatDP}}
#' @param obj A "CatDP" object.
#' @param n integer, number of samples.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return integer, the categorical samples.
#' @export
#' @examples
#' x <- sample(1L:10L,size = 40,replace = TRUE)
#' obj <- CatDP()
#' ss <- sufficientStatistics(obj=obj,x=x)
#' posterior(obj = obj,ss = ss)
#' rPosteriorPredictive(obj = obj,n=200L)
#' @import stats
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
rPosteriorPredictive.CatDP <- function(obj,n=1L,...){
    vapply(1L:n,function(i){
        if(runif(1)<obj$gamma$pH0) return(obj$gamma$newLabel)
        else return(sample.int(n =length(obj$gamma$nk),size = 1,prob = obj$gamma$prop))
    },FUN.VALUE = integer(1),USE.NAMES = FALSE)
}

#' @title Create objects of type "CatHDP".
#' @description
#' Create an object of type "CatHDP", which represents the Categorical-Hierarchical-Dirichlet-Process(Multinomial-Hierarchical-Dirichlet-Process) conjugate structure on positive integers:
#'      \deqn{G_j|gamma \sim DP(gamma,U), j = 1:J}
#'      \deqn{pi_j|G_j,alpha \sim DP(alpha,G_j)}
#'      \deqn{z|pi_j \sim Categorical(pi_j)}
#'      \deqn{k|z,G_j \sim Categorical(G_j), if z is a sample from the base measure G}
#' where DP(gamma,U) is a Dirichlet Process on positive integers, gamma is the "concentration parameter", U is the "base measure" of this Dirichlet process, U is an uniform distribution on all positive integers.  DP(gamma,G_j) is a Dirichlet Process on integers with concentration parameter alpha and base measure G_j. Categorical() is the Categorical distribution. See \code{dCategorical} for the definition of the Categorical distribution. \cr
#' In the case of CatHDP, z and k can only be positive integers. \cr
#' This object will be used as a place for recording and accumulating information in the related inference/sampling functions such as posterior(), posteriorDiscard(), rPosteriorPredictive() and so on.
#' @seealso \code{\link{posterior.CatHDP}},\code{\link{posteriorDiscard.CatHDP}} ...
#' @param objCopy an object of type "CatHDP". If "objCopy" is not NULL, the function create a new "CatHDP" object by copying the content from objCopy, otherwise this new object will be created by using "ENV" and "gamma". Default NULL.
#' @param ENV environment, specify where the object will be created.
#' @param gamma list, a named list of parameters, gamma=list(gamma,alpha,j). Where gamma$gamma is a numeric value specifying the concentration parameter of DP(gamma,U), gamma$alpha is a numeric value specifying the concentration parameter of DP(alpha,G_j), gamma$j is the number of groups J.
#' @return An object of class "CatHDP".
#' @export
#' @examples
#' obj <- CatHDP(gamma=list(gamma=1,alpha=1,j=2))
#' obj #print the content
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
CatHDP <- function(objCopy=NULL,ENV=parent.frame(),gamma=list(
                                                       gamma=1, #concentration parameter for prior DP1
                                                       alpha=1,  #concentration parameter for DP2
                                                       j=2     #number of groups
                                                   )){
    object <- new.env(parent=ENV)
    if(!is.null(objCopy)){
        if(!.is(objCopy,"CatHDP")) stop("'objCopy' must be of class 'CatHDP'")
        object$gamma <- objCopy$gamma
        object$Z1 <- objCopy$Z1
        object$Z12map <- objCopy$Z12map
        object$Z2 <- objCopy$Z2
    }else{
        if(!missing(gamma))
            if((!is.list(gamma)) |
               (!all(names(gamma) %in% c("gamma","alpha","j"))))
                stop("gamma must be of list(gamma,j,alpha)")
        gamma$j <- as.integer(gamma$j)
        object$gamma <- gamma
        object$Z1 <- CatDP(ENV = object,gamma = list(alpha=gamma$gamma)) #a CatDP object for partition distribution
        object$Z12map <- replicate(n = gamma$j,expr = integer(0),simplify = FALSE)
        object$Z2 <- lapply(1L:gamma$j,function(J){CatDP(ENV = object,gamma = list(alpha=gamma$alpha))})
        
    }
    class(object) <- c('CatHDP',"BayesianBrick")
    return(object)
}

#' Print the content of an CatHDP object
#' @param x An object of type "CatHDP".
#' @param ... Additional arguments to be passed to other inherited types.
#' @return None.
#' @export
print.CatHDP <- function(x,...){
    cat("The partition distribution is governed by a 'CatDP' class with following parameters:\n")
    print(x$Z1)
    cat("\n\n\n")
    cat("The partition mapping between Z1 and Z2 is:\n")
    print(x$Z12map)
    cat("\n\n\n")
    cat("The partition distribution for each group is governed by following 'CatDP' objects:\n")
    print(x$Z2)
}

#' @title Update a "CatHDP" object with sample sufficient statistics
#' @description
#' For the model structure:
#'      \deqn{G_j|gamma \sim DP(gamma,U), j = 1:J}
#'      \deqn{pi_j|G_j,alpha \sim DP(alpha,G_j)}
#'      \deqn{z|pi_j \sim Categorical(pi_j)}
#'      \deqn{k|z,G_j \sim Categorical(G_j), if z is a sample from the base measure G}
#' where DP(gamma,U) is a Dirichlet Process on positive integers, gamma is the "concentration parameter", U is the "base measure" of this Dirichlet process, U is an uniform distribution on all positive integers.  DP(gamma,G_j) is a Dirichlet Process on integers with concentration parameter alpha and base measure G_j. Categorical() is the Categorical distribution. See \code{dCategorical} for the definition of the Categorical distribution. \cr
#' In the case of CatHDP, z and k can only be positive integers. \cr
#' Update the prior knowledge by adding the information of newly observed samples z and k. The model structure and prior parameters are stored in a "CatHDP" object, the prior parameters in this object will be updated after running this function.
#' @seealso \code{\link{CatHDP}},\code{\link{posteriorDiscard.CatHDP}}
#' @param obj A "CatHDP" object.
#' @param ss1 Sufficient statistics of k. In CatHDP case the sufficient statistic of sample k is k itself(if k is a integer vector with all positive values).
#' @param ss2 Sufficient statistics of z. In CatHDP case the sufficient statistic of sample z is z itself(if z is a integer vector with all positive values).
#' @param j integer, group label.
#' @param w Sample weights, default NULL.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return None. the model stored in "obj" will be updated based on "ss1" and "ss2".
#' @export
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
posterior.CatHDP <- function(obj,ss1,ss2,j,w=NULL,...){
    if(missing(ss2)|missing(ss1)) stop("'ss2' and 'ss1' must be specified")
    if(length(ss2)>1L) stop("posterior.CatHDP can only update from observations one at a time, for now.")
    if(!is.integer(ss2) | ss2<=0L) stop("'ss2' must be a positive integer.")
    if(!is.integer(ss1) | ss1<=0L) stop("'ss1' must be a positive integer.")
    if(!is.integer(j) | j<=0L) stop("'j' must be a positive integer.")

    if(length(obj$Z2)<j){                #if there's new group label observed
        obj$Z2 <- c(obj$Z2,
                    replicate(n = j-length(obj$Z2), expr = do.call("CatDP",list(ENV=obj,gamma=list(alpha=obj$gamma$alpha))),simplify = FALSE))
        obj$Z12map <- c(obj$Z12map,
                        replicate(n = j-length(obj$Z2),expr = integer(0),simplify = FALSE))
    }
    posterior(obj = obj$Z2[[j]],ss=ss2,w=w)
    if(is.null(w)) w <- 1
    if(obj$Z2[[j]]$gamma$nk[ss2]==w){
        posterior(obj = obj$Z1,ss=ss1)
        obj$Z12map[[j]][ss2] <- ss1
    }
}

#' @title Update a "CatHDP" object with sample sufficient statistics
#' @description
#' For the model structure:
#'      \deqn{G_j|gamma \sim DP(gamma,U), j = 1:J}
#'      \deqn{pi_j|G_j,alpha \sim DP(alpha,G_j)}
#'      \deqn{z|pi_j \sim Categorical(pi_j)}
#'      \deqn{k|z,G_j \sim Categorical(G_j), if z is a sample from the base measure G}
#' where DP(gamma,U) is a Dirichlet Process on positive integers, gamma is the "concentration parameter", U is the "base measure" of this Dirichlet process, U is an uniform distribution on all positive integers.  DP(gamma,G_j) is a Dirichlet Process on integers with concentration parameter alpha and base measure G_j. Categorical() is the Categorical distribution. See \code{dCategorical} for the definition of the Categorical distribution. \cr
#' In the case of CatHDP, z and k can only be positive integers. \cr
#' Contrary to posterior(), this function will update the prior knowledge by removing the information of observed samples z and k. The model structure and prior parameters are stored in a "CatDP" object, the prior parameters in this object will be updated after running this function.
#' @seealso \code{\link{CatHDP}},\code{\link{posteriorDiscard.CatHDP}}
#' @param obj A "CatHDP" object.
#' @param ss1 Sufficient statistics of k. In CatHDP case the sufficient statistic of sample k is k itself(if k is a integer vector with all positive values).
#' @param ss2 Sufficient statistics of z. In CatHDP case the sufficient statistic of sample z is z itself(if z is a integer vector with all positive values).
#' @param j integer, group label.
#' @param w Sample weights, default NULL.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return None. the model stored in "obj" will be updated based on "ss1" and "ss2".
#' @export
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
posteriorDiscard.CatHDP <- function(obj,ss1,ss2,j,w=NULL,...){
    if(missing(ss1)|missing(ss2)|missing(j)) stop("'ss1','ss2' and 'j'  must all be specified")
    posteriorDiscard(obj = obj$Z2[[j]],ss=ss2,w=w)
    if(obj$Z2[[j]]$gamma$nk[ss2]==0){
        posteriorDiscard(obj = obj$Z1,ss=ss1)
        obj$Z12map[[j]][ss2] <- as.integer(NA) #very important step
    }
}

#' @title Posterior predictive density function of a "CatHDP" object
#' @description
#' Generate the the density value of the posterior predictive distribution of the following structure:
#'      \deqn{G_j|gamma \sim DP(gamma,U), j = 1:J}
#'      \deqn{pi_j|G_j,alpha \sim DP(alpha,G_j)}
#'      \deqn{z|pi_j \sim Categorical(pi_j)}
#'      \deqn{k|z,G_j \sim Categorical(G_j), if z is a sample from the base measure G}
#' where DP(gamma,U) is a Dirichlet Process on positive integers, gamma is the "concentration parameter", U is the "base measure" of this Dirichlet process, U is an uniform distribution on all positive integers.  DP(gamma,G_j) is a Dirichlet Process on integers with concentration parameter alpha and base measure G_j. Categorical() is the Categorical distribution. See \code{dCategorical} for the definition of the Categorical distribution. \cr
#' In the case of CatHDP, z and k can only be positive integers. \cr
#' The model structure and prior parameters are stored in a "CatHDP" object. \cr
#' Posterior predictive density = p(z,k|alpha,gamma,U,j)
#'
#' @seealso \code{\link{CatHDP}}, \code{\link{dPosteriorPredictive.CatHDP}}
#' @param obj A "CatHDP" object.
#' @param z integer, the elements of the vector must all greater than 0, the samples of a Categorical distribution.
#' @param k integer, the elements of the vector must all greater than 0, the samples of a Categorical distribution.
#' @param j integer, group label.
#' @param LOG Return the log density if set to "TRUE".
#' @param ... Additional arguments to be passed to other inherited types.
#' @return A numeric vector, the posterior predictive density.
#' @export
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
dPosteriorPredictive.CatHDP <- function(obj,z,k,j,LOG=TRUE,...){
    if(missing(z)|missing(j)|missing(k)) stop("'z','j' and 'k' must all be specified")
    if(length(j)>1L) stop("'j' must be of length 1")
    if(length(z)!=length(k)) stop("length of 'z' 'k' don't match")
    if(!is.integer(z)|!is.integer(j)|!is.integer(k)) stop("'z','j','k' must be a positive integer.")
    if(length(z)==1){
        p <- {if(is.na(obj$Z12map[[j]][z])){dPosteriorPredictive(obj = obj$Z1,x = k,LOG = FALSE)*obj$Z2[[j]]$gamma$pH0}else{dPosteriorPredictive(obj = obj$Z2[[j]],x = z,LOG = FALSE)}}
    }else{
        p <- vapply(seq_along(z),function(i){
            if(is.na(obj$Z12map[[j]][z[i]])){
                dPosteriorPredictive(obj = obj$Z1,x = k[i],LOG = FALSE)*obj$Z2[[j]]$gamma$pH0
            }else{
                dPosteriorPredictive(obj = obj$Z2[[j]],x = z[i],LOG = FALSE)
            }
        },FUN.VALUE = numeric(1),USE.NAMES = FALSE)
    }
    if(LOG) p <- log(p)
    p
}

#' @title Generate random samples from the posterior predictive distribution of a "CatHDP" object
#' @description
#' Generate random samples from the posterior predictive distribution of the following structure:
#'      \deqn{G_j|gamma \sim DP(gamma,U), j = 1:J}
#'      \deqn{pi_j|G_j,alpha \sim DP(alpha,G_j)}
#'      \deqn{z|pi_j \sim Categorical(pi_j)}
#'      \deqn{k|z,G_j \sim Categorical(G_j), if z is a sample from the base measure G}
#' where DP(gamma,U) is a Dirichlet Process on positive integers, gamma is the "concentration parameter", U is the "base measure" of this Dirichlet process, U is an uniform distribution on all positive integers.  DP(gamma,G_j) is a Dirichlet Process on integers with concentration parameter alpha and base measure G_j. Categorical() is the Categorical distribution. See \code{dCategorical} for the definition of the Categorical distribution. \cr
#' In the case of CatHDP, z and k can only be positive integers. \cr
#' The model structure and prior parameters are stored in a "CatHDP" object. \cr
#' Posterior predictive is a distribution of z,k|alpha,gamm,U,G_j.
#' 
#' @seealso \code{\link{CatHDP}}, \code{\link{dPosteriorPredictive.CatHDP}}
#' @param obj A "CatHDP" object.
#' @param n integer, number of samples.
#' @param j integer, group label.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return integer, the categorical samples.
#' @export
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
rPosteriorPredictive.CatHDP <- function(obj,n=1L,j,...){
    if(missing(j)) stop("'j'must be specified")
    if(n>1L) stop("for now only support n=1L")
    allz <- which(obj$Z2[[j]]$gamma$nk>0)
    allk <- which(obj$Z1$gamma$nk>0)
    zs <- c(allz,rep(obj$Z2[[j]]$gamma$newLabel,length(allk)+1L))
    ks <- c(obj$Z12map[[j]][allz],allk,obj$Z1$gamma$newLabel)
    probs <- dPosteriorPredictive(obj = obj,z=zs,k=ks,j=j,LOG = FALSE)
    idx <- sample.int(length(zs),size = n,replace = TRUE,prob = probs)
    c(z=zs[idx],k=ks[idx])
}

#' @title Create objects of type "CatHDP2".
#' @description
#' Create an object of type "CatHDP2" that represents the Categorical-Hierarchical-Dirichlet-Process of two Dirichlet Process hierarchies, which is basically CatHDP with an additional layer of Dirichlet Process:
#'      \deqn{G_m |eta \sim DP(eta,U), m = 1:M}
#'      \deqn{G_mj|gamma \sim DP(gamma,G_m), j = 1:J_m}
#'      \deqn{pi_mj|G_mj,alpha \sim DP(alpha,G_mj)}
#'      \deqn{z|pi_mj \sim Categorical(pi_mj)}
#'      \deqn{k|z,G_mj \sim Categorical(G_mj), if z is a sample from the base measure G_mj}
#'      \deqn{u|k,G_m \sim Categorical(G_m), if k is a sample from the base measure G_m}
#' where DP(eta,U) is a Dirichlet Process on positive integers, eta is the "concentration parameter", U is the "base measure" of this Dirichlet process, U is an uniform distribution on all positive integers.  DP(gamma,G_m) is a Dirichlet Process on integers with concentration parameter gamma and base measure G_m. DP(alpha,G_mj) is a Dirichlet Process on integers with concentration parameter alpha and base measure G_mj. Categorical() is the Categorical distribution. See \code{dCategorical} for the definition of the Categorical distribution. \cr
#' In the case of CatHDP2, u, z and k can only be positive integers. \cr
#' This object will be used as a place for recording and accumulating information in the related inference/sampling functions such as posterior(), posteriorDiscard(), dPosteriorPredictive(), rPosteriorPredictive() and so on.
#'
#' @seealso \code{\link{posterior.CatHDP2}},\code{\link{posteriorDiscard.CatHDP2}} ...
#' @param objCopy an object of type "CatHDP2". If "objCopy" is not NULL, the function create a new "CatHDP2" object by copying the content from objCopy, otherwise this new object will be created by using "ENV" and "gamma". Default NULL.
#' @param ENV environment, specify where the object will be created.
#' @param gamma list, a named list of parameters, gamma=list(eta,gamma,alpha,m,j). Where gamma$eta is a numeric value specifying the concentration parameter of DP(eta,U), gamma$gamma is a numeric value specifying the concentration parameter of DP(gamma,G_m), gamma$alpha is a numeric value specifying the concentration parameter of DP(alpha,G_mj), gamma$m is the number of groups M, gamma$j is the number of subgroups in each group, must satisfy length(gamma$j)=gamma$m.
#' @return An object of class "CatHDP2".
#' @export
#' @examples
#' obj <- CatHDP2(gamma=list(eta=1,gamma=1,alpha=1,m=2,j=c(2,3)))
#' obj #print the content
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
CatHDP2 <- function(objCopy=NULL,ENV=parent.frame(),gamma=list(
                                                        eta=1,
                                                        gamma=1, #concentration parameter for prior DP1
                                                        alpha=1,      #concentration parameter for DP2
                                                        m=3,
                                                        j=c(2,3,4)     #number of groups per bigger group
                                                    )){
    object <- new.env(parent=ENV)
    if(!is.null(objCopy)){
        if(!.is(objCopy,"CatHDP2")) stop("'objCopy' must be of class 'CatHDP2'")
        object$gamma <- objCopy$gamma
        object$Z1 <- objCopy$Z1
        object$Z12map <- objCopy$Z12map
        object$Z2 <- objCopy$Z2
    }else{
        if(!missing(gamma))
            if((!is.list(gamma)) |
               (!all(names(gamma) %in% c("eta","gamma","alpha","m","j"))))
                stop("gamma must be of list(eta,gamma,alpha, m,j)")
        gamma$j <- as.integer(gamma$j); gamma$m <- as.integer(gamma$m)
        if(length(gamma$j)==1) gamma$j <- rep(gamma$j,gamma$m)
        if(length(gamma$j)!=gamma$m) stop("j should be a integer vector of number of inner groups satisfying length(j)=m, m is the number of outer groups")
        object$gamma <- gamma
        object$Z1 <- CatDP(ENV = object,gamma = list(alpha=gamma$eta)) #a CatDP object for partition distribution
        object$Z12map <- replicate(n = gamma$m,expr = integer(0),simplify = FALSE)
        object$Z2 <- lapply(1L:gamma$m,function(M){CatHDP(ENV = object,gamma = list(gamma=gamma$gamma,alpha=gamma$alpha,j=gamma$j[M]))})
    }
    class(object) <- c('CatHDP2',"BayesianBrick")
    return(object)

}

#' @title Print the content of an CatHDP2 object
#' @param x An object of type "CatHDP2".
#' @param ... Additional arguments to be passed to other inherited types.
#' @return None.
#' @export
print.CatHDP2 <- function(x,...){
    cat("The partition distribution is governed by a 'CatDP' class with following parameters:\n")
    print(x$Z1)
    cat("\n\n\n")
    cat("The partition mapping between Z1 and Z2 is:\n")
    print(x$Z12map)
    cat("\n\n\n")
    cat("The partition distribution for each group is governed by following 'CatDP' objects:\n")
    print(x$Z2)
}



#' @title Update a "CatHDP2" object with sample sufficient statistics
#' @description
#' For the model structure:
#'      \deqn{G_m |eta \sim DP(eta,U), m = 1:M}
#'      \deqn{G_mj|gamma \sim DP(gamma,G_m), j = 1:J_m}
#'      \deqn{pi_mj|G_mj,alpha \sim DP(alpha,G_mj)}
#'      \deqn{z|pi_mj \sim Categorical(pi_mj)}
#'      \deqn{k|z,G_mj \sim Categorical(G_mj), if z is a sample from the base measure G_mj}
#'      \deqn{u|k,G_m \sim Categorical(G_m), if k is a sample from the base measure G_m}
#' where DP(eta,U) is a Dirichlet Process on positive integers, eta is the "concentration parameter", U is the "base measure" of this Dirichlet process, U is an uniform distribution on all positive integers.  DP(gamma,G_m) is a Dirichlet Process on integers with concentration parameter gamma and base measure G_m. DP(alpha,G_mj) is a Dirichlet Process on integers with concentration parameter alpha and base measure G_mj. Categorical() is the Categorical distribution. See \code{dCategorical} for the definition of the Categorical distribution. \cr
#' In the case of CatHDP2, u, z and k can only be positive integers. \cr
#' Update the prior knowledge by adding the information of newly observed samples u, z and k. The model structure and prior parameters are stored in a "CatHDP2" object, the prior parameters in this object will be updated after running this function.
#' @seealso \code{\link{CatHDP2}},\code{\link{posteriorDiscard.CatHDP2}}
#' @param obj A "CatHDP2" object.
#' @param ss1 Sufficient statistics of u. In CatHDP2 case the sufficient statistic of sample u is u itself(if u is a integer vector with all positive values).
#' @param ss2 Sufficient statistics of k. In CatHDP2 case the sufficient statistic of sample k is k itself(if k is a integer vector with all positive values).
#' @param ss3 Sufficient statistics of z. In CatHDP2 case the sufficient statistic of sample z is z itself(if z is a integer vector with all positive values).
#' @param m integer, group label.
#' @param j integer, subgroup label.
#' @param w Sample weights, default NULL.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return None. the model stored in "obj" will be updated based on "ss1" and "ss2".
#' @export
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
posterior.CatHDP2 <- function(obj,ss1,ss2,ss3,m,j,w=NULL,...){
    if(missing(ss3)|missing(ss2)|missing(ss1)|missing(m)|missing(j)) stop("ss3, ss2, ss1, m and j must be specified")
    if(length(ss3)>1L) stop("posterior.CatHDP2 can only update from observations one at a time, for now.")
    if(!is.integer(ss1) | ss1<=0L) stop("'ss1' must be a positive integer.")

    if(length(obj$Z2)<m){                #if there's new group label observed
        obj$Z2 <- c(obj$Z2,
                    lapply(1L:(m-length(obj$Z2)),function(M){CatHDP(ENV = obj,gamma = list(gamma=obj$gamma$gamma,alpha=obj$gamma$alpha,j=1L))}))
        obj$Z12map <- c(obj$Z12map,
                        replicate(n = m-length(obj$Z2),expr = integer(0),simplify = FALSE))
    }
    posterior.CatHDP(obj = obj$Z2[[m]],ss1 = ss2,ss2=ss3,w=w,j=j)

    if(obj$Z2[[m]]$Z1$gamma$nk[ss2]==1){
        posterior(obj = obj$Z1,ss=ss1)
        obj$Z12map[[m]][ss2] <- ss1
    }
}

#' @title Update a "CatHDP2" object with sample sufficient statistics
#' @description
#' For the model structure:
#'      \deqn{G_m |eta \sim DP(eta,U), m = 1:M}
#'      \deqn{G_mj|gamma \sim DP(gamma,G_m), j = 1:J_m}
#'      \deqn{pi_mj|G_mj,alpha \sim DP(alpha,G_mj)}
#'      \deqn{z|pi_mj \sim Categorical(pi_mj)}
#'      \deqn{k|z,G_mj \sim Categorical(G_mj), if z is a sample from the base measure G_mj}
#'      \deqn{u|k,G_m \sim Categorical(G_m), if k is a sample from the base measure G_m}
#' where DP(eta,U) is a Dirichlet Process on positive integers, eta is the "concentration parameter", U is the "base measure" of this Dirichlet process, U is an uniform distribution on all positive integers.  DP(gamma,G_m) is a Dirichlet Process on integers with concentration parameter gamma and base measure G_m. DP(alpha,G_mj) is a Dirichlet Process on integers with concentration parameter alpha and base measure G_mj. Categorical() is the Categorical distribution. See \code{dCategorical} for the definition of the Categorical distribution. \cr
#' In the case of CatHDP2, u, z and k can only be positive integers. \cr
#' Contrary to posterior(), this function will update the prior knowledge by removing the information of observed samples u, z and k. The model structure and prior parameters are stored in a "CatDP" object, the prior parameters in this object will be updated after running this function.
#'
#' @seealso \code{\link{CatHDP2}},\code{\link{posteriorDiscard.CatHDP2}}
#' @param obj A "CatHDP2" object.
#' @param ss1 Sufficient statistics of u. In CatHDP2 case the sufficient statistic of sample u is u itself(if u is a integer vector with all positive values).
#' @param ss2 Sufficient statistics of k. In CatHDP2 case the sufficient statistic of sample k is k itself(if k is a integer vector with all positive values).
#' @param ss3 Sufficient statistics of z. In CatHDP2 case the sufficient statistic of sample z is z itself(if z is a integer vector with all positive values).
#' @param m integer, group label.
#' @param j integer, subgroup label.
#' @param w Sample weights, default NULL.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return None. the model stored in "obj" will be updated based on "ss1" and "ss2".
#' @export
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
posteriorDiscard.CatHDP2 <- function(obj,ss1,ss2,ss3,m,j,w=NULL,...){
    if(missing(ss1)|missing(ss2)|missing(ss3)|missing(m)|missing(j)) stop("'ss1','ss2', 'ss3', 'm' and 'j'  must all be specified")
    posteriorDiscard.CatHDP(obj = obj$Z2[[m]],ss1=ss2,ss2=ss3,j=j,w=w)
    if(obj$Z2[[m]]$Z1$gamma$nk[ss2]==0){
        posteriorDiscard(obj = obj$Z1,ss=ss1)
        obj$Z12map[[m]][ss2] <- as.integer(NA) #very important step
    }
}

#' @title Posterior predictive density function of a "CatHDP" object
#' @description
#' Generate the the density value of the posterior predictive distribution of the following structure:
#'      \deqn{G_m |eta \sim DP(eta,U), m = 1:M}
#'      \deqn{G_mj|gamma \sim DP(gamma,G_m), j = 1:J_m}
#'      \deqn{pi_mj|G_mj,alpha \sim DP(alpha,G_mj)}
#'      \deqn{z|pi_mj \sim Categorical(pi_mj)}
#'      \deqn{k|z,G_mj \sim Categorical(G_mj), if z is a sample from the base measure G_mj}
#'      \deqn{u|k,G_m \sim Categorical(G_m), if k is a sample from the base measure G_m}
#' where DP(eta,U) is a Dirichlet Process on positive integers, eta is the "concentration parameter", U is the "base measure" of this Dirichlet process, U is an uniform distribution on all positive integers.  DP(gamma,G_m) is a Dirichlet Process on integers with concentration parameter gamma and base measure G_m. DP(alpha,G_mj) is a Dirichlet Process on integers with concentration parameter alpha and base measure G_mj. Categorical() is the Categorical distribution. See \code{dCategorical} for the definition of the Categorical distribution. \cr
#' In the case of CatHDP2, u, z and k can only be positive integers. \cr
#' The model structure and prior parameters are stored in a "CatHDP" object. \cr
#' Posterior predictive density =  p(u,z,k|alpha,gamm,eta,U,G_mj,G_m).
#'
#' @seealso \code{\link{CatHDP}}, \code{\link{dPosteriorPredictive.CatHDP}}
#' @param obj A "CatHDP" object.
#' @param u integer, the elements of the vector must all greater than 0, the samples of a Categorical distribution.
#' @param z integer, the elements of the vector must all greater than 0, the samples of a Categorical distribution.
#' @param k integer, the elements of the vector must all greater than 0, the samples of a Categorical distribution.
#' @param m integer, group label.
#' @param j integer, subgroup label.
#' @param LOG Return the log density if set to "TRUE".
#' @param ... Additional arguments to be passed to other inherited types.
#' @return A numeric vector, the posterior predictive density.
#' @export
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
dPosteriorPredictive.CatHDP2 <- function(obj,u,k,z,m,j,LOG=TRUE,...){
    if(missing(z)|missing(j)|missing(k)) stop("'z','j' and 'k' must all be specified")
    if(length(j)>1L|length(m)>1L) stop("'j' and 'm' must be of length 1")
    if(length(z)!=length(k)|length(z)!=length(u)) stop("length of 'z' 'k' 'u' don't match")
    if(!is.integer(z)|!is.integer(j)|!is.integer(k)|!is.integer(u)|!is.integer(m)) stop("'z','j','k','u','m' must be a positive integer.")
    
    if(length(k)==1){
        p <- {if(is.na(obj$Z12map[[m]][k])){dPosteriorPredictive.CatDP(obj = obj$Z1,x = u,LOG = FALSE)*obj$Z2[[m]]$Z1$gamma$pH0*obj$Z2[[m]]$Z2[[j]]$gamma$pH0}else{dPosteriorPredictive.CatHDP(obj = obj$Z2[[m]],z=z, k=k,j=j,LOG = FALSE)}}
    }else{
        p <- vapply(seq_along(k),function(i){
            if(is.na(obj$Z12map[[m]][k[i]])){
                dPosteriorPredictive.CatDP(obj = obj$Z1,x = u[i],LOG = FALSE)*obj$Z2[[m]]$Z1$gamma$pH0
            }else{
                dPosteriorPredictive.CatHDP(obj = obj$Z2[[m]],z = z[i],k=k[i],j=j,LOG = FALSE)
            }
        },FUN.VALUE = numeric(1),USE.NAMES = FALSE)
    }
    if(LOG) p <- log(p)
    p
}

#' @title Generate random samples from the posterior predictive distribution of a "CatHDP2" object
#' @description
#' Generate random samples from the posterior predictive distribution of the following structure:
#'      \deqn{G_m |eta \sim DP(eta,U), m = 1:M}
#'      \deqn{G_mj|gamma \sim DP(gamma,G_m), j = 1:J_m}
#'      \deqn{pi_mj|G_mj,alpha \sim DP(alpha,G_mj)}
#'      \deqn{z|pi_mj \sim Categorical(pi_mj)}
#'      \deqn{k|z,G_mj \sim Categorical(G_mj), if z is a sample from the base measure G_mj}
#'      \deqn{u|k,G_m \sim Categorical(G_m), if k is a sample from the base measure G_m}
#' where DP(eta,U) is a Dirichlet Process on positive integers, eta is the "concentration parameter", U is the "base measure" of this Dirichlet process, U is an uniform distribution on all positive integers.  DP(gamma,G_m) is a Dirichlet Process on integers with concentration parameter gamma and base measure G_m. DP(alpha,G_mj) is a Dirichlet Process on integers with concentration parameter alpha and base measure G_mj. Categorical() is the Categorical distribution. See \code{dCategorical} for the definition of the Categorical distribution. \cr
#' In the case of CatHDP2, u, z and k can only be positive integers. \cr
#' The model structure and prior parameters are stored in a "CatHDP2" object. \cr
#' Posterior predictive is a distribution of u,z,k|alpha,gamm,eta,U,G_mj,G_m.
#' 
#' @seealso \code{\link{CatHDP2}}, \code{\link{dPosteriorPredictive.CatHDP2}}
#' @param obj A "CatHDP2" object.
#' @param n integer, number of samples.
#' @param m integer, group label.
#' @param j integer, subgroup label.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return integer, the categorical samples.
#' @export
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
rPosteriorPredictive.CatHDP2 <- function(obj,n=1L,m,j,...){
    if(missing(j)) stop("'j'must be specified")
    if(n>1L) stop("for now only support n=1L")
    allz <- which(obj$Z2[[m]]$Z2[[j]]$gamma$nk>0)
    allk <- which(obj$Z2[[m]]$Z1$gamma$nk>0)
    allu <- which(obj$Z1$gamma$nk>0)

    zs <- c(allz, rep(obj$Z2[[m]]$Z2[[j]]$gamma$newLabel,length(allk)+length(allu)+1L) )
    tmp <- c(obj$Z2[[m]]$Z12map[[j]][allz],allk)
    ks <- c(tmp, rep(obj$Z2[[m]]$Z1$gamma$newLabel,length(allu)+1L))
    us <- c(obj$Z12map[[m]][tmp],allu,obj$Z1$gamma$newLabel)
    
    probs <- vapply(seq_along(zs),function(i){
        dPosteriorPredictive.CatHDP2(obj = obj,u=us[i],k=ks[i],z=zs[i],m=m,j=j,LOG = FALSE)
    },FUN.VALUE = numeric(1),USE.NAMES = FALSE)
    idx <- sample.int(length(zs),size = n,replace = TRUE,prob = probs)
    c(u=us[idx],k=ks[idx],z=zs[idx])
}

#' @title Create objects of type "DP".
#' @description
#' Create an object of type "DP", which represents the Dirichlet-Process model structure:
#'      \deqn{pi|alpha \sim DP(alpha,U)}
#'      \deqn{z|pi \sim Categorical(pi)}
#'      \deqn{theta_z|psi \sim H0(psi)}
#'      \deqn{x|theta_z,z \sim F(theta_z)}
#' where DP(alpha,U) is a Dirichlet Process on positive integers, alpha is the "concentration parameter" of the Dirichlet Process, U is the "base measure" of this Dirichlet process. The choice of F() and H0() can be described by an arbitrary "BasicBayesian" object such as "GaussianGaussian","GaussianInvWishart","GaussianNIW", "GaussianNIG", "CatDirichlet", and "CatDP". See \code{?BasicBayesian} for definition of "BasicBayesian" objects, and see for example \code{?GaussianGaussian} for specific "BasicBayesian" instances. As a summary, An "DP" object is simply a combination of a "CatDP" object (see \code{?CatDP}) and an object of any "BasicBayesian" type.\cr
#' The "DP" object will be used as a place for recording and accumulating information in the related inference/sampling functions such as posterior(), posteriorDiscard(), dPosteriorPredictive(), rPosteriorPredictive() and so on.
#' @seealso \code{\link{BasicBayesian}},\code{\link{GaussianNIW}},\code{\link{GaussianNIG}},\code{\link{CatDirichlet}},\code{\link{CatDP}},\code{\link{posterior.DP}},\code{\link{posteriorDiscard.DP}},\code{\link{marginalLikelihood.DP}} ...
#' @param objCopy an object of type "DP". If "objCopy" is not NULL, the function create a new "DP" object by copying the content from objCopy, otherwise this new object will be created by using "ENV" and "gamma". Default NULL.
#' @param ENV environment, specify where the object will be created.
#' @param gamma list, a named list of parameters, gamma=list(alpha,H0aF,parH0). Where gamma$alpha is a numeric value specifying the concentration parameter of the Dirichlet Process. gamma$H0aF is the name of the "BasicBayesian" object which specifies the structure of H0() and F(). gamma$parH0 is the parameters passed to the selected H0aF. For example, if gamma$H0aF="GaussianNIW", then gamma$parH0 should be a named list of NIW parameters: gamma$parH0=list(m,k,v,S), where gamma$parH0$m is a numeric "location" parameter; gamma$parH0$S is a symmetric positive definite matrix representing the "scale" parameters; gamma$parH0$k and gamma$parH0$v are numeric values.
#' @return An object of class "DP".
#' @export
#' @examples
#' \donttest{
#' 
#' ## This is an example of Gibbs sampling on Gaussian mixture models, using Dirichlet Process.
#' 
#' ## generate some Gaussian mixture samples for demo purpose
#' x <- rbind(
#'     rGaussian(50,mu = c(-1.5,1.5),Sigma = matrix(c(0.1,0.03,0.03,0.1),2,2)),
#'     rGaussian(60,mu = c(-1.5,-1.5),Sigma = matrix(c(0.8,0.5,0.5,0.8),2,2)),
#'     rGaussian(70,mu = c(1.5,1.5),Sigma = matrix(c(0.3,0.05,0.05,0.3),2,2)),
#'     rGaussian(50,mu = c(1.5,-1.5),Sigma = matrix(c(0.5,-0.08,-0.08,0.5),2,2))
#' )
#' 
#' ## Step1: Initialize----------------------------------------------
#' maxit <- 100                            #iterative for maxit times
#' burnin <- 50                            #number of burnin samples
#' z <- matrix(1L,nrow(x),maxit-burnin)    #labels
#' ## create an "GaussianNIW" object to track all the changes.
#' obj <- DP(gamma = list(alpha=1,H0aF="GaussianNIW",parH0=list(m=c(0,0),k=0.001,v=2,S=diag(2))))
#' ss <- sufficientStatistics(obj,x=x,foreach = TRUE) #sufficient statistics of each x
#' N <- nrow(x)
#' for(i in 1L:N){ # initialize labels before Gibbs sampling
#'     z[i,1] <- rPosteriorPredictive(obj = obj,n=1,x=x[i,,drop=FALSE])
#'     posterior(obj = obj,ss = ss[[i]], z = z[i,1])
#' }
#' 
#' ## Step2: Main Gibbs sampling loop--------------------------------
#' it <- 0                                 #iteration tracker
#' pb <- txtProgressBar(min = 0,max = maxit,style = 3)
#' while(TRUE){
#'     if(it>burnin) colIdx <- it-burnin
#'     else colIdx <- 1
#'     for(i in 1L:N){
#'         ## remove the sample information from the posterior
#'         posteriorDiscard(obj = obj,ss = ss[[i]],z=z[i,colIdx])
#'         ## get a new sample
#'         z[i,colIdx] <- rPosteriorPredictive(obj = obj,n=1,x=x[i,,drop=FALSE])
#'         ## add the new sample information to the posterior
#'         posterior(obj = obj,ss = ss[[i]],z=z[i,colIdx])
#'     }
#'     if(it>burnin & colIdx<ncol(z)) z[,colIdx+1] <- z[,colIdx] #copy result of previous iteration
#'     it <- it+1
#'     setTxtProgressBar(pb,it)
#'     if(it>=maxit){cat("\n");break}
#'     plot(x=x[,1],y=x[,2],col=z[,colIdx]) #to see how the labels change
#' }
#' 
#' ## Step3: Estimate group labels of each observation---------------
#' col <- apply(z,1,function(l){
#'     tmp <- table(l)
#'     ## pick the most frequent group label in the samples to be the best estimate
#'     names(tmp)[which.max(tmp)]
#' })
#' plot(x=x[,1],y=x[,2],col=col)
#' }
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
#' @references Li, Yuelin, Elizabeth Schofield, and Mithat Gönen. "A tutorial on Dirichlet process mixture modeling." Journal of Mathematical Psychology 91 (2019): 128-144.
DP <- function(objCopy=NULL,ENV=parent.frame(),gamma=list(
                                                 alpha=1, #concentration parameter
                                                 H0aF="GaussianNIW", #name of the base measure and the observation distribution
                                                 parH0=list(m=0,k=1,v=2,S=1) #pars for the basemeasure
                                                 )){
    object <- new.env(parent=ENV)
    if(!is.null(objCopy)){
        if(!.is(objCopy,"DP")) stop("'objCopy' must be of class 'DP'")
        object$H0aF <- objCopy$H0aF
        object$Z <- objCopy$Z
        object$H <- objCopy$H
        object$X <- objCopy$X
    }else{
        if(!missing(gamma))
            if((!is.list(gamma)) |
               (!all(names(gamma) %in% c("alpha","H0aF","parH0"))))
                stop("gamma must be of list(alpha,H0aF,parH0)")
        if(!is.character(gamma$H0aF) | length(gamma$H0aF)!=1) stop("'H0aF' must be a length 1 character specifying the name of the BasicBayesian class")
        object$H0aF <- gamma$H0aF
        object$Z <- CatDP(ENV = object,gamma = list(alpha=gamma$alpha)) #a CatDP object for partition distribution
        object$H <- do.call(gamma$H0aF,list(ENV = object,gamma=gamma$parH0)) #a object of class specified by H0aF for observation distribution, new() will create the object according to class representation of setClass(), so use do.call() instead
        object$X <- list()                                      #the list of observation distributions for each partition
    }
    class(object) <- c('DP',"BayesianBrick")
    return(object)
}

#' @title print the content of a "DP" object
#' @param x An object of type "DP".
#' @param ... Additional arguments to be passed to other inherited types.
#' @return None.
#' @export
print.DP <- function(x,...){
    cat("The partition distribution is governed by a 'CatDP' class with following parameters:\n")
    print(x$Z)
    cat("\n\n\n")
    cat("The observation distribution is governed by a '",x$H0aF,"' object with following parameters::\n",sep = "")
    print(x$H)
    cat("\n\n\n")
    cat("The observation distributions of different partitions is governed by following '",x$H0aF,"' objects:\n",sep = "")
    print(x$X)
}

#' @title Sufficient statistics of a "DP" object
#' @description
#' For following model structure:
#'      \deqn{pi|alpha \sim DP(alpha,U)}
#'      \deqn{z|pi \sim Categorical(pi)}
#'      \deqn{theta_z|psi \sim H0(psi)}
#'      \deqn{x|theta_z,z \sim F(theta_z)}
#' where DP(alpha,U) is a Dirichlet Process on positive integers, alpha is the "concentration parameter" of the Dirichlet Process, U is the "base measure" of this Dirichlet process. The choice of F() and H0() can be described by an arbitrary "BasicBayesian" object such as "GaussianGaussian","GaussianInvWishart","GaussianNIW", "GaussianNIG", "CatDirichlet", and "CatDP". See \code{?BasicBayesian} for definition of "BasicBayesian" objects, and see for example \code{?GaussianGaussian} for specific "BasicBayesian" instances. As a summary, An "DP" object is simply a combination of a "CatDP" object (see \code{?CatDP}) and an object of any "BasicBayesian" type.\cr
#' The sufficient statistics of a set of samples x in a "DP" object is the same sufficient statistics of the "BasicBayesian" inside the "DP", see examples.
#' @seealso \code{\link{DP}}, \code{\link{sufficientStatistics_Weighted.DP}} 
#' @param obj A "DP" object.
#' @param x Random samples of the "BasicBayesian" object.
#' @param ... further arguments passed to the corresponding sufficientStatistics method of the "BasicBayesian" object.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return Return the sufficient statistics of the corresponding BasicBayesian type, see examples.
#' @export
#' @examples
#' obj1 <- DP(gamma=list(alpha=1,H0aF="GaussianNIW",parH0=list(m=1,k=1,v=1,S=1)))
#' obj2 <- DP(gamma=list(alpha=1,H0aF="CatDirichlet",parH0=list(alpha=1,uniqueLabels=letters)))
#' x1 <- rnorm(100)
#' x2 <- sample(letters,100,replace = TRUE)
#' sufficientStatistics(obj = obj1,x=x1)
#' sufficientStatistics(obj = obj2,x=x2)
#' sufficientStatistics(obj = obj1,x=x1,foreach = TRUE)
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
sufficientStatistics.DP <- function(obj,x,...){
    sufficientStatistics(obj = obj$H,x=x,...)
}

#' @title Weighted sufficient statistics of a "DP" object
#' @description
#' For following model structure:
#'      \deqn{pi|alpha \sim DP(alpha,U)}
#'      \deqn{z|pi \sim Categorical(pi)}
#'      \deqn{theta_z|psi \sim H0(psi)}
#'      \deqn{x|theta_z,z \sim F(theta_z)}
#' where DP(alpha,U) is a Dirichlet Process on positive integers, alpha is the "concentration parameter" of the Dirichlet Process, U is the "base measure" of this Dirichlet process. The choice of F() and H0() can be described by an arbitrary "BasicBayesian" object such as "GaussianGaussian","GaussianInvWishart","GaussianNIW", "GaussianNIG", "CatDirichlet", and "CatDP". See \code{?BasicBayesian} for definition of "BasicBayesian" objects, and see for example \code{?GaussianGaussian} for specific "BasicBayesian" instances. As a summary, An "DP" object is simply a combination of a "CatDP" object (see \code{?CatDP}) and an object of any "BasicBayesian" type.\cr
#' The sufficient statistics of a set of samples x in a "DP" object is the same sufficient statistics of the "BasicBayesian" inside the "DP", see examples.
#' @seealso \code{\link{DP}}, \code{\link{sufficientStatistics.DP}} 
#' @param obj A "DP" object.
#' @param x Random samples of the "BasicBayesian" object.
#' @param w numeric, sample weights.
#' @param ... further arguments passed to the corresponding sufficientStatistics method of the "BasicBayesian" object.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return Return the sufficient statistics of the corresponding BasicBayesian type, see examples.
#' @export
#' @examples
#' obj1 <- DP(gamma=list(alpha=1,H0aF="GaussianNIW",parH0=list(m=1,k=1,v=1,S=1)))
#' obj2 <- DP(gamma=list(alpha=1,H0aF="CatDirichlet",parH0=list(alpha=1,uniqueLabels=letters)))
#' x1 <- rnorm(100)
#' x2 <- sample(letters,100,replace = TRUE)
#' w <- runif(100)
#' sufficientStatistics_Weighted(obj = obj1,x=x1,w=w)
#' sufficientStatistics_Weighted(obj = obj2,x=x2,w=w)
#' sufficientStatistics_Weighted(obj = obj1,x=x1,w=w,foreach = TRUE)
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
sufficientStatistics_Weighted.DP <- function(obj,x,w,...){
    sufficientStatistics_Weighted(obj = obj$H,x=x,w=w,...)
}

#' @title Update a "DP" object with sample sufficient statistics
#' @description
#' For the model structure:
#'      \deqn{pi|alpha \sim DP(alpha,U)}
#'      \deqn{z|pi \sim Categorical(pi)}
#'      \deqn{theta_z|psi \sim H0(psi)}
#'      \deqn{x|theta_z,z \sim F(theta_z)}
#' where DP(alpha,U) is a Dirichlet Process on positive integers, alpha is the "concentration parameter" of the Dirichlet Process, U is the "base measure" of this Dirichlet process. The choice of F() and H0() can be described by an arbitrary "BasicBayesian" object such as "GaussianGaussian","GaussianInvWishart","GaussianNIW", "GaussianNIG", "CatDirichlet", and "CatDP". See \code{?BasicBayesian} for definition of "BasicBayesian" objects, and see for example \code{?GaussianGaussian} for specific "BasicBayesian" instances. As a summary, An "DP" object is simply a combination of a "CatDP" object (see \code{?CatDP}) and an object of any "BasicBayesian" type.\cr
#' This function will update the prior knowledge by adding the information of newly observed samples x. The model structure and prior parameters are stored in a "DP" object, the prior parameters in this object will be updated after running this function.
#' @seealso \code{\link{DP}},\code{\link{posteriorDiscard.DP}},\code{\link{sufficientStatistics.DP}}
#' @param obj A "DP" object.
#' @param ss Sufficient statistics of x of the "BasicBayesian" object, must be a list of sufficient statistics for each of the observations. Use sufficientStatistics(...,foreach=TRUE) to generate ss. See examples.
#' @param z integer, the partition label of the parameter space where the observation x is drawn from.
#' @param w Sample weights, default NULL.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return None. the model stored in "obj" will be updated based on "ss".
#' @export
#' @examples
#' x <- rnorm(40)
#' z <- sample(1L:10L,size = 40,replace = TRUE)
#' obj <- DP()
#' ss <- sufficientStatistics(obj = obj,x=x,foreach = TRUE) #must use foreach=TRUE
#' for(i in 1L:length(z)) posterior(obj = obj,ss = ss[[i]],z=z[i])
#' obj
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
posterior.DP <- function(obj,ss,z,w=NULL,...){
    if(missing(ss)|missing(z)) stop("'ss' and 'z' must be specified")
    if(length(z)>1L) stop("posterior.DP can only update from observations one at a time")
    if(!is.integer(z) | z<=0L) stop("'z' must be a positive integer.")
    if(length(obj$X)<z)
        obj$X <- c(obj$X,
                   replicate(n = z-length(obj$X), expr = do.call(obj$H0aF,list(objCopy=obj$H,ENV=obj)),simplify = FALSE))
    posterior(obj = obj$Z,ss=z,w=w)
    posterior(obj = obj$X[[z]],ss = ss,w=w)
}

#' @title Update a "DP" object with sample sufficient statistics
#' @description
#' For the model structure:
#'      \deqn{pi|alpha \sim DP(alpha,U)}
#'      \deqn{z|pi \sim Categorical(pi)}
#'      \deqn{theta_z|psi \sim H0(psi)}
#'      \deqn{x|theta_z,z \sim F(theta_z)}
#' where DP(alpha,U) is a Dirichlet Process on positive integers, alpha is the "concentration parameter" of the Dirichlet Process, U is the "base measure" of this Dirichlet process. The choice of F() and H0() can be described by an arbitrary "BasicBayesian" object such as "GaussianGaussian","GaussianInvWishart","GaussianNIW", "GaussianNIG", "CatDirichlet", and "CatDP". See \code{?BasicBayesian} for definition of "BasicBayesian" objects, and see for example \code{?GaussianGaussian} for specific "BasicBayesian" instances. As a summary, An "DP" object is simply a combination of a "CatDP" object (see \code{?CatDP}) and an object of any "BasicBayesian" type.\cr
#' Contrary to posterior(), this function will update the prior knowledge by removing the information of observed samples x. The model structure and prior parameters are stored in a "CatDP" object, the prior parameters in this object will be updated after running this function.
#' @seealso \code{\link{DP}},\code{\link{posteriorDiscard.DP}},\code{\link{sufficientStatistics.DP}}
#' @param obj A "DP" object.
#' @param ss Sufficient statistics of x of the "BasicBayesian" object, must be a list of sufficient statistics for each of the observations. Use sufficientStatistics(...,foreach=TRUE) to generate ss. See examples.
#' @param z integer, the partition label of the parameter space where the observation x is drawn from.
#' @param w Sample weights, default NULL.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return None. the model stored in "obj" will be updated based on "ss".
#' @export
#' @examples
#' x <- rnorm(40)
#' z <- sample(1L:10L,size = 40,replace = TRUE)
#' obj <- DP()
#' ss <- sufficientStatistics(obj = obj,x=x,foreach = TRUE) #must use foreach=TRUE
#' for(i in 1L:length(z)) posterior(obj = obj,ss = ss[[i]],z=z[i])
#' obj
#' for(i in 1L:length(z)) posteriorDiscard(obj = obj,ss = ss[[i]],z=z[i])
#' obj
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
posteriorDiscard.DP <- function(obj,ss,z,w=NULL,...){
    if(missing(ss)|missing(z)) stop("'ss' and 'z' must be specified")
    posteriorDiscard(obj = obj$Z,ss=z,w=w)
    posteriorDiscard(obj = obj$X[[z]],ss = ss,w=w)
}

#' @title Marginal likelihood for Dirichlet Process(DP)
#' @param obj an object.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return numeric, the marginal likelihood.
#' @export
marginalLikelihood.DP <- function(obj,...){
    stop("marginalLikelihood for this type not implemented yet")
}

#' @title Posterior predictive density function of a Dirichlet Process object
#' @description
#' Generate the the density value of the posterior predictive distribution of the following structure:
#'      \deqn{pi|alpha \sim DP(alpha,U)}
#'      \deqn{z|pi \sim Categorical(pi)}
#'      \deqn{theta_z|psi \sim H0(psi)}
#'      \deqn{x|theta_z,z \sim F(theta_z)}
#' where DP(alpha,U) is a Dirichlet Process on positive integers, alpha is the "concentration parameter" of the Dirichlet Process, U is the "base measure" of this Dirichlet process. The choice of F() and H0() can be described by an arbitrary "BasicBayesian" object such as "GaussianGaussian","GaussianInvWishart","GaussianNIW", "GaussianNIG", "CatDirichlet", and "CatDP". See \code{?BasicBayesian} for definition of "BasicBayesian" objects, and see for example \code{?GaussianGaussian} for specific "BasicBayesian" instances. As a summary, An "DP" object is simply a combination of a "CatDP" object (see \code{?CatDP}) and an object of any "BasicBayesian" type.\cr
#' The model structure and prior parameters are stored in a "DP" object. \cr
#' Posterior predictive density = p(x,z|alpha,psi).
#'
#' @seealso \code{\link{DP}}, \code{\link{dPosteriorPredictive.DP}}, \code{\link{marginalLikelihood.DP}}
#' @param obj A "DP" object.
#' @param x Random samples of the "BasicBayesian" object.
#' @param z integer, the partition label of the parameter space where the observation x is drawn from.
#' @param LOG Return the log density if set to "TRUE".
#' @param ... Additional arguments to be passed to other inherited types.
#' @return A numeric vector, the posterior predictive density.
#' @export
#' @examples
#' x <- rnorm(4)
#' z <- sample(1L:10L,size = 4,replace = TRUE)
#' obj <- DP()
#' ss <- sufficientStatistics(obj = obj,x=x,foreach = TRUE)
#' for(i in 1L:length(x)) posterior(obj = obj,ss=ss[[i]],z=z[i])
#' xnew <- rnorm(10)
#' znew <- sample(1L:10L,size = 10,replace = TRUE)
#' dPosteriorPredictive(obj = obj,x=xnew,z=znew)
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
dPosteriorPredictive.DP <- function(obj,x,z,LOG=TRUE,...){
    if(missing(x)|missing(z)) stop("'x' and 'z' must be specified")
    if(is.vector(x)){
        x <- matrix(x, ncol = 1)
    }else if(!.is(x,"matrix")){
        stop("'x' must be a vector(univariate) or a matrix(multivariate)!")
    }
    if(!is.vector(z)) z <- as.vector(z)
    if(!is.integer(z)) stop("'z' must be a integer vector!")
    if(length(z)!=nrow(x)) stop("number of samples in 'z' and 'x' don't match")

    if(length(z)==1L){
        logp <- dPosteriorPredictive(obj = obj$Z,x = z,LOG = TRUE) + {if(length(obj$X)<z) dPosteriorPredictive(obj = obj$H,x = x,LOG = TRUE) else dPosteriorPredictive(obj = obj$X[[z]],x = x,LOG = TRUE)}
    }else{
        logp <- vapply(1L:length(z),function(i){
            dPosteriorPredictive(obj = obj$Z,x = z[i],LOG = TRUE) + {if(length(obj$X)<z[i]) dPosteriorPredictive(obj = obj$H,x = x[i,,drop=FALSE],LOG = TRUE) else dPosteriorPredictive(obj = obj$X[[z[i]]],x = x[i,,drop=FALSE],LOG = TRUE)}
        },FUN.VALUE = numeric(1),USE.NAMES = FALSE)
    }
    if(!LOG) logp <- exp(logp)
    logp
}

#' @title Generate random samples from the posterior predictive distribution of a "DP" object
#' @description
#' Generate random samples from the posterior predictive distribution of the following structure:
#'      \deqn{pi|alpha \sim DP(alpha,U)}
#'      \deqn{z|pi \sim Categorical(pi)}
#'      \deqn{theta_z|psi \sim H0(psi)}
#'      \deqn{x|theta_z,z \sim F(theta_z)}
#' where DP(alpha,U) is a Dirichlet Process on positive integers, alpha is the "concentration parameter" of the Dirichlet Process, U is the "base measure" of this Dirichlet process. The choice of F() and H0() can be described by an arbitrary "BasicBayesian" object such as "GaussianGaussian","GaussianInvWishart","GaussianNIW", "GaussianNIG", "CatDirichlet", and "CatDP". See \code{?BasicBayesian} for definition of "BasicBayesian" objects, and see for example \code{?GaussianGaussian} for specific "BasicBayesian" instances. As a summary, An "DP" object is simply a combination of a "CatDP" object (see \code{?CatDP}) and an object of any "BasicBayesian" type.\cr
#' The model structure and prior parameters are stored in a "DP" object. \cr
#' This function will generate random samples from the distribution  z|alpha,psi,x.
#' 
#' @seealso \code{\link{DP}}, \code{\link{dPosteriorPredictive.DP}}
#' @param obj A "DP" object.
#' @param n integer, number of samples.
#' @param x Random samples of the "BasicBayesian" object.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return integer, the categorical samples.
#' @export
#' @examples
#' x <- rnorm(4)
#' z <- sample(1L:10L,size = 4,replace = TRUE)
#' obj <- DP()
#' ss <- sufficientStatistics(obj = obj,x=x,foreach = TRUE)
#' for(i in 1L:length(x)) posterior(obj = obj,ss=ss[[i]],z=z[i])
#' xnew <- rnorm(10)
#' znew <- sample(1L:10L,size = 10,replace = TRUE)
#' rPosteriorPredictive(obj = obj,n=1,x=xnew[5])
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
rPosteriorPredictive.DP <- function(obj,n=1,x,...){
    if(missing(x)) stop("'x' must be specified")
    if(is.vector(x)){
        x <- matrix(x, ncol = 1)
    }else if(!.is(x,"matrix")){
        stop("'x' must be a vector(univariate) or a matrix(multivariate)!")
    }
    if(nrow(x)>1L) stop("'x' should be a matrix of only 1 row or a vector of length 1")

    zs <- which(obj$Z$gamma$nk!=0)
    if(length(zs)>0){
        probs <- log(obj$Z$gamma$prop[zs]*obj$Z$gamma$pFreq)+
            vapply(zs,function(z){
                dPosteriorPredictive(obj = obj$X[[z]],x=x,LOG = TRUE)
            },FUN.VALUE = numeric(1),USE.NAMES = FALSE)
    }else{
        probs <- numeric()
    }
    zs <- c(zs,obj$Z$gamma$newLabel)
    probs <- c(probs,log(obj$Z$gamma$pH0)+dPosteriorPredictive(obj = obj$H,x = x,LOG = TRUE))
    probs <- exp(probs-max(probs))
    probs <- probs/sum(probs)
    zs[sample.int(length(zs),size = n,replace = TRUE,prob = probs)]
}

#' @title Create objects of type "HDP".
#' @description
#' Create an object of type "HDP", which represents the Hierarchical-Dirichlet-Process conjugate structure:
#'      \deqn{G_j|gamma \sim DP(gamma,U), j = 1:J}
#'      \deqn{pi_j|G_j,alpha \sim DP(alpha,G_j)}
#'      \deqn{z|pi_j \sim Categorical(pi_j)}
#'      \deqn{k|z,G_j \sim Categorical(G_j), if z is a sample from the base measure G}
#'      \deqn{theta_k|psi \sim H0(psi)}
#'      \deqn{x|theta_k,k \sim F(theta_k)}
#' where DP(gamma,U) is a Dirichlet Process on positive integers, gamma is the "concentration parameter", U is the "base measure" of this Dirichlet process, U is an uniform distribution on all positive integers.  DP(gamma,G_j) is a Dirichlet Process on integers with concentration parameter alpha and base measure G_j. The choice of F() and H0() can be described by an arbitrary "BasicBayesian" object such as "GaussianGaussian","GaussianInvWishart","GaussianNIW", "GaussianNIG", "CatDirichlet", and "CatDP". See \code{?BasicBayesian} for definition of "BasicBayesian" objects, and see for example \code{?GaussianGaussian} for specific "BasicBayesian" instances. As a summary, An "HDP" object is simply a combination of a "CatHDP" object (see \code{?CatHDP}) and an object of any "BasicBayesian" type.\cr
#' In the case of HDP, z and k can only be positive integers. \cr
#' This object will be used as a place for recording and accumulating information in the related inference/sampling functions such as posterior(), posteriorDiscard(), dPosteriorPredictive(), rPosteriorPredictive() and so on.
#' @seealso \code{\link{BasicBayesian}},\code{\link{GaussianNIW}},\code{\link{GaussianNIG}},\code{\link{CatDirichlet}},\code{\link{CatHDP}},\code{\link{posterior.HDP}},\code{\link{posteriorDiscard.HDP}},\code{\link{marginalLikelihood.HDP}} ...
#' @param objCopy an object of type "HDP". If "objCopy" is not NULL, the function create a new "HDP" object by copying the content from objCopy, otherwise this new object will be created by using "ENV" and "gamma". Default NULL.
#' @param ENV environment, specify where the object will be created.
#' @param gamma list, a named list of parameters, gamma=list(gamma,alpha,j,H0aF,parH0). Where gamma$gamma is a numeric value specifying the concentration parameter of DP(gamma,U), gamma$alpha is a numeric value specifying the concentration parameter of DP(alpha,G_j), gamma$j is the number of groups J. gamma$H0aF is the name of the "BasicBayesian" object which specifies the structure of H0() and F(). gamma$parH0 is the parameters passed to the selected H0aF. For example, if gamma$H0aF="GaussianNIW", then gamma$parH0 should be a named list of NIW parameters: gamma$parH0=list(m,k,v,S), where gamma$parH0$m is a numeric "location" parameter; gamma$parH0$S is a symmetric positive definite matrix representing the "scale" parameters; gamma$parH0$k and gamma$parH0$v are numeric values.
#' @return An object of class "HDP".
#' @export
#' @examples
#' \donttest{
#'     
#' ## This is an example of Gibbs sampling on an hierarchical mixture model,
#' ## using Hierarchical Dirichlet Process.
#'
#' ## load some hierarchical mixture data, check ?mmhData for details.
#' data(mmhData)
#' x <- mmhData$x
#' js <- mmhData$groupLabel
#' 
#' 
#' ## Step1: initialize--------------------------------------------------
#' z <- rep(1L,nrow(x))
#' k <- rep(1L,nrow(x))
#' ## create a HDP object to track all the changes
#' obj <- HDP(gamma = list(gamma=1,j=max(js),alpha=1,
#'            H0aF="GaussianNIW",
#'            parH0=list(m=c(0,0),k=0.001,v=2,S=diag(2)*0.01)))
#' ss <- sufficientStatistics(obj$H,x=x,foreach = TRUE) #sufficient statistics
#' N <- length(ss)
#' for(i in 1L:N){# initialize k and z
#'     tmp <- rPosteriorPredictive(obj = obj,n=1,x=x[i,,drop=FALSE],j=js[i])
#'     z[i] <- tmp["z"]
#'     k[i] <- tmp["k"]
#'     posterior.HDP(obj = obj,ss = ss[[i]],ss1 = k[i],ss2 = z[i],j = js[i])
#' }
#' 
#' ## Step2: main Gibbs loop---------------------------------------------
#' maxit <- 20                             #iterative for maxit times
#' it <- 0                                 #iteration tracker
#' pb <- txtProgressBar(min = 0,max = maxit,style = 3)
#' while(TRUE){
#'     for(i in 1L:N){
#'         ## remove the sample from the posterior info
#'         posteriorDiscard(obj = obj,ss = ss[[i]],ss1=k[i],ss2=z[i],j=js[i])
#'         ## resample a new partition
#'         tmp <- rPosteriorPredictive(obj = obj,n=1,x=x[i,,drop=FALSE],j=js[i])
#'         z[i] <- tmp["z"]
#'         k[i] <- tmp["k"]
#'         posterior(obj = obj,ss = ss[[i]], ss1=k[i],ss2 = z[i],j=js[i])
#'     }
#'     plot(x=x[,1],y=x[,2],col=k)         #to visualize the group label dynamics
#'     it <- it+1
#'     setTxtProgressBar(pb,it)
#'     if(it>=maxit){cat("\n");break}
#' }
#' 
#' }
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
HDP <- function(objCopy=NULL,ENV=parent.frame(),gamma=list(
                                                    gamma=1, #concentration parameter for prior DP1
                                                    j=2L,     #number of groups
                                                    alpha=1,     #concentration parameter for DP2
                                                    H0aF="GaussianNIW", #name of the base measure and the observation distribution
                                                    parH0=list(m=0,k=1,v=2,S=1) #pars for the basemeasure
                                                )){
    object <- new.env(parent=ENV)
    if(!is.null(objCopy)){
        if(!.is(objCopy,"HDP")) stop("'objCopy' must be of class 'HDP'")
        object$gamma <- objCopy$gamma
        object$Z <- objCopy$Z
        object$H <- objCopy$H
        object$X <- objCopy$X
    }else{
        if(!missing(gamma))
            if((!is.list(gamma)) |
               (!all(names(gamma) %in% c("gamma","alpha","j","H0aF","parH0"))))
                stop("gamma must be of list(gamma,j,alpha,H0aF,parH0)")
        if(!is.character(gamma$H0aF) | length(gamma$H0aF)!=1) stop("'H0aF' must be a length 1 character specifying the name of the BasicBayesian class")
        if(!is.integer(gamma$j)) gamma$j <- as.integer(gamma$j)
        object$gamma <- gamma
        object$Z <- CatHDP(ENV = object,gamma = list(gamma=gamma$gamma,j=gamma$j,alpha=gamma$alpha))
        object$H <- do.call(gamma$H0aF,list(ENV = object,gamma=gamma$parH0)) #a object of class specified by H0aF for observation distribution, new() will create the object according to class representation of setClass(), so use do.call() instead
        object$X <- list()                                      #the list of observation distributions for each partition
    }
    class(object) <- c('HDP',"BayesianBrick")
    return(object)

}

#' @title print the content of a "HDP" object
#' @param x An object of type "HDP".
#' @param ... Additional arguments to be passed to other inherited types.
#' @return None.
#' @export
print.HDP <- function(x,...){
    cat("The partition distribution is governed by a 'CatHDP' class with following parameters:\n")
    print(x$Z)
    cat("\n\n\n")
    cat("The observation distribution is governed by a '",x$gamma$H0aF,"' object with following parameters::\n",sep = "")
    print(x$H)
    cat("\n\n\n")
    cat("The observation distributions of different partitions is governed by following '",x$gamma$H0aF,"' objects:\n",sep = "")
    print(x$X)
}

#' @title Sufficient statistics of a "HDP" object
#' @description
#' For following model structure:
#'      \deqn{G_j|gamma \sim DP(gamma,U), j = 1:J}
#'      \deqn{pi_j|G_j,alpha \sim DP(alpha,G_j)}
#'      \deqn{z|pi_j \sim Categorical(pi_j)}
#'      \deqn{k|z,G_j \sim Categorical(G_j), if z is a sample from the base measure G}
#'      \deqn{theta_k|psi \sim H0(psi)}
#'      \deqn{x|theta_k,k \sim F(theta_k)}
#' where DP(gamma,U) is a Dirichlet Process on positive integers, gamma is the "concentration parameter", U is the "base measure" of this Dirichlet process, U is an uniform distribution on all positive integers.  DP(gamma,G_j) is a Dirichlet Process on integers with concentration parameter alpha and base measure G_j. The choice of F() and H0() can be described by an arbitrary "BasicBayesian" object such as "GaussianGaussian","GaussianInvWishart","GaussianNIW", "GaussianNIG", "CatDirichlet", and "CatDP". See \code{?BasicBayesian} for definition of "BasicBayesian" objects, and see for example \code{?GaussianGaussian} for specific "BasicBayesian" instances. As a summary, An "HDP" object is simply a combination of a "CatHDP" object (see \code{?CatHDP}) and an object of any "BasicBayesian" type.\cr
#' In the case of HDP, z and k can only be positive integers. \cr
#' The sufficient statistics of a set of samples x in a "HDP" object is the same sufficient statistics of the "BasicBayesian" inside the "HDP", see examples.
#' @seealso \code{\link{HDP}}, \code{\link{sufficientStatistics_Weighted.HDP}} 
#' @param obj A "HDP" object.
#' @param x Random samples of the "BasicBayesian" object.
#' @param ... further arguments passed to the corresponding sufficientStatistics method of the "BasicBayesian" object.
#' @return Return the sufficient statistics of the corresponding BasicBayesian type, see examples.
#' @export
#' @examples
#' ## a HDP with Gaussian NIW observations
#' obj1 <- HDP(gamma=list(gamma=1,alpha=1,j=2,
#'                        H0aF="GaussianNIW",
#'                        parH0=list(m=0,k=1,v=2,S=1)))
#' ## a HDP with Categorical-Dirichlet observations
#' obj2 <- HDP(gamma=list(gamma=1,alpha=1,j=2,
#'             H0aF="CatDirichlet",
#'             parH0=list(alpha=1,uniqueLabels=letters[1:3])))
#' x1 <- rnorm(100)
#' x2 <- sample(letters[1:3],100,replace = TRUE)
#' sufficientStatistics(obj = obj1,x=x1,foreach = TRUE)
#' sufficientStatistics(obj = obj1,x=x1,foreach = FALSE)
#' sufficientStatistics(obj = obj2,x=x2,foreach = FALSE)
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
sufficientStatistics.HDP <- function(obj,x,...){
    sufficientStatistics(obj = obj$H,x=x,...)
}

#' @title Weighted sufficient statistics of a "HDP" object
#' @description
#' For following model structure:
#'      \deqn{G_j|gamma \sim DP(gamma,U), j = 1:J}
#'      \deqn{pi_j|G_j,alpha \sim DP(alpha,G_j)}
#'      \deqn{z|pi_j \sim Categorical(pi_j)}
#'      \deqn{k|z,G_j \sim Categorical(G_j), if z is a sample from the base measure G}
#'      \deqn{theta_k|psi \sim H0(psi)}
#'      \deqn{x|theta_k,k \sim F(theta_k)}
#' where DP(gamma,U) is a Dirichlet Process on positive integers, gamma is the "concentration parameter", U is the "base measure" of this Dirichlet process, U is an uniform distribution on all positive integers.  DP(gamma,G_j) is a Dirichlet Process on integers with concentration parameter alpha and base measure G_j. The choice of F() and H0() can be described by an arbitrary "BasicBayesian" object such as "GaussianGaussian","GaussianInvWishart","GaussianNIW", "GaussianNIG", "CatDirichlet", and "CatDP". See \code{?BasicBayesian} for definition of "BasicBayesian" objects, and see for example \code{?GaussianGaussian} for specific "BasicBayesian" instances. As a summary, An "HDP" object is simply a combination of a "CatHDP" object (see \code{?CatHDP}) and an object of any "BasicBayesian" type.\cr
#' In the case of HDP, z and k can only be positive integers. \cr
#' The sufficient statistics of a set of samples x in a "HDP" object is the same sufficient statistics of the "BasicBayesian" inside the "HDP", see examples.
#'
#' @seealso \code{\link{HDP}}, \code{\link{sufficientStatistics.HDP}} 
#' @param obj A "HDP" object.
#' @param x Random samples of the "BasicBayesian" object.
#' @param w numeric, sample weights.
#' @param ... further arguments passed to the corresponding sufficientStatistics method of the "BasicBayesian" object.
#' @return Return the sufficient statistics of the corresponding BasicBayesian type, see examples.
#' @export
#' @examples
#' ## a HDP with Gaussian NIW observations
#' obj1 <- HDP(gamma=list(gamma=1,alpha=1,j=2,
#'                        H0aF="GaussianNIW",
#'                        parH0=list(m=0,k=1,v=2,S=1)))
#' ## a HDP with Categorical-Dirichlet observations
#' obj2 <- HDP(gamma=list(gamma=1,alpha=1,j=2,
#'             H0aF="CatDirichlet",
#'             parH0=list(alpha=1,uniqueLabels=letters[1:3])))
#' x1 <- rnorm(100)
#' x2 <- sample(letters[1:3],100,replace = TRUE)
#' w <- runif(100)
#' sufficientStatistics_Weighted(obj = obj1,x=x1,w=w,foreach = FALSE)
#' sufficientStatistics_Weighted(obj = obj2,x=x2,w=w,foreach = FALSE)
#' sufficientStatistics_Weighted(obj = obj1,x=x1,w=w,foreach = TRUE)
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
sufficientStatistics_Weighted.HDP <- function(obj,x,w,...){
    sufficientStatistics_Weighted(obj = obj$H,x=x,w=w,...)
}

#' @title Update a "HDP" object with sample sufficient statistics
#' @description
#' For the model structure:
#'      \deqn{G_j|gamma \sim DP(gamma,U), j = 1:J}
#'      \deqn{pi_j|G_j,alpha \sim DP(alpha,G_j)}
#'      \deqn{z|pi_j \sim Categorical(pi_j)}
#'      \deqn{k|z,G_j \sim Categorical(G_j), if z is a sample from the base measure G}
#'      \deqn{theta_k|psi \sim H0(psi)}
#'      \deqn{x|theta_k,k \sim F(theta_k)}
#' where DP(gamma,U) is a Dirichlet Process on positive integers, gamma is the "concentration parameter", U is the "base measure" of this Dirichlet process, U is an uniform distribution on all positive integers.  DP(gamma,G_j) is a Dirichlet Process on integers with concentration parameter alpha and base measure G_j. The choice of F() and H0() can be described by an arbitrary "BasicBayesian" object such as "GaussianGaussian","GaussianInvWishart","GaussianNIW", "GaussianNIG", "CatDirichlet", and "CatDP". See \code{?BasicBayesian} for definition of "BasicBayesian" objects, and see for example \code{?GaussianGaussian} for specific "BasicBayesian" instances. As a summary, An "HDP" object is simply a combination of a "CatHDP" object (see \code{?CatHDP}) and an object of any "BasicBayesian" type.\cr
#' In the case of HDP, z and k can only be positive integers. \cr
#' This function will update the prior knowledge by adding the information of newly observed samples x, z and k. The model structure and prior parameters are stored in a "HDP" object, the prior parameters in this object will be updated after running this function.
#'
#' @seealso \code{\link{HDP}},\code{\link{posteriorDiscard.HDP}},\code{\link{sufficientStatistics.HDP}}
#' @param obj A "HDP" object.
#' @param ss Sufficient statistics of x of the "BasicBayesian" object, must be a list of sufficient statistics for each of the observations. Use sufficientStatistics(...,foreach=TRUE) to generate ss.
#' @param ss1 Sufficient statistics of k. In HDP case the sufficient statistic of sample k is k itself(if k is a integer vector with all positive values).
#' @param ss2 Sufficient statistics of z. In HDP case the sufficient statistic of sample z is z itself(if z is a integer vector with all positive values).
#' @param j integer, group label.
#' @param w Sample weights, default NULL.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return None. the model stored in "obj" will be updated based on "ss", "ss1" and "ss2".
#' @export
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
posterior.HDP <- function(obj,ss,ss1,ss2,j,w=NULL,...){
    if(missing(ss)|missing(ss2)|missing(j)|missing(ss1)) stop("'ss','ss2','j' and 'ss1' must all be specified")
    if(length(ss2)>1L) stop("posterior.HDP can only update from observations one at a time, for now")
    if(!is.integer(ss2) | ss2<=0L) stop("'ss2' must be a positive integer.")
    if(!is.integer(j) | j<=0L) stop("'j' must be a positive integer.")
    if(!is.integer(ss1) | ss1<=0L) stop("'ss1' must be a positive integer.")
    if(length(obj$X)<ss1)
        obj$X <- c(obj$X,
                   replicate(n = ss1-length(obj$X), expr = do.call(obj$gamma$H0aF,list(objCopy=obj$H,ENV=obj)),simplify = FALSE))
    posterior(obj = obj$Z,ss1 = ss1,ss2 = ss2,j = j,w = w)
    posterior(obj = obj$X[[ss1]],ss = ss,w=w)
}

#' @title Update a "HDP" object with sample sufficient statistics
#' @description
#' For the model structure:
#'      \deqn{G_j|gamma \sim DP(gamma,U), j = 1:J}
#'      \deqn{pi_j|G_j,alpha \sim DP(alpha,G_j)}
#'      \deqn{z|pi_j \sim Categorical(pi_j)}
#'      \deqn{k|z,G_j \sim Categorical(G_j), if z is a sample from the base measure G}
#'      \deqn{theta_k|psi \sim H0(psi)}
#'      \deqn{x|theta_k,k \sim F(theta_k)}
#' where DP(gamma,U) is a Dirichlet Process on positive integers, gamma is the "concentration parameter", U is the "base measure" of this Dirichlet process, U is an uniform distribution on all positive integers.  DP(gamma,G_j) is a Dirichlet Process on integers with concentration parameter alpha and base measure G_j. The choice of F() and H0() can be described by an arbitrary "BasicBayesian" object such as "GaussianGaussian","GaussianInvWishart","GaussianNIW", "GaussianNIG", "CatDirichlet", and "CatDP". See \code{?BasicBayesian} for definition of "BasicBayesian" objects, and see for example \code{?GaussianGaussian} for specific "BasicBayesian" instances. As a summary, An "HDP" object is simply a combination of a "CatHDP" object (see \code{?CatHDP}) and an object of any "BasicBayesian" type.\cr
#' In the case of HDP, z and k can only be positive integers. \cr
#' Contrary to posterior(), this function will update the prior knowledge by removing the information of observed samples x. The model structure and prior parameters are stored in a "CatDP" object, the prior parameters in this object will be updated after running this function.
#'
#' @seealso \code{\link{HDP}},\code{\link{posteriorDiscard.HDP}},\code{\link{sufficientStatistics.HDP}}
#' @param obj A "HDP" object.
#' @param ss Sufficient statistics of x of the "BasicBayesian" object, must be a list of sufficient statistics for each of the observations. Use sufficientStatistics(...,foreach=TRUE) to generate ss.
#' @param ss1 Sufficient statistics of k. In HDP case the sufficient statistic of sample k is k itself(if k is a integer vector with all positive values).
#' @param ss2 Sufficient statistics of z. In HDP case the sufficient statistic of sample z is z itself(if z is a integer vector with all positive values).
#' @param j integer, group label.
#' @param w Sample weights, default NULL.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return None. the model stored in "obj" will be updated based on "ss1" and "ss2".
#' @export
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
posteriorDiscard.HDP <- function(obj,ss,ss1,ss2,j,w=NULL,...){
    if(missing(ss)|missing(ss2)|missing(j)|missing(ss1)) stop("'ss','ss2','j' and 'ss1' must all be specified")
    posteriorDiscard(obj = obj$Z,ss1 = ss1,ss2 = ss2,j = j,w = w)
    posteriorDiscard(obj = obj$X[[ss1]],ss = ss,w=w)
}

#' @title Marginal likelihood for HDP
#' @param obj an object.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return numeric, the marginal likelihood.
#' @export
marginalLikelihood.HDP <- function(obj,...){
    stop("marginalLikelihood for this type not implemented yet")
}

#' @title Posterior predictive density function of a "HDP" object
#' @description
#' Generate the the density value of the posterior predictive distribution of the following structure:
#'      \deqn{G_j|gamma \sim DP(gamma,U), j = 1:J}
#'      \deqn{pi_j|G_j,alpha \sim DP(alpha,G_j)}
#'      \deqn{z|pi_j \sim Categorical(pi_j)}
#'      \deqn{k|z,G_j \sim Categorical(G_j), if z is a sample from the base measure G}
#'      \deqn{theta_k|psi \sim H0(psi)}
#'      \deqn{x|theta_k,k \sim F(theta_k)}
#' where DP(gamma,U) is a Dirichlet Process on positive integers, gamma is the "concentration parameter", U is the "base measure" of this Dirichlet process, U is an uniform distribution on all positive integers.  DP(gamma,G_j) is a Dirichlet Process on integers with concentration parameter alpha and base measure G_j. The choice of F() and H0() can be described by an arbitrary "BasicBayesian" object such as "GaussianGaussian","GaussianInvWishart","GaussianNIW", "GaussianNIG", "CatDirichlet", and "CatDP". See \code{?BasicBayesian} for definition of "BasicBayesian" objects, and see for example \code{?GaussianGaussian} for specific "BasicBayesian" instances. As a summary, An "HDP" object is simply a combination of a "CatHDP" object (see \code{?CatHDP}) and an object of any "BasicBayesian" type.\cr
#' In the case of HDP, z and k can only be positive integers. \cr
#' The model structure and prior parameters are stored in a "HDP" object. \cr
#' Posterior predictive density = p(x,z,k|gamma,alpha,psi) when x is not NULL, or p(z,k|gamma,alpha,psi) when x is NULL.
#'
#' @seealso \code{\link{HDP}}, \code{\link{dPosteriorPredictive.HDP}}, \code{\link{marginalLikelihood.HDP}}
#' @param obj A "HDP" object.
#' @param x Random samples of the "BasicBayesian" object.
#' @param z integer.
#' @param k integer, the partition label of the parameter space where the observation x is drawn from.
#' @param j integer, group label.
#' @param LOG Return the log density if set to "TRUE".
#' @param ... Additional arguments to be passed to other inherited types.
#' @return A numeric vector, the posterior predictive density.
#' @export
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
dPosteriorPredictive.HDP <- function(obj,x=NULL,z,k,j,LOG=TRUE,...){
    pzk <- dPosteriorPredictive(obj = obj$Z,z=z,k=k,j=j,LOG = FALSE)
    if(!is.null(x)){
        if(is.vector(x)){
            x <- matrix(x, ncol = 1)
        }else if(!.is(x,"matrix")){
            stop("'x' must be a vector(univariate) or a matrix(multivariate)!")
        }
        if(nrow(x)>1L) stop("There can only be one observation in 'x', so 'x' must be a matrix of only 1 row.")
        if(length(obj$Z$Z1$gamma$nk)>0){
            lpx <- vapply(seq_along(obj$Z$Z1$gamma$nk),function(ki){
                if(ki>0) dPosteriorPredictive(obj = obj$X[[ki]],x = x,LOG = TRUE)
                else as.numeric(NA)
            },FUN.VALUE = numeric(1),USE.NAMES = FALSE)
            lpx <- lpx[k]
            lpx[is.na(lpx)] <- dPosteriorPredictive(obj = obj$H,x = x,LOG = TRUE)
        }else{
            lpx <- rep(dPosteriorPredictive(obj = obj$H,x = x,LOG = TRUE),length(k))
        }
        out <- log(pzk)+lpx             #log probability
        if(LOG) return(out)
        else return(exp(out))
    }else{
        if(LOG) return(log(pzk))
        else return(pzk)
    }
}

#' @title Generate random samples from the posterior predictive distribution of a "HDP" object
#' @description
#' Generate random samples from the posterior predictive distribution of the following structure:
#'      \deqn{G_j|gamma \sim DP(gamma,U), j = 1:J}
#'      \deqn{pi_j|G_j,alpha \sim DP(alpha,G_j)}
#'      \deqn{z|pi_j \sim Categorical(pi_j)}
#'      \deqn{k|z,G_j \sim Categorical(G_j), if z is a sample from the base measure G}
#'      \deqn{theta_k|psi \sim H0(psi)}
#'      \deqn{x|theta_k,k \sim F(theta_k)}
#' where DP(gamma,U) is a Dirichlet Process on positive integers, gamma is the "concentration parameter", U is the "base measure" of this Dirichlet process, U is an uniform distribution on all positive integers.  DP(gamma,G_j) is a Dirichlet Process on integers with concentration parameter alpha and base measure G_j. The choice of F() and H0() can be described by an arbitrary "BasicBayesian" object such as "GaussianGaussian","GaussianInvWishart","GaussianNIW", "GaussianNIG", "CatDirichlet", and "CatDP". See \code{?BasicBayesian} for definition of "BasicBayesian" objects, and see for example \code{?GaussianGaussian} for specific "BasicBayesian" instances. As a summary, An "HDP" object is simply a combination of a "CatHDP" object (see \code{?CatHDP}) and an object of any "BasicBayesian" type.\cr
#' In the case of HDP, z and k can only be positive integers. \cr
#' The model structure and prior parameters are stored in a "HDP" object. \cr
#' This function will generate random samples from the distribution z,k|gamma,alpha,psi,x.
#' 
#' @seealso \code{\link{HDP}}, \code{\link{dPosteriorPredictive.HDP}}
#' @param obj A "HDP" object.
#' @param n integer, number of samples.
#' @param x Random samples of the "BasicBayesian" object.
#' @param j integer, group label.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return integer, the categorical samples.
#' @export
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
rPosteriorPredictive.HDP <- function(obj,n=1,x,j,...){
    if(missing(x)|missing(j)) stop("'x' and 'j'must be specified")
    if(n>1) stop("for now only support n=1")
    allz <- which(obj$Z$Z2[[j]]$gamma$nk>0)
    allk <- which(obj$Z$Z1$gamma$nk>0)
    zs <- c(allz,rep(obj$Z$Z2[[j]]$gamma$newLabel,length(allk)+1L))
    ks <- c(obj$Z$Z12map[[j]][allz],allk,obj$Z$Z1$gamma$newLabel)
    probs <- dPosteriorPredictive(obj = obj,x=x,z=zs,k=ks,j=j,LOG = TRUE)
    probs <- exp(probs-max(probs))
    idx <- sample.int(length(zs),size = n,replace = TRUE,prob = probs)
    c(z=zs[idx],k=ks[idx])
}

#' @title Create objects of type "HDP2".
#' @description
#' Create an object of type "HDP2", which represents the Hierarchical-Dirichlet-Process with two Dirichlet-Process hierarchies:
#'      \deqn{G_m |eta \sim DP(eta,U), m = 1:M}
#'      \deqn{G_mj|gamma,G_m \sim DP(gamma,G_m), j = 1:J_m}
#'      \deqn{pi_mj|G_mj,alpha \sim DP(alpha,G_mj)}
#'      \deqn{z|pi_mj \sim Categorical(pi_mj)}
#'      \deqn{k|z,G_mj \sim Categorical(G_mj), if z is a sample from the base measure G_mj}
#'      \deqn{u|k,G_m \sim Categorical(G_m), if k is a sample from the base measure G_m}
#'      \deqn{theta_u|psi \sim H0(psi)}
#'      \deqn{x|theta_u,u \sim F(theta_u)}
#' where DP(eta,U) is a Dirichlet Process on positive integers, eta is the "concentration parameter", U is the "base measure" of this Dirichlet process, U is an uniform distribution on all positive integers.  DP(gamma,G_m) is a Dirichlet Process on integers with concentration parameter gamma and base measure G_m. DP(alpha,G_mj) is a Dirichlet Process on integers with concentration parameter alpha and base measure G_mj. The choice of F() and H0() can be described by an arbitrary "BasicBayesian" object such as "GaussianGaussian","GaussianInvWishart","GaussianNIW", "GaussianNIG", "CatDirichlet", and "CatDP". See \code{?BasicBayesian} for definition of "BasicBayesian" objects, and see for example \code{?GaussianGaussian} for specific "BasicBayesian" instances. As a summary, An "HDP2" object is simply a combination of a "CatHDP2" object (see \code{?CatHDP2}) and an object of any "BasicBayesian" type.\cr
#' In the case of HDP2, u, z and k can only be positive integers. \cr
#' This object will be used as a place for recording and accumulating information in the related inference/sampling functions such as posterior(), posteriorDiscard(), MAP() and so on.
#' 
#' @seealso \code{\link{BasicBayesian}},\code{\link{GaussianNIW}},\code{\link{GaussianNIG}},\code{\link{CatDirichlet}},\code{\link{CatHDP2}},\code{\link{posterior.HDP2}},\code{\link{posteriorDiscard.HDP2}},\code{\link{marginalLikelihood.HDP2}} ...
#' @param objCopy an object of type "HDP2". If "objCopy" is not NULL, the function create a new "HDP2" object by copying the content from objCopy, otherwise this new object will be created by using "ENV" and "gamma". Default NULL.
#' @param ENV environment, specify where the object will be created.
#' @param gamma list, a named list of parameters, gamma=list(eta,gamma,alpha,m,j,H0aF,parH0), where gamma$eta is a numeric value specifying the concentration parameter of DP(eta,U), gamma$gamma is a numeric value specifying the concentration parameter of DP(gamma,G_m), gamma$alpha is a numeric value specifying the concentration parameter of DP(alpha,G_mj), gamma$m is the number of groups M, gamma$j is the number of subgroups in each group, must satisfy length(gamma$j)=gamma$m. gamma$H0aF is the name of the "BasicBayesian" object which specifies the structure of H0() and F(). gamma$parH0 is the parameters passed to the selected H0aF. For example, if gamma$H0aF="GaussianNIW", then gamma$parH0 should be a named list of NIW parameters: gamma$parH0=list(m,k,v,S), where gamma$parH0$m is a numeric "location" parameter; gamma$parH0$S is a symmetric positive definite matrix representing the "scale" parameters; gamma$parH0$k and gamma$parH0$v are numeric values.
#' @return An object of class "HDP2".
#' @export
#' @examples
#' \donttest{
#'
#' ## This is an example of Gibbs sampling on a hierarchical mixture model, using HDP2.
#'
#' ## load some hierarchical mixture data, check ?mmhhData for details.
#' data(mmhhData)
#' x <- mmhhData$x
#' ms <- mmhhData$groupLabel
#' js <- mmhhData$subGroupLabel
#'
#' ## Step1: initialize--------------------------------------------------
#' maxit <- 50                            #iterative for maxit times
#' z <- rep(1L,nrow(x))
#' k <- rep(1L,nrow(x))
#' u <- rep(1L,nrow(x))
#' obj <- HDP2(gamma = list(eta=1,gamma=1,alpha=1,m=2L,j=c(10L,20L),
#'             H0aF="GaussianNIW",
#'             parH0=list(m=c(0,0),k=0.001,v=2,S=diag(2)*0.001)))
#' ss <- sufficientStatistics(obj$H,x=x,foreach = TRUE) #sufficient statistics
#' N <- length(ss)
#' for(i in 1L:N){                         #initialize z k and u
#'    tmp <- rPosteriorPredictive(obj = obj,n=1,x=x[i,,drop=FALSE],m=ms[i],j=js[i])
#'    z[i] <- tmp["z"]
#'    k[i] <- tmp["k"]
#'    u[i] <- tmp["u"]
#'    posterior(obj = obj,ss = ss[[i]],ss1 = u[i],ss2 = k[i],ss3 = z[i],m=ms[i],j = js[i])
#' }
#'
#' ## Step2: main Gibbs loop---------------------------------------------
#' it <- 0                                 #iteration tracker
#' pb <- txtProgressBar(min = 0,max = maxit,style = 3)
#' while(TRUE){
#'    for(i in 1L:N){
#'        ##remove the sample from the posterior info
#'        posteriorDiscard(obj = obj,ss = ss[[i]],ss1=u[i],ss2=k[i],ss3 = z[i],m=ms[i],j=js[i])
#'        ##resample a new partition
#'        tmp <- rPosteriorPredictive(obj = obj,n=1L,x=x[i,,drop=FALSE],m=ms[i],j=js[i])
#'        z[i] <- tmp["z"]
#'        k[i] <- tmp["k"]
#'        u[i] <- tmp["u"]
#'        posterior(obj = obj,ss = ss[[i]], ss1=u[i],ss2 = k[i],ss3 = z[i],m=ms[i],j=js[i])
#'    }
#'    plot(x=x[,1],y=x[,2],col=u)
#'    it <- it+1
#'    setTxtProgressBar(pb,it)
#'    if(it>=maxit){cat("\n");break}
#' }
#'
#' }
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
HDP2 <- function(objCopy=NULL,ENV=parent.frame(),gamma=list(
                                                     eta=1,
                                                     gamma=1, #concentration parameter for prior DP1
                                                     alpha=1,     #concentration parameter for DP2
                                                     m=3L,
                                                     j=c(2,3,4),     #number of groups
                                                     
                                                     H0aF="GaussianNIW", #name of the base measure and the observation distribution
                                                     parH0=list(m=0,k=1,v=2,S=1) #pars for the basemeasure
                                                 )){
    object <- new.env(parent=ENV)
    if(!is.null(objCopy)){
        if(!.is(objCopy,"HDP2")) stop("'objCopy' must be of class 'HDP2'")
        object$gamma <- objCopy$gamma
        object$Z <- objCopy$Z
        object$H <- objCopy$H
        object$X <- objCopy$X
    }else{
        if(!missing(gamma))
            if((!is.list(gamma)) |
               (!all(names(gamma) %in% c("eta","gamma","alpha","m","j","H0aF","parH0"))))
                stop("gamma must be of list(eta,gamma,alpha,m,j,H0aF,parH0)")
        if(!is.character(gamma$H0aF) | length(gamma$H0aF)!=1) stop("'H0aF' must be a length 1 character specifying the name of the BasicBayesian class")
        if(!is.integer(gamma$j)) gamma$j <- as.integer(gamma$j)
        if(!is.integer(gamma$m)) gamma$m <- as.integer(gamma$m)
        object$gamma <- gamma
        object$Z <- CatHDP2(ENV = object,gamma = list(eta=gamma$eta,gamma=gamma$gamma,alpha=gamma$alpha,m=gamma$m,j=gamma$j))
        object$H <- do.call(gamma$H0aF,list(ENV = object,gamma=gamma$parH0)) #a object of class specified by H0aF for observation distribution, new() will create the object according to class representation of setClass(), so use do.call() instead
        object$X <- list()                                      #the list of observation distributions for each partition
    }
    class(object) <- c('HDP2',"BayesianBrick")
    return(object)
}

#' print the content of a "HDP2" object
#' @param x An object of type "HDP2".
#' @param ... Additional arguments to be passed to other inherited types.
#' @return None.
#' @export
print.HDP2 <- function(x,...){
    cat("The partition distribution is governed by a 'CatHDP2' class with following parameters:\n")
    print(x$Z)
    cat("\n\n\n")
    cat("The observation distribution is governed by a '",x$gamma$H0aF,"' object with following parameters::\n",sep = "")
    print(x$H)
    cat("\n\n\n")
    cat("The observation distributions of different partitions is governed by following '",x$gamma$H0aF,"' objects:\n",sep = "")
    print(x$X)
}

#' @title Sufficient statistics of a "HDP2" object
#' @description
#' For following model structure:
#'      \deqn{G_m |eta \sim DP(eta,U), m = 1:M}
#'      \deqn{G_mj|gamma,G_m \sim DP(gamma,G_m), j = 1:J_m}
#'      \deqn{pi_mj|G_mj,alpha \sim DP(alpha,G_mj)}
#'      \deqn{z|pi_mj \sim Categorical(pi_mj)}
#'      \deqn{k|z,G_mj \sim Categorical(G_mj), if z is a sample from the base measure G_mj}
#'      \deqn{u|k,G_m \sim Categorical(G_m), if k is a sample from the base measure G_m}
#'      \deqn{theta_u|psi \sim H0(psi)}
#'      \deqn{x|theta_u,u \sim F(theta_u)}
#' where DP(eta,U) is a Dirichlet Process on positive integers, eta is the "concentration parameter", U is the "base measure" of this Dirichlet process, U is an uniform distribution on all positive integers.  DP(gamma,G_m) is a Dirichlet Process on integers with concentration parameter gamma and base measure G_m. DP(alpha,G_mj) is a Dirichlet Process on integers with concentration parameter alpha and base measure G_mj. The choice of F() and H0() can be described by an arbitrary "BasicBayesian" object such as "GaussianGaussian","GaussianInvWishart","GaussianNIW", "GaussianNIG", "CatDirichlet", and "CatDP". See \code{?BasicBayesian} for definition of "BasicBayesian" objects, and see for example \code{?GaussianGaussian} for specific "BasicBayesian" instances. As a summary, An "HDP2" object is simply a combination of a "CatHDP2" object (see \code{?CatHDP2}) and an object of any "BasicBayesian" type.\cr
#' In the case of HDP2, u, z and k can only be positive integers. \cr
#' The sufficient statistics of a set of samples x in a "HDP2" object is the same sufficient statistics of the "BasicBayesian" inside the "HDP2", see examples.
#' @seealso \code{\link{HDP2}}, \code{\link{sufficientStatistics_Weighted.HDP2}} 
#' @param obj A "HDP2" object.
#' @param x Random samples of the "BasicBayesian" object.
#' @param ... further arguments passed to the corresponding sufficientStatistics method of the "BasicBayesian" object.
#' @return Return the sufficient statistics of the corresponding BasicBayesian type, see examples.
#' @export
#' @examples
#' ## a HDP2 with Gaussian NIW observations
#' obj1 <- HDP2(gamma=list(gamma=1,alpha=1,j=2,m=2,
#'                         H0aF="GaussianNIW",
#'                         parH0=list(m=0,k=1,v=2,S=1)))
#' ## a HDP2 with Categorical-Dirichlet observations
#' obj2 <- HDP2(gamma=list(gamma=1,alpha=1,j=2,m=2,
#'                         H0aF="CatDirichlet",
#'                         parH0=list(alpha=1,uniqueLabels=letters[1:3])))
#' x1 <- rnorm(100)
#' x2 <- sample(letters[1:3],100,replace = TRUE)
#' sufficientStatistics(obj = obj1,x=x1,foreach = FALSE)
#' sufficientStatistics(obj = obj2,x=x2,foreach = FALSE)
#' sufficientStatistics(obj = obj1,x=x1,foreach = TRUE)
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
sufficientStatistics.HDP2 <- function(obj,x,...){
    sufficientStatistics(obj = obj$H,x=x,...)
}

#' @title Weighted sufficient statistics of a "HDP2" object
#' @description
#' For following model structure:
#'      \deqn{G_m |eta \sim DP(eta,U), m = 1:M}
#'      \deqn{G_mj|gamma,G_m \sim DP(gamma,G_m), j = 1:J_m}
#'      \deqn{pi_mj|G_mj,alpha \sim DP(alpha,G_mj)}
#'      \deqn{z|pi_mj \sim Categorical(pi_mj)}
#'      \deqn{k|z,G_mj \sim Categorical(G_mj), if z is a sample from the base measure G_mj}
#'      \deqn{u|k,G_m \sim Categorical(G_m), if k is a sample from the base measure G_m}
#'      \deqn{theta_u|psi \sim H0(psi)}
#'      \deqn{x|theta_u,u \sim F(theta_u)}
#' where DP(eta,U) is a Dirichlet Process on positive integers, eta is the "concentration parameter", U is the "base measure" of this Dirichlet process, U is an uniform distribution on all positive integers.  DP(gamma,G_m) is a Dirichlet Process on integers with concentration parameter gamma and base measure G_m. DP(alpha,G_mj) is a Dirichlet Process on integers with concentration parameter alpha and base measure G_mj. The choice of F() and H0() can be described by an arbitrary "BasicBayesian" object such as "GaussianGaussian","GaussianInvWishart","GaussianNIW", "GaussianNIG", "CatDirichlet", and "CatDP". See \code{?BasicBayesian} for definition of "BasicBayesian" objects, and see for example \code{?GaussianGaussian} for specific "BasicBayesian" instances. As a summary, An "HDP2" object is simply a combination of a "CatHDP2" object (see \code{?CatHDP2}) and an object of any "BasicBayesian" type.\cr
#' In the case of HDP2, u, z and k can only be positive integers. \cr
#' The sufficient statistics of a set of samples x in a "HDP2" object is the same sufficient statistics of the "BasicBayesian" inside the "HDP2", see examples.
#' @seealso \code{\link{HDP2}}, \code{\link{sufficientStatistics.HDP2}} 
#' @param obj A "HDP2" object.
#' @param x Random samples of the "BasicBayesian" object.
#' @param w numeric, sample weights.
#' @param ... further arguments passed to the corresponding sufficientStatistics method of the "BasicBayesian" object.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return Return the sufficient statistics of the corresponding BasicBayesian type, see examples.
#' @export
#' @examples
#' ## a HDP2 with Gaussian NIW observations
#' obj1 <- HDP2(gamma=list(gamma=1,alpha=1,j=2,m=2,
#'                         H0aF="GaussianNIW",
#'                         parH0=list(m=0,k=1,v=2,S=1)))
#' ## a HDP2 with Categorical-Dirichlet observations
#' obj2 <- HDP2(gamma=list(gamma=1,alpha=1,j=2,m=2,
#'                         H0aF="CatDirichlet",
#'                         parH0=list(alpha=1,uniqueLabels=letters[1:3])))
#' x1 <- rnorm(100)
#' x2 <- sample(letters[1:3],100,replace = TRUE)
#' w <- runif(100)
#' sufficientStatistics_Weighted(obj = obj1,x=x1,w=w,foreach = FALSE)
#' sufficientStatistics_Weighted(obj = obj2,x=x2,w=w,foreach = FALSE)
#' sufficientStatistics_Weighted(obj = obj1,x=x1,w=w,foreach = TRUE)
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
sufficientStatistics_Weighted.HDP2 <- function(obj,x,w,...){
    sufficientStatistics_Weighted(obj = obj$H,x=x,w=w,...)
}

#' @title Update a "HDP2" object with sample sufficient statistics
#' @description
#' For the model structure:
#'      \deqn{G_m |eta \sim DP(eta,U), m = 1:M}
#'      \deqn{G_mj|gamma,G_m \sim DP(gamma,G_m), j = 1:J_m}
#'      \deqn{pi_mj|G_mj,alpha \sim DP(alpha,G_mj)}
#'      \deqn{z|pi_mj \sim Categorical(pi_mj)}
#'      \deqn{k|z,G_mj \sim Categorical(G_mj), if z is a sample from the base measure G_mj}
#'      \deqn{u|k,G_m \sim Categorical(G_m), if k is a sample from the base measure G_m}
#'      \deqn{theta_u|psi \sim H0(psi)}
#'      \deqn{x|theta_u,u \sim F(theta_u)}
#' where DP(eta,U) is a Dirichlet Process on positive integers, eta is the "concentration parameter", U is the "base measure" of this Dirichlet process, U is an uniform distribution on all positive integers.  DP(gamma,G_m) is a Dirichlet Process on integers with concentration parameter gamma and base measure G_m. DP(alpha,G_mj) is a Dirichlet Process on integers with concentration parameter alpha and base measure G_mj. The choice of F() and H0() can be described by an arbitrary "BasicBayesian" object such as "GaussianGaussian","GaussianInvWishart","GaussianNIW", "GaussianNIG", "CatDirichlet", and "CatDP". See \code{?BasicBayesian} for definition of "BasicBayesian" objects, and see for example \code{?GaussianGaussian} for specific "BasicBayesian" instances. As a summary, An "HDP2" object is simply a combination of a "CatHDP2" object (see \code{?CatHDP2}) and an object of any "BasicBayesian" type.\cr
#' In the case of HDP2, u, z and k can only be positive integers. \cr
#' This function will update the prior knowledge by adding the information of newly observed samples x, z and k. The model structure and prior parameters are stored in a "HDP2" object, the prior parameters in this object will be updated after running this function.
#'
#' @seealso \code{\link{HDP2}},\code{\link{posteriorDiscard.HDP2}},\code{\link{sufficientStatistics.HDP2}}
#' @param obj A "HDP2" object.
#' @param ss Sufficient statistics of x of the "BasicBayesian" object, must be a list of sufficient statistics for each of the observations. Use sufficientStatistics(...,foreach=TRUE) to generate ss.
#' @param ss1 Sufficient statistics of u. In HDP2 case the sufficient statistic of sample u is u itself(if u is a integer vector with all positive values).
#' @param ss2 Sufficient statistics of k. In HDP2 case the sufficient statistic of sample k is k itself(if k is a integer vector with all positive values).
#' @param ss3 Sufficient statistics of z. In HDP2 case the sufficient statistic of sample z is z itself(if z is a integer vector with all positive values).
#' @param m integer, group label.
#' @param j integer, subgroup label.
#' @param w Sample weights, default NULL.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return None. the model stored in "obj" will be updated based on "ss", "ss1", "ss2"and "ss3".
#' @export
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
posterior.HDP2 <- function(obj,ss,ss1,ss2,ss3,m,j,w=NULL,...){
    if(length(ss2)>1L) stop("posterior.HDP2 can only update from observations one at a time, for now")
    if(length(obj$X)<ss1)
        obj$X <- c(obj$X,
                   replicate(n = ss1-length(obj$X), expr = do.call(obj$gamma$H0aF,list(objCopy=obj$H,ENV=obj)),simplify = FALSE))
    posterior(obj = obj$Z,ss1 = ss1,ss2 = ss2,ss3 = ss3,m=m,j = j,w = w)
    posterior(obj = obj$X[[ss1]],ss = ss,w=w)
}

#' @title Update a "HDP2" object with sample sufficient statistics
#' @description
#' For the model structure:
#'      \deqn{G_m |eta \sim DP(eta,U), m = 1:M}
#'      \deqn{G_mj|gamma,G_m \sim DP(gamma,G_m), j = 1:J_m}
#'      \deqn{pi_mj|G_mj,alpha \sim DP(alpha,G_mj)}
#'      \deqn{z|pi_mj \sim Categorical(pi_mj)}
#'      \deqn{k|z,G_mj \sim Categorical(G_mj), if z is a sample from the base measure G_mj}
#'      \deqn{u|k,G_m \sim Categorical(G_m), if k is a sample from the base measure G_m}
#'      \deqn{theta_u|psi \sim H0(psi)}
#'      \deqn{x|theta_u,u \sim F(theta_u)}
#' where DP(eta,U) is a Dirichlet Process on positive integers, eta is the "concentration parameter", U is the "base measure" of this Dirichlet process, U is an uniform distribution on all positive integers.  DP(gamma,G_m) is a Dirichlet Process on integers with concentration parameter gamma and base measure G_m. DP(alpha,G_mj) is a Dirichlet Process on integers with concentration parameter alpha and base measure G_mj. The choice of F() and H0() can be described by an arbitrary "BasicBayesian" object such as "GaussianGaussian","GaussianInvWishart","GaussianNIW", "GaussianNIG", "CatDirichlet", and "CatDP". See \code{?BasicBayesian} for definition of "BasicBayesian" objects, and see for example \code{?GaussianGaussian} for specific "BasicBayesian" instances. As a summary, An "HDP2" object is simply a combination of a "CatHDP2" object (see \code{?CatHDP2}) and an object of any "BasicBayesian" type.\cr
#' In the case of HDP2, u, z and k can only be positive integers. \cr
#' Contrary to posterior(), this function will update the prior knowledge by removing the information of observed samples x. The model structure and prior parameters are stored in a "CatDP" object, the prior parameters in this object will be updated after running this function.
#'
#' @seealso \code{\link{HDP2}},\code{\link{posteriorDiscard.HDP2}},\code{\link{sufficientStatistics.HDP2}}
#' @param obj A "HDP2" object.
#' @param ss Sufficient statistics of x of the "BasicBayesian" object, must be a list of sufficient statistics for each of the observations. Use sufficientStatistics(...,foreach=TRUE) to generate ss.
#' @param ss1 Sufficient statistics of u. In HDP2 case the sufficient statistic of sample u is u itself(if u is a integer vector with all positive values).
#' @param ss2 Sufficient statistics of k. In HDP2 case the sufficient statistic of sample k is k itself(if k is a integer vector with all positive values).
#' @param ss3 Sufficient statistics of z. In HDP2 case the sufficient statistic of sample z is z itself(if z is a integer vector with all positive values).
#' @param m integer, group label.
#' @param j integer, subgroup label.
#' @param w Sample weights, default NULL.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return None. the model stored in "obj" will be updated based on "ss", "ss1", "ss2"and "ss3".
#' @export
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
posteriorDiscard.HDP2 <- function(obj,ss,ss1,ss2,ss3,m,j,w=NULL,...){
    if(missing(ss)|missing(j)|missing(m)) stop("'ss', 'j' and 'm' must all be specified")
    posteriorDiscard(obj = obj$Z,ss1 = ss1,ss2 = ss2,ss3 = ss3,m = m,j = j,w = w)
    posteriorDiscard(obj = obj$X[[ss1]],ss = ss,w=w)
}

#' @title Marginal likelihood for HDP2
#' @param obj an object.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return numeric, the marginal likelihood.
#' @export
marginalLikelihood.HDP2 <- function(obj,...){
    stop("marginalLikelihood for this type not implemented yet")
}

#' @title Posterior predictive density function of a "HDP2" object
#' @description
#' Generate the the density value of the posterior predictive distribution of the following structure:
#'      \deqn{G_m |eta \sim DP(eta,U), m = 1:M}
#'      \deqn{G_mj|gamma,G_m \sim DP(gamma,G_m), j = 1:J_m}
#'      \deqn{pi_mj|G_mj,alpha \sim DP(alpha,G_mj)}
#'      \deqn{z|pi_mj \sim Categorical(pi_mj)}
#'      \deqn{k|z,G_mj \sim Categorical(G_mj), if z is a sample from the base measure G_mj}
#'      \deqn{u|k,G_m \sim Categorical(G_m), if k is a sample from the base measure G_m}
#'      \deqn{theta_u|psi \sim H0(psi)}
#'      \deqn{x|theta_u,u \sim F(theta_u)}
#' where DP(eta,U) is a Dirichlet Process on positive integers, eta is the "concentration parameter", U is the "base measure" of this Dirichlet process, U is an uniform distribution on all positive integers.  DP(gamma,G_m) is a Dirichlet Process on integers with concentration parameter gamma and base measure G_m. DP(alpha,G_mj) is a Dirichlet Process on integers with concentration parameter alpha and base measure G_mj. The choice of F() and H0() can be described by an arbitrary "BasicBayesian" object such as "GaussianGaussian","GaussianInvWishart","GaussianNIW", "GaussianNIG", "CatDirichlet", and "CatDP". See \code{?BasicBayesian} for definition of "BasicBayesian" objects, and see for example \code{?GaussianGaussian} for specific "BasicBayesian" instances. As a summary, An "HDP2" object is simply a combination of a "CatHDP2" object (see \code{?CatHDP2}) and an object of any "BasicBayesian" type.\cr
#' In the case of HDP2, u, z and k can only be positive integers. \cr
#' The model structure and prior parameters are stored in a "HDP2" object. \cr
#' Posterior predictive density = p(u,z,k,x|eta,gamma,alpha,psi) when x is not NULL, or p(u,z,k|eta,gamma,alpha,psi) when x is NULL.
#'
#' @seealso \code{\link{HDP2}}, \code{\link{dPosteriorPredictive.HDP2}}, \code{\link{marginalLikelihood.HDP2}}
#' @param obj A "HDP2" object.
#' @param x Random samples of the "BasicBayesian" object.
#' @param z integer.
#' @param k integer.
#' @param u integer, the partition label of the parameter space where the observation x is drawn from.
#' @param m integer, group label.
#' @param j integer, subgroup label.
#' @param LOG Return the log density if set to "TRUE".
#' @param ... Additional arguments to be passed to other inherited types.
#' @return A numeric vector, the posterior predictive density.
#' @export
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
dPosteriorPredictive.HDP2 <- function(obj,x=NULL,u,k,z,m,j,LOG=TRUE,...){
    pzk <- dPosteriorPredictive(obj = obj$Z,u=u,z=z,k=k,m=m,j=j,LOG = FALSE)
    if(!is.null(x)){
        if(is.vector(x)){
            x <- matrix(x, ncol = 1)
        }else if(!.is(x,"matrix")){
            stop("'x' must be a vector(univariate) or a matrix(multivariate)!")
        }
        if(nrow(x)>1L) stop("There can only be one observation in 'x', so 'x' must be a matrix of only 1 row.")
        if(length(obj$Z$Z1$gamma$nk)>0L){
            lpx <- vapply(seq_along(obj$Z$Z1$gamma$nk),function(ui){
                if(ui>0) dPosteriorPredictive(obj = obj$X[[ui]],x = x,LOG = TRUE)
                else as.numeric(NA)
            },FUN.VALUE = numeric(1),USE.NAMES = FALSE)
            lpx <- lpx[u]
            lpx[is.na(lpx)] <- dPosteriorPredictive(obj = obj$H,x = x,LOG = TRUE)
        }else{
            lpx <- rep(dPosteriorPredictive(obj = obj$H,x = x,LOG = TRUE),length(k))
        }
        out <- log(pzk)+lpx             #log probability
        if(LOG) return(out)
        else return(exp(out))
    }else{
        if(LOG) return(log(pzk))
        else return(pzk)
    }
}

#' @title Generate random samples from the posterior predictive distribution of a "HDP2" object
#' @description
#' Generate random samples from the posterior predictive distribution of the following structure:
#'      \deqn{G_m |eta \sim DP(eta,U), m = 1:M}
#'      \deqn{G_mj|gamma,G_m \sim DP(gamma,G_m), j = 1:J_m}
#'      \deqn{pi_mj|G_mj,alpha \sim DP(alpha,G_mj)}
#'      \deqn{z|pi_mj \sim Categorical(pi_mj)}
#'      \deqn{k|z,G_mj \sim Categorical(G_mj), if z is a sample from the base measure G_mj}
#'      \deqn{u|k,G_m \sim Categorical(G_m), if k is a sample from the base measure G_m}
#'      \deqn{theta_u|psi \sim H0(psi)}
#'      \deqn{x|theta_u,u \sim F(theta_u)}
#' where DP(eta,U) is a Dirichlet Process on positive integers, eta is the "concentration parameter", U is the "base measure" of this Dirichlet process, U is an uniform distribution on all positive integers.  DP(gamma,G_m) is a Dirichlet Process on integers with concentration parameter gamma and base measure G_m. DP(alpha,G_mj) is a Dirichlet Process on integers with concentration parameter alpha and base measure G_mj. The choice of F() and H0() can be described by an arbitrary "BasicBayesian" object such as "GaussianGaussian","GaussianInvWishart","GaussianNIW", "GaussianNIG", "CatDirichlet", and "CatDP". See \code{?BasicBayesian} for definition of "BasicBayesian" objects, and see for example \code{?GaussianGaussian} for specific "BasicBayesian" instances. As a summary, An "HDP2" object is simply a combination of a "CatHDP2" object (see \code{?CatHDP2}) and an object of any "BasicBayesian" type.\cr
#' In the case of HDP2, u, z and k can only be positive integers. \cr
#' The model structure and prior parameters are stored in a "HDP2" object. \cr
#' This function will generate random samples from the distribution u,z,k|eta,gamma,alpha,psi,x.
#' 
#' @seealso \code{\link{HDP2}}, \code{\link{dPosteriorPredictive.HDP2}}
#' @param obj A "HDP2" object.
#' @param n integer, number of samples.
#' @param x Random samples of the "BasicBayesian" object.
#' @param m integer, group label.
#' @param j integer, subgroup label.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return integer, the categorical samples.
#' @export
#' @references Teh, Yee W., et al. "Sharing clusters among related groups: Hierarchical Dirichlet processes." Advances in neural information processing systems. 2005.
rPosteriorPredictive.HDP2 <- function(obj,n=1,x,m,j,...){
    if(missing(x)|missing(j)|missing(m)) stop("'x', 'm'and 'j'must be specified")
    if(n>1) stop("for now only support n=1")

    allz <- which(obj$Z$Z2[[m]]$Z2[[j]]$gamma$nk>0)
    allk <- which(obj$Z$Z2[[m]]$Z1$gamma$nk>0)
    allu <- which(obj$Z$Z1$gamma$nk>0)

    zs <- c(allz, rep(obj$Z$Z2[[m]]$Z2[[j]]$gamma$newLabel,length(allk)+length(allu)+1L) )
    tmp <- c(obj$Z$Z2[[m]]$Z12map[[j]][allz],allk)
    ks <- c(tmp, rep(obj$Z$Z2[[m]]$Z1$gamma$newLabel,length(allu)+1L))
    us <- c(obj$Z$Z12map[[m]][tmp],allu,rep(obj$Z$Z1$gamma$newLabel))

    probs <- dPosteriorPredictive.HDP2(obj = obj,x=x,u=us,k=ks,z=zs,m=m,j=j,LOG = TRUE)
    probs <- exp(probs-max(probs))
    idx <- sample.int(length(zs),size = n,replace = TRUE,prob = probs)
    c(u=us[idx],k=ks[idx],z=zs[idx])
}
