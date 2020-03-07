
#' @include Bayesian_Bricks.r







#' @title Inverse of a positive definite symmetric matrix
#' @description Use Cholesky decomposition to calculate the inverse where S = A'A, A is a upper diagonal matrix then inv(S) = inv(A)inv(A)'.
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

#' Inference in joint Gaussian distribution
#'
#' For the model structure \cr
#'    x1,x2|mu,Sigma ~ Gaussian(mu,Sigma) \cr
#'    x1|x2,mu,Sigma ~ Gaussian(mu12,Sigma12)
#'
#' @seealso \code{\link{linearGaussian}}
#' @param x2 numeric, an sample of X2, satisfying length(x2)<D, D is the dimension of the joint distribution.
#' @param mu numeric, length D mean vector. mu=c(mu_X1,mu_X2)/.
#' @param Sigma DxD covariance matrix. At least one of Sigma and Precision should be non-NULL.
#' @param Precision DxD precision matrix, satisfying Precision = inverse(Sigma). At least one of Sigma and Precision should be non-NULL.
#' @return A named list containing the conditional mean and covariance matrix.
#' @export
#' @examples
#' tmp <- matrix(runif(100),20,5)
#' S <- crossprod(tmp)                 #some synthetic covariance matrix
#' P <- solve(S)
#' m <- runif(5)
#' x2 <- runif(3)
#' inferenceJointGaussian(x2 = x2,mu = m,Precision = P)
inferenceJointGaussian <- function(x2,mu,Sigma=NULL,Precision=NULL){
    if(missing(x2)|missing(mu)) stop("'x2' or 'mu' not specified!")
    if(is.null(Sigma) & is.null(Precision)) stop("Error in inferenceJointGaussian(): At least one of Sigma and Precision should be non-NULL")
    if(is.null(Precision)) Precision <- pdsInverse(Sigma)
    D <- nrow(Precision)
    if(D!=length(mu)) stop("Error in inferenceJointGaussian(): Dimension doesn't mathc! nrow(Precision) != length(mu)")
    D1 <- D-length(x2)                  #dimension of X1
    if(D1<=0) stop("Error in inferenceJointGaussian(): length(mu2) should be strictly smaller than D")
    
    
    P11 <- Precision[1:D1,1:D1,drop=FALSE]
    P12 <- Precision[1:D1,(D1+1):D,drop=FALSE]

    Sigma12 <- pdsInverse(P11)
    mu12 <- as.vector(mu[1:D1] - Sigma12%*%P12%*%drop(x2-mu[(D1+1):D])) #use drop() to remove dimension of x, incase x is matrix of 1 row

    list(mu12=mu12,Sigma12=Sigma12)
}

#' Linear Gaussian systems
#'
#' For the model structure \cr
#'    x1 ~ Gaussian(mu1,Sigma1) \cr
#'    x2 ~ Gaussian(Ax1+b,Sigma21) \cr
#'    x1|x2,mu1,Sigma1,A,b,Sigma21 ~ Gaussian(mu12,Sigma12)
#'
#' @seealso \code{\link{inferenceJointGaussian}} 
#' @param x2 numberic, an sample of x2, length(x2)=d2.
#' @param mu1 numeric, marginal mean of x1, length(m1)=d1.
#' @param Sigma1 matrix, marginal covariance matrix of x1. At least one of Sigma1 and Precision1 should be non-NULL.
#' @param Precision1 marginal precision matrix of x1,satisfying Precision1 = inverse(Sigma1). At least one of Sigma1 and Precision1 should be non-NULL.
#' @param A matrix, dim(A) = d1 x d2.
#' @param b numeric, length(b)=d1.
#' @param Sigma21 matrix, dim(Sigma21) = d1 x d1, either Precision21 or Sigma21 should be non-NULL.
#' @param Precision21 matrix, satisfying Precision21 = inverse(Sigma21), either Precision21 or Sigma21 should be non-NULL.
#' @return A named list containing the conditional mean and covariance matrix.
#' @export
#' @examples
#' x2 <- rGaussian(1,mu = c(0,0.2),Sigma = diag(2))
#' mu1 <- c(0,0); Sigma1 <- diag(2); A <- matrix(c(2,1,1,1),2,2); b <- runif(2); Sigma21 <- diag(2)*0.2
#' linearGaussian(x2=x2,mu1 = mu1,Sigma1 = Sigma1,A=A,b=b,Sigma21 = Sigma21)
linearGaussian <- function(x2,mu1,Sigma1=NULL,Precision1=NULL,A,b,Sigma21=NULL,Precision21=NULL){
    if(missing(x2)|missing(mu1)) stop("'x2' or 'mu1' not specified!")
    if(is.null(Sigma1) & is.null(Precision1)) stop("Error in linearGaussian(): At least one of Sigma1 and Precision1 should be non-NULL")
    if(is.null(Precision1)) Precision1 <- pdsInverse(Sigma1)
    if(is.null(Sigma1)) Sigma1 <- pdsInverse(Precision1)
    if(is.null(Sigma21) & is.null(Precision21)) stop("Error in linearGaussian(): At least one of Sigma21 and Precision21 should be non-NULL")
    if(is.null(Precision21)) Precision21 <- pdsInverse(Sigma21)

    Precision12 <- Precision1 + t(A)%*%Precision21%*%A
    Sigma12 <- pdsInverse(Precision12)
    mu12 <- Sigma12 %*% (t(A)%*%Precision21%*%drop(x2-b) + Precision1%*%drop(mu1))

    list(
        mu12=mu12,
        Sigma12=Sigma12## ,
    )
}

#' Random generation for Gaussian distribution
#'
#' Generate random samples from a Gaussian distribution.
#'
#' @seealso \code{\link{dGaussian}}
#' @param n integer, number of samples.
#' @param mu numeric, mean vector.
#' @param Sigma matrix, covariance matrix, one of Sigma and A should be non-NULL.
#' @param A matrix, the Cholesky decomposition of Sigma, an upper triangular matrix, one of Sigma and A should be non-NULL.
#' @return A matrix of n rows and length(mu) columns.
#' @export
#' @examples
#' \donttest{
#' x <- rGaussian(1000,mu = c(1,1),Sigma = matrix(c(1,0.5,0.5,3),2,2))
#' plot(x)
#' }
#' @import stats
rGaussian <- function(n,mu,Sigma=NULL,A=NULL){
    if(missing(n)|missing(mu)) stop("'n' or 'mu' not specified!")
    d <- length(mu)                     #dimension
    if(is.null(A) & !is.null(Sigma)) A <- chol(Sigma)
    if(is.null(A) & is.null(Sigma)) stop("Error in rGaussian(): at least one of 'Sigma' and 'A' should be non-NULL!")
    X <- matrix(rnorm(n*d,mean = 0,sd = 1),nrow = d,ncol = n)
    X <- A%*%X + mu
    
    t(X)
}

#' Density function of Gaussian distribution
#'
#' Get the density of a set of  samples from a Gaussian distribution
#'
#' @seealso \code{\link{rGaussian}}
#' @param x matrix, when x is a numeric vector, it will be converted to a matrix with 1 column!
#' @param mu numeric, mean vector.
#' @param Sigma matrix, covariance matrix, one of Sigma and A should be non-NULL.
#' @param A matrix, the Cholesky decomposition of Sigma, an upper triangular matrix, one of Sigma and A should be non-NULL.
#' @param LOG logical, return log density of LOG=TRUE, default TRUE.
#' @return A numeric vector.
#' @export
#' @examples
#' \donttest{
#' plot(
#'    dGaussian(x=seq(-5,5,length.out = 1000),mu = 0,Sigma = 1,LOG = FALSE)
#'    ,type = "l"
#' )
#' }
dGaussian <- function(x,mu,Sigma=NULL,A=NULL,LOG=TRUE){
    if(missing(x)|missing(mu)) stop("'x' or 'mu' not specified!")

    if(is.vector(x)){
        x <- matrix(x,ncol = 1)
    }else if(!.is(x,"matrix")){
        stop("'x' must be a vector(for univariate t) or matrix(for multivariate t)!")
    }

    if(is.null(A) & !is.null(Sigma)) A <- chol(Sigma)
    if(is.null(A) & is.null(Sigma)) stop("Error in dGaussian(): at least one of 'Sigma' and 'A' should be non-NULL!")
    
    p <- ncol(x)
    b <- backsolve(A,diag(p))
    logp <- -(p/2)*log(2*pi) + sum(log(diag(b))) - 0.5*colSums((crossprod(b,(t(x)-mu)))^2)
    if(!LOG) logp <- exp(logp)
    logp
}

#' Random Generation for (multivariate) t distribution
#'
#' Generate random samples from a (multivariate) t distribution
#'
#' @seealso \code{\link{dT}}
#' @param n integer, number of samples.
#' @param mu numeric, mean vector.
#' @param Sigma matrix, covariance matrix, one of Sigma and A should be non-NULL.
#' @param A matrix, the Cholesky decomposition of Sigma, an upper triangular matrix, one of Sigma and A should be non-NULL.
#' @param df numeric, degrees of freedom.
#' @return A matrix of n rows and length(mu) columns, each row is a sample.
#' @export
#' @examples
#' \donttest{
#' x <- rT(1000,mu = c(1,1),Sigma = matrix(c(1,0.5,0.5,3),2,2))
#' plot(x)
#' }
#' @import stats
rT <- function(n,mu,Sigma=NULL,A=NULL,df=1){
    if(missing(n)|missing(mu)) stop("'n' and 'mu' must all be specified!")

    if(is.null(A) & !is.null(Sigma)) A <- chol(Sigma)
    if(is.null(A) & is.null(Sigma)) stop("Error in rT(): at least one of 'Sigma' and 'A' should be non-NULL!")
    
    if(df == 0 |is.infinite(df))
        return(rGaussian(n=n,mu=mu,A=A))
    if(df<0) stop("'df' can not be smaller than 0")
    
    X <- rGaussian(n=n, mu=rep(0,nrow(A)),A=A)/sqrt(rchisq(n,df)/df)
    return(sweep(X,2,mu,"+"))
}

#' Density function for (multivariate) t distribution
#'
#' Get the density of a set of samples from a t distribution
#'
#' @seealso \code{\link{rT}}
#' @param x matrix, when x is a numeric vector, it will be converted to a matrix with 1 column!
#' @param mu numeric, mean vector.
#' @param Sigma matrix, covariance matrix, one of Sigma and A should be non-NULL.
#' @param A matrix, the Cholesky decomposition of Sigma, an upper triangular matrix, one of Sigma and A should be non-NULL.
#' @param df numeric, degrees of freedom.
#' @param LOG logical, return log density of LOG=TRUE, default TRUE.
#' @return A numeric vector, the probability densities.
#' @export
#' @examples
#' \donttest{
#' plot(
#'    dT(x=seq(-5,5,length.out = 1000),mu = 0,Sigma = 1,LOG = FALSE)
#'    ,type = "l"
#' )
#' }
dT <- function(x,mu,Sigma=NULL,A=NULL,df=1,LOG=TRUE){
    if(missing(x)|missing(mu)) stop("'x' and 'mu' must all be specified!")
    
    if(is.vector(x)){
        x <- matrix(x,ncol = 1)
    }else if(!.is(x,"matrix")){
        stop("'x' must be a vector(for univariate t) or matrix(for multivariate t)!")
    }
    
    d <- ncol(x)                        #dimension
    if(df==0|is.infinite(df))
        return(dGaussian(x=x,mu=mu,Sigma = Sigma,A=A,LOG = LOG))
    if(df<0) stop("'df' can not be smaller than 0")

    if(is.null(A) & !is.null(Sigma)) A <- chol(Sigma)
    if(is.null(A) & is.null(Sigma)) stop("Error in dT(): at least one of 'Sigma' and 'A' should be non-NULL!")
    
    b <- backsolve(A, t(x)-mu,transpose = TRUE)
    b <- colSums(b^2)
    out <- lgamma((d + df)/2) - (lgamma(df/2)+sum(log(diag(A)))+d/2*log(pi*df)) - 0.5 * (df + d) * log1p(b/df)

    if(!LOG)  out <- exp(out)

    out
}

#' Create objects of type "GaussianNIW".
#'
#' Create an object of type "GaussianNIW", which represents the Gaussian-Normal-Inverse-Wishart (Gaussian-NIW) conjugate structure: \cr
#'      theta|gamma ~ NIW(gamma) \cr
#'      x|theta ~ Gaussian(theta) \cr
#' where theta = (mu,Sigma) is the Gaussian parameter, gamma = (m,k,v,S) is the Normal-Inverse-Wishart(NIW) parameter. \cr
#' This object will be used as a place for recording and accumulating information in the related inference/sampling functions such as posterior(), posteriorDiscard(), MAP() and so on.
#'
#' @seealso \code{\link{posterior.GaussianNIW}},\code{\link{posteriorDiscard.GaussianNIW}},\code{\link{MAP.GaussianNIW}},\code{\link{marginalLikelihood.GaussianNIW}} ...
#' @param objCopy An object of type "GaussianNIW". If "objCopy" is not NULL, the function create a new "GaussianNIW" object by copying the content from objCopy, otherwise this new object will be created by using "ENV" and "gamma". Default NULL.
#' @param ENV environment, specify in which environment the object will be created.
#' @param gamma list, a named list of NIW parameters, gamma=list(m,k,v,S). Where gamma$m is a numeric "location" parameter; gamma$S is a symmetric positive definite matrix representing the "scale" parameters; gamma$k and gamma$v are numeric values.
#' @return An object of class "GaussianNIW".
#' @export
#' @examples
#' obj <- GaussianNIW(gamma=list(m=c(0,1),k=0.0001,v=2,S=diag(2)))
#' obj #print the content
#' @references Murphy, Kevin P. "Conjugate Bayesian analysis of the Gaussian distribution." def 1.22 (2007): 16.
#' @references Gelman, Andrew, et al. "Bayesian Data Analysis Chapman & Hall." CRC Texts in Statistical Science (2004).
GaussianNIW <- function(objCopy=NULL,ENV=parent.frame(),gamma=list(m=0,k=1,v=2,S=1)){
    object <- BasicBayesian(ENV = ENV)
    if(!is.null(objCopy)){
        if(!.is(objCopy,"GaussianNIW")) stop("'objCopy' must be of class 'GaussianNIW'")
        object$gamma <- objCopy$gamma
        object$H <- objCopy$H
        object$F <- objCopy$F
    }else{
        if(!missing(gamma))
            if((!is.list(gamma)) |
               (!all(names(gamma) %in% c("m","k","v","S"))))
                stop("gamma must be of list(m,k,v,S)")
        object$gamma <- gamma
        object$H <- "NIW"
        object$F <- "Gaussian"
    }
    class(object) <- c("GaussianNIW",class(object))
    return(object)
}

#' Sufficient statistics of a "GaussianNIW" object
#'
#' For following Gaussian-NIW model structure: \cr
#'      theta|gamma ~ NIW(gamma) \cr
#'      x|theta ~ Gaussian(theta) \cr
#' where theta = (mu,Sigma) is the Gaussian parameter, gamma = (m,k,v,S) is the Normal-Inverse-Wishart(NIW) parameter. \cr
#' The sufficient statistics of a set of samples x (each row of x is a sample) are: \cr
#'  the effective number of samples N=nrow(x); \cr
#'  the sample sum xsum = colSums(x); \cr
#'  the uncentered scatter matrix S = t(x)%*%x.
#'
#' @seealso \code{\link{GaussianNIW}}, \code{\link{sufficientStatistics_Weighted.GaussianNIW}}
#' @param obj A "GaussianNIW" object.
#' @param x matrix, Gaussian samples, when x is a matrix, each row is a sample of dimension ncol(x). when x is a vector, x is length(x) samples of dimension 1.
#' @param foreach logical, if foreach=TRUE, will return a list of sufficient statistics for each row of x, otherwise will return the sufficient statistics of x as a whole.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return If foreach=TRUE, will return a list of sufficient statistics for each row of x, otherwise will return the sufficient statistics of x as a whole.
#' @export
#' @examples
#' x <- rGaussian(10,mu = c(-1.5,1.5),Sigma = matrix(c(0.1,0.03,0.03,0.1),2,2))
#' obj <- GaussianNIW()                    #an GaussianNIW object
#' sufficientStatistics(obj=obj,x=x,foreach = FALSE)
#' sufficientStatistics(obj=obj,x=x,foreach = TRUE)
#' @references Murphy, Kevin P. "Conjugate Bayesian analysis of the Gaussian distribution." def 1.22 (2007): 16.
#' @references Gelman, Andrew, et al. "Bayesian Data Analysis Chapman & Hall." CRC Texts in Statistical Science (2004).
sufficientStatistics.GaussianNIW <- function(obj,x,foreach=FALSE,...){
    if(missing(x)) stop("'x' must be specified")
    if(is.vector(x)){
        x <- matrix(x, ncol = 1)
    }else if(!.is(x,"matrix")){
        stop("'x' must be a vector(for univariate t) or matrix(for multivariate t)!")
    }
    if(foreach){
        sapply(1:nrow(x),function(i){
            ss <- list(N=1,
                       xsum=x[i,,drop=TRUE],
                       S=t(x[i,,drop=FALSE])%*%x[i,,drop=FALSE])
            class(ss) <- "ssGaussian"
            ss
        },simplify = FALSE,USE.NAMES = FALSE)
    }else{
        ss <- list(N=nrow(x),
                   xsum=colSums(x),
                   S = t(x)%*%x)
        class(ss) <- "ssGaussian"
        ss
    }
}

#' Weighted sufficient statistics for a "GaussianNIW" object
#'
#' For following Gaussian-NIW model structure: \cr
#'      theta|gamma ~ NIW(gamma) \cr
#'      x|theta ~ Gaussian(theta) \cr
#' where theta = (mu,Sigma) is the Gaussian parameter, gamma = (m,k,v,S) is the Normal-Inverse-Wishart(NIW) parameter. \cr
#' The sufficient statistics of a set of samples x (each row of x is a sample) and weights w are: \cr
#'  the effective number of samples N=sum(w); \cr
#'  the sample sum xsum = colSums(x*w); \cr
#'  the uncentered scatter matrix S = t(w*x)%*%x.
#'
#' @seealso \code{\link{GaussianNIW}}, \code{\link{sufficientStatistics.GaussianNIW}}
#' @param obj A "GaussianNIW" object.
#' @param x, matrix, Gaussian samples, when x is a matrix, each row is a sample of dimension ncol(x). when x is a vector, x is length(x) samples of dimension 1.
#' @param w numeric, sample weights.
#' @param foreach logical, if foreach=TRUE, will return a list of sufficient statistics for each row of x, otherwise will return the sufficient statistics of x as a whole.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return If foreach=TRUE, will return a list of sufficient statistics for each row of x, otherwise will return the sufficient statistics of x as a whole.
#' @export
#' @examples
#' x <- rGaussian(10,mu = c(-1.5,1.5),Sigma = matrix(c(0.1,0.03,0.03,0.1),2,2))
#' obj <- GaussianNIW()                    #an GaussianNIW object
#' w <- runif(10)
#' sufficientStatistics_Weighted(obj=obj,x=x,w=w,foreach = FALSE)
#' sufficientStatistics_Weighted(obj=obj,x=x,w=w,foreach = TRUE)
#' @references Murphy, Kevin P. "Conjugate Bayesian analysis of the Gaussian distribution." def 1.22 (2007): 16.
#' @references Gelman, Andrew, et al. "Bayesian Data Analysis Chapman & Hall." CRC Texts in Statistical Science (2004).
sufficientStatistics_Weighted.GaussianNIW<- function(obj,x,w,foreach=FALSE,...){
    if(missing(x)|missing(x)) stop("'x' or 'w' not specified!")
    if(is.vector(x)){
        x <- matrix(x, ncol = 1)
    }else if(!.is(x,"matrix")){
        stop("'x' must be a vector(for univariate t) or matrix(for multivariate t)!")
    }
    if(length(w)!=nrow(x)) stop("Error in sufficientStatisticsGaussian_Weighted(): number of weights and observations don't match")
    if(foreach){
        sapply(1:nrow(x),function(i){
            ss <- list(N=w[i],
                 xsum=x[i,,drop=TRUE]*w[i],
                 S=t(w[i]*x[i,,drop=FALSE])%*%x[i,,drop=FALSE])
            class(ss) <- "ssGaussian"
            ss
        },simplify = FALSE,USE.NAMES = FALSE)
    }else{
        ss <- list(N=sum(w),
                   xsum=colSums(x*w),
                   S = t(w*x)%*%x)
        class(ss) <- "ssGaussian"
        ss
    }
}

#' Update a "GaussianNIW" object with sample sufficient statistics
#'
#' For the model structure: \cr
#'      theta|gamma ~ NIW(gamma) \cr
#'      x|theta ~ Gaussian(theta) \cr
#' where theta = (mu,Sigma) is the Gaussian parameter, gamma = (m,k,v,S) is the Normal-Inverse-Wishart(NIW) parameter. \cr
#' Update gamma by adding the information of newly observed samples x. The model structure and prior parameters are stored in a "GaussianNIW" object, the prior parameters in this object will be updated after running this function.
#'
#' @seealso \code{\link{GaussianNIW}},\code{\link{posteriorDiscard.GaussianNIW}},\code{\link{sufficientStatistics.GaussianNIW}}
#' @param obj A "GaussianNIW" object.
#' @param ss Sufficient statistics of x. In Gaussian-NIW case the sufficient statistic of sample x is a object of type "ssGaussian", it can be  generated by the function sufficientStatistics().
#' @param w Sample weights, default NULL.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return None. the gamma stored in "obj" will be updated based on "ss".
#' @export
#' @examples
#' x <- rGaussian(1000,mu = c(1,1),Sigma = matrix(c(1,0.5,0.5,3),2,2))
#' w <- runif(1000)
#' obj <- GaussianNIW(gamma=list(m=c(0,0),k=1,v=2,S=diag(2)))
#' obj
#' ss <- sufficientStatistics_Weighted(obj = obj,x=x,w=w,foreach = TRUE)
#' for(i in 1L:length(ss)) posterior(obj = obj,ss = ss[[i]])
#' obj
#' @references Murphy, Kevin P. "Conjugate Bayesian analysis of the Gaussian distribution." def 1.22 (2007): 16.
#' @references Gelman, Andrew, et al. "Bayesian Data Analysis Chapman & Hall." CRC Texts in Statistical Science (2004).
posterior.GaussianNIW <- function(obj,ss,w=NULL,...){
    
    if(missing(ss)) stop("'ss' not specified!")
    if(!.is(ss,"ssGaussian")) stop("'ss' must be of class 'ssGaussian', you need to use sufficientStatistics() to generate 'ssGaussian' objects")
    k0 <- obj$gamma$k
    m0 <- obj$gamma$m


    obj$gamma$k <- obj$gamma$k+ss$N
    obj$gamma$m <- (k0*m0+ss$xsum)/obj$gamma$k
    obj$gamma$v <- obj$gamma$v+ss$N
    obj$gamma$S <- obj$gamma$S+ss$S+(k0*m0)%*%t(m0)-(obj$gamma$k*obj$gamma$m)%*%t(obj$gamma$m)
}

#' Update a "GaussianNIW" object with sample sufficient statistics
#'
#' For the model structure: \cr
#'      theta|gamma ~ NIW(gamma) \cr
#'      x|theta ~ Gaussian(theta) \cr
#' where theta = (mu,Sigma) is the Gaussian parameter, gamma = (m,k,v,S) is the Normal-Inverse-Wishart(NIW) parameter. \cr
#' Contrary to posterior(), this function will update gamma by removing the information of observed samples x. The model structure and prior parameters are stored in a "GaussianNIW" object, the prior parameters in this object will be updated after running this function.
#'
#' @seealso \code{\link{GaussianNIW}},\code{\link{posterior.GaussianNIW}}
#' @param obj A "GaussianNIW" object.
#' @param ss Sufficient statistics of x. In Gaussian-NIW case the sufficient statistic of sample x is a object of type "ssGaussian", it can be  generated by the function sufficientStatistics().
#' @param w Sample weights,default NULL.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return None. the gamma stored in "obj" will be updated with the information in "ss".
#' @export
#' @examples
#' ## generate some random Gaussian samples
#' x <- rGaussian(1000,mu = c(1,1),Sigma = matrix(c(1,0.5,0.5,3),2,2))
#' w <- runif(1000)
#' ## add information to 'obj' then remove them one by one
#' obj <- GaussianNIW(gamma=list(m=c(0.2,3),k=1,v=2,S=diag(2)))
#' ss <- sufficientStatistics_Weighted(obj = obj,x=x,w=w,foreach = TRUE)
#' for(i in 1L:length(ss)) posterior(obj = obj,ss=ss[[i]])
#' obj
#' for(i in 1L:length(ss)) posteriorDiscard(obj = obj,ss=ss[[i]])
#' obj
#' ## add information to 'obj' then remove them as a whole
#' obj <- GaussianNIW(gamma=list(m=c(0.2,3),k=1,v=2,S=diag(2)))
#' ssAll <- sufficientStatistics_Weighted(obj = obj,x=x,w=w,foreach = FALSE)
#' posterior(obj = obj,ss = ssAll)
#' obj
#' posteriorDiscard(obj = obj,ss = ssAll)
#' obj
#' @references Murphy, Kevin P. "Conjugate Bayesian analysis of the Gaussian distribution." def 1.22 (2007): 16.
#' @references Gelman, Andrew, et al. "Bayesian Data Analysis Chapman & Hall." CRC Texts in Statistical Science (2004).
posteriorDiscard.GaussianNIW <- function(obj,ss,w=NULL,...){
    
    if(missing(ss)) stop("'ss' not specified!")
    if(!.is(ss,"ssGaussian")) stop("'ss' must be of class 'ssGaussian', you need to use sufficientStatistics() to generate 'ssGaussian' objects")
    kN <- obj$gamma$k
    mN <- obj$gamma$m
    
    obj$gamma$k <- kN-ss$N
    obj$gamma$m <- (mN*kN-ss$xsum)/obj$gamma$k
    obj$gamma$v <- obj$gamma$v-ss$N
    obj$gamma$S <- obj$gamma$S-ss$S-(obj$gamma$k*obj$gamma$m)%*%t(obj$gamma$m)+(kN*mN)%*%t(mN)
}

#' Maximum A Posteriori(MAP) estimate of a "GaussianNIW" object
#'
#' Generate the MAP estimate of "theta" in following Gaussian-NIW structure: \cr
#'      theta|gamma ~ NIW(gamma) \cr
#'      x|theta ~ Gaussian(theta) \cr
#' where theta = (mu,Sigma) is the Gaussian parameter, gamma = (m,k,v,S) is the Normal-Inverse-Wishart(NIW) parameter. \cr
#' The model structure and prior parameters are stored in a "GaussianNIW" object. \cr
#' MAP is theta_MAP = argmax_theta p(theta|gamma,x).
#'
#' @seealso \code{\link{GaussianNIW}}
#' @param obj A "GaussianNIW" object.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return A named list, the MAP estimate of "theta".
#' @export
#' @examples
#' ## update the piror with new observations then calculate the MAP estimate
#' x <- rGaussian(1000,mu = c(1,1),Sigma = matrix(c(1,0.5,0.5,3),2,2))
#' w <- runif(1000)
#' obj <- GaussianNIW(gamma=list(m=c(0,0),k=1,v=2,S=diag(2)))
#' ss <- sufficientStatistics_Weighted(obj = obj,x=x,w=w,foreach = TRUE)
#' for(i in 1L:length(ss)) posterior(obj = obj,ss=ss[[i]])
#' MAP(obj)
#' @references Murphy, Kevin P. "Conjugate Bayesian analysis of the Gaussian distribution." def 1.22 (2007): 16.
#' @references Gelman, Andrew, et al. "Bayesian Data Analysis Chapman & Hall." CRC Texts in Statistical Science (2004).
MAP.GaussianNIW <- function(obj,...){
    D <- length(obj$gamma$m)                      #dimension
    list(muMAP=obj$gamma$m,
         sigmaMAP=obj$gamma$S/(obj$gamma$v+D+2))
}

#' Mean Posterior Estimate(MPE) of a "GaussianNIW" object
#'
#' Generate the MPE of "theta" in following GaussianNIW structure: \cr
#'      theta|gamma ~ NIW(gamma) \cr
#'      x|theta ~ Gaussian(theta) \cr
#' where theta = (mu,Sigma) is the Gaussian parameter, gamma = (m,k,v,S) is the Normal-Inverse-Wishart(NIW) parameter. \cr
#' The model structure and prior parameters are stored in a "GaussianNIW" object. \cr
#' MPE is theta_MPE = mean_theta(theta|gamma,x)
#'
#' @seealso \code{\link{GaussianNIW}}
#' @param obj A "GaussianNIW" object.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return A named list, the MPE estimate of "theta".
#' @export
#' @references Murphy, Kevin P. "Conjugate Bayesian analysis of the Gaussian distribution." def 1.22 (2007): 16.
#' @references Gelman, Andrew, et al. "Bayesian Data Analysis Chapman & Hall." CRC Texts in Statistical Science (2004).
MPE.GaussianNIW <- function(obj,...){
    stop("MPE method for class 'GaussianNIW' is not implemented yet")
}

#' Marginal likelihood of a "GaussianNIW" object
#'
#' Generate the marginal likelihood of the following model structure: \cr
#'      theta|gamma ~ NIW(gamma) \cr
#'      x|theta ~ Gaussian(theta) \cr
#' where theta = (mu,Sigma) is the Gaussian parameter, gamma = (m,k,v,S) is the Normal-Inverse-Wishart(NIW) parameter. \cr
#' The model structure and prior parameters are stored in a "GaussianNIW" object. \cr
#' Marginal likelihood is the likelihood of x|gamma.
#'
#' @seealso @seealso \code{\link{GaussianNIW}}, \code{\link{marginalLikelihood_bySufficientStatistics.GaussianNIW}}
#' @param obj A "GaussianNIW" object.
#' @param x matrix, or the ones that can be converted to matrix. each row of x is an observationobservation matrix, or the ones that can be converted to matrix. each row of x is an observation
#' @param LOG Return the log density if set to "TRUE".
#' @param ... Additional arguments to be passed to other inherited types.
#' @return numeric, the marginal likelihood.
#' @export
#' @examples
#' x <- rGaussian(1000,mu = c(1,1),Sigma = matrix(c(1,0.5,0.5,3),2,2))
#' obj <- GaussianNIW(gamma=list(m=c(0,0),k=1,v=2,S=diag(2)))
#' marginalLikelihood(obj = obj,x=x)
#' ## or...
#' ss <- sufficientStatistics(obj=obj,x=x,foreach = FALSE)
#' marginalLikelihood_bySufficientStatistics(obj = obj,ss=ss)
#' @references Murphy, Kevin P. "Conjugate Bayesian analysis of the Gaussian distribution." def 1.22 (2007): 16.
#' @references Gelman, Andrew, et al. "Bayesian Data Analysis Chapman & Hall." CRC Texts in Statistical Science (2004).
marginalLikelihood.GaussianNIW <- function(obj,x,LOG=TRUE,...){
    if(missing(x)) stop("'x' not specified!")
    if(is.vector(x)){
        x <- matrix(x,ncol = 1)
    }else if(!.is(x,"matrix")){
        stop("'x' must be a vector(for univariate t) or matrix(for multivariate t)!")
    }
    ss <- sufficientStatistics.GaussianNIW(obj=obj,x=x,foreach = FALSE)
    marginalLikelihood_bySufficientStatistics.GaussianNIW(obj = obj,ss=ss,LOG = LOG)
}

#' Marginal likelihood of a "GaussianNIW" object, usnig sufficient statistics
#'
#' Generate the marginal likelihood of a set of observations of the following model structure: \cr
#'      theta|gamma ~ NIW(gamma) \cr
#'      x|theta ~ Gaussian(theta) \cr
#' where theta = (mu,Sigma) is the Gaussian parameter, gamma = (m,k,v,S) is the Normal-Inverse-Wishart(NIW) parameter. \cr
#' The model structure and prior parameters are stored in a "GaussianNIW" object. \cr
#' Marginal likelihood is the likelihood of x|gamma.
#'
#' @seealso @seealso \code{\link{GaussianNIW}}, \code{\link{marginalLikelihood.GaussianNIW}}
#' @param obj A "GaussianNIW" object.
#' @param ss Sufficient statistics of x. In Gaussian-NIW case the sufficient statistic of sample x is a object of type "ssGaussian", it can be  generated by the function sufficientStatistics().
#' @param LOG Return the log density if set to "TRUE".
#' @param ... Additional arguments to be passed to other inherited types.
#' @return numeric, the marginal likelihood.
#' @export
#' @examples
#' x <- rGaussian(1000,mu = c(1,1),Sigma = matrix(c(1,0.5,0.5,3),2,2))
#' obj <- GaussianNIW(gamma=list(m=c(0,0),k=1,v=2,S=diag(2)))
#' marginalLikelihood(obj = obj,x=x)
#' ## or...
#' ss <- sufficientStatistics(obj=obj,x=x,foreach = FALSE)
#' marginalLikelihood_bySufficientStatistics(obj = obj,ss=ss)
#' @references Murphy, Kevin P. "Conjugate Bayesian analysis of the Gaussian distribution." def 1.22 (2007): 16.
#' @references Gelman, Andrew, et al. "Bayesian Data Analysis Chapman & Hall." CRC Texts in Statistical Science (2004).
marginalLikelihood_bySufficientStatistics.GaussianNIW <- function(obj,ss,LOG=TRUE,...){
    if(missing(ss)) stop("'ss' not specified!")
    if(!.is(ss,"ssGaussian")) stop("'ss' must be of class 'ssGaussian', you need to use sufficientStatistics() to generate 'ssGaussian' objects")
    
    lmvgamma <- function(a,p){
        sapply(a,function(ai){
            p*(p-1)/4*log(pi)+ sum(lgamma(ai+(1-(1:p))/2))
        },simplify = TRUE)
    }
    det <- function(m){                 #det doesn't support scalar values, so write a wrapper around it
        if(is.matrix(m)){
            base::det(m)
        }else if(is.vector(m)){
            base::prod(m)
        }
    }
    
    
    obj2 <- GaussianNIW(objCopy = obj) #copy obj to obj2
    posterior.GaussianNIW(obj = obj2,ss=ss) #update the posteriors
    
    D <- length(obj2$gamma$m)                     #dimension
    
    logp <- -ss$N*D/2*log(pi) + D/2*log(obj$gamma$k/obj2$gamma$k) + obj$gamma$v/2*log(det(obj$gamma$S)) - obj2$gamma$v/2*log(det(obj2$gamma$S)) + lmvgamma(obj2$gamma$v/2,D) - lmvgamma(obj$gamma$v/2,D)
    if(!LOG) logp <- exp(logp)
    
    logp
}

#' Posterior predictive density function of a "GaussianNIW" object
#'
#' Generate the the density value of the posterior predictive distribution of the following structure: \cr
#'      theta|gamma ~ NIW(gamma) \cr
#'      x|theta ~ Gaussian(theta) \cr
#' where theta = (mu,Sigma) is the Gaussian parameter, gamma = (m,k,v,S) is the Normal-Inverse-Wishart(NIW) parameter. \cr
#' The model structure and prior parameters are stored in a "GaussianNIW" object. \cr
#' Posterior predictive is a distribution of x|gamma.
#'
#' @seealso @seealso \code{\link{GaussianNIW}}, \code{\link{dPosteriorPredictive.GaussianNIW}}, \code{\link{marginalLikelihood.GaussianNIW}}
#' @param obj A "GaussianNIW" object.
#' @param x numeric/integer/character vector, observed Categorical samples.
#' @param LOG Return the log density if set to "TRUE".
#' @param ... Additional arguments to be passed to other inherited types.
#' @return A numeric vector, the posterior predictive density.
#' @export
#' @examples
#' x <- rGaussian(1000,mu = c(1,1),Sigma = matrix(c(1,0.5,0.5,3),2,2))
#' obj <- GaussianNIW(gamma=list(m=c(0,0),k=1,v=2,S=diag(2)))
#' ## out1 and out2 it should have the same values:
#' out1 <- dPosteriorPredictive(obj = obj, x = x,LOG = TRUE)
#' out2 <- numeric(nrow(x))
#' for(i in 1:nrow(x))
#' out2[i] <- marginalLikelihood(obj,x=x[i,,drop=FALSE],LOG = TRUE)
#' max(abs(out1-out2))
#' @references Murphy, Kevin P. "Conjugate Bayesian analysis of the Gaussian distribution." def 1.22 (2007): 16.
#' @references Gelman, Andrew, et al. "Bayesian Data Analysis Chapman & Hall." CRC Texts in Statistical Science (2004).
dPosteriorPredictive.GaussianNIW <- function(obj,x,LOG=TRUE,...){
    if(missing(x)) stop("'x' not specified!")
    if(is.vector(x)){
        x <- matrix(x,ncol = 1)
    }else if(!.is(x,"matrix")){
        stop("'x' must be a vector(for univariate t) or matrix(for multivariate t)!")
    }



    if(obj$gamma$v-length(obj$gamma$m)+1 < 0) stop("In the parameters of NIW, 'v' must be greater than p-1, where p is the dimension of the Gaussian variable. This error can be resolved by setting a larger 'v' when initializing the GaussianNIW object.")
    d <- length(obj$gamma$m)
    dT(x=x,mu=obj$gamma$m,Sigma = (obj$gamma$k+1)/obj$gamma$k/(obj$gamma$v-d+1)*obj$gamma$S,df = obj$gamma$v-d+1,LOG = LOG)
}

#' Posterior predictive random generation of a "GaussianNIW" object
#'
#' Generate random samples from the posterior predictive distribution of the following structure: \cr
#'      theta|gamma ~ NIW(gamma) \cr
#'      x|theta ~ Gaussian(theta) \cr
#' where theta = (mu,Sigma) is the Gaussian parameter, gamma = (m,k,v,S) is the Normal-Inverse-Wishart(NIW) parameter. \cr
#' The model structure and prior parameters are stored in a "GaussianNIW" object. \cr
#' Posterior predictive is a distribution of x|gamma.
#' 
#' @seealso @seealso \code{\link{GaussianNIW}}, \code{\link{dPosteriorPredictive.GaussianNIW}}
#' @param obj A "GaussianNIW" object.
#' @param n integer, number of samples.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return A matrix of n rows, each row is a sample.
#' @export
#' @examples
#' obj <- GaussianNIW(gamma=list(m=c(0,0),k=1,v=2,S=diag(2)))
#' rPosteriorPredictive(obj=obj,20)
#' @references Murphy, Kevin P. "Conjugate Bayesian analysis of the Gaussian distribution." def 1.22 (2007): 16.
#' @references Gelman, Andrew, et al. "Bayesian Data Analysis Chapman & Hall." CRC Texts in Statistical Science (2004).
rPosteriorPredictive.GaussianNIW <- function(obj,n,...){
    if(obj$gamma$v-length(obj$gamma$m)+1 < 0) stop("In the parameters of NIW, 'v' must be greater than p-1, where p is the dimension of the Gaussian variable. This error can be resolved by setting a larger 'v' when initializing the GaussianNIW object.")
    d <- length(obj$gamma$m)
    rT(n=n,mu=obj$gamma$m,Sigma = (obj$gamma$k+1)/obj$gamma$k/(obj$gamma$v-d+1)*obj$gamma$S,df = obj$gamma$v-d+1)
}

#' Create objects of type "GaussianNIG".
#'
#' Create an object of type "GaussianNIG", which represents the Gaussian-Normal-Inverse-Gamma (Gaussian-NIG) conjugate structure: \cr
#'      beta,sigma^2|gamma ~ NIG(gamma) \cr
#'      x|beta,sigma^2,X ~ Gaussian(X * beta,sigma^2) \cr
#' where gamma = (m,V,a,b) is the Normal-Inverse-Gamma(NIG) parameter, "m" is a numeric "location" parameter; "V" is a symmetric positive definite matrix representing the "scale" parameters; "a" and "b" are the "shape" and "rate" parameter of the Inverse Gamma distribution. \cr
#' This object will be used as a place for recording and accumulating information in the related inference/sampling functions such as posterior(), posteriorDiscard(), MAP() and so on.
#'
#' @seealso \code{\link{posterior.GaussianNIG}},\code{\link{posteriorDiscard.GaussianNIG}},\code{\link{MAP.GaussianNIG}},\code{\link{marginalLikelihood.GaussianNIG}} ...
#' @param objCopy An object of type "GaussianNIG". If "objCopy" is not NULL, the function create a new "GaussianNIG" object by copying the content from objCopy, otherwise this new object will be created by using "ENV" and "gamma". Default NULL.
#' @param ENV environment, specify in which environment the object will be created
#' @param gamma list, a named list of NIG parameters, gamma=list(m,V,a,b). Where gamma$m is a numeric "location" parameter; gamma$V is a symmetric positive definite matrix representing the "scale" parameters; gamma$a and gamma$b are the "shape" and "rate" parameter of the Inverse Gamma distribution.
#' @return An object of class "GaussianNIG".
#' @export
#' @examples
#' \donttest{
#' X <- 1:20                               #generate some linear data
#' x <- rnorm(20)+ X*0.3                   #generate some linear data
#' obj <- GaussianNIG(gamma=list(m=0,V=1,a=1,b=0)) #create a GaussianNIG object
#' ss <- sufficientStatistics(obj = obj,X=X,x=x)   #the sufficient statistics of X and x
#' posterior(obj = obj,ss = ss)                    #add the infomation to the posterior
#' MAP(obj)                                        #get the MAP estimate of beta and sigma^2
#' ## print the whole content, "invV" and "mVm" in the output are temporary variables.
#' obj
#' }
#' @references Banerjee, Sudipto. "Bayesian Linear Model: Gory Details." Dowloaded from http://www. biostat. umn. edu/~ ph7440 (2008).
GaussianNIG <- function(objCopy=NULL,ENV=parent.frame(),gamma=list(m=0,V=1,a=1,b=1)){
    object <- BasicBayesian(ENV = ENV)
    if(!is.null(objCopy)){
        if(!.is(objCopy,"GaussianNIG")) stop("'objCopy' must be of class 'GaussianNIG'")
        object$gamma <- objCopy$gamma
        object$H <- objCopy$H
        object$F <- objCopy$F
    }else{
        if(!missing(gamma))
            if((!is.list(gamma)) |
               (!all(names(gamma) %in% c("m","V","a","b"))))
                stop("gamma must be of list(m,V,a,b)")
        object$gamma <- gamma
        object$gamma$invV <- pdsInverse(S=gamma$V,returnUpper = FALSE)                 #intermediate variable
        object$gamma$mVm <- t(object$gamma$m)%*%object$gamma$invV%*%object$gamma$m     #intermediate variable
        object$H <- "NIG"
        object$F <- "Gaussian"
    }
    class(object) <- c("GaussianNIG",class(object))
    return(object)
}

#' Sufficient statistics of a "GaussianNIG" object
#'
#' For following Gaussian-NIG model structure: \cr
#'      beta,sigma^2|gamma ~ NIG(gamma) \cr
#'      x|beta,sigma^2,X ~ Gaussian(X%*%beta,sigma^2) \cr
#' where gamma = (m,V,a,b) is the Normal-Inverse-Gamma(NIG) parameter, "m" is a numeric "location" parameter; "V" is a symmetric positive definite matrix representing the "scale" parameters; "a" and "b" are the "shape" and "rate" parameter of the Inverse Gamma distribution. \cr
#' The model structure and prior parameters are stored in a "GaussianNIG" object. \cr
#' The sufficient statistics of a set of samples (x,X) are: \cr
#'  the effective number of samples N=nrow(X) or length(x); \cr
#'  the covariance of X and x SXx=t(X)%*%x \cr
#'  the covariance of X SX=t(X)%*%X \cr
#'  the covariance of x Sx=t(x)%*%x
#'
#' @seealso \code{\link{GaussianNIG}}, \code{\link{sufficientStatistics_Weighted.GaussianNIG}}
#' @param obj A "GaussianNIG" object.
#' @param x numeric, must satisfy length(x) = nrow(X)
#' @param X matrix, must satisfy length(x) = nrow(X)
#' @param foreach logical, if foreach=TRUE, will return a list of sufficient statistics for each (x,X), otherwise will return the sufficient statistics as a whole.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return If foreach=TRUE, will return a list of sufficient statistics for each row of (x,X), otherwise will return the sufficient statistics of (x,X) as a whole.
#' @export
#' @examples
#' obj <- GaussianNIG(gamma=list(m=0,V=1,a=1,b=0))
#' X <- 1:20
#' x <- rnorm(20)+ X*0.3
#' sufficientStatistics(obj = obj,X=X,x=x)
#' sufficientStatistics(obj = obj,X=X,x=x,foreach = TRUE)
#' @references Banerjee, Sudipto. "Bayesian Linear Model: Gory Details." Dowloaded from http://www. biostat. umn. edu/~ ph7440 (2008).
sufficientStatistics.GaussianNIG <- function(obj,x,X,foreach=FALSE,...){
    if(missing(x)|missing(X)) stop("'x' and 'X' must be specified")
    if(!is.vector(x)) x <- as.vector(x)
    if(is.vector(X)){
        X <- matrix(X, ncol = 1)
    }else if(!.is(X,"matrix")){
        stop("'X' must be a vector or matrix!")
    }
    if(length(x)!=nrow(X)) stop("number of observations don't match")
    if(ncol(X)!=length(obj$gamma$m)) stop("number of columns in X doesn't match number of elements in beta")
    if(foreach){
        sapply(1:nrow(X),function(i){
            ss <- list(N=1,
                       SXy=crossprod(X[i,,drop=FALSE],x[i]),
                       SX=crossprod(X[i,,drop=FALSE]),
                       Sy=crossprod(x[i]))
            class(ss) <- "ssGaussianLinear"
            ss
        },simplify = FALSE,USE.NAMES = FALSE)
    }else{
        ss <- list(N=nrow(X),
                   SXy=crossprod(X,x),
                   SX=crossprod(X),
                   Sy=crossprod(x)
                   )
        class(ss) <- "ssGaussianLinear"
        ss
    }
}

#' Weighted sufficient statistics of a "GaussianNIG" object
#'
#' For following Gaussian-NIG model structure: \cr
#'      beta,sigma^2|gamma ~ NIG(gamma) \cr
#'      x|beta,sigma^2,X ~ Gaussian(X%*%beta,sigma^2) \cr
#' where gamma = (m,V,a,b) is the Normal-Inverse-Gamma(NIG) parameter, "m" is a numeric "location" parameter; "V" is a symmetric positive definite matrix representing the "scale" parameters; "a" and "b" are the "shape" and "rate" parameter of the Inverse Gamma distribution. \cr
#' The model structure and prior parameters are stored in a "GaussianNIG" object. \cr
#' This object will be used as a place for recording and accumulating information in the related inference/sampling functions such as posterior(), posteriorDiscard(), MAP() and so on. \cr
#' The sufficient statistics of a set of samples (x,X) and weights ware: \cr
#'  the effective number of samples N=sum(w); \cr
#'  the covariance of X and x SXx=t(w*X)%*%x \cr
#'  the covariance of X SX=t(w*X)%*%X \cr
#'  the covariance of x Sx=t(w*x)%*%x
#'
#' @seealso \code{\link{GaussianNIG}}, \code{\link{sufficientStatistics.GaussianNIG}}
#' @param obj A "GaussianNIG" object.
#' @param x numeric, must satisfy length(x) = nrow(X).
#' @param w numeric, sample weights.
#' @param X matrix, must satisfy length(x) = nrow(X).
#' @param foreach logical, if foreach=TRUE, will return a list of sufficient statistics for each (x,X), otherwise will return the sufficient statistics as a whole.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return If foreach=TRUE, will return a list of sufficient statistics for each row of (x,X), otherwise will return the sufficient statistics of (x,X) as a whole.
#' @export
#' @examples
#' obj <- GaussianNIG(gamma=list(m=0,V=1,a=1,b=0))
#' X <- 1:20
#' x <- rnorm(20)+ X*0.3
#' w <- runif(20)
#' sufficientStatistics_Weighted(obj = obj,X=X,x=x,w=w)
#' sufficientStatistics_Weighted(obj = obj,X=X,x=x,w=w,foreach = TRUE)
#' @references Banerjee, Sudipto. "Bayesian Linear Model: Gory Details." Dowloaded from http://www. biostat. umn. edu/~ ph7440 (2008).
sufficientStatistics_Weighted.GaussianNIG<- function(obj,x,w,X,foreach=FALSE,...){
    if(missing(x)|missing(w)|missing(X)) stop("'x', 'w' and 'X' must be specified")
    if(!is.vector(w)) w <- as.vector(w)
    if(!is.vector(x)) x <- as.vector(x)
    if(is.vector(X)){
        X <- matrix(X, ncol = 1)
    }else if(!.is(X,"matrix")){
        stop("'X' must be a vector or matrix!")
    }
    if(length(x)!=nrow(X)) stop("number of observations between x and X don't match")
    if(length(x)!=length(w)) stop("number of observations between x and w don't match")
    if(ncol(X)!=length(obj$gamma$m)) stop("number of columns in X doesn't match number of elements in beta")
    if(foreach){
        sapply(1:nrow(X),function(i){
            ss <- list(N=w[i],
                       SXy=crossprod(X[i,,drop=FALSE],x[i]*w[i]),
                       SX=crossprod(w[i]*X[i,,drop=FALSE],X[i,,drop=FALSE]),
                       Sy=crossprod(x[i]*w[i],x[i]))
            class(ss) <- "ssGaussianLinear"
            ss
        },simplify = FALSE,USE.NAMES = FALSE)
    }else{
        ss <- list(N=sum(w),
                   SXy=crossprod(X,w*x),
                   SX=crossprod(w*X,X),
                   Sy=crossprod(x*w,x)
                   )
        class(ss) <- "ssGaussianLinear"
        ss
    }

}

#' Update a "GaussianNIG" object with sample sufficient statistics
#'
#' For the model structure: \cr
#'      beta,sigma^2|gamma ~ NIG(gamma) \cr
#'      x|beta,sigma^2,X ~ Gaussian(X%*%beta,sigma^2) \cr
#' where gamma = (m,V,a,b) is the Normal-Inverse-Gamma(NIG) parameter, "m" is a numeric "location" parameter; "V" is a symmetric positive definite matrix representing the "scale" parameters; "a" and "b" are the "shape" and "rate" parameter of the Inverse Gamma distribution. \cr
#' The model structure and prior parameters are stored in a "GaussianNIG" object. \cr
#' Update gamma by adding the information of newly observed samples (x,X). The model structure and prior parameters are stored in a "GaussianNIG" object, the prior parameters in this object will be updated after running this function.
#'
#' @seealso \code{\link{GaussianNIG}},\code{\link{posteriorDiscard.GaussianNIG}},\code{\link{sufficientStatistics.GaussianNIG}}
#' @param obj A "GaussianNIG" object.
#' @param ss Sufficient statistics of (x,X). In Gaussian-NIG case the sufficient statistic of sample (x,X) is a object of type "ssGaussianLinear", it can be  generated by the function sufficientStatistics().
#' @param w Sample weights, default NULL.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return None. the gamma stored in "obj" will be updated based on "ss".
#' @export
#' @examples
#' obj <- GaussianNIG(gamma=list(m=0,V=1,a=1,b=0))
#' X <- 1:20
#' x <- rnorm(20)+ X*0.3
#' ss <- sufficientStatistics(obj = obj,X=X,x=x)
#' posterior(obj = obj,ss = ss)
#' obj
#' @references Banerjee, Sudipto. "Bayesian Linear Model: Gory Details." Dowloaded from http://www. biostat. umn. edu/~ ph7440 (2008).
posterior.GaussianNIG <- function(obj,ss,w=NULL,...){
    
    if(missing(ss)) stop("'ss' not specified!")
    if(!.is(ss,"ssGaussianLinear")) stop("'ss' must be of class 'ssGaussianLinear', you need to use sufficientStatistics() to generate 'ssGaussianLinear' objects")
    invV0 <- obj$gamma$invV
    mVm0 <- obj$gamma$mVm
    
    obj$gamma$invV <- obj$gamma$invV+ss$SX
    obj$gamma$V <- pdsInverse(S=obj$gamma$invV,returnUpper = FALSE)
    obj$gamma$m <- obj$gamma$V %*% (invV0%*%obj$gamma$m+ss$SXy)
    obj$gamma$mVm <- t(obj$gamma$m)%*%obj$gamma$invV%*%obj$gamma$m
    obj$gamma$a <- obj$gamma$a+ss$N/2
    obj$gamma$b <- drop(obj$gamma$b+0.5*(mVm0 + ss$Sy - obj$gamma$mVm)) #remove dimension, otherwise there will be an error in rPosteriorPredictive

}

#' Update a "GaussianNIG" object with sample sufficient statistics
#'
#' For the model structure: \cr
#'      beta,sigma^2|gamma ~ NIG(gamma) \cr
#'      x|beta,sigma^2,X ~ Gaussian(X%*%beta,sigma^2) \cr
#' where gamma = (m,V,a,b) is the Normal-Inverse-Gamma(NIG) parameter, "m" is a numeric "location" parameter; "V" is a symmetric positive definite matrix representing the "scale" parameters; "a" and "b" are the "shape" and "rate" parameter of the Inverse Gamma distribution. \cr
#' The model structure and prior parameters are stored in a "GaussianNIG" object. \cr
#' Contrary to posterior(), this function will update gamma by removing the information of observed samples (x,X). The model structure and prior parameters are stored in a "GaussianNIG" object, the prior parameters in this object will be updated after running this function.
#'
#' @seealso \code{\link{GaussianNIG}},\code{\link{posterior.GaussianNIG}}
#' @param obj A "GaussianNIG" object.
#' @param ss Sufficient statistics of (x,X). In Gaussian-NIG case the sufficient statistic of sample (x,X) is a object of type "ssGaussianLinear", it can be  generated by the function sufficientStatistics().
#' @param w Sample weights,default NULL.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return None. the gamma stored in "obj" will be updated with the information in "ss".
#' @export
#' @examples
#' X <- 1:20
#' x <- rnorm(20)+ X*0.3
#' ## add information, then remove them
#' obj <- GaussianNIG(gamma=list(m=0,V=1,a=1,b=1))
#' ss <- sufficientStatistics(obj = obj,X=X,x=x)
#' posterior(obj = obj,ss = ss)
#' obj
#' posteriorDiscard(obj=obj,ss=ss)
#' obj
#' ## or, add information, then remove them one by one
#' obj <- GaussianNIG(gamma=list(m=0,V=1,a=1,b=1))
#' ssEach <- sufficientStatistics(obj = obj,X=X,x=x,foreach = TRUE)
#' for(sss in ssEach) posterior(obj = obj,ss = sss)
#' obj
#' for(sss in ssEach) posteriorDiscard(obj = obj,ss = sss)
#' obj
#' @references Banerjee, Sudipto. "Bayesian Linear Model: Gory Details." Dowloaded from http://www. biostat. umn. edu/~ ph7440 (2008).
posteriorDiscard.GaussianNIG <- function(obj,ss,w=NULL,...){
    
    if(missing(ss)) stop("'ss' not specified!")
    if(!.is(ss,"ssGaussianLinear")) stop("'ss' must be of class 'ssGaussianLinear', you need to use sufficientStatistics() to generate 'ssGaussianLinear' objects")
    invVN <- obj$gamma$invV
    mVmN <- obj$gamma$mVm

    obj$gamma$invV <- obj$gamma$invV-ss$SX
    obj$gamma$V <- pdsInverse(S=obj$gamma$invV,returnUpper = FALSE)
    obj$gamma$m <- obj$gamma$V %*% (invVN%*%obj$gamma$m-ss$SXy)
    obj$gamma$mVm <- t(obj$gamma$m)%*%obj$gamma$invV%*%obj$gamma$m
    obj$gamma$a <- obj$gamma$a-ss$N/2
    obj$gamma$b <- drop(obj$gamma$b-0.5*(obj$gamma$mVm + ss$Sy - mVmN)) #remove dimension, otherwise there will be an error in rPosteriorPredictive

}

#' MAP estimate of a "GaussianNIG" object
#'
#' Generate the MAP estimate of "theta" in following Gaussian-NIG structure: \cr
#'      beta,sigma^2|gamma ~ NIG(gamma) \cr
#'      x|beta,sigma^2,X ~ Gaussian(X%*%beta,sigma^2) \cr
#' where gamma = (m,V,a,b) is the Normal-Inverse-Gamma(NIG) parameter, "m" is a numeric "location" parameter; "V" is a symmetric positive definite matrix representing the "scale" parameters; "a" and "b" are the "shape" and "rate" parameter of the Inverse Gamma distribution. \cr
#' The model structure and prior parameters are stored in a "GaussianNIG" object. \cr
#' MAP is {beta,sigma^2}_MAP = argmax_{beta,sigma^2} p(beta,sigma^2|gamma,x,X).
#'
#' @seealso \code{\link{GaussianNIG}}
#' @param obj A "GaussianNIG" object.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return A named list, the MAP estimate of beta and sigma^2.
#' @export
#' @examples
#' obj <- GaussianNIG(gamma=list(m=0,V=1,a=1,b=1))
#' X <- 1:20
#' x <- rnorm(20)+ X*0.3
#' ss <- sufficientStatistics(obj = obj,X=X,x=x)
#' posterior(obj = obj,ss = ss)
#' MAP(obj)
#' @references Banerjee, Sudipto. "Bayesian Linear Model: Gory Details." Dowloaded from http://www. biostat. umn. edu/~ ph7440 (2008).
MAP.GaussianNIG <- function(obj,...){
    D <- length(obj$gamma$m)                      #dimension
    list(betaMAP=obj$gamma$m,
         sigmaMAP=obj$gamma$b/(obj$gamma$a+1+D/2))
}

#' MPE estimate of a "GaussianNIG" object
#'
#' Generate the MPE estimate of "theta" in following Gaussian-NIG structure: \cr
#'      beta,sigma^2|gamma ~ NIG(gamma) \cr
#'      x|beta,sigma^2,X ~ Gaussian(X%*%beta,sigma^2) \cr
#' where gamma = (m,V,a,b) is the Normal-Inverse-Gamma(NIG) parameter, "m" is a numeric "location" parameter; "V" is a symmetric positive definite matrix representing the "scale" parameters; "a" and "b" are the "shape" and "rate" parameter of the Inverse Gamma distribution. \cr
#' The model structure and prior parameters are stored in a "GaussianNIG" object. \cr
#' MPE is {beta,sigma^2}_MPE = E(beta,sigma^2|gamma,x,X), E() is the expectation function.
#'
#' @seealso \code{\link{GaussianNIG}}
#' @param obj A "GaussianNIG" object.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return A named list, the MPE estimate of beta and sigma^2.
#' @export
#' @references Banerjee, Sudipto. "Bayesian Linear Model: Gory Details." Dowloaded from http://www. biostat. umn. edu/~ ph7440 (2008).
MPE.GaussianNIG <- function(obj,...){
    stop("MPE method for class 'GaussianNIG' is not implemented yet")
}

#' Marginal likelihood of a "GaussianNIG" object
#'
#' Generate the marginal likelihood of the following model structure: \cr
#'      beta,sigma^2|gamma ~ NIG(gamma) \cr
#'      x|beta,sigma^2,X ~ Gaussian(X%*%beta,sigma^2) \cr
#' where gamma = (m,V,a,b) is the Normal-Inverse-Gamma(NIG) parameter, "m" is a numeric "location" parameter; "V" is a symmetric positive definite matrix representing the "scale" parameters; "a" and "b" are the "shape" and "rate" parameter of the Inverse Gamma distribution. \cr
#' The model structure and prior parameters are stored in a "GaussianNIG" object. \cr
#' Marginal likelihood is the likelihood of x|gamma,X.
#'
#' @seealso @seealso \code{\link{GaussianNIG}}, \code{\link{marginalLikelihood_bySufficientStatistics.GaussianNIG}}
#' @param obj A "GaussianNIG" object.
#' @param x numeric, must satisfy length(x) = nrow(X).
#' @param X matrix, must satisfy length(x) = nrow(X).
#' @param LOG Return the log density if set to "TRUE".
#' @param ... Additional arguments to be passed to other inherited types.
#' @return numeric, the marginal likelihood.
#' @export
#' @examples
#' obj <- GaussianNIG(gamma=list(m=0,V=1,a=1,b=1))
#' X <- 1:20
#' x <- rnorm(20)+ X*0.3
#' marginalLikelihood(obj = obj,x = x, X = X)
#' marginalLikelihood(obj = obj,x = x, X = X,LOG = FALSE)
#' @references Banerjee, Sudipto. "Bayesian Linear Model: Gory Details." Dowloaded from http://www. biostat. umn. edu/~ ph7440 (2008).
marginalLikelihood.GaussianNIG <- function(obj,x,X,LOG=TRUE,...){
    if(missing(x)|missing(X)) stop("'x' and 'X' not specified!")
    if(!is.vector(x)) x <- as.vector(x)
    if(is.vector(X)){
        X <- matrix(X,ncol = 1)
    }else if(!.is(X,"matrix")){
        stop("'X' must be a vector or matrix!")
    }
    if(ncol(X)!=length(obj$gamma$m)) stop("number of columns in X doesn't match number of elements in beta")
    ss <- sufficientStatistics(obj=obj,x=x,X=X,foreach = FALSE)
    marginalLikelihood_bySufficientStatistics(obj = obj,ss=ss,LOG = LOG)
}

#' Marginal likelihood of a "GaussianNIG" object, usnig sufficient statistics
#'
#' Generate the marginal likelihood of a set of observations of the following model structure: \cr
#'      beta,sigma^2|gamma ~ NIG(gamma) \cr
#'      x|beta,sigma^2,X ~ Gaussian(X%*%beta,sigma^2) \cr
#' where gamma = (m,V,a,b) is the Normal-Inverse-Gamma(NIG) parameter, "m" is a numeric "location" parameter; "V" is a symmetric positive definite matrix representing the "scale" parameters; "a" and "b" are the "shape" and "rate" parameter of the Inverse Gamma distribution. \cr
#' The model structure and prior parameters are stored in a "GaussianNIG" object. \cr
#' Marginal likelihood is the likelihood of x|gamma,X.
#'
#' @seealso @seealso \code{\link{GaussianNIG}}, \code{\link{marginalLikelihood.GaussianNIG}}
#' @param obj A "GaussianNIG" object.
#' @param ss Sufficient statistics of (x,X). In Gaussian-NIG case the sufficient statistic of sample (x,X) is a object of type "ssGaussianLinear", it can be  generated by the function sufficientStatistics().
#' @param LOG Return the log density if set to "TRUE".
#' @param ... Additional arguments to be passed to other inherited types.
#' @return numeric, the marginal likelihood.
#' @export
#' @examples
#' obj <- GaussianNIG(gamma=list(m=0,V=1,a=1,b=1))
#' X <- 1:20
#' x <- rnorm(20)+ X*0.3
#' ss <- sufficientStatistics(obj=obj,x=x,X=X,foreach=FALSE)
#' marginalLikelihood_bySufficientStatistics(obj = obj,ss = ss)
#' marginalLikelihood_bySufficientStatistics(obj = obj,ss = ss,LOG = FALSE)
#' @references Banerjee, Sudipto. "Bayesian Linear Model: Gory Details." Dowloaded from http://www. biostat. umn. edu/~ ph7440 (2008).
marginalLikelihood_bySufficientStatistics.GaussianNIG <- function(obj,ss,LOG=TRUE,...){
    if(missing(ss)) stop("'ss' not specified!")
    if(!.is(ss,"ssGaussianLinear")) stop("'ss' must be of class 'ssGaussianLinear', you need to use sufficientStatistics() to generate 'ssGaussianLinear' objects")
    det <- function(m){                 #det doesn't support scalar values, so write a wrapper around it
        if(is.matrix(m)){
            base::det(m)
        }else if(is.vector(m)){
            base::prod(m)
        }
    }
    
    obj2 <- GaussianNIG(objCopy = obj) #copy obj to obj2
    posterior.GaussianNIG(obj = obj2,ss=ss) #update the posteriors

    logp <- log(det(obj2$gamma$V))/2 + lgamma(obj2$gamma$a) - obj2$gamma$a*log(obj2$gamma$b) -
        (log(det(obj$gamma$V))/2 + lgamma(obj$gamma$a) - obj$gamma$a*log(obj$gamma$b))-
        ss$N/2*log(2*pi)

    if(!LOG) logp <- exp(logp)
    
    drop(logp)                          #remove dimensions
}

#' Posterior predictive density function of a "GaussianNIG" object
#'
#' Generate the the density value of the posterior predictive distribution of the following structure: \cr
#'      beta,sigma^2|gamma ~ NIG(gamma) \cr
#'      x|beta,sigma^2,X ~ Gaussian(X%*%beta,sigma^2) \cr
#' where gamma = (m,V,a,b) is the Normal-Inverse-Gamma(NIG) parameter, "m" is a numeric "location" parameter; "V" is a symmetric positive definite matrix representing the "scale" parameters; "a" and "b" are the "shape" and "rate" parameter of the Inverse Gamma distribution. \cr
#' The model structure and prior parameters are stored in a "GaussianNIG" object. \cr
#' Posterior predictive is the distribution of x|gamma,X.
#'
#' @seealso @seealso \code{\link{GaussianNIG}}, \code{\link{dPosteriorPredictive.GaussianNIG}}, \code{\link{marginalLikelihood.GaussianNIG}}
#' @param obj A "GaussianNIG" object.
#' @param x numeric, must satisfy length(x) = nrow(X).
#' @param X matrix, must satisfy length(x) = nrow(X).
#' @param LOG Return the log density if set to "TRUE".
#' @param ... Additional arguments to be passed to other inherited types.
#' @return A numeric vector, the posterior predictive density.
#' @export
#' @examples
#' obj <- GaussianNIG(gamma=list(m=0,V=1,a=1,b=1))
#' X <- 1:20
#' x <- rnorm(20)+ X*0.3
#' ## out1 and out2 it should have the same values:
#' out1 <- dPosteriorPredictive(obj = obj, x = x,X=X,LOG = TRUE)
#' out2 <- numeric(length(x))
#' for(i in 1:length(x))
#' out2[i] <- marginalLikelihood(obj,x=x[i],X=X[i],LOG = TRUE)
#' max(abs(out1-out2))
#' @references Banerjee, Sudipto. "Bayesian Linear Model: Gory Details." Dowloaded from http://www. biostat. umn. edu/~ ph7440 (2008).
dPosteriorPredictive.GaussianNIG <- function(obj,x,X,LOG=TRUE,...){
    if(missing(x)|missing(X)) stop("'x' and 'X' not specified!")
    if(!is.vector(x)) x <- as.vector(x)
    if(is.vector(X)){
        X <- matrix(X,ncol = 1)
    }else if(!.is(X,"matrix")){
        stop("'X' must be a vector or matrix!")
    }
    if(nrow(X)!=length(x)) stop("nrow(X) and length(x) don't match!")
    if(ncol(X)!=length(obj$gamma$m)) stop("number of columns in X doesn't match number of elements in beta")
    sapply(1L:nrow(X),function(i){
        dT(x=x[i],mu=X[i,,drop=FALSE]%*%obj$gamma$m,Sigma = obj$gamma$b/obj$gamma$a*(1+X[i,,drop=FALSE]%*%obj$gamma$V%*%t(X[i,,drop=FALSE])),df = 2*obj$gamma$a,LOG = LOG)
    })
}

#' Posterior predictive random generation of a "GaussianNIG" object
#'
#' Generate random samples from the posterior predictive distribution of the following structure: \cr
#'      beta,sigma^2|gamma ~ NIG(gamma) \cr
#'      x|beta,sigma^2,X ~ Gaussian(X%*%beta,sigma^2) \cr
#' where gamma = (m,V,a,b) is the Normal-Inverse-Gamma(NIG) parameter, "m" is a numeric "location" parameter; "V" is a symmetric positive definite matrix representing the "scale" parameters; "a" and "b" are the "shape" and "rate" parameter of the Inverse Gamma distribution. \cr
#' The model structure and prior parameters are stored in a "GaussianNIG" object. \cr
#' Posterior predictive is a distribution of x|gamma,X.
#' 
#' @seealso @seealso \code{\link{GaussianNIG}}, \code{\link{dPosteriorPredictive.GaussianNIG}}
#' @param obj A "GaussianNIG" object.
#' @param n integer, number of samples.
#' @param X matrix, the location of the prediction, each row is a location.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return A matrix of n rows and nrow(X) columns, each row is a sample.
#' @export
#' @examples
#' obj <- GaussianNIG(gamma=list(m=c(1,1),V=diag(2),a=1,b=1))
#' X <- matrix(runif(20),ncol=2)
#' rPosteriorPredictive(obj=obj,n=3,X=X)
#' @references Banerjee, Sudipto. "Bayesian Linear Model: Gory Details." Dowloaded from http://www. biostat. umn. edu/~ ph7440 (2008).
rPosteriorPredictive.GaussianNIG <- function(obj,n,X,...){
    if(missing(X)) stop("'X' not specified!")
    if(!.is(X,"matrix")){
        stop("'X' must be a matrix!")
    }
    if(ncol(X)!=length(obj$gamma$m)) stop("number of columns in X doesn't match number of elements in beta")
    N <- nrow(X)
    rT(n=n,mu=X%*%obj$gamma$m,Sigma = obj$gamma$b/obj$gamma$a*(diag(N)+X%*%obj$gamma$V%*%t(X)),df = 2*obj$gamma$a)
}
