
#' @include Bayesian_Bricks.r Gamma_Inference.r






#' @title Inference in joint Gaussian distribution
#' @description
#' For the model structure
#'    \deqn{x1,x2|mu,Sigma ~ Gaussian(mu,Sigma)}
#'    \deqn{x1|x2,mu,Sigma ~ Gaussian(mu12,Sigma12)}
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

#' @title Linear Gaussian systems
#' @description
#' For the model structure:
#'    \deqn{x1 ~ Gaussian(mu1,Sigma1)}
#'    \deqn{x2 ~ Gaussian(Ax1+b,Sigma21)}
#'    \deqn{x1|x2,mu1,Sigma1,A,b,Sigma21 ~ Gaussian(mu12,Sigma12)}
#' @seealso \code{\link{inferenceJointGaussian}} 
#' @param x2 numeric, an sample of x2, length(x2)=d2.
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

#' @title Random generation for Gaussian distribution
#' @description
#' Generate random samples from a Gaussian distribution. For a random vector x, the density function of a (multivariate) Gaussian distribution is defined as:
#'    \deqn{sqrt(2 pi^p |Sigma|)^{-1} exp(-1/2 (x-mu )^T Sigma^{-1} (x-mu))}
#' where p is the dimension of x.
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
    X <- t(A)%*%X + mu
    
    t(X)
}

#' @title Density function of Gaussian distribution
#' @description
#' Get the density of a set of samples from a (multivariate) Gaussian distribution. For a random vector x, the density function is defined as:
#'    \deqn{sqrt(2 pi^p |Sigma|)^{-1} exp(-1/2 (x-mu )^T Sigma^{-1} (x-mu))}
#' where p is the dimension of x.
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
    if(ncol(x)!=length(mu)) stop("dimension of 'x' and 'mu' don't match")

    if(is.null(A) & !is.null(Sigma)) A <- chol(Sigma)
    if(is.null(A) & is.null(Sigma)) stop("Error in dGaussian(): at least one of 'Sigma' and 'A' should be non-NULL!")
    
    p <- ncol(x)
    b <- backsolve(A,diag(p))
    logp <- -(p/2)*log(2*pi) + sum(log(diag(b))) - 0.5*colSums((crossprod(b,(t(x)-mu)))^2)
    if(!LOG) logp <- exp(logp)
    logp
}

#' @title Random Generation for (multivariate) t distribution
#' @description
#' Generate random samples from a (multivariate) t distribution. For a random vector x, the density function is defined as:
#'     \deqn{Gamma((df + p)/2) / (Gamma(df/2)df^{p/2} pi ^{p/2} |Sigma|^{1/2}) [1+1/df (x-df)^T Sigma^{-1} (x-df)]^{-(df +p)/2}}
#' Where p is the dimension of x.
#' @seealso \code{\link{dT}}
#' @param n integer, number of samples.
#' @param mu numeric, mean vector.
#' @param Sigma matrix, Sigma is proportional to the covariance matrix of x, one of Sigma and A should be non-NULL.
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

#' @title Density function for (multivariate) t distribution
#' @description
#' Get the density of a set of samples from a t distribution. For a random vector x, the density function is defined as:
#'    \deqn{Gamma((df + p)/2) / (Gamma(df/2)df^{p/2} pi ^{p/2} |Sigma|^{1/2}) [1+1/df (x-df)^T Sigma^{-1} (x-df)]^{-(df +p)/2}}
#' Where p is the dimension of x.
#' @seealso \code{\link{rT}}
#' @param x matrix, when x is a numeric vector, it will be converted to a matrix with 1 column!
#' @param mu numeric, mean vector.
#' @param Sigma matrix, Sigma is proportional to the covariance matrix of x, one of Sigma and A should be non-NULL.
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

#' @title Density function for Normal-Inverse-Wishart (NIW) distribution.
#' @description
#' Get the density of a NIW sample. For a random vector mu, and a random matrix Sigma, the density function is defined as:
#'    \deqn{sqrt(2 pi^p |Sigma/k|)^{-1} exp(-1/2 (mu-m )^T (Sigma/k)^{-1} (mu-m)) (2^{(v p)/2} Gamma_p(v/2) |S|^{-v/2})^{-1} |Sigma|^{(-v-p-1)/2} exp(-1/2 tr(Sigma^{-1} S))}
#' Where p is the dimension of mu and Sigma.
#' @seealso \code{\link{rNIW}}
#' @param mu numeric, the Gaussian sample.
#' @param Sigma matrix, a symmetric positive definite matrix, the Inverse-Wishart sample.
#' @param m numeric, mean of mu.
#' @param k numeric, precision of mu.
#' @param v numeric, degree of freedom of Sigma.
#' @param S numeric, a symmetric positive definite scale matrix of Sigma, S is proportional to E(Sigma).
#' @param LOG logical, return log density of LOG=TRUE, default TRUE.
#' @return A numeric vector, the probability density of (mu,Sigma).
#' @export
#' @examples
#' S <- crossprod(matrix(rnorm(15),5,3))
#' Sigma <- crossprod(matrix(rnorm(15),5,3))
#' mu <- runif(3)
#' m <- runif(3)
#' dNIW(mu=mu,Sigma=Sigma,m=m,k=2,v=4,S=S,LOG = TRUE)
#' @references O'Hagan, Anthony, and Jonathan J. Forster. Kendall's advanced theory of statistics, volume 2B: Bayesian inference. Vol. 2. Arnold, 2004.
#' @references MARolA, K. V., JT KBNT, and J. M. Bibly. Multivariate analysis. AcadeInic Press, Londres, 1979.
dNIW <- function(mu,Sigma,m,k,v,S,LOG=TRUE){
    if(is.vector(mu)){
        mu <- matrix(mu,nrow=1)
    }
    if(ncol(mu)!=nrow(Sigma)) stop("dimension of 'mu' and 'Sigma' don't match")
    logp <- dGaussian(x=mu,mu = m,Sigma = Sigma/k,LOG = TRUE)+
        dInvWishart(x=Sigma,df=v,scale=S,LOG = TRUE)
    if(!LOG) logp <- exp(logp)
    logp
}

#' @title Random number generation for Normal-Inverse-Wishart (NIW) distribution.
#' @description
#' Generate a NIW sample.For a random vector mu, and a random matrix Sigma, the density function is defined as:
#'    \deqn{sqrt(2 pi^p |Sigma/k|)^{-1} exp(-1/2 (mu-m )^T (Sigma/k)^{-1} (mu-m)) (2^{(v p)/2} Gamma_p(v/2) |S|^{-v/2})^{-1} |Sigma|^{(-v-p-1)/2} exp(-1/2 tr(Sigma^{-1} S))}
#' Where p is the dimension of mu and Sigma.
#' @seealso \code{\link{dNIW}}
#' @param m numeric, mean of mu.
#' @param k numeric, precision of mu.
#' @param v numeric, degree of freedom of Sigma.
#' @param S numeric, a symmetric positive definite scale matrix of Sigma, S is proportional to E(Sigma).
#' @return A list of two list(mu,Sigma), where 'mu' is a numeric vector, the Gaussian sample; 'Sigma' is a symmetric positive definite matrix, the Inverse-Wishart sample.
#' @export
#' @examples
#' rNIW(m=runif(3),k=0.001,v=5,S=crossprod(matrix(rnorm(15),5,3)))
#' @references O'Hagan, Anthony, and Jonathan J. Forster. Kendall's advanced theory of statistics, volume 2B: Bayesian inference. Vol. 2. Arnold, 2004.
#' @references MARolA, K. V., JT KBNT, and J. M. Bibly. Multivariate analysis. AcadeInic Press, Londres, 1979.
rNIW <- function(m,k,v,S){
    Sigma <- rInvWishart(df = v,scale = S)
    mu <- drop(rGaussian(n=1,mu = m,Sigma = Sigma/k))
    list(mu=mu,Sigma=Sigma)
}

#' @title Create objects of type "GaussianGaussian".
#' @description
#' Create an object of type "GaussianGaussian", which represents the Gaussian and Gaussian conjugate structure:
#'      \deqn{x ~ Gaussian(mu,Sigma)}
#'      \deqn{mu ~ Gaussian(m,S)}
#' Where Sigma is known. Gaussian() is the Gaussian distribution. See \code{?dGaussian} for the definition of Gaussian distribution.\cr
#' The created object will be used as a place for recording and accumulating information in the related inference/sampling functions such as posterior(), posteriorDiscard(), MAP(), marginalLikelihood(), dPosteriorPredictive(), rPosteriorPredictive() and so on.
#' @seealso \code{\link{posterior.GaussianGaussian}},\code{\link{posteriorDiscard.GaussianGaussian}},\code{\link{MAP.GaussianGaussian}},\code{\link{MPE.GaussianGaussian}},\code{\link{marginalLikelihood.GaussianGaussian}},\code{\link{rPosteriorPredictive.GaussianGaussian}},\code{\link{dPosteriorPredictive.GaussianGaussian}}.
#' @param objCopy an object of type "GaussianGaussian". If "objCopy" is not NULL, the function create a new "GaussianGaussian" object by copying the content from objCopy, otherwise this new object will be created by using "ENV" and "gamma". Default NULL.
#' @param ENV environment, specify where the object will be created.
#' @param gamma list, a named list of parameters, gamma=list(Sigma,m,S). Where gamma$Sigma is the known covariance matrix of x, gamma$m and gamma$S are the prior mean and covariance matrix of mu.
#' @return An object of class "GaussianGaussian".
#' @export
#' @examples
#' obj <- GaussianGaussian(gamma=list(Sigma=matrix(c(2,1,1,2),2,2),m=c(0.2,0.5),S=diag(2)))
#' obj #print the content
#' @references Gelman, Andrew, et al. Bayesian data analysis. CRC press, 2013.
GaussianGaussian <- function(objCopy=NULL,ENV=parent.frame(),gamma=list(Sigma=1,m=0,S=1)){
    object <- BasicBayesian(ENV = ENV)
    if(!is.null(objCopy)){
        if(!.is(objCopy,"GaussianGaussian")) stop("'objCopy' must be of class 'GaussianGaussian'")
        object$gamma <- objCopy$gamma
        object$H <- objCopy$H
        object$F <- objCopy$F
    }else{
        if(!missing(gamma))
            if((!is.list(gamma)) |
               (!all(names(gamma) %in% c("Sigma","m","S"))))
                stop("gamma must be of list(Sigma,m,S)")
        object$gamma <- gamma
        object$H <- "Gaussian"
        object$F <- "Gaussian"
        object$gamma$priorPrecision <- pdsInverse(gamma$S)
        object$gamma$observedPrecision <- pdsInverse(gamma$Sigma)
    }
    class(object) <- c("GaussianGaussian",class(object))
    return(object)
}

#' @title Sufficient statistics of a "GaussianGaussian" object
#' @description
#' For following model structure:
#'      \deqn{x ~ Gaussian(mu,Sigma)}
#'      \deqn{mu ~ Gaussian(m,S)}
#' Where Sigma is known. Gaussian() is the Gaussian distribution. See \code{?dGaussian} for the definition of Gaussian distribution.\cr
#' The sufficient statistics of a set of samples x (each row of x is a sample) are: \cr
#' \itemize{
#'  \item the effective number of samples N=nrow(x)
#'  \item the sample sum xsum = colSums(x)
#' }
#' @seealso \code{\link{GaussianGaussian}}, \code{\link{sufficientStatistics_Weighted.GaussianGaussian}} 
#' @param obj A "GaussianGaussian" object.
#' @param x matrix, Gaussian samples, when x is a matrix, each row is a sample of dimension ncol(x). when x is a vector, x is length(x) samples of dimension 1.
#' @param foreach logical, specifying whether to return the sufficient statistics for each observation. Default FALSE.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return If foreach=TRUE, will return a list of sufficient statistics for each row of x, otherwise will return the sufficient statistics of x as a whole.
#' @export
#' @examples
#' obj <- GaussianGaussian(gamma=list(Sigma=matrix(c(2,1,1,2),2,2),m=c(0.2,0.5),S=diag(2)))
#' x <- rGaussian(100,c(0,0),Sigma = matrix(c(2,1,1,2),2,2))
#' sufficientStatistics(obj = obj,x=x)
#' sufficientStatistics(obj = obj,x=x,foreach=TRUE)
#' @references Gelman, Andrew, et al. Bayesian data analysis. CRC press, 2013.
sufficientStatistics.GaussianGaussian <- function(obj,x,foreach=FALSE,...){
    if(missing(x)) stop("'x' must be specified")
    if(is.vector(x)){
        x <- matrix(x, ncol = 1)
    }else if(!.is(x,"matrix")){
        stop("'x' must be a vector(for univariate t) or matrix(for multivariate t)!")
    }
    if(foreach){
        sapply(1:nrow(x),function(i){
            ss <- list(N=1,
                       xsum=x[i,,drop=TRUE])
            class(ss) <- "ssGaussianMean"
            ss
        },simplify = FALSE,USE.NAMES = FALSE)
    }else{
        ss <- list(N=nrow(x),
                   xsum=colSums(x))
        class(ss) <- "ssGaussianMean"
        ss
    }
}

#' @title Weighted sufficient statistics of a "GaussianGaussian" object
#' @description
#' For following model structure:
#'      \deqn{x ~ Gaussian(mu,Sigma)}
#'      \deqn{mu ~ Gaussian(m,S)}
#' Where Sigma is known. Gaussian() is the Gaussian distribution. See \code{?dGaussian} for the definition of Gaussian distribution.\cr
#' The sufficient statistics of a set of samples x (each row of x is a sample) and weights w are: \cr
#' \itemize{
#'  \item the effective number of samples N=sum(w)
#'  \item the sample sum xsum = colSums(x*w)
#' }
#' @seealso \code{\link{GaussianGaussian}}, \code{\link{sufficientStatistics.GaussianGaussian}} 
#' @param obj A "GaussianGaussian" object.
#' @param x matrix, Gaussian samples, when x is a matrix, each row is a sample of dimension ncol(x). when x is a vector, x is length(x) samples of dimension 1.
#' @param w numeric, sample weights.
#' @param foreach logical, specifying whether to return the sufficient statistics for each observation. Default FALSE.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return If foreach=TRUE, will return a list of sufficient statistics for each row of x, otherwise will return the sufficient statistics of x as a whole.
#' @export
#' @examples
#' obj <- GaussianGaussian(gamma=list(Sigma=matrix(c(2,1,1,2),2,2),m=c(0.2,0.5),S=diag(2)))
#' x <- rGaussian(100,c(0,0),Sigma = matrix(c(2,1,1,2),2,2))
#' w <- runif(100)
#' sufficientStatistics_Weighted(obj=obj,x=x,w=w,foreach = FALSE)
#' sufficientStatistics_Weighted(obj=obj,x=x,w=w,foreach = TRUE)
#' @references Gelman, Andrew, et al. Bayesian data analysis. CRC press, 2013.
sufficientStatistics_Weighted.GaussianGaussian<- function(obj,x,w,foreach=FALSE,...){
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
                       xsum=x[i,,drop=TRUE]*w[i])
            class(ss) <- "ssGaussianMean"
            ss
        },simplify = FALSE,USE.NAMES = FALSE)
    }else{
        ss <- list(N=sum(w),
                   xsum=colSums(x*w))
        class(ss) <- "ssGaussianMean"
        ss
    }
}

#' @title Update a "GaussianGaussian" object with sample sufficient statistics
#' @description
#' For the model structure:
#'      \deqn{x ~ Gaussian(mu,Sigma)}
#'      \deqn{mu ~ Gaussian(m,S)}
#' Where Sigma is known. Gaussian() is the Gaussian distribution. See \code{?dGaussian} for the definition of Gaussian distribution.\cr
#' Update (m,S) by adding the information of newly observed samples x. The model structure and prior parameters are stored in a "GaussianGaussian" object, the prior parameters in this object will be updated after running this function.
#' @seealso \code{\link{GaussianGaussian}},\code{\link{posteriorDiscard.GaussianGaussian}},\code{\link{sufficientStatistics.GaussianGaussian}}
#' @param obj A "GaussianGaussian" object.
#' @param ss Sufficient statistics of x. In Gaussian-Gaussian case the sufficient statistic of sample x is a object of type "ssGaussianMean", it can be  generated by the function sufficientStatistics().
#' @param w Sample weights, default NULL.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return None. the gamma stored in "obj" will be updated based on "ss".
#' @export
#' @examples
#' obj <- GaussianGaussian(gamma=list(Sigma=matrix(c(2,1,1,2),2,2),m=c(0.2,0.5),S=diag(2)))
#' obj
#' x <- rGaussian(100,c(0,0),Sigma = matrix(c(2,1,1,2),2,2))
#' ss <- sufficientStatistics(obj=obj,x=x,foreach = FALSE)
#' posterior(obj = obj,ss = ss)
#' obj
#' @references Gelman, Andrew, et al. Bayesian data analysis. CRC press, 2013.
posterior.GaussianGaussian <- function(obj,ss,w=NULL,...){
    
    if(missing(ss)) stop("'ss' not specified!")
    if(!.is(ss,"ssGaussianMean")) stop("'ss' must be of class 'ssGaussianMean', you need to use sufficientStatistics() to generate 'ssGaussianMean' objects")

    tmp <- obj$gamma$priorPrecision + obj$gamma$observedPrecision*ss$N
    obj$gamma$S <- pdsInverse(S = tmp,returnUpper = FALSE)
    obj$gamma$m <- obj$gamma$S %*% (obj$gamma$priorPrecision%*%obj$gamma$m+obj$gamma$observedPrecision%*%ss$xsum)
    obj$gamma$priorPrecision <- tmp
}

#' @title Update a "GaussianGaussian" object with sample sufficient statistics
#' @description
#' For the model structure:
#'      \deqn{x ~ Gaussian(mu,Sigma)}
#'      \deqn{mu ~ Gaussian(m,S)}
#' Where Sigma is known. Gaussian() is the Gaussian distribution. See \code{?dGaussian} for the definition of Gaussian distribution.\cr
#' Contrary to posterior(), this function will update (m,S) by removing the information of observed samples x. The model structure and prior parameters are stored in a "GaussianGaussian" object, the prior parameters in this object will be updated after running this function.
#' @seealso \code{\link{GaussianGaussian}},\code{\link{posterior.GaussianGaussian}},\code{\link{sufficientStatistics.GaussianGaussian}}
#' @param obj A "GaussianGaussian" object.
#' @param ss Sufficient statistics of x. In Gaussian-Gaussian case the sufficient statistic of sample x is a object of type "ssGaussianMean", it can be generated by the function sufficientStatistics().
#' @param w Sample weights, default NULL.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return None. the gamma stored in "obj" will be updated based on "ss".
#' @export
#' @examples
#' obj <- GaussianGaussian(gamma=list(Sigma=matrix(c(2,1,1,2),2,2),m=c(0.2,0.5),S=diag(2)))
#' obj
#' x <- rGaussian(100,c(0,0),Sigma = matrix(c(2,1,1,2),2,2))
#' ss <- sufficientStatistics(obj=obj,x=x,foreach = FALSE)
#' ## update prior into posterior
#' posterior(obj = obj,ss = ss)
#' obj
#' ## remove the information, back to prior
#' posteriorDiscard(obj = obj,ss = ss)
#' obj
#' @references Gelman, Andrew, et al. Bayesian data analysis. CRC press, 2013.
posteriorDiscard.GaussianGaussian <- function(obj,ss,w=NULL,...){
    
    if(missing(ss)) stop("'ss' not specified!")
    if(!.is(ss,"ssGaussianMean")) stop("'ss' must be of class 'ssGaussianMean', you need to use sufficientStatistics() to generate 'ssGaussianMean' objects")

    tmp <- obj$gamma$priorPrecision
    obj$gamma$priorPrecision <- obj$gamma$priorPrecision - obj$gamma$observedPrecision*ss$N
    obj$gamma$S <- pdsInverse(S = obj$gamma$priorPrecision,returnUpper = FALSE)
    obj$gamma$m <- obj$gamma$S %*% (tmp %*% obj$gamma$m - obj$gamma$observedPrecision %*% ss$xsum)
}

#' @title Maximum A Posteriori (MAP) estimate of a "GaussianGaussian" object
#' @description
#' Generate the MAP estimate of mu in following model structure:
#'      \deqn{x ~ Gaussian(mu,Sigma)}
#'      \deqn{mu ~ Gaussian(m,S)}
#' Where Sigma is known. Gaussian() is the Gaussian distribution. See \code{?dGaussian} for the definition of Gaussian distribution.\cr
#' The model structure and prior parameters are stored in a "GaussianGaussian" object. \cr
#' The MAP estimates are:
#' \itemize{
#'   \item (mu_MAP) = argmax p(mu|m,S,x,Sigma)
#' }
#' @seealso \code{\link{GaussianGaussian}}
#' @param obj A "GaussianGaussian" object.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return numeric vector, the MAP estimate of "mu".
#' @export
#' @examples
#' obj <- GaussianGaussian(gamma=list(Sigma=matrix(c(2,1,1,2),2,2),m=c(0.2,0.5),S=diag(2)))
#' x <- rGaussian(100,c(0,0),Sigma = matrix(c(2,1,1,2),2,2))
#' ss <- sufficientStatistics(obj=obj,x=x,foreach = FALSE)
#' ## update prior into posterior
#' posterior(obj = obj,ss = ss)
#' ## get the MAP estimate of mu
#' MAP(obj)
#' @references Gelman, Andrew, et al. Bayesian data analysis. CRC press, 2013.
MAP.GaussianGaussian <- function(obj,...){
    obj$gamma$m
}

#' @title Mean Posterior Estimate (MPE) of a "GaussianGaussian" object
#' @description
#' Generate the MAP estimate of mu in following model structure:
#'      \deqn{x ~ Gaussian(mu,Sigma)}
#'      \deqn{mu ~ Gaussian(m,S)}
#' Where Sigma is known. Gaussian() is the Gaussian distribution. See \code{?dGaussian} for the definition of Gaussian distribution.\cr
#' The model structure and prior parameters are stored in a "GaussianGaussian" object. \cr
#' The MPE estimates is:
#' \itemize{
#'   \item mu_MPE = E(mu|m,S,x,Sigma)
#' }
#' @seealso \code{\link{GaussianGaussian}}
#' @param obj A "GaussianGaussian" object.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return numeric vector, the MPE estimate of "mu".
#' @export
#' @examples
#' obj <- GaussianGaussian(gamma=list(Sigma=matrix(c(2,1,1,2),2,2),m=c(0.2,0.5),S=diag(2)))
#' x <- rGaussian(100,c(0,0),Sigma = matrix(c(2,1,1,2),2,2))
#' ss <- sufficientStatistics(obj=obj,x=x,foreach = FALSE)
#' ## update prior into posterior
#' posterior(obj = obj,ss = ss)
#' ## get the MPE estimate of mu
#' MPE(obj)
#' @references Gelman, Andrew, et al. Bayesian data analysis. CRC press, 2013.
MPE.GaussianGaussian <- function(obj,...){
    obj$gamma$m
}

#' @title Marginal likelihood of a "GaussianGaussian" object
#' @description
#' Generate the marginal likelihood of the following model structure:
#'      \deqn{x ~ Gaussian(mu,Sigma)}
#'      \deqn{mu ~ Gaussian(m,S)}
#' Where Sigma is known. Gaussian() is the Gaussian distribution. See \code{?dGaussian} for the definition of Gaussian distribution.\cr
#' The model structure and prior parameters are stored in a "GaussianGaussian" object. \cr
#' Marginal likelihood = p(x|m,S,Sigma)
#' @seealso \code{\link{GaussianGaussian}}, \code{\link{marginalLikelihood_bySufficientStatistics.GaussianGaussian}}
#' @param obj A "GaussianGaussian" object.
#' @param x matrix, or the ones that can be converted to matrix. each row of x is an observation matrix, or the ones that can be converted to matrix.
#' @param LOG Return the log density if set to "TRUE".
#' @param ... Additional arguments to be passed to other inherited types.
#' @return numeric, the marginal likelihood.
#' @export
#' @examples
#' obj <- GaussianGaussian(gamma=list(Sigma=matrix(c(2,1,1,2),2,2),m=c(0.2,0.5),S=diag(2)))
#' x <- rGaussian(100,c(0,0),Sigma = matrix(c(2,1,1,2),2,2))
#' marginalLikelihood(obj = obj,x=x,LOG = TRUE)
#' marginalLikelihood(obj = obj,x=x,LOG = FALSE)
#' @references Gelman, Andrew, et al. Bayesian data analysis. CRC press, 2013.
marginalLikelihood.GaussianGaussian <- function(obj,x,LOG=TRUE,...){

    if(missing(x)) stop("'x' not specified!")
    if(is.vector(x)){
        x <- matrix(x,ncol = 1)
    }else if(!.is(x,"matrix")){
        stop("'x' must be a vector(for univariate t) or matrix(for multivariate t)!")
    }
    ss <- sufficientStatistics.GaussianGaussian(obj=obj,x=x,foreach = TRUE)
    obj2 <- GaussianGaussian(objCopy = obj) #copy obj to obj2
    logp <- numeric(nrow(x))
    for(i in 1:nrow(x)){
        logp[i] <- dPosteriorPredictive.GaussianGaussian(obj = obj2,x=x[i,,drop=FALSE],LOG = TRUE)
        posterior.GaussianGaussian(obj = obj2,ss=ss[[i]])
    }

    if(!LOG) prod(exp(logp))
    else sum(logp)
}

#' @title Marginal likelihood of a "GaussianGaussian" object, using sufficient statistics
#' @description
#' Generate the marginal likelihood of the following model structure:
#'      \deqn{x ~ Gaussian(mu,Sigma)}
#'      \deqn{mu ~ Gaussian(m,S)}
#' Where Sigma is known. Gaussian() is the Gaussian distribution. See \code{?dGaussian} for the definition of Gaussian distribution.\cr
#' The model structure and prior parameters are stored in a "GaussianGaussian" object. \cr
#' Marginal likelihood = p(x|m,S,Sigma)
#' @seealso \code{\link{GaussianGaussian}}, \code{\link{marginalLikelihood.GaussianGaussian}}
#' @param obj A "GaussianGaussian" object.
#' @param ss Sufficient statistics of x. In Gaussian-Gaussian case the sufficient statistic of sample x is a object of type "ssGaussianMean", it can be  generated by the function sufficientStatistics().
#' @param LOG Return the log density if set to "TRUE".
#' @param ... Additional arguments to be passed to other inherited types.
#' @return numeric, the marginal likelihood.
#' @export
marginalLikelihood_bySufficientStatistics.GaussianGaussian <- function(obj,ss,LOG=TRUE,...){
    stop("marginalLikelihood_bySufficientStatistics() for class 'GaussianGaussian' is not implemented yet, please use marginalLikelihood() instead.")
    if(missing(ss)) stop("'ss' not specified!")
    if(!.is(ss,"ssGaussianMean")) stop("'ss' must be of class 'ssGaussianMean', you need to use sufficientStatistics() to generate 'ssGaussianMean' objects")
}

#' @title Posterior predictive density function of a "GaussianGaussian" object
#' @description
#' Generate the the density value of the posterior predictive distribution of the following structure:
#'      \deqn{x ~ Gaussian(mu,Sigma)}
#'      \deqn{mu ~ Gaussian(m,S)}
#' Where Sigma is known. Gaussian() is the Gaussian distribution. See \code{?dGaussian} for the definition of Gaussian distribution.\cr
#' The model structure and prior parameters are stored in a "GaussianGaussian" object. \cr
#' Posterior predictive density is p(x|m,S,Sigma).
#' @seealso \code{\link{GaussianGaussian}}, \code{\link{dPosteriorPredictive.GaussianGaussian}}, \code{\link{marginalLikelihood.GaussianGaussian}}
#' @param obj A "GaussianGaussian" object.
#' @param x  matrix, or the ones that can be converted to matrix. each row of x is an observation matrix, or the ones that can be converted to matrix.
#' @param LOG Return the log density if set to "TRUE".
#' @param ... Additional arguments to be passed to other inherited types.
#' @return A numeric vector of the same length as nrow(x), the posterior predictive density.
#' @export
#' @examples
#' obj <- GaussianGaussian(gamma=list(Sigma=matrix(c(2,1,1,2),2,2),m=c(0.2,0.5),S=diag(2)))
#' x <- rGaussian(100,c(0,0),Sigma = matrix(c(2,1,1,2),2,2))
#' dPosteriorPredictive(obj = obj,x=x,LOG = TRUE)
#' dPosteriorPredictive(obj = obj,x=x,LOG = FALSE)
#' @references Gelman, Andrew, et al. Bayesian data analysis. CRC press, 2013.
dPosteriorPredictive.GaussianGaussian <- function(obj,x,LOG=TRUE,...){
    if(missing(x)) stop("'x' not specified!")
    if(is.vector(x)){
        x <- matrix(x,ncol = 1)
    }else if(!.is(x,"matrix")){
        stop("'x' must be a vector(for univariate t) or matrix(for multivariate t)!")
    }
    dGaussian(x=x,mu=obj$gamma$m,Sigma = obj$gamma$S+obj$gamma$Sigma,LOG = LOG)
}

#' @title Posterior predictive random generation of a "GaussianGaussian" object
#' @description
#' Generate random samples from the posterior predictive distribution of the following structure:
#'      \deqn{x ~ Gaussian(mu,Sigma)}
#'      \deqn{mu ~ Gaussian(m,S)}
#' Where Sigma is known. Gaussian() is the Gaussian distribution. See \code{?dGaussian} for the definition of Gaussian distribution.\cr
#' The model structure and prior parameters are stored in a "GaussianGaussian" object. \cr
#' Posterior predictive is a distribution of x|m,S,Sigma.
#' @seealso \code{\link{GaussianGaussian}}, \code{\link{dPosteriorPredictive.GaussianGaussian}}
#' @param obj A "GaussianGaussian" object.
#' @param n integer, number of samples.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return A matrix of n rows, each row is a sample.
#' @export
#' @examples
#' obj <- GaussianGaussian(gamma=list(Sigma=matrix(c(2,1,1,2),2,2),m=c(0.2,0.5),S=diag(2)))
#' rPosteriorPredictive(obj=obj,20)
#' @references Gelman, Andrew, et al. Bayesian data analysis. CRC press, 2013.
rPosteriorPredictive.GaussianGaussian <- function(obj,n=1,...){
    rGaussian(n=n,mu=obj$gamma$m,Sigma = obj$gamma$S+obj$gamma$Sigma)
}

#' @title Create objects of type "GaussianInvWishart".
#' @description
#' Create an object of type "GaussianInvWishart", which represents the Gaussian and Inverse-Wishart conjugate structure:
#'     \deqn{x ~ Gaussian(mu,Sigma)}
#'     \deqn{Sigma ~ InvWishart(v,S)}
#' mu is known. Gaussian() is the Gaussian distribution. See \code{?dGaussian} and \code{?dInvWishart} for the definition of the distributions.\cr
#' The created object will be used as a place for recording and accumulating information in the related inference/sampling functions such as posterior(), posteriorDiscard(), MAP(), marginalLikelihood(), dPosteriorPredictive(), rPosteriorPredictive() and so on.
#' @seealso \code{\link{posterior.GaussianInvWishart}},\code{\link{posteriorDiscard.GaussianInvWishart}},\code{\link{MAP.GaussianInvWishart}},\code{\link{MPE.GaussianInvWishart}},\code{\link{marginalLikelihood.GaussianInvWishart}},\code{\link{rPosteriorPredictive.GaussianInvWishart}},\code{\link{dPosteriorPredictive.GaussianInvWishart}} ...
#' @param objCopy an object of type "GaussianInvWishart". If "objCopy" is not NULL, the function create a new "GaussianInvWishart" object by copying the content from objCopy, otherwise this new object will be created by using "ENV" and "gamma". Default NULL.
#' @param ENV environment, specify where the object will be created.
#' @param gamma list, a named list of parameters, gamma=list(mu,v,S). Where gamma$mu is the known mean vector of x, gamma$v and gamma$S are the prior degree of freedom and scale of Sigma.
#' @return An object of class "GaussianInvWishart".
#' @export
#' @examples
#' obj <- GaussianInvWishart(gamma=list(mu=c(-1.5,1.5),v=3,S=diag(2)))
#' obj #print the content
#' @references Gelman, Andrew, et al. Bayesian data analysis. CRC press, 2013.
GaussianInvWishart <- function(objCopy=NULL,ENV=parent.frame(),gamma=list(mu=0,v=3,S=1)){
    object <- BasicBayesian(ENV = ENV)
    if(!is.null(objCopy)){
        if(!.is(objCopy,"GaussianInvWishart")) stop("'objCopy' must be of class 'GaussianInvWishart'")
        object$gamma <- objCopy$gamma
        object$H <- objCopy$H
        object$F <- objCopy$F
    }else{
        if(!missing(gamma))
            if((!is.list(gamma)) |
               (!all(names(gamma) %in% c("mu","v","S"))))
                stop("gamma must be of list(mu,v,S)")
        gamma$mu <- as.vector(gamma$mu)
        object$gamma <- gamma
        object$H <- "InvWishart"
        object$F <- "Gaussian"
    }
    class(object) <- c("GaussianInvWishart",class(object))
    return(object)
}

#' @title Sufficient statistics of a "GaussianInvWishart" object
#' @description
#' For following model structure:
#'     \deqn{x ~ Gaussian(mu,Sigma)}
#'     \deqn{Sigma ~ InvWishart(v,S)}
#' mu is known. Gaussian() is the Gaussian distribution. See \code{?dGaussian} and \code{?dInvWishart} for the definition of the distributions.\cr
#' The sufficient statistics of a set of samples x (each row of x is a sample) are: \cr
#' \itemize{
#'  \item the effective number of samples N=nrow(x)
#'  \item the centered sample scatter matrix S = (t(x)-mu)^T %*% (t(x)-mu)
#' }
#' @seealso \code{\link{GaussianInvWishart}}, \code{\link{sufficientStatistics_Weighted.GaussianInvWishart}} 
#' @param obj A "GaussianInvWishart" object.
#' @param x matrix, Gaussian samples, when x is a matrix, each row is a sample of dimension ncol(x). when x is a vector, x is length(x) samples of dimension 1.
#' @param foreach logical, specifying whether to return the sufficient statistics for each observation. Default FALSE.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return If foreach=TRUE, will return a list of sufficient statistics for each row of x, otherwise will return the sufficient statistics of x as a whole.
#' @export
#' @examples
#' obj <- GaussianInvWishart(gamma=list(mu=c(-1.5,1.5),v=3,S=diag(2)))
#' x <- rGaussian(10,mu = c(-1.5,1.5),Sigma = matrix(c(0.1,0.03,0.03,0.1),2,2))
#' sufficientStatistics(obj=obj,x=x,foreach = FALSE)
#' sufficientStatistics(obj=obj,x=x,foreach = TRUE)
#' @references Gelman, Andrew, et al. Bayesian data analysis. CRC press, 2013.
#' @references MARolA, K. V., JT KBNT, and J. M. Bibly. Multivariate analysis. AcadeInic Press, Londres, 1979.
sufficientStatistics.GaussianInvWishart <- function(obj,x,foreach=FALSE,...){
    if(missing(x)) stop("'x' must be specified")
    if(is.vector(x)){
        x <- matrix(x, ncol = 1)
    }else if(!.is(x,"matrix")){
        stop("'x' must be a vector(for univariate t) or matrix(for multivariate t)!")
    }
    if(foreach){
        sapply(1:nrow(x),function(i){
            ss <- list(N=1,
                       S=crossprod(t(x[i,,drop=TRUE]-obj$gamma$mu)))
            class(ss) <- "ssGaussianVar"
            ss
        },simplify = FALSE,USE.NAMES = FALSE)
    }else{
        ss <- list(N=nrow(x),
                   S = crossprod(t(t(x)-obj$gamma$mu)))
        class(ss) <- "ssGaussianVar"
        ss
    }
}

#' @title Weighted sufficient statistics of a "GaussianInvWishart" object
#' @description
#' For following model structure:
#'     \deqn{x ~ Gaussian(mu,Sigma)}
#'     \deqn{Sigma ~ InvWishart(v,S)}
#' mu is known. Gaussian() is the Gaussian distribution. See \code{?dGaussian} and \code{?dInvWishart} for the definition of the distributions.\cr
#' The sufficient statistics of a set of samples x (each row of x is a sample)  and weights w are: \cr
#' \itemize{
#'  \item the effective number of samples N=sum(w)
#'  \item the centered sample scatter matrix S = (w*(t(x)-mu))^T %*% (t(x)-mu)
#' }
#' @seealso \code{\link{GaussianInvWishart}}, \code{\link{sufficientStatistics.GaussianInvWishart}} 
#' @param obj A "GaussianInvWishart" object.
#' @param x matrix, Gaussian samples, when x is a matrix, each row is a sample of dimension ncol(x). when x is a vector, x is length(x) samples of dimension 1.
#' @param w numeric, sample weights.
#' @param foreach logical, specifying whether to return the sufficient statistics for each observation. Default FALSE.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return If foreach=TRUE, will return a list of sufficient statistics for each row of x, otherwise will return the sufficient statistics of x as a whole.
#' @export
#' @examples
#' obj <- GaussianInvWishart(gamma=list(mu=c(-1.5,1.5),v=4,S=diag(2)))
#' x <- rGaussian(10,mu = c(-1.5,1.5),Sigma = matrix(c(0.1,0.03,0.03,0.1),2,2))
#' w <- runif(10)
#' sufficientStatistics_Weighted(obj=obj,x=x,w=w,foreach = FALSE)
#' sufficientStatistics_Weighted(obj=obj,x=x,w=w,foreach = TRUE)
#' @references Gelman, Andrew, et al. Bayesian data analysis. CRC press, 2013.
#' @references MARolA, K. V., JT KBNT, and J. M. Bibly. Multivariate analysis. AcadeInic Press, Londres, 1979.
sufficientStatistics_Weighted.GaussianInvWishart<- function(obj,x,w,foreach=FALSE,...){
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
                       S=w[i]*crossprod(t(x[i,,drop=TRUE]-obj$gamma$mu)))
            class(ss) <- "ssGaussianVar"
            ss
        },simplify = FALSE,USE.NAMES = FALSE)
    }else{
        x <- t(t(x)-obj$gamma$mu)
        ss <- list(N=sum(w),
                   S = t(w*x)%*%x)
        class(ss) <- "ssGaussianVar"
        ss
    }
}

#' @title Update a "GaussianInvWishart" object with sample sufficient statistics
#' @description
#' For the model structure:
#'     \deqn{x ~ Gaussian(mu,Sigma)}
#'     \deqn{Sigma ~ InvWishart(v,S)}
#' mu is known. Gaussian() is the Gaussian distribution. See \code{?dGaussian} and \code{?dInvWishart} for the definition of the distributions.\cr
#' Update (v,S) by adding the information of newly observed samples x. The model structure and prior parameters are stored in a "GaussianInvWishart" object, the prior parameters in this object will be updated after running this function.
#' @seealso \code{\link{GaussianInvWishart}},\code{\link{posteriorDiscard.GaussianInvWishart}},\code{\link{sufficientStatistics.GaussianInvWishart}}
#' @param obj A "GaussianInvWishart" object.
#' @param ss Sufficient statistics of x. In Gaussian and Inverse-Wishart case the sufficient statistic of sample x is a object of type "ssGaussianVar", it can be  generated by the function sufficientStatistics().
#' @param w Sample weights, default NULL.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return None. the gamma stored in "obj" will be updated based on "ss".
#' @export
#' @examples
#' obj <- GaussianInvWishart(gamma=list(mu=c(-1.5,1.5),v=3,S=diag(2)))
#' x <- rGaussian(10,mu = c(-1.5,1.5),Sigma = matrix(c(0.1,0.03,0.03,0.1),2,2))
#' ss <- sufficientStatistics(obj=obj,x=x,foreach = FALSE)
#' obj
#' posterior(obj,ss = ss)
#' obj
#' @references Gelman, Andrew, et al. Bayesian data analysis. CRC press, 2013.
#' @references MARolA, K. V., JT KBNT, and J. M. Bibly. Multivariate analysis. AcadeInic Press, Londres, 1979.
posterior.GaussianInvWishart <- function(obj,ss,w=NULL,...){
    if(missing(ss)) stop("'ss' not specified!")
    if(!.is(ss,"ssGaussianVar")) stop("'ss' must be of class 'ssGaussianVar', you need to use sufficientStatistics() to generate 'ssGaussianVar' objects")
    obj$gamma$v <- obj$gamma$v+ss$N
    obj$gamma$S <- obj$gamma$S+ss$S
}

#' @title Update a "GaussianInvWishart" object with sample sufficient statistics
#' @description
#' For the model structure:
#'     \deqn{x ~ Gaussian(mu,Sigma)}
#'     \deqn{Sigma ~ InvWishart(v,S)}
#' mu is known. Gaussian() is the Gaussian distribution. See \code{?dGaussian} and \code{?dInvWishart} for the definition of the distributions.\cr
#' Contrary to posterior(), this function will update (v,S) by removing the information of observed samples x. The model structure and prior parameters are stored in a "GaussianInvWishart" object, the prior parameters in this object will be updated after running this function.
#' @seealso \code{\link{GaussianInvWishart}},\code{\link{posterior.GaussianInvWishart}},\code{\link{sufficientStatistics.GaussianInvWishart}}
#' @param obj A "GaussianInvWishart" object.
#' @param ss Sufficient statistics of x. In Gaussian and Inverse-Wishart case the sufficient statistic of sample x is a object of type "ssGaussianVar", it can be generated by the function sufficientStatistics().
#' @param w Sample weights, default NULL.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return None. the gamma stored in "obj" will be updated based on "ss".
#' @export
#' @examples
#' obj <- GaussianInvWishart(gamma=list(mu=c(-1.5,1.5),v=3,S=diag(2)))
#' x <- rGaussian(100,mu = c(-1.5,1.5),Sigma = matrix(c(0.1,0.03,0.03,0.1),2,2))
#' ss <- sufficientStatistics(obj=obj,x=x,foreach = FALSE)
#' obj
#' posterior(obj=obj,ss = ss)
#' obj
#' posteriorDiscard(obj=obj,ss=ss)
#' obj
#' @references Gelman, Andrew, et al. Bayesian data analysis. CRC press, 2013.
#' @references MARolA, K. V., JT KBNT, and J. M. Bibly. Multivariate analysis. AcadeInic Press, Londres, 1979.
posteriorDiscard.GaussianInvWishart <- function(obj,ss,w=NULL,...){
    if(missing(ss)) stop("'ss' not specified!")
    if(!.is(ss,"ssGaussianVar")) stop("'ss' must be of class 'ssGaussianVar', you need to use sufficientStatistics() to generate 'ssGaussianVar' objects")
    obj$gamma$v <- obj$gamma$v-ss$N
    obj$gamma$S <- obj$gamma$S-ss$S
}

#' @title Maximum A Posteriori (MAP) estimate of a "GaussianInvWishart" object
#' @description
#' Generate the MAP estimate of Sigma in following model structure:
#'     \deqn{x ~ Gaussian(mu,Sigma)}
#'     \deqn{Sigma ~ InvWishart(v,S)}
#' mu is known. Gaussian() is the Gaussian distribution. See \code{?dGaussian} and \code{?dInvWishart} for the definition of the distributions.\cr
#' The model structure and prior parameters are stored in a "GaussianInvWishart" object. \cr
#' The MAP estimates are:
#' \itemize{
#'   \item (Sigma_MAP) = argmax p(Sigma|v,S,x,mu)
#' }
#' @seealso \code{\link{GaussianInvWishart}}
#' @param obj A "GaussianInvWishart" object.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return matrix, the MAP estimate of "Sigma".
#' @export
#' @examples
#' obj <- GaussianInvWishart(gamma=list(mu=c(-1.5,1.5),v=3,S=diag(2)))
#' x <- rGaussian(100,mu = c(-1.5,1.5),Sigma = matrix(c(0.1,0.03,0.03,0.1),2,2))
#' ss <- sufficientStatistics(obj=obj,x=x,foreach = FALSE)
#' posterior(obj=obj,ss = ss)
#' MAP(obj)
#' @references Gelman, Andrew, et al. Bayesian data analysis. CRC press, 2013.
#' @references MARolA, K. V., JT KBNT, and J. M. Bibly. Multivariate analysis. AcadeInic Press, Londres, 1979.
MAP.GaussianInvWishart <- function(obj,...){
    D <- length(obj$gamma$mu)                      #dimension
    obj$gamma$S/(obj$gamma$v+D+1)
}

#' @title Mean Posterior Estimate (MPE) of a "GaussianInvWishart" object
#' @description
#' Generate the MPE estimate of Sigma in following model structure:
#'     \deqn{x ~ Gaussian(mu,Sigma)}
#'     \deqn{Sigma ~ InvWishart(v,S)}
#' mu is known. Gaussian() is the Gaussian distribution. See \code{?dGaussian} and \code{?dInvWishart} for the definition of the distributions.\cr
#' The model structure and prior parameters are stored in a "GaussianInvWishart" object. \cr
#' The MPE estimates are:
#' \itemize{
#'   \item (Sigma_MPE) = E(Sigma|v,S,x,mu)
#' }
#' @seealso \code{\link{GaussianInvWishart}}
#' @param obj A "GaussianInvWishart" object.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return matrix, the MPE estimate of "Sigma".
#' @export
#' @examples
#' obj <- GaussianInvWishart(gamma=list(mu=c(-1.5,1.5),v=3,S=diag(2)))
#' x <- rGaussian(100,mu = c(-1.5,1.5),Sigma = matrix(c(0.1,0.03,0.03,0.1),2,2))
#' ss <- sufficientStatistics(obj=obj,x=x,foreach = FALSE)
#' posterior(obj=obj,ss = ss)
#' MPE(obj)
#' @references Gelman, Andrew, et al. Bayesian data analysis. CRC press, 2013.
#' @references MARolA, K. V., JT KBNT, and J. M. Bibly. Multivariate analysis. AcadeInic Press, Londres, 1979.
MPE.GaussianInvWishart <- function(obj,...){
    D <- length(obj$gamma$mu)                      #dimension
    if(obj$gamma$v-D-1<=0) stop("obj$gamma$v must be greater than p+1, 'p' is the dimension of the Gaussian observations.")
    obj$gamma$S/(obj$gamma$v-D-1)
}

#' @title Marginal likelihood of a "GaussianInvWishart" object
#' @description
#' Generate the marginal likelihood of the following model structure:
#'     \deqn{x ~ Gaussian(mu,Sigma)}
#'     \deqn{Sigma ~ InvWishart(v,S)}
#' mu is known. Gaussian() is the Gaussian distribution. See \code{?dGaussian} and \code{?dInvWishart} for the definition of the distributions.\cr
#' The model structure and prior parameters are stored in a "GaussianInvWishart" object. \cr
#' Marginal likelihood = p(x|v,S,mu)
#' @seealso \code{\link{GaussianInvWishart}}, \code{\link{marginalLikelihood_bySufficientStatistics.GaussianInvWishart}}
#' @param obj A "GaussianInvWishart" object.
#' @param x matrix, or the ones that can be converted to matrix. each row of x is an observation matrix, or the ones that can be converted to matrix.
#' @param LOG Return the log density if set to "TRUE".
#' @param ... Additional arguments to be passed to other inherited types.
#' @return numeric, the marginal likelihood.
#' @export
#' @examples
#' obj <- GaussianInvWishart(gamma=list(mu=c(-1.5,1.5),v=3,S=diag(2)))
#' x <- rGaussian(100,mu = c(-1.5,1.5),Sigma = matrix(c(0.1,0.03,0.03,0.1),2,2))
#' xNew <- rGaussian(100,mu = c(-1.5,1.5),Sigma = matrix(c(0.1,0.03,0.03,0.1),2,2))
#' ss <- sufficientStatistics(obj=obj,x=x,foreach = FALSE)
#' ## update piror with x
#' posterior(obj=obj,ss = ss)
#' ## use the posterior to calculate the likelihood of xNew
#' marginalLikelihood(obj = obj,x = xNew,LOG = TRUE)
#' @references Gelman, Andrew, et al. Bayesian data analysis. CRC press, 2013.
#' @references MARolA, K. V., JT KBNT, and J. M. Bibly. Multivariate analysis. AcadeInic Press, Londres, 1979.
marginalLikelihood.GaussianInvWishart <- function(obj,x,LOG=TRUE,...){
    if(missing(x)) stop("'x' not specified!")
    if(is.vector(x)){
        x <- matrix(x,ncol = 1)
    }else if(!.is(x,"matrix")){
        stop("'x' must be a vector(for univariate t) or matrix(for multivariate t)!")
    }
    if(ncol(x)!=length(obj$gamma$mu)) stop("dimensions of 'x' and obj$gamma$mu don't match!")
    ss <- sufficientStatistics.GaussianInvWishart(obj=obj,x=x,foreach = TRUE)
    obj2 <- GaussianInvWishart(objCopy = obj) #copy obj to obj2
    logp <- numeric(nrow(x))
    for(i in 1:nrow(x)){
        logp[i] <- dPosteriorPredictive.GaussianInvWishart(obj = obj2,x=x[i,,drop=FALSE],LOG = TRUE)
        posterior.GaussianInvWishart(obj = obj2,ss=ss[[i]])
    }

    if(!LOG) prod(exp(logp))
    else sum(logp)

}

#' @title Marginal likelihood of a "GaussianInvWishart" object, using sufficient statistics
#' @description
#' Generate the marginal likelihood of the following model structure:
#'     \deqn{x ~ Gaussian(mu,Sigma)}
#'     \deqn{Sigma ~ InvWishart(v,S)}
#' mu is known. Gaussian() is the Gaussian distribution. See \code{?dGaussian} and \code{?dInvWishart} for the definition of the distributions.\cr
#' The model structure and prior parameters are stored in a "GaussianInvWishart" object. \cr
#' Marginal likelihood = p(x|v,S,mu)
#' @seealso \code{\link{GaussianInvWishart}}, \code{\link{marginalLikelihood.GaussianInvWishart}}
#' @param obj A "GaussianInvWishart" object.
#' @param ss Sufficient statistics of x. In Gaussian and Inverse-Wishart case the sufficient statistic of sample x is a object of type "ssGaussianVar", it can be  generated by the function sufficientStatistics().
#' @param LOG Return the log density if set to "TRUE".
#' @param ... Additional arguments to be passed to other inherited types.
#' @return numeric, the marginal likelihood.
#' @export
marginalLikelihood_bySufficientStatistics.GaussianInvWishart <- function(obj,ss,LOG=TRUE,...){
    stop("marginalLikelihood_bySufficientStatistics() for class 'GaussianInvWishart' is not implemented yet, please use marginalLikelihood() instead.")
    if(missing(ss)) stop("'ss' not specified!")
    if(!.is(ss,"ssGaussianVar")) stop("'ss' must be of class 'ssGaussianVar', you need to use sufficientStatistics() to generate 'ssGaussianMean' objects")
}

#' @title Posterior predictive density function of a "GaussianInvWishart" object
#' @description
#' Generate the the density value of the posterior predictive distribution of the following structure:
#'     \deqn{x ~ Gaussian(mu,Sigma)}
#'     \deqn{Sigma ~ InvWishart(v,S)}
#' mu is known. Gaussian() is the Gaussian distribution. See \code{?dGaussian} and \code{?dInvWishart} for the definition of the distributions.\cr
#' The model structure and prior parameters are stored in a "GaussianInvWishart" object. \cr
#' Posterior predictive density is p(x|v,S,mu).
#' @seealso \code{\link{GaussianInvWishart}}, \code{\link{dPosteriorPredictive.GaussianInvWishart}}, \code{\link{marginalLikelihood.GaussianInvWishart}}
#' @param obj A "GaussianInvWishart" object.
#' @param x  matrix, or the ones that can be converted to matrix. each row of x is an observation matrix, or the ones that can be converted to matrix.
#' @param LOG Return the log density if set to "TRUE".
#' @param ... Additional arguments to be passed to other inherited types.
#' @return A numeric vector of the same length as nrow(x), the posterior predictive density.
#' @export
#' @examples
#' obj <- GaussianInvWishart(gamma=list(mu=c(-1.5,1.5),v=3,S=diag(2)))
#' x <- rGaussian(100,mu = c(-1.5,1.5),Sigma = matrix(c(0.1,0.03,0.03,0.1),2,2))
#' xNew <- rGaussian(100,mu = c(-1.5,1.5),Sigma = matrix(c(0.1,0.03,0.03,0.1),2,2))
#' ss <- sufficientStatistics(obj=obj,x=x,foreach = FALSE)
#' ## update piror with x
#' posterior(obj=obj,ss = ss)
#' ## use the posterior to calculate the probability of observing each xNew
#' dPosteriorPredictive(obj = obj,x = xNew,LOG = TRUE)
#' @references Gelman, Andrew, et al. Bayesian data analysis. CRC press, 2013.
#' @references MARolA, K. V., JT KBNT, and J. M. Bibly. Multivariate analysis. AcadeInic Press, Londres, 1979.
dPosteriorPredictive.GaussianInvWishart <- function(obj,x,LOG=TRUE,...){
    if(missing(x)) stop("'x' not specified!")
    if(is.vector(x)){
        x <- matrix(x,ncol = 1)
    }else if(!.is(x,"matrix")){
        stop("'x' must be a vector(for univariate t) or matrix(for multivariate t)!")
    }
    d <- length(obj$gamma$mu)           #dimension
    if(ncol(x)!=d) stop("dimensions of 'x' and obj$gamma$mu don't match!")
    if(obj$gamma$v-d+1 < 0) stop("In the parameters of NIW, 'v' must be greater than p-1, where p is the dimension of the Gaussian variable. This error can be resolved by setting a larger 'v' when initializing the GaussianInvWishart object.")
    dT(x=x,mu=obj$gamma$mu,Sigma = obj$gamma$S/(obj$gamma$v-d+1),df = obj$gamma$v-d+1,LOG = LOG)
}

#' @title Posterior predictive random generation of a "GaussianInvWishart" object
#' @description
#' Generate random samples from the posterior predictive distribution of the following structure:
#'     \deqn{x ~ Gaussian(mu,Sigma)}
#'     \deqn{Sigma ~ InvWishart(v,S)}
#' mu is known. Gaussian() is the Gaussian distribution. See \code{?dGaussian} and \code{?dInvWishart} for the definition of the distributions.\cr
#' The model structure and prior parameters are stored in a "GaussianInvWishart" object. \cr
#' Posterior predictive is a distribution of x|v,S,mu.
#' @seealso \code{\link{GaussianInvWishart}}, \code{\link{dPosteriorPredictive.GaussianInvWishart}}
#' @param obj A "GaussianInvWishart" object.
#' @param n integer, number of samples.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return A matrix of n rows, each row is a sample.
#' @export
#' @examples
#' obj <- GaussianInvWishart(gamma=list(mu=c(-1.5,1.5),v=3,S=diag(2)))
#' x <- rGaussian(100,mu = c(-1.5,1.5),Sigma = matrix(c(0.1,0.03,0.03,0.1),2,2))
#' ss <- sufficientStatistics(obj=obj,x=x,foreach = FALSE)
#' ## use x to update the prior informatoin
#' posterior(obj=obj,ss = ss)
#' ## use the posterior to generate new samples
#' rPosteriorPredictive(obj = obj,n=20)
#' @references Gelman, Andrew, et al. Bayesian data analysis. CRC press, 2013.
#' @references MARolA, K. V., JT KBNT, and J. M. Bibly. Multivariate analysis. AcadeInic Press, Londres, 1979.
rPosteriorPredictive.GaussianInvWishart <- function(obj,n,...){
    d <- length(obj$gamma$mu)            #dimension
    if(obj$gamma$v-d+1 < 0) stop("In the parameters of NIW, 'v' must be greater than p-1, where p is the dimension of the Gaussian variable. This error can be resolved by setting a larger 'v' when initializing the GaussianInvWishart object.")
    rT(n=n,mu=obj$gamma$m,Sigma = obj$gamma$S/(obj$gamma$v-d+1),df = obj$gamma$v-d+1)
}

#' @title Create objects of type "GaussianNIW".
#' @description
#' Create an object of type "GaussianNIW", which represents the Gaussian-Normal-Inverse-Wishart (Gaussian-NIW) conjugate structure:
#'      \deqn{mu,Sigma|m,k,v,S ~ NIW(m,k,v,S)}
#'      \deqn{x|mu,Sigma ~ Gaussian(mu,Sigma)}
#' Where NIW() is the Normal-Inverse-Wishart distribution, Gaussian() is the Gaussian distribution. See \code{?dNIW} and \code{dGaussian} for the definitions of these distribution.\cr
#' This object will be used as a place for recording and accumulating information in the related inference/sampling functions such as posterior(), posteriorDiscard(), MAP(), marginalLikelihood(), dPosteriorPredictive(), rPosteriorPredictive() and so on.
#' @seealso \code{\link{posterior.GaussianNIW}},\code{\link{posteriorDiscard.GaussianNIW}},\code{\link{MAP.GaussianNIW}},\code{\link{MPE.GaussianNIW}},\code{\link{marginalLikelihood.GaussianNIW}},\code{\link{rPosteriorPredictive.GaussianNIW}},\code{\link{dPosteriorPredictive.GaussianNIW}} ...
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

#' @title Sufficient statistics of a "GaussianNIW" object
#' @description
#' For following Gaussian-NIW model structure:
#'      \deqn{mu,Sigma|m,k,v,S ~ NIW(m,k,v,S)}
#'      \deqn{x|mu,Sigma ~ Gaussian(mu,Sigma)}
#' Where NIW() is the Normal-Inverse-Wishart distribution, Gaussian() is the Gaussian distribution. See \code{?dNIW} and \code{dGaussian} for the definitions of these distribution.\cr
#' The sufficient statistics of a set of samples x (each row of x is a sample) are: \cr
#' \itemize{
#'  \item the effective number of samples N=nrow(x)
#'  \item the sample sum xsum = colSums(x)
#'  \item the uncentered scatter matrix S = t(x)%*%x
#' }
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

#' @title Weighted sufficient statistics for a "GaussianNIW" object
#' @description
#' For following Gaussian-NIW model structure:
#'      \deqn{mu,Sigma|m,k,v,S ~ NIW(m,k,v,S)}
#'      \deqn{x|mu,Sigma ~ Gaussian(mu,Sigma)}
#' Where NIW() is the Normal-Inverse-Wishart distribution, Gaussian() is the Gaussian distribution. See \code{?dNIW} and \code{dGaussian} for the definitions of these distribution.\cr
#' The sufficient statistics of a set of samples x (each row of x is a sample) and weights w are:
#' \itemize{
#'  \item the effective number of samples N=sum(w)
#'  \item the sample sum xsum = colSums(x*w)
#'  \item the uncentered scatter matrix S = t(w*x)%*%x
#' }
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

#' @title Update a "GaussianNIW" object with sample sufficient statistics
#' @description
#' For the model structure:
#'      \deqn{mu,Sigma|m,k,v,S ~ NIW(m,k,v,S)}
#'      \deqn{x|mu,Sigma ~ Gaussian(mu,Sigma)}
#' Where NIW() is the Normal-Inverse-Wishart distribution, Gaussian() is the Gaussian distribution. See \code{?dNIW} and \code{dGaussian} for the definitions of these distribution.\cr
#' Update (m,k,v,S) by adding the information of newly observed samples x. The model structure and prior parameters are stored in a "GaussianNIW" object, the prior parameters in this object will be updated after running this function.
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

#' @title Update a "GaussianNIW" object with sample sufficient statistics
#' @description
#' For the model structure:
#'      \deqn{mu,Sigma|m,k,v,S ~ NIW(m,k,v,S)}
#'      \deqn{x|mu,Sigma ~ Gaussian(mu,Sigma)}
#' Where NIW() is the Normal-Inverse-Wishart distribution, Gaussian() is the Gaussian distribution. See \code{?dNIW} and \code{dGaussian} for the definitions of these distribution.\cr
#' Contrary to posterior(), this function will update (m,k,v,S) by removing the information of observed samples x. The model structure and prior parameters are stored in a "GaussianNIW" object, the prior parameters in this object will be updated after running this function.
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

#' @title Maximum A Posteriori (MAP) estimate of a "GaussianNIW" object
#' @description
#' Generate the MAP estimate of (mu,Sigma) in following Gaussian-NIW structure:
#'      \deqn{mu,Sigma|m,k,v,S ~ NIW(m,k,v,S)}
#'      \deqn{x|mu,Sigma ~ Gaussian(mu,Sigma)}
#' Where NIW() is the Normal-Inverse-Wishart distribution, Gaussian() is the Gaussian distribution. See \code{?dNIW} and \code{dGaussian} for the definitions of these distribution.\cr
#' The model structure and prior parameters are stored in a "GaussianNIW" object. \cr
#' The MAP estimates are:
#' \itemize{
#'   \item (mu_MAP,Sigma_MAP) = argmax p(mu,Sigma|m,k,v,S,x)
#' }
#'
#' @seealso \code{\link{GaussianNIW}}
#' @param obj A "GaussianNIW" object.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return A named list, the MAP estimate of mu and Sigma.
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

#' @title Mean Posterior Estimate (MPE) of a "GaussianNIW" object
#' @description
#' Generate the MPE of (mu,Sigma) in following GaussianNIW structure:
#'      \deqn{mu,Sigma|m,k,v,S ~ NIW(m,k,v,S)}
#'      \deqn{x|mu,Sigma ~ Gaussian(mu,Sigma)}
#' Where NIW() is the Normal-Inverse-Wishart distribution, Gaussian() is the Gaussian distribution. See \code{?dNIW} and \code{dGaussian} for the definitions of these distribution.\cr
#' The model structure and prior parameters are stored in a "GaussianNIW" object. \cr
#' The MPE estimates are:
#' \itemize{
#'   \item (mu_MPE,Sigma_MPE) = E(mu,Sigma|m,k,v,S,x)
#' }
#' @seealso \code{\link{GaussianNIW}}
#' @param obj A "GaussianNIW" object.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return A named list, the MPE estimate of mu and Sigma.
#' @export
#' @references Murphy, Kevin P. "Conjugate Bayesian analysis of the Gaussian distribution." def 1.22 (2007): 16.
#' @references Gelman, Andrew, et al. "Bayesian Data Analysis Chapman & Hall." CRC Texts in Statistical Science (2004).
MPE.GaussianNIW <- function(obj,...){
    if(obj$gamma$v-D-1<0) stop("'v' must be greater than D+1 for the expactation of the covariance matrix to exist. D is the dimension of the Gaussian random variable.")
    D <- length(obj$gamma$m)                      #dimension
    list(muMAP=obj$gamma$m,
         sigmaMAP=obj$gamma$S/(obj$gamma$v-D-1))
}

#' @title Marginal likelihood of a "GaussianNIW" object
#' @description
#' Generate the marginal likelihood of the following model structure:
#'      \deqn{mu,Sigma|m,k,v,S ~ NIW(m,k,v,S)}
#'      \deqn{x|mu,Sigma ~ Gaussian(mu,Sigma)}
#' Where NIW() is the Normal-Inverse-Wishart distribution, Gaussian() is the Gaussian distribution. See \code{?dNIW} and \code{dGaussian} for the definitions of these distribution.\cr
#' The model structure and prior parameters are stored in a "GaussianNIW" object. \cr
#' Marginal likelihood = p(x|m,k,v,S)
#'
#' @seealso \code{\link{GaussianNIW}}, \code{\link{marginalLikelihood_bySufficientStatistics.GaussianNIW}}
#' @param obj A "GaussianNIW" object.
#' @param x matrix, or the ones that can be converted to matrix. each row of x is an observation matrix, or the ones that can be converted to matrix. each row of x is an observation
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

#' @title Marginal likelihood of a "GaussianNIW" object, using sufficient statistics
#' @description
#' Generate the marginal likelihood of a set of observations of the following model structure:
#'      \deqn{mu,Sigma|m,k,v,S ~ NIW(m,k,v,S)}
#'      \deqn{x|mu,Sigma ~ Gaussian(mu,Sigma)}
#' Where NIW() is the Normal-Inverse-Wishart distribution, Gaussian() is the Gaussian distribution. See \code{?dNIW} and \code{dGaussian} for the definitions of these distribution.\cr
#' The model structure and prior parameters are stored in a "GaussianNIW" object. \cr
#' Marginal likelihood = p(x|m,k,v,S)
#' @seealso \code{\link{GaussianNIW}}, \code{\link{marginalLikelihood.GaussianNIW}}
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
    
    
    
    obj2 <- GaussianNIW(objCopy = obj) #copy obj to obj2
    posterior.GaussianNIW(obj = obj2,ss=ss) #update the posteriors
    
    D <- length(obj2$gamma$m)                     #dimension
    
    logp <- -ss$N*D/2*log(pi) + D/2*log(obj$gamma$k/obj2$gamma$k) + obj$gamma$v/2*pdsDeterminant(S=obj$gamma$S,LOG=TRUE) - obj2$gamma$v/2*pdsDeterminant(S=obj2$gamma$S,LOG=TRUE) + .lmvgamma(obj2$gamma$v/2,D) - .lmvgamma(obj$gamma$v/2,D)
    if(!LOG) logp <- exp(logp)
    
    logp
}

#' @title Posterior predictive density function of a "GaussianNIW" object
#' @description
#' Generate the the density value of the posterior predictive distribution of the following structure:
#'      \deqn{mu,Sigma|m,k,v,S ~ NIW(m,k,v,S)}
#'      \deqn{x|mu,Sigma ~ Gaussian(mu,Sigma)}
#' Where NIW() is the Normal-Inverse-Wishart distribution, Gaussian() is the Gaussian distribution. See \code{?dNIW} and \code{dGaussian} for the definitions of these distribution.\cr
#' The model structure and prior parameters are stored in a "GaussianNIW" object. \cr
#' Posterior predictive density is p(x|m,k,v,S).
#' @seealso \code{\link{GaussianNIW}}, \code{\link{dPosteriorPredictive.GaussianNIW}}, \code{\link{marginalLikelihood.GaussianNIW}}
#' @param obj A "GaussianNIW" object.
#' @param x  matrix, or the ones that can be converted to matrix. each row of x is an observation matrix, or the ones that can be converted to matrix.
#' @param LOG Return the log density if set to "TRUE".
#' @param ... Additional arguments to be passed to other inherited types.
#' @return A numeric vector of the same length as nrow(x), the posterior predictive density.
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

#' @title Posterior predictive random generation of a "GaussianNIW" object
#' @description
#' Generate random samples from the posterior predictive distribution of the following structure:
#'      \deqn{mu,Sigma|m,k,v,S ~ NIW(m,k,v,S)}
#'      \deqn{x|mu,Sigma ~ Gaussian(mu,Sigma)}
#' Where NIW() is the Normal-Inverse-Wishart distribution, Gaussian() is the Gaussian distribution. See \code{?dNIW} and \code{dGaussian} for the definitions of these distribution.\cr
#' The model structure and prior parameters are stored in a "GaussianNIW" object. \cr
#' Posterior predictive is a distribution of x|m,k,v,S.
#' @seealso \code{\link{GaussianNIW}}, \code{\link{dPosteriorPredictive.GaussianNIW}}
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

#' @title Create objects of type "GaussianNIG".
#' @description
#' Create an object of type "GaussianNIG", which represents the Gaussian and Normal-Inverse-Gamma (Gaussian-NIG) conjugate structure:
#'    \deqn{x ~ Gaussian(X^T beta,sigma^2)}
#'    \deqn{sigma^2 ~ InvGamma(a,b)}
#'    \deqn{beta ~ Gaussian(m,sigma^2 V)}
#' Where InvGamma() is the Inverse-Gamma distribution, Gaussian() is the Gaussian distribution. See \code{?dInvGamma} and \code{dGaussian} for the definitions of these distribution.\cr
#' This object will be used as a place for recording and accumulating information in the related inference/sampling functions such as posterior(), posteriorDiscard(), MAP(), marginalLikelihood(), dPosteriorPredictive(), rPosteriorPredictive() and so on.
#' @seealso \code{\link{posterior.GaussianNIG}},\code{\link{posteriorDiscard.GaussianNIG}},\code{\link{MAP.GaussianNIG}},\code{\link{MPE.GaussianNIG}},\code{\link{marginalLikelihood.GaussianNIG}},\code{\link{dPosteriorPredictive.GaussianNIG}}, \code{\link{rPosteriorPredictive.GaussianNIG}} ...
#' @param objCopy An object of type "GaussianNIG". If "objCopy" is not NULL, the function create a new "GaussianNIG" object by copying the content from objCopy, otherwise this new object will be created by using "ENV" and "gamma". Default NULL.
#' @param ENV environment, specify in which environment the object will be created
#' @param gamma list, a named list of NIG parameters, gamma=list(m,V,a,b). Where gamma$m is a numeric "location" parameter; gamma$V is a symmetric positive definite matrix representing the "scale" parameters; gamma$a and gamma$b are the "shape" and "scale" parameter of the Inverse Gamma distribution.
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
#' @references Banerjee, Sudipto. "Bayesian Linear Model: Gory Details." Downloaded from http://www. biostat. umn. edu/~ ph7440 (2008).
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

#' @title Sufficient statistics of a "GaussianNIG" object
#' @description
#' For following Gaussian-NIG model structure:
#'    \deqn{x ~ Gaussian(X^T beta,sigma^2)}
#'    \deqn{sigma^2 ~ InvGamma(a,b)}
#'    \deqn{beta ~ Gaussian(m,sigma^2 V)}
#' Where InvGamma() is the Inverse-Gamma distribution, Gaussian() is the Gaussian distribution. See \code{?dInvGamma} and \code{dGaussian} for the definitions of these distribution.\cr
#' The model structure and prior parameters are stored in a "GaussianNIG" object. \cr
#' The sufficient statistics of a set of samples (x,X) are:
#' \itemize{
#'  \item the effective number of samples N=nrow(X) or length(x)
#'  \item the covariance of X and x SXx=t(X)%*%x
#'  \item the covariance of X SX=t(X)%*%X
#'  \item the covariance of x Sx=t(x)%*%x
#' }
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
#' @references Banerjee, Sudipto. "Bayesian Linear Model: Gory Details." Downloaded from http://www. biostat. umn. edu/~ ph7440 (2008).
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

#' @title Weighted sufficient statistics of a "GaussianNIG" object
#' @description
#' For following Gaussian-NIG model structure:
#'    \deqn{x ~ Gaussian(X^T beta,sigma^2)}
#'    \deqn{sigma^2 ~ InvGamma(a,b)}
#'    \deqn{beta ~ Gaussian(m,sigma^2 V)}
#' Where InvGamma() is the Inverse-Gamma distribution, Gaussian() is the Gaussian distribution. See \code{?dInvGamma} and \code{dGaussian} for the definitions of these distribution.\cr
#' The model structure and prior parameters are stored in a "GaussianNIG" object. \cr
#' This object will be used as a place for recording and accumulating information in the related inference/sampling functions such as posterior(), posteriorDiscard(), MAP(), marginalLikelihood(), dPosteriorPredictive(), rPosteriorPredictive() and so on. \cr
#' The sufficient statistics of a set of samples (x,X) and weights ware:
#' \itemize{
#'  \item the effective number of samples N=sum(w); \cr
#'  \item the covariance of X and x SXx=t(w*X)%*%x \cr
#'  \item the covariance of X SX=t(w*X)%*%X \cr
#'  \item the covariance of x Sx=t(w*x)%*%x
#' }
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
#' @references Banerjee, Sudipto. "Bayesian Linear Model: Gory Details." Downloaded from http://www. biostat. umn. edu/~ ph7440 (2008).
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

#' @title Update a "GaussianNIG" object with sample sufficient statistics
#' @description
#' For the model structure:
#'    \deqn{x ~ Gaussian(X^T beta,sigma^2)}
#'    \deqn{sigma^2 ~ InvGamma(a,b)}
#'    \deqn{beta ~ Gaussian(m,sigma^2 V)}
#' Where InvGamma() is the Inverse-Gamma distribution, Gaussian() is the Gaussian distribution. See \code{?dInvGamma} and \code{dGaussian} for the definitions of these distribution.\cr
#' The model structure and prior parameters are stored in a "GaussianNIG" object. \cr
#' Update (m,V,a,b) by adding the information of newly observed samples (x,X). The model structure and prior parameters are stored in a "GaussianNIG" object, the prior parameters in this object will be updated after running this function.
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
#' @references Banerjee, Sudipto. "Bayesian Linear Model: Gory Details." Downloaded from http://www. biostat. umn. edu/~ ph7440 (2008).
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

#' @title Update a "GaussianNIG" object with sample sufficient statistics
#' @description
#' For the model structure:
#'    \deqn{x ~ Gaussian(X^T beta,sigma^2)}
#'    \deqn{sigma^2 ~ InvGamma(a,b)}
#'    \deqn{beta ~ Gaussian(m,sigma^2 V)}
#' Where InvGamma() is the Inverse-Gamma distribution, Gaussian() is the Gaussian distribution. See \code{?dInvGamma} and \code{dGaussian} for the definitions of these distribution.\cr
#' The model structure and prior parameters are stored in a "GaussianNIG" object. \cr
#' Contrary to posterior(), this function will update (m,V,a,b) by removing the information of observed samples (x,X). The model structure and prior parameters are stored in a "GaussianNIG" object, the prior parameters in this object will be updated after running this function.
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
#' @references Banerjee, Sudipto. "Bayesian Linear Model: Gory Details." Downloaded from http://www. biostat. umn. edu/~ ph7440 (2008).
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

#' @title Maximum A Posteriori (MAP) estimate of a "GaussianNIG" object
#' @description
#' Generate the MAP estimate of (beta,sigma^2) in following Gaussian-NIG structure:
#'    \deqn{x ~ Gaussian(X^T beta,sigma^2)}
#'    \deqn{sigma^2 ~ InvGamma(a,b)}
#'    \deqn{beta ~ Gaussian(m,sigma^2 V)}
#' Where InvGamma() is the Inverse-Gamma distribution, Gaussian() is the Gaussian distribution. See \code{?dInvGamma} and \code{dGaussian} for the definitions of these distribution.\cr
#' The model structure and prior parameters are stored in a "GaussianNIG" object. \cr
#' The MAP estimates are:
#' \itemize{
#'  \item (beta,sigma^2)_MAP = argmax p(beta,sigma^2|m,V,a,b,x,X)
#' }
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
#' @references Banerjee, Sudipto. "Bayesian Linear Model: Gory Details." Downloaded from http://www. biostat. umn. edu/~ ph7440 (2008).
MAP.GaussianNIG <- function(obj,...){
    D <- length(obj$gamma$m)                      #dimension
    list(betaMAP=obj$gamma$m,
         sigmaMAP=obj$gamma$b/(obj$gamma$a+1+D/2))
}

#' @title Mean Posterior Estimate (MPE) of a "GaussianNIG" object
#' @description
#' Generate the MPE estimate of (beta,sigma^2) in following Gaussian-NIG structure:
#'    \deqn{x ~ Gaussian(X^T beta,sigma^2)}
#'    \deqn{sigma^2 ~ InvGamma(a,b)}
#'    \deqn{beta ~ Gaussian(m,sigma^2 V)}
#' Where InvGamma() is the Inverse-Gamma distribution, Gaussian() is the Gaussian distribution. See \code{?dInvGamma} and \code{dGaussian} for the definitions of these distribution.\cr
#' The model structure and prior parameters are stored in a "GaussianNIG" object. \cr
#' The MPEs are E(beta,sigma^2|m,V,a,b,X,x)
#' @seealso \code{\link{GaussianNIG}}
#' @param obj A "GaussianNIG" object.
#' @param ... Additional arguments to be passed to other inherited types.
#' @return A named list, the MPE estimate of beta and sigma^2.
#' @export
#' @references Banerjee, Sudipto. "Bayesian Linear Model: Gory Details." Downloaded from http://www. biostat. umn. edu/~ ph7440 (2008).
MPE.GaussianNIG <- function(obj,...){
    stop("MPE method for class 'GaussianNIG' is not implemented yet")
}

#' @title Marginal likelihood of a "GaussianNIG" object
#' @description
#' Generate the marginal likelihood of the following model structure:
#'    \deqn{x ~ Gaussian(X^T beta,sigma^2)}
#'    \deqn{sigma^2 ~ InvGamma(a,b)}
#'    \deqn{beta ~ Gaussian(m,sigma^2 V)}
#' Where InvGamma() is the Inverse-Gamma distribution, Gaussian() is the Gaussian distribution. See \code{?dInvGamma} and \code{dGaussian} for the definitions of these distribution.\cr
#' The model structure and prior parameters are stored in a "GaussianNIG" object. \cr
#' Marginal likelihood = p(x|m,V,a,b,X).
#' @seealso \code{\link{GaussianNIG}}, \code{\link{marginalLikelihood_bySufficientStatistics.GaussianNIG}}
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
#' @references Banerjee, Sudipto. "Bayesian Linear Model: Gory Details." Downloaded from http://www. biostat. umn. edu/~ ph7440 (2008).
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

#' @title Marginal likelihood of a "GaussianNIG" object, using sufficient statistics
#' @description
#' Generate the marginal likelihood of a set of observations of the following model structure:
#'    \deqn{x ~ Gaussian(X^T beta,sigma^2)}
#'    \deqn{sigma^2 ~ InvGamma(a,b)}
#'    \deqn{beta ~ Gaussian(m,sigma^2 V)}
#' Where InvGamma() is the Inverse-Gamma distribution, Gaussian() is the Gaussian distribution. See \code{?dInvGamma} and \code{dGaussian} for the definitions of these distribution.\cr
#' The model structure and prior parameters are stored in a "GaussianNIG" object. \cr
#' Marginal likelihood = p(x|m,V,a,b,X)
#' @seealso \code{\link{GaussianNIG}}, \code{\link{marginalLikelihood.GaussianNIG}}
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
#' @references Banerjee, Sudipto. "Bayesian Linear Model: Gory Details." Downloaded from http://www. biostat. umn. edu/~ ph7440 (2008).
marginalLikelihood_bySufficientStatistics.GaussianNIG <- function(obj,ss,LOG=TRUE,...){
    if(missing(ss)) stop("'ss' not specified!")
    if(!.is(ss,"ssGaussianLinear")) stop("'ss' must be of class 'ssGaussianLinear', you need to use sufficientStatistics() to generate 'ssGaussianLinear' objects")
    
    obj2 <- GaussianNIG(objCopy = obj) #copy obj to obj2
    posterior.GaussianNIG(obj = obj2,ss=ss) #update the posteriors

    logp <- pdsDeterminant(S=obj2$gamma$V,LOG=TRUE)/2 + lgamma(obj2$gamma$a) - obj2$gamma$a*log(obj2$gamma$b) -
        (pdsDeterminant(S=obj$gamma$V,LOG=TRUE)/2 + lgamma(obj$gamma$a) - obj$gamma$a*log(obj$gamma$b))-
        ss$N/2*log(2*pi)

    if(!LOG) logp <- exp(logp)
    
    drop(logp)                          #remove dimensions
}

#' @title Posterior predictive density function of a "GaussianNIG" object
#' @description
#' Generate the the density value of the posterior predictive distribution of the following structure:
#'    \deqn{x ~ Gaussian(X^T beta,sigma^2)}
#'    \deqn{sigma^2 ~ InvGamma(a,b)}
#'    \deqn{beta ~ Gaussian(m,sigma^2 V)}
#' Where InvGamma() is the Inverse-Gamma distribution, Gaussian() is the Gaussian distribution. See \code{?dInvGamma} and \code{dGaussian} for the definitions of these distribution.\cr
#' The model structure and prior parameters are stored in a "GaussianNIG" object. \cr
#' Posterior predictive density is p(x|m,V,a,b,X).
#' @seealso \code{\link{GaussianNIG}}, \code{\link{dPosteriorPredictive.GaussianNIG}}, \code{\link{marginalLikelihood.GaussianNIG}}
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
#' @references Banerjee, Sudipto. "Bayesian Linear Model: Gory Details." Downloaded from http://www. biostat. umn. edu/~ ph7440 (2008).
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

#' @title Posterior predictive random generation of a "GaussianNIG" object
#' @description
#' Generate random samples from the posterior predictive distribution of the following structure:
#'    \deqn{x ~ Gaussian(X^T beta,sigma^2)}
#'    \deqn{sigma^2 ~ InvGamma(a,b)}
#'    \deqn{beta ~ Gaussian(m,sigma^2 V)}
#' Where InvGamma() is the Inverse-Gamma distribution, Gaussian() is the Gaussian distribution. See \code{?dInvGamma} and \code{dGaussian} for the definitions of these distribution.\cr
#' The model structure and prior parameters are stored in a "GaussianNIG" object. \cr
#' Posterior predictive is a distribution of x|m,V,a,b,X
#' @seealso \code{\link{GaussianNIG}}, \code{\link{dPosteriorPredictive.GaussianNIG}}
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
#' @references Banerjee, Sudipto. "Bayesian Linear Model: Gory Details." Downloaded from http://www. biostat. umn. edu/~ ph7440 (2008).
rPosteriorPredictive.GaussianNIG <- function(obj,n,X,...){
    if(missing(X)) stop("'X' not specified!")
    if(!.is(X,"matrix")){
        stop("'X' must be a matrix!")
    }
    if(ncol(X)!=length(obj$gamma$m)) stop("number of columns in X doesn't match number of elements in beta")
    N <- nrow(X)
    rT(n=n,mu=X%*%obj$gamma$m,Sigma = obj$gamma$b/obj$gamma$a*(diag(N)+X%*%obj$gamma$V%*%t(X)),df = 2*obj$gamma$a)
}
