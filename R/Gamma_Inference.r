
#' @include Bayesian_Bricks.r


#' @title Density function of Inverse-Gamma distribution
#' @description
#' For a random variable x, the density function of Inverse-Gamma distribution is defined as:
#'    \deqn{(rate^shape)/Gamma(shape) x^{-shape-1} exp(-rate/x)}
#' @param x numeric, positive numeric values.
#' @param shape numeric, the shape parameter of gamma distribution.
#' @param scale numeric, the scale, or inverse-scale parameter of gamma distribution. The 'rate' parameter in Gamma is the 'scale' parameter in InvGamma
#' @param LOG logical, return log density of LOG=TRUE, default TRUE.
#' @return A numeric vector, the density values.
#' @export
dInvGamma <- function(x,shape,scale,LOG=TRUE){
    logp <- shape*log(scale) - lgamma(shape)+(-shape-1)*log(x) - scale/x
    if(!LOG) logp <- exp(logp)
    logp
}

#' @title Random number generation of Inverse-Gamma distribution
#' @description
#' Generation random samples from Inverse-Gamma distribution. For a random variable x, the density function is defined as:
#'    \deqn{(rate^shape)/Gamma(shape) x^{-shape-1} exp(-rate/x)}
#' Where Gamma() is the Gamma function.
#' @param n integer, number of samples to be generated.
#' @param shape numeric, the shape parameter of gamma distribution.
#' @param scale numeric, the scale, or inverse-scale parameter of gamma distribution. The 'rate' parameter in Gamma is the 'scale' parameter in InvGamma.
#' @return A numeric vector, samples of Inverse-Gamma distribution.
#' @export
#' @import stats
rInvGamma <- function(n,shape,scale){
    1/rgamma(n=n,shape = shape,rate = scale)
}

#' @title Random generation for Wishart distribution
#' @description
#' Generate random samples from Wishart distribution. For a random matrix x, the density function of Wishart distribution is defined as:
#'     \deqn{(2^{(df p)/2} Gamma_p(df/2) |rate|^{-df/2})^{-1} |x|^{(df-p-1)/2} exp(-1/2 tr(x rate))}
#' Where x is a pxp symmetric positive definite matrix, Gamma_p() is the multivariate Gamma function of dimension p.
#' @param df numeric, the degree of freedom.
#' @param rate matrix, a symmetric positive-definite matrix, the 'rate', or 'inverse-scale' parameter. The 'rate' parameter in Wishart is the 'scale' parameter in InvWishart
#' @param scale, matrix, the inverse of rate. Only one of 'rate' and 'scale' should be non-NULL.
#' @return A symmetric positive-definite matrix.
#' @export
#' @examples
#' rate <- crossprod(matrix(rnorm(15),5,3)) #the prior inverse-scale
#' m <- matrix(0,3,3)
#' ## get 1000 samples and calculate the sample mean
#' for(i in 1:100){
#'     m <- m+rWishart(df=5,rate=rate)/100
#' }
#' ## m should roughly equal to df*inverse(rate):
#' m
#' pdsInverse(rate)*5
#' ## try generating samples with 'rate' parameter:
#' scale <- pdsInverse(rate)
#' m2 <- matrix(0,3,3)
#' for(i in 1:100){
#'     m2 <- m2+rWishart(df=5,scale=scale)/100
#' }
#' ## m2 should roughly equal df*scale:
#' m2
#' 5*scale
#' @references Smith, W. B., and R. R. Hocking. "Algorithm as 53: Wishart variate generator." Journal of the Royal Statistical Society. Series C (Applied Statistics) 21.3 (1972): 341-345.
#' @import stats
rWishart <- function(df,rate=NULL,scale=NULL){
    if(missing(df)) stop("'df' not specified!")
    if(is.null(rate) & is.null(scale)) stop("one of 'rate' or 'scale' should be non-NULL.")
    if(!is.null(rate)){
        if(!.is(rate,"matrix")){
            stop("'rate' must be a matrix!")
        }else if(nrow(rate)!=ncol(rate)){
            stop("'rate' must be a square matrix!")
        }
        if(df<nrow(rate)) stop("'df' must be greater than nrow(rate)!")
        D <- nrow(rate)                     #dimension
    }else if(!is.null(scale)){
        if(!.is(scale,"matrix")){
            stop("'scale' must be a matrix!")
        }else if(nrow(scale)!=ncol(scale)){
            stop("'scale' must be a square matrix!")
        }
        if(df<nrow(scale)) stop("'df' must be greater than nrow(rate)!")
        D <- nrow(scale)                     #dimension
    }

    
    A <- matrix(0,D,D)
    diag(A) <- sqrt(vapply((df+1)-1:D,function(dfchisq){
        rchisq(n=1,df = dfchisq)
    },FUN.VALUE = numeric(1)))
    if(D>1L){
        A[upper.tri(A)] <- rnorm(D*(D-1L)/2L)
    }

    if(!is.null(scale)){
        return(crossprod(A%*%chol(scale)))
    }else{
        return(crossprod(A%*%t(pdsInverse(S=rate,returnUpper = TRUE))))
    }
}

#' @title Density function of Wishart distribution
#' @description
#' For a random matrix x, the density function of Wishart distribution is defined as:
#'     \deqn{(2^{(df p)/2} Gamma_p(df/2) |rate|^{-df/2})^{-1} |x|^{(df-p-1)/2} exp(-1/2 tr(x rate))}
#' Where x is a pxp symmetric positive definite matrix, Gamma_p() is the multivariate Gamma function of dimension p.
#' @param x matrix, a symmetric positive-definite matrix.
#' @param df numeric, the degree of freedom.
#' @param rate matrix, a symmetric positive-definite matrix, the 'rate', or 'inverse-scale' parameter. The 'rate' parameter in Wishart is the 'scale' parameter in InvWishart
#' @param LOG logical, return log density of LOG=TRUE, default TRUE.
#' @return A numeric vector, the density values.
#' @export
#' @examples
#' ##generate a symmetric positive-definite matrix
#' x <- crossprod(matrix(rnorm(15),5,3))
#' rate <- crossprod(matrix(rnorm(15),5,3)) #the prior inverse-scale of x
#' dWishart(x,df = 5,rate = rate,LOG = TRUE)
#' dWishart(x,df = 5,rate = rate,LOG = FALSE)
#' @references Wishart, John. "The generalized product moment distribution in samples from a normal multivariate population." Biometrika (1928): 32-52.
#' @references MARolA, K. V., JT KBNT, and J. M. Bibly. Multivariate analysis. AcadeInic Press, Londres, 1979.
dWishart <- function(x,df,rate,LOG=TRUE){
    if(missing(df)|missing(rate)) stop("'n' or 'rate' not specified!")
    if(!.is(x,"matrix")){
        stop("'x' must be a matrix!")
    }else if(nrow(x)!=ncol(x)){
        stop("'x' must be a square matrix!")
    }
    if(!.is(rate,"matrix")){
        stop("'rate' must be a matrix!")
    }else if(nrow(rate)!=ncol(rate)){
        stop("'rate' must be a square matrix!")
    }
    if(df<nrow(rate)) stop("'df' must be greater than nrow(rate)!")
    p <- ncol(x)                        #dimension
    logp <- (df-p-1)/2 * pdsDeterminant(S=x,LOG=TRUE) - sum(diag(x%*%rate))/2 -
        ( df*p/2*log(2) + .lmvgamma(a=df/2,p=p) - df/2*pdsDeterminant(S=rate,LOG = TRUE))
    if(!LOG) logp <- exp(logp)
    logp
}

#' @title Random generation for Inverse-Wishart distribution
#' @description
#' Generate random samples from Inverse-Wishart distribution. For a random matrix x, the density function of Inverse-Wishart is defined as:
#'   \deqn{(2^{(df p)/2} Gamma_p(df/2) |scale|^{-df/2})^{-1} |x|^{(-df-p-1)/2} exp(-1/2 tr(x^{-1} scale))}
#' Where x is a pxp symmetric positive definite matrix, Gamma_p() is the multivariate Gamma function of dimension p.
#' @param df numeric, the degree of freedom.
#' @param scale matrix, a symmetric positive-definite matrix, the 'scale' parameter. The 'rate' parameter in Wishart is the 'scale' parameter in InvWishart.
#' @return A symmetric positive-definite matrix.
#' @export
#' @examples
#' scale <- crossprod(matrix(rnorm(15),5,3)) # the prior scale
#' m <- matrix(0,3,3)
#' ## get 1000 samples and calculate the sample mean
#' for(i in 1:1000){
#'     m <- m+rInvWishart(df=5,scale=scale)/1000
#' }
#' ## m should roughly equal scale/(df-p-1), p is the dimension.
#' m
#' scale/(5-3-1)
#' @references Hoff, Peter D. A first course in Bayesian statistical methods. Vol. 580. New York: Springer, 2009.
rInvWishart <- function(df,scale){
    if(missing(df)|missing(scale)) stop("'df' or 'scale' not specified!")
    if(df<nrow(scale)) stop("'df' must be greater than nrow(scale)!")
    pdsInverse(rWishart(df = df,rate = scale))      #rate in Wishart is scale in InverseWishart
}

#' @title Density function of Inverse-Wishart distribution
#' @description
#' For a random matrix x, The density function of Inverse-Wishart distribution is defined as:
#'   \deqn{(2^{(df p)/2} Gamma_p(df/2) |scale|^{-df/2})^{-1} |x|^{(-df-p-1)/2} exp(-1/2 tr(x^{-1} scale))}
#' Where x is a pxp symmetric positive definite matrix, Gamma_p() is the multivariate Gamma function of dimension p.
#' @param x matrix, a symmetric positive-definite matrix.
#' @param df numeric, the degree of freedom.
#' @param scale matrix, a symmetric positive-definite matrix, the 'scale' parameter. The 'rate' parameter in Wishart is the 'scale' parameter in InvWishart.
#' @param LOG logical, return log density of LOG=TRUE, default TRUE.
#' @return A numeric vector, the density values.
#' @export
#' @examples
#' x <- crossprod(matrix(rnorm(15),5,3)) #generate a symmetric positive-definite matrix
#' scale <- crossprod(matrix(rnorm(15),5,3)) #the prior scale of x
#' dInvWishart(x,df = 5,scale = scale,LOG = TRUE)
#' dInvWishart(x,df = 5,scale = scale,LOG = FALSE)
#' @references Wishart, John. "The generalized product moment distribution in samples from a normal multivariate population." Biometrika (1928): 32-52.
#' @references MARolA, K. V., JT KBNT, and J. M. Bibly. Multivariate analysis. AcadeInic Press, Londres, 1979.
dInvWishart <- function(x,df,scale,LOG=TRUE){
    if(!.is(x,"matrix")){
        stop("'x' must be a matrix!")
    }else if(nrow(x)!=ncol(x)){
        stop("'x' must be a square matrix!")
    }
    if(df<nrow(scale)) stop("'df' must be greater than nrow(scale)!")
    p <- ncol(x)                        #dimension
    logp <- -(df+p+1)/2 * pdsDeterminant(S=x,LOG=TRUE) - sum(diag(solve(x)%*%scale))/2 -
        ( df*p/2*log(2) + .lmvgamma(a=df/2,p=p) - df/2*pdsDeterminant(S=scale,LOG = TRUE))
    if(!LOG) logp <- exp(logp)
    logp
    
}
