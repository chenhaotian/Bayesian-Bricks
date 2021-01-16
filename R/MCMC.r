
#' Metropolis-Hastings sampler
#' 
#' sample xhat from the target distribution p(xhat), with given proposal distribution q(xhat|x), the acceptance rate is: \cr
#'   min( 1 , dp(xhat)/dp(x) * dq(x|xhat)dq(xhat|x) ) \cr
#' where dp() is the density function of the target distribution, dq() the the density function of the proposal distribution. A new sample xhat is drawn from the sampler of the proposal distribution rq(). See examples.
#' 
#' @param nsamples integer, number of samples to draw
#' @param xini initial sample, the chain of samples starts from here. xini must be a matrix of one row, or a numeric vector that will be converted to a matrix of one row.
#' @param dp function(x), the LOG density function of the target distribution log dp(x), DON'T FORGET THE LOG.
#' @param dq function(xhat,x), the LOG density function of the proposal distribution log dq(xhat | x), DON'T FORGET THE LOG.
#' @param rq function(x), the generator of the proposal distribution rq(xhat | x).
#' @return a matrix of nsamples rows.
#' @export
#' @examples
#' \donttest{
#' 
#' ## example1: independent Metropolis-Hastings algorithm, get 5000 samples from Beta(2.7,6.3)
#' ## with independent uniform proposal U(0,1), and independent normal proposal N(0.5,1).
#' 
#' ## step1: define p() and q()
#' dp <- function(x) if(x>0&x<1) dbeta(x,2.7,6.3,log = TRUE) else -Inf
#' dq1 <- function(xnew,x) 0               #uniform proposal log density
#' dq2 <- function(xnew,x) dnorm(xnew,0.5,1) #normal proposal log density
#' rq1 <- function(x) runif(1,0,1)           #uniform proposal sampler
#' rq2 <- function(x) rnorm(1,0.5,1)         #normal proposal sampler
#' ## step2: get 5000 samples, with two different proposals
#' X1 <- MetropolisHastings(nsamples = 5000,xini = runif(1,0,1),dp=dp,dq=dq1,rq=rq1)
#' X2 <- MetropolisHastings(nsamples = 5000,xini = runif(1,0,1),dp=dp,dq=dq2,rq=rq2)
#' ## step3: plot the result, calculate acceptance rate
#' sum(diff(X1)!=0)/nrow(X1)                 #the acceptance rate of uniform proposal
#' sum(diff(X2)!=0)/nrow(X2)                 #the acceptance rate of normal proposal
#' ## Clearly Uniform, compare to Normal, can better resemble Beta, so the acceptance rate is higher
#' ## plot the results
#' hist(X1)
#' hist(X2)
#' hist(rbeta(5000,2.7,6.3))
#' 
#' ## example2: independent Metropolis-Hastings algorithm, sample from an improper distribution
#' ## p(x) = -|x|+1, where -1<x<1, with independent uniform proposal U(-1,1)
#' 
#' ## step1: define p() and q()
#' dp <- function(x) log(-abs(x)+1)        #log dp
#' dq <- function(xnew,x) 1
#' rq <- function(x) runif(1,-1,1)         #make sure -1<x<1
#' ## step2: get 5000 samples
#' X <- MetropolisHastings(nsamples = 5000,xini = runif(1,-1,1),dp=dp,dq=dq,rq=rq)
#' ## step3: plot the result, calculate acceptance rate
#' hist(X)
#' sum(diff(X)!=0)/nrow(X)                 #the acceptance rate
#' 
#' ## example3: random walk Metropolis-Hastings algorithm, sample from a
#' ## normal mixture 0.2*N(1,1)+0.8*N(-5,1), with symmetric proposal xhat ~ U(x-l,x+l),
#' ## compare different values of l.
#' 
#' ## step1: define p() and q()
#' dp <- function(x) log(dnorm(x,1,1)*0.2+dnorm(x,-5,1)*0.8)
#' ## a symmetric proposal has no influence to the acceptance rate, so a constant function
#' ## would suffice.
#' dq <- function(xnew,x) 1
#' rq1 <- function(x) runif(1,x-0.01,x+0.01)
#' rq2 <- function(x) runif(1,x-2,x+2)
#' ## step2: get 5000 samples
#' X1 <- MetropolisHastings(nsamples = 5000,xini = rnorm(1),dp=dp,dq=dq,rq=rq1)
#' X2 <- MetropolisHastings(nsamples = 50000,xini = rnorm(1),dp=dp,dq=dq,rq=rq2)
#' ## step3: plot the result, calculate acceptance rate
#' sum(diff(X1)!=0)/nrow(X1)
#' sum(diff(X2)!=0)/nrow(X2)
#' ## plot the results
#' hist(X1,xlim = c(-10,5))
#' hist(X2,xlim = c(-10,5))
#' hist(c(rnorm(1000,1,1),rnorm(4000,-5,1)),xlim = c(-10,5))
#'
#' ## note that X1 has a higher acceptance rate comparing to X2, though it performs poorer.
#' ## So we use Kolmogorov-Smirnov to test the real performance of X1 and X2:
#' ks.test(jitter(X1),c(rnorm(1000,1,1),rnorm(4000,-5,1))) #ks.test() assumes continuous
#' ## samples doesn't contain equal values, otherwise there will be a warning.so use jitter() to
#' ## remove the equals
#' ks.test(jitter(X2),c(rnorm(1000,1,1),rnorm(4000,-5,1)))
#' ## it turns out that even though X2 looks better from the histogram, it still doesn't
#' ## pass the KS test.
#' 
#' ## example4: hybrid Metropolis-Hastings algorithm, questions same as previous example,
#' ## but use a mixutre proposal instead.
#' 
#' ## we use mixture proposal to capture both the local and global areas of the target distribution
#' ## step1: define p() and q()
#' dp <- function(x) log(dnorm(x,1,1)*0.2+dnorm(x,-5,1)*0.8)
#' ## a symmetric proposal has no influence to the acceptance rate,
#' ## so a constant function would suffice.
#' dq <- function(xnew,x) 1
#' ##70% local, 30% global
#' rq <- function(x) if(runif(1)<0.7) runif(1,x-0.1,x+0.1) else runif(1,x-3,x+3)
#' ## step2: get 5000 samples
#' X <- MetropolisHastings(nsamples = 5000,xini = rnorm(1),dp=dp,dq=dq,rq=rq)
#' ## step3: plot the result, calculate acceptance rate
#' sum(diff(X)!=0)/nrow(X)
#' ## plot the results
#' hist(X,xlim = c(-10,5))
#' hist(c(rnorm(1000,1,1),rnorm(4000,-5,1)),xlim = c(-10,5))
#' 
#' ## perform the KS test again, this time it says there's no significance difference between the MH
#' ## and the real samples. ks.test() assumes continuous samples doesn't contain equal values,
#' ## otherwise there will be a warning.so use jitter() to remove the equals
#' ks.test(jitter(X),c(rnorm(1000,1,1),rnorm(4000,-5,1))) 
#' }
MetropolisHastings <- function(nsamples,xini,dp,dq,rq){
    if(is.vector(xini)) xini <- matrix(xini,nrow = 1L)
    else if(!is.matrix(xini)) stop("xini must be a numeric vector or a numeric matrix of only one row")
    if(!is.integer(nsamples)) nsamples <- as.integer(nsamples)
    if(nsamples<0L) stop("Error in MetropolisHastings(): nsamples must be greater than 0!")
    out <- matrix(0,nrow = nsamples+1L,ncol = length(xini))
    out[1L,] <- xini
    if(!is.function(dp) | !is.function(dq) | !is.function(rq)) stop("dp, dq and rq must be functions!")
    xnew <- numeric(ncol(out))          #initialize place holder
    x <- xnew                           #initialize place holder
    for(i in 2L:(nsamples+1L)){
        x <- out[i-1L,]
        xnew <- rq(x)
        if(runif(1) < exp(dp(xnew) - dp(x)-
                          dq(x,xnew) + dq(xnew,x))){
            out[i,] <- xnew
        }else{
            out[i,] <- out[i-1L,]
        }
    }
    out[-1,,drop=FALSE]
}

