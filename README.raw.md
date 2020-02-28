# bbricks

[![Build Status](https://travis-ci.com/chenhaotian/Bayesian-Bricks.svg?token=hByNmnrjfd4L3sAdyCVy&branch=master)](https://travis-ci.com/chenhaotian/Bayesian-Bricks)

**bbricks** provides a collection of inference tools and conditional probability distributions(CPDs) to facilitate Bayesian modeling. This is a package designed for statisticians, and the ones who want to learn the basic Bayesian mindsets.

See [Mindset](#mindset) and [Examples](#examples) to get started.

**Installation:**

```R
# install development version from GitHub:
# install.packages("devtools")
devtools::install_github("chenhaotian/Bayesian-Bricks")
```

**Contents:**

[Mindset](#Mindset)

[Examples](#examples)

+ [Mixture of Gaussian](#mixture-of-gaussian)
+ [Dirichlet Process Mixture Model](#dirichlet-process-mixture-model)
+ [Mixture Model with Partially Observed Cluster Labels](#mixture-model-with-partially-observed-cluster-labels)
+ [Hierarchical Mixture Models](#hierarchical-mixture-models)
+ [Topic Modeling with HDP](topic-modeling-with-hdp)
+ [Hierarchical Topic Modeling with HDP2](#hierarchical-topic-modeling-with-hdp2)
+ [Infinite State Hidden Markov Model](#infinite-state-hidden-markov-model)
+ [Bayesian Linear Regression](#bayesian-linear-regression)
+ [Hierarchical Bayesian: Estimate Cancer Mortality Rates](hierarchical-bayesian-estimate-cancer-mortality-rates)



## Mindset

The idea of **bbricks** came from the fact that modeling in Bayesian statistics is nothing more than applying set of **tasks** on a specific **model structure**.

Where the most frequently applied **tasks** are:
+ Update prior info into posterior when new samples are observed.
+ Sample from the posterior distribution.
+ Calculate marginal likelihood of the posterior distribution.
+ Calculate posterior predictive of the posterior distribution.
+ ...

And the **model structure**s are just generalizations of $3$ basic Bayesian modeling structures:
![](./notes_pictures/basicStructures.png)
Where

+ $(a)$ is the most basic "parameter-observation" structure. Models like Gaussian, Gamma and Exponential are in this category.
+ $(b)$ is the "prior-posterior" structure. Models like Gaussian-NIW(Gaussian observations with NIW prior), Categorical-Dirichlet(Categorical observations with Dirichlet prior) are in this category.
+ $(c)$ is the "hierarchical-Bayesian" structure, and $(d)$ is the same hierarchical structure but with more hierarchies. Models like Hierarchical Dirichlet Process(HDP) and HDP with additional hidden layers are in this category.

**bbricks** tries to provide a `type/class` for each basic **model structure** and a `function/method` for each **task**.

See [Examples](#examples) for details.

## Examples

Here's a list of examples:

[Mixture of Gaussian](#mixture-of-gaussian)

[Dirichlet Process Mixture Model](#dirichlet-process-mixture-model)

[Mixture Model with Partially Observed Cluster Labels](#mixture-model-with-partially-observed-cluster-labels)

[Hierarchical Mixture Models](#hierarchical-mixture-models)

[Topic Modeling with HDP](topic-modeling-with-hdp)

[Hierarchical Topic Modeling with HDP2](#hierarchical-topic-modeling-with-hdp2)

[Infinite State Hidden Markov Model](#infinite-state-hidden-markov-model)

[Bayesian Linear Regression](#bayesian-linear-regression)

[Estimate Cancer Mortality Rates](estimate-cancer-mortality-rates)



### Mixture of Gaussian

A mixture of Gaussian has following graph structure:

![](./notes_pictures/mixtureModel.png)

The CPDs are:
$$
\begin{align}
\pi |\alpha &\sim Dirichlet(\alpha) \\
z|\pi & \sim Categorical(\pi) \\
\theta_z | \gamma &\sim NIW(\gamma) \\
x| z,\theta_z &\sim  Gaussian(\theta_z)
\end{align}
$$

Where $NIW(\gamma)$ is the Normal-Inverse-Wishart distribution with parameter $\gamma = (m,k,v,S)$. $m$ is a numeric vector representing the "location parameter", $S$ is a symmetric positive definitive matrix representing the "scale parameter", $k$ and $v$ are degree of freedoms.

A mixture model can be see as a combination of two "prior-posterior" structures(As shown in [Mindset](#mindset) graph $(b)$): One Categorical-Dirichlet structure $\alpha \rightarrow \pi \rightarrow z$ for the hidden cluster labels. and one Gaussian-NIW structure $\gamma \rightarrow \theta_z \rightarrow x$ for the observation distribution.

In **bbricks** these two structures are initialized with a `CatDirichlet` object and a `GaussianNIW` object. See R example for details:

```R
## Get the MAP estimate of pi and theta using EM algorithm.

library(bbricks)

## load some mixture of Gaussian samples.
## mmData is a numeric matrix with 2 columns, each row is a sample
## see ?mmData for details
data(mmData)
K <- 4L                                 #number of clusters(mixtures components)
z <- matrix(runif(nrow(mmData)*K),nrow(mmData),K) #the expected cluster label of each observation
allK <- 1L:K    #temp variable, all component labels
allZ <- rep(allK,each=nrow(mmData))     #temp variable, all possible cluster labels for all observations
## z, pi and alpha are distributed as a Categorical-Dirichlet sturcture:
mc <- CatDirichlet(gamma = list(alpha=0.5,uniqueLabels=allK)) # create a CatDirichlet object to track the posterior info, see ?CatDirichlet for details
## each component distribution is a Gaussian-NIW structure:
ec <- replicate(K,GaussianNIW(gamma = list(m=c(0,0),k=0.00001,v=2,S=diag(2)))) # create a GaussianNIW object to track the posterior info of each mixture component, see ?GaussianNIW for details
mcMAP <- MAP(mc)                        #initialize the MAP estimate of pi
ecMAP <- replicate(K,list(muMAP=runif(2),sigmaMAP=diag(2)),simplify = FALSE) #initialize the MAP estimate of theta
## The main EM loop
maxit <- 100                            #number of EM loops
it <- 1
while(it<=maxit){
    ## E-step---------------------------------------------------------
    ## calculate the expected cluster labels: p(z|pi,theta)
    for(k in allK) z[,k] <- dGaussian(x=mmData,mu = ecMAP[[k]]$muMAP,Sigma=ecMAP[[k]]$sigmaMAP)+log(mcMAP[k])
    z <- exp(z-logsumexp(z))            #use logsumexp() to avoid numerical underflow
    ## M-step---------------------------------------------------------
    ## calculate the weighted sufficient statistics, based on results of the E-step:
    ssComponents <- lapply(allK,function(k){
        sufficientStatistics_Weighted(obj = ec[[k]],x=mmData,w=z[,k])
    })                                  #the weighted sufficient statistics of each Gaussian component
    ssPi <- sufficientStatistics_Weighted(obj = mc,x=allZ,w=as.vector(z)) #the weighted sufficient statistics of the cluster label distribution
    ## use the sufficient statistics to update the prior distributions:
    for(k in allK) posterior(obj = ec[[k]],ss=ssComponents[[k]]) #update component distributions
    posterior(obj = mc,ss = ssPi)                                #update cluster label distribution
    ## calculate the MAP estimates from posterior:
    mcMAP <- MAP(mc)
    ecMAP <- lapply(ec,MAP)
    ## Reset the priors for next EM loop-----------------------------------------
    ## to prepare for the next EM iteration, discard the sufficient statistics info from the posteriors:
    for(k in allK) posteriorDiscard(obj = ec[[k]],ss=ssComponents[[k]])
    posteriorDiscard(obj = mc,ss = ssPi)
    ## increase the iteration counter
    it <- it+1
}

plot(mmData,col=apply(z,1,which.max)) #plot the best estimates
mcMAP                                 #the MAP estimate of pi
ecMAP                                 #the MAP estimate of theta_z
```



### Dirichlet Process Mixture Model

The graph structure of Dirichlet Process Mixture Model(DP-MM) is exactly the same as a standard mixture model, except that the number of mixture components is not predetermined:

![](./notes_pictures/mixtureModelDP.png)

The CPDs of a DP-MM is similar to the ones shown in [Mixture of Gaussian](#mixture-of-gaussian), the only difference is the distribution of $\pi|\alpha$ is a Dirichlet process rather than a Dirichlet distribution, for example if the observations are Gaussian distributed, the CPDs will be:
$$
\begin{align}
\pi |\alpha &\sim DP(\alpha) \\
z|\pi & \sim Categorical(\pi) \\
\theta_z | \gamma &\sim NIW(\gamma) \\
x| z,\theta_z &\sim  Gaussian(\theta_z)
\end{align}
$$
Where $DP(\alpha)$ is a Dirichlet process on positive integers with "concentration parameter" $\alpha$, the "base measure", which is an uniform distribution on positive integers, is omitted from the formula.  $NIW(\gamma)$ is the Normal-Inverse-Wishart distribution with parameter $\gamma = (m,k,v,S)$. $m$ is a numeric vector representing the "location parameter", $S$ is a symmetric positive definitive matrix representing the "scale parameter", $k$ and $v$ are degree of freedoms.

A DP-MM can be see as a combination of two "prior-posterior" structures(As shown in [Mindset](#mindset) graph $(b)$): One Categorical-DirichletProcess structure for the hidden cluster label distribution $\alpha \rightarrow \pi \rightarrow z$, which we call it a "DP on positive integers". And one structure for the observation distribution $\gamma \rightarrow \theta_z \rightarrow x$.

To simplify the calculations, **bbricks** provides an `"DP"` type to represent all Dirichlet process structures. An object of type  `"DP"` is in essence a combination of a `"CatDP"` object, which encodes the $\alpha \rightarrow \pi \rightarrow z$ structure, i.e. a Dirichlet process on positive integers, and an arbitrary `"BasicBayesian"` object, which encodes the $\gamma \rightarrow \theta_z \rightarrow x$ structure. (in **bbricks**, all models with same structure as [Mindset](#mindset) graph $(b)$ are `"BasicBayesian" `s, such as `"GaussianNIW"`, `"CatDirichlet"` and even `"CatDP"`) 

 See R example for details:

```R
## Sample cluster labels z from DP-MM using Gibbs sampling

library(bbricks)

## load some mixture of Gaussian samples.
## mmData is a numeric matrix with 2 columns, each row is a sample
## see ?mmData for details
data(mmData)
maxit <- 100                            #number of total samples
burnin <- 50                            #number of burnin samples
## Step1: Initialization -----------------------------------------
obj <- DP(gamma = list(alpha=10,H0aF="GaussianNIW",parH0=list(m=c(0,0),k=0.001,v=2,S=diag(2)))) #create a DP object to track all the changes, the DP object in this case is a combination of a CatDP object and a GaussianNIW object
z <- matrix(1L,nrow(mmData),maxit-burnin)    #place-holder for the sampled z
ss <- sufficientStatistics(obj,x=mmData,foreach = TRUE) #sufficient statistics of each observed sample
N <- nrow(mmData)
for(i in 1L:N){ # initialize labels before Gibbs sampling
    z[i,1] <- rPosteriorPredictive(obj = obj,n=1,x=mmData[i,,drop=FALSE])
    posterior(obj = obj,ss = ss[[i]], z = z[i,1])
}
## Step2: Main Gibbs sampling loop--------------------------------
it <- 1                                 #iteration tracker
pb <- txtProgressBar(min = 0,max = maxit,style = 3)
while(it<=maxit){
    if(it>burnin) colIdx <- it-burnin
    else colIdx <- 1
    for(i in 1L:N){
        ## remove the sample information from the posterior
        posteriorDiscard(obj = obj,ss = ss[[i]],z=z[i,colIdx])
        ## get a new sample
        z[i,colIdx] <- rPosteriorPredictive(obj = obj,n=1,x=mmData[i,,drop=FALSE])
        ## add the new sample information to the posterior
        posterior(obj = obj,ss = ss[[i]],z=z[i,colIdx])
    }
    if(it>burnin & colIdx<ncol(z)) z[,colIdx+1] <- z[,colIdx] #copy result of previous iteration
    it <- it+1
    setTxtProgressBar(pb,it)
    if(it>=maxit){cat("\n");break}
    plot(x=mmData[,1],y=mmData[,2],col=z[,colIdx]) #to see how the labels change in each iteration
}
## Step3: Estimate group labels of each observation---------------
## pick the most frequent z as the best estimate
zBest <- apply(z,1,function(l){
    tmp <- table(l)
    names(tmp)[which.max(tmp)]
})
plot(x=mmData[,1],y=mmData[,2],col=zBest)
```



### Mixture Model with Partially Observed Cluster Labels

In the dataset `mmData` of the previous example, what if we know the 50, 100, 150 and 200th samples belong to 4 different clusters(they are shown as different color and shapes in the graph below), how should we incorporate this information in the model?

![](./notes_pictures/mixtureModelPO.png)

With DP-MM, one only need to **1.** update the DP prior (as defined in previous R example) with the information of the 4 observed samples, and **2.** use the updated prior as the prior of the Gibbs sampling procedure.  These 2 steps can be achieved by adding following code after `obj <- DP(...)` in the previous R example:

```R
## 1. add the information of the 4 observed samples to the DP object
ssObserved <- sufficientStatistics(obj=obj,x=mmData[c(50,100,150,200),,drop=FALSE],foreach = TRUE)
for(i in 1L:4L) posterior(obj = obj,ss = ssObserved[[i]], z = i) # the choice of cluster label 'z' for the 4 observed samples are arbitrary, as long as they are different from each other. In this example I simply use z=1L:4L.
## 2. remove the 4 samples from the upcoming Gibbs sampling procedure
mmData <- mmData[-c(50,100,150,200),]
```

Run the code, and the result will be:

![](./notes_pictures/mixtureModelPO2.png)



### Hierarchical Mixture Models

In a hierarchical mixture model, the observation $x$ are generated by some unknown mixture components and are split into $J$ groups, all $J$ groups share the same set of mixture components but with different mixture weights.

Hierarchical Dirichlet Process(HDP) is a natural representation of a hierarchical mixture model, It has following graph structure:

![](./notes_pictures/hierarchicalMixtureModel.png)

If the observations are Gaussian distributed, the CPDs will be:
$$
\begin{align}
G_j|\gamma & \sim DP(\gamma) \\
\pi |\alpha,G_j & \sim DP(\alpha,G_j) \\
z | \pi & \sim Categorical(\pi) \\
k| z,G_j & \sim Categorical(G_j) \text{, if }z\text{ is a sample from } G_j\\
\theta_k | \psi & \sim NIW(\psi) \\
x|\theta_k & \sim Gaussian(\theta_k)
\end{align}
$$
 $DP(\gamma)$ is a Dirichlet process on positive integers with "concentration parameter" $\gamma$, the "base measure", which is an uniform distribution on positive integers, is omitted from the formula.  $DP(\alpha,G_j)$ is a Dirichlet process with concentration parameter $\alpha$ and base measure $G_j$. $NIW(\psi)$ is the Normal-Inverse-Wishart distribution with parameter $\psi = (m,k,v,S)$. $m$ is a numeric vector representing the "location parameter", $S$ is a symmetric positive definitive matrix representing the "scale parameter", $k$ and $v$ are degree of freedoms.

The distribution of $(\gamma, \alpha, G_j, \pi_j , z, k)$ is a "HDP on positive integers". HDP on positive integers are usually represented in a much simpler and compact way(though not easier to understand) in most literatures:

![](./notes_pictures/HDP.png)

From the compact representation we can see that HDP on positive integers is following the "Hierarchical Bayesian" structure shown in [Mindset](#mindset) graph $(c)$.

To simplify the calculations, **bbricks** provides an `"HDP"` type to represent all Hierarchical Dirichlet process structures. An object of type  `"HDP"` is in essence a combination of a `"CatHDP"` object, which encodes the distribution of $(\gamma, \alpha, G_j, \pi_j , z, k)$, i.e. a HDP on positive integers; and an arbitrary `"BasicBayesian"` object, which encodes the $\gamma \rightarrow \theta_z \rightarrow x$ structure. (in **bbricks**, all models with same structure as [Mindset](#mindset) graph $(b)$ are `"BasicBayesian" `s, such as `"GaussianNIW"`, `"CatDirichlet"` and even `"CatDP"`) 

 See R example for details:

```R
## Sample cluster labels k from HDP-MM using Gibbs sampling

library(bbricks)

## load some mixture of Gaussian samples.
## mmhData is a list of two elements. mmhData$x is a matrix of Gaussian observations, each row is an observation; mmhData$groupLabel is the group label of each observations.
## see ?mmhData for details
data(mmhData)
x <- mmhData$x
js <- mmhData$groupLabel
## Step1: Initialization------------------------------------------
maxit <- 50                             #iterative for maxit times
burnin <- 30                            #number of burn in samples
## create a HDP object to track all the changes, the HDP object in this case is a combination of a CatHDP object and a GaussianNIW object:
obj <- HDP(gamma = list(gamma=1,j=max(js),alpha=1,
                        H0aF="GaussianNIW",
                        parH0=list(m=c(0,0),k=0.001,v=2,S=diag(2)*0.01)))
ss <- sufficientStatistics(obj$H,x=x,foreach = TRUE) #sufficient statistics
set.seed(1)
z <- rep(1L,nrow(x))
k <- matrix(1L,nrow(x),maxit-burnin)    #place-holder for the sampled k
N <- length(ss)
for(i in 1L:N){# initialize k and z
    tmp <- rPosteriorPredictive(obj = obj,n=1,x=x[i,,drop=FALSE],j=js[i])
    z[i] <- tmp["z"]
    k[i,1] <- tmp["k"]
    posterior.HDP(obj = obj,ss = ss[[i]],ss1 = k[i],ss2 = z[i],j = js[i])
}
## Step2: main Gibbs loop---------------------------------------------
it <- 1                                 #iteration tracker
pb <- txtProgressBar(min = 0,max = maxit,style = 3)
while(it<=maxit){
    if(it>burnin) colIdx <- it-burnin
    else colIdx <- 1
    for(i in 1L:N){
        ## remove the sample from the posterior info
        posteriorDiscard(obj = obj,ss = ss[[i]],ss1=k[i,colIdx],ss2=z[i],j=js[i])
        ## resample a new partition
        tmp <- rPosteriorPredictive(obj = obj,n=1,x=x[i,,drop=FALSE],j=js[i])
        z[i] <- tmp["z"]
        k[i,colIdx] <- tmp["k"]
        ## add the information of the new sample
        posterior(obj = obj,ss = ss[[i]], ss1=k[i,colIdx],ss2 = z[i],j=js[i])
    }
    if(it>burnin & colIdx<ncol(k)) k[,colIdx+1] <- k[,colIdx] #copy result of previous iteration
    it <- it+1
    plot(x=x[,1],y=x[,2],col=k[,colIdx])         #to visualize the group label dynamics
    setTxtProgressBar(pb,it)
}
## Step3: Estimate group labels of each observation---------------
## pick the most frequent k as the best estimate
kBest <- apply(k,1,function(l){
    tmp <- table(l)
    names(tmp)[which.max(tmp)]
})
plot(x=x[,1],y=x[,2],col=kBest)
```



### Topic Modeling with HDP

A topic model is a hierarchical mixture model(See [Hierarchical Mixture Models](#hierarchical-mixture-models)) with categorical observations:
$$
\begin{align}
G_j|\gamma & \sim DP(\gamma) \\
\pi |\alpha,G_j & \sim DP(\alpha,G_j) \\
z | \pi & \sim Categorical(\pi) \\
k| z,G_j & \sim Categorical(G_j) \text{, if }z\text{ is a sample from } G_j\\
\pi_k | \psi & \sim Dirichlet(\psi) \\
x|\pi_k & \sim Categorical(\pi_k)
\end{align}
$$




### Hierarchical Topic Modeling with HDP2



### Infinite State Hidden Markov Model



### Bayesian Linear Regression

A Bayesian linear regression model has following graph structure:

![](./notes_pictures/bayesianLinearRegression.png)

The CPDs are:
$$
\begin{align}
\beta,\sigma^2 | \gamma & \sim NIG(\gamma) \\
x| \beta, \sigma^2,X & \sim Gaussian(X\beta,\sigma^2)
\end{align}
$$
Where $NIG(\gamma)$ is the Normal-Inverse-Gamma distribution with parameter $\gamma=(m,V,a,b)$, $m$ and $V$ are the "location" and "scale" parameters, $a$ and $b$ are the "shape" and "rate" parameters.

The distribution of $\gamma \rightarrow (\beta,\sigma^2) \rightarrow x$ is a basic prior-posterior structure as shown in [Mindset](#mindset) graph $(b)$. **bbricks** provides an object type "GaussianNIG" to represent such a structures. See the R example below:

```R
## Bayesian linear regression

library(bbricks)

## lrData is a list of two elements. lrData$x is the sample set of the dependent variable; lrData$X is the sample set of the independent variable
## see ?lrData for details
data(lrData)
X <- lrData$X                           #a matrix of 1 column
x <- lrData$x                           #a numeric vector
obj <- GaussianNIG(gamma=list(m=0,V=1,a=1,b=0)) #create a GaussianNIG object
ss <- sufficientStatistics(obj = obj,X=X,x=x)   #the sufficient statistics of X and x
posterior(obj = obj,ss = ss)                    #add the infomation to the posterior
bsMAP <- MAP(obj)                               #get the MAP estimate of beta and sigma^2
bsMAP                                           #print the MAP estimate
## plot the MAP estimate of the regression line
plot(X,X%*%bsMAP$betaMAP,type = "l")
points(X,x,pch=20)
```



### Hierarchical Bayesian: Estimate Cancer Mortality Rates