# bbricks 0.1.4 (2020-05-07)

## New features

Add method `dAllIndicators` to all DP family distributions. `dAllIndicators` will return the probabilities of all the possible values of the hidden indicator variables.

## Bug fixes

+ Fixed bug: DP, HDP and HDP2 objects don't initialize observation distribution when created.
+ Fixed bug: `rPosteriorPredictive.HDP` and `rPosteriorPredictive.HDP2` return samples with unexpected name such as "z.z".
+ Fixed bug: fix the error `Error in obj$Z$Z2[[j]] : subscript out of bounds` when `rPosteriorPredictive.HDP` try to draw samples from an unrecorded group. 

## Improvements

+ Fix graphical model representation error in package vignette.
+ Fix wrong conditional probability distributions in DP related documentations.
+ Add HMM and HDP-HMM (iHMM) examples to the vignette

# bbricks 0.1.3 (2020-04-26)

Enhance documentations.

## Improvements
+ fix some display errors in function documentations.
+ add example "hierarchical Bayesian linear regression" to the package vignette.

# bbricks 0.1.2 (2020-04-14)

Add features and enhance documentations. Fix minor bugs.

## New features

**New Model Structures:**
+ Linear Gaussian system.
+ Gaussian and Gaussian conjugate structure.
+ Gaussian and Inverse-Wishart conjugate structure.

**New Inference Tasks:**
+ Sample from the posterior distribution.
+ Calculate density values from the posterior distribution.

**New Distributions:**
+ Inverse Gamma distribution.
+ Wishart and Inverse-Wishart distributions.

## Bug fixes

Fixed bug: rGaussian returns "non-conformable arrays" error when sample mean "mu" is a 1 column matrix.

## Improvements
Reformat all the function documentations, now they are more readable.

# bbricks 0.1.1 (2020-03-07)

Fix bug.

## Bug fixes
Fixed bug: posterior.CatHDP, posteriorDiscard.CatHDP, posterior.CatHDP2,posteriorDiscard.CatHDP2 update posterior improperly when sample weight "w" is not NULL. 


# bbricks 0.1.0 (2020-02-28)
This is the first release of bbricks.

## New features
**New Model Structures:**
+ Gaussian and NIW conjugate. 
+ Gaussian and NIG conjugate.
+ Categorical and Dirichlet conjugate.
+ Categorical and Dirichlet process.
+ Categorical and 
+ Dirichlet process.
+ Hierarchical Dirichlet process.
+ Hierarchical Dirichlet process with two layers of hierarchies.

**New Inference Tasks:**
+ Update posterior info.
+ Calculate sample sufficient statistics.
+ Calculate MAP estimate.
+ Calculate marginal likelihood.
+ Calculate posterior predictive density.
+ Sample from posterior predictive distribution.

## Improvements
None

## Bug fixes
None
