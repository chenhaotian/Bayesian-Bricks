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
