# gp: Gaussian Process in R

This R package offers an engine for fitting Gaussian Process (GP) models, with the following features:

1. The Gaussian process can act as a latent variable in the model, as a part of arbitrary user-defined likelihood.

2. User defines the Gaussian process using a formula, as a sum of covariance functions. The GP thus can have multiple components and make predictions, including standard errors, for the whole model, 
or for each component separately. 

3. User defines the likelihood the same way as in the RTMB package, which is used for automatic differentiation (not for model fitting though). 

4. The package has its own fitting engine, which uses Laplace approximation.


## Install instructions

Install the package directly from github:

```r
library(devtools)

install_github("https://github.com/telenskyt/gp/")
```

You can also clone the github repository like this:

```
git clone https://github.com/telenskyt/gp/
```

## Start using the package

See the vignette [A simple example](articles/a-simple-example.html) for introduction to the package. 
