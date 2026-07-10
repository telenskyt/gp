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

install_github("https://github.com/telenskyt/gp/", force = TRUE, INSTALL_opts = c("--with-keep.source"))
```

You can also clone the github repository like this:

```
git clone https://github.com/telenskyt/gp/
```

## Install Intel OneMKL Math Kernel Library

The matrix operations will be sped up significantly using this library. It is using special CPU instructions designed for matrix operations. 
It works not only on Intel CPUs, but also on AMD and others, since they have similar instruction sets.

1. Install OneMKL from <https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html>

2. Benchmark your R installation before integration into R: download and run [R-benchmark-25.R](https://mac.r-project.org/benchmarks/R-benchmark-25.R) from <https://mac.r-project.org/benchmarks/>

3. [Install OneMKL into R](https://www.intel.com/content/www/us/en/developer/articles/technical/using-onemkl-with-r.html) and perhaps [learn more about it](https://www.intel.com/content/www/us/en/developer/articles/technical/accelerate-r-code-with-math-kernel-library.html)

4. Repeat the benchmark from point 2; you should see significant speed improvement!

## Start using the package

See the vignette [A simple example](articles/a-simple-example.html) for introduction to the package. 
