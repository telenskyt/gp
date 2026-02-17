# version 0.1.2 (Feb 17, 2026)

- added support for likelihood hyperparameters! With full AD (automatic differentiation) support from RTMB!

- solved the numerical problems in gpFitLaplace() that appeared with the nest survival likelihood in Lapwings
	- solved numerical problems when optimise() couldn't find optimum step size because of infinite psi(). Introduced a remedy that will try to find a region of step_size that would have finite psi(). This solves the problem at the moment (problem documented in d:\tomas\cejky\cejky_clean\errors\02-numericka-obj_Inf,optimise\) (c1a8148)
	- added parameter num.correct.W.tol  and introduced small numeric corrections of negative W values, along with some new warning() (19186dd)
	R/gpFitLaplace.R

- added NLL and null model NLL to the CV stats (bfe895f)

- new mechanism for predict(type = "response") - possible to define predictor.fun parameter to gp(), which can be more complex than just a simple link function. This will be usually needed for likelihoods with extra hyperparameters, like e.g. the nest survival in Lapwing (a268168)

- added new method gpGetCVModel() to get one of the GP models used for cross-validation


## Other changes
- added new parameter tr.max.lines for gpFitCV() and parallelJobWrapper()
- added `...` argument to gpFit() for passing arguments to fitting methods (now, just gpFitLaplace) (98bd8e2)
- d0() shouldn't check for infinite values - there is already code in gpFitLaplace to take care of that, which allows fallback from use_f_start (6f139cf)
R/d0.R
- passing the parameter `parallel` to parallelJobWrapper(), so now in case of parallel = FALSE it will behave differently. (I was hoping to be able to achieve reasonable debugging of gpFitCV(), but I failed at the moment, because foreach()'s eval breaks normal options(error = recover) style debugging) (7a9b0ca)
R/gpFitCV.R

- by default, stay interactive when parallel == FALSE (6333865)
- changed options(warn = 2) to options(warn = 1)! warn = 2 is too crazy.
- added some help pages
- allow matrices in the input data - important for model.matrix() in the data (5a250c0)
R/gpData.R

# version 0.1.1 (Feb 10, 2026)

- added cross-validation gpFitCV()
- added cross-validation evaluation - cv_eval_bern()
	- with likelihood ratio statistics - made p.null a mandatory parameter (doesn't make sense to try to derive it in this function, as I did in the Atlas!)
	- doc with references
	- added a cross-validation example to the simple_example.R
- completed "A simple example" vignette
- significant updates to the pkgdown website (https://telenskyt.github.io/gp)
- added a new covariance function cov.I.factor.sigma2()
- gp() interface change: user should now specify NEGATIVE log likelihood instead of log likelihood, to comply with the RTMB/TMB convention (6f07453)

## Fixes

- gpDetermineSize(): handle a special case when GP_factor is not a factor by which any table would be indexed (e.g. for formula f <- ~ i:1 + habitat:(I|habitat) in cejka). Anyway, this change wasn't enough to make this special case work. (c33625b)
- fixes in gpFitCV(): 
	- set a clear criteria of how fold.fact must relate to gp$GP_factor - introduced reindexing of the fold.fact to gp$GP_factor when needed! (ac35bec)
	- fixed a problem in start.from.model (model list is not a named list, so it needs different approach) (5ab6565)
- parallelJobWrapper() now reports traceback()! Finally figured how to do it!!!
	- Changed the default for log.fn and dump.fn to NULL in parallelJobWrapper() and foreach2(). (4a5fb70)
- added some help pages 
- and many more


# version 0.1.0 (Dec 18, 2025)

First release.
