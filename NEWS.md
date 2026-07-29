# version 0.2.0 (Jul 29, 2026)

## Major new features

- Introduced a completely new class-based likelihood template interface with `likTempPhased`. This, for the first time, allows automatic integration of user-defined likelihoods seamlessly into fitting, prediction, cross-validation, evaluation of predictive performance, and simulation. A likelihood can now define separate terms, link, processing and likelihood stages. Built-in Bernoulli and Poisson likelihoods are included.

- Introduced the new `"Laplace-Fisher"` fitting method for likelihoods whose Hessian is not positive semidefinite. This overcomes the limitation of the widespread Laplace Approximation. Fisher information matrix is calculated automatically for the supported likelihoods, or can estimated by sampling for any supplied likelihood template.

- Added approximations for leave-one-out and leave-group-out cross-validation with `gpLGO()` and `gpLGOpred()`. It calculates LGO-CV predictions without the need to fit the model N times for every fold. It supports general non-diagonal hessian matrix `W`, produces predictions and likelihood-specific performance statistics, and stores the results in `gp$fit$LGOCV`. An optional dense block implementation can accelerate diagonal-`W` models with many folds. Unlike `gpFitCV()`, this approximation operates at a fixed set of hyperparameters.

- Generalized the Laplace and Laplace-Fisher approximations in fitting beyond diagonal likelihood Hessians. Now even non-diagonal Hessian matrices are supported. The higher-order  derivative calculations and likelihood-hyperparameter gradients, prediction support, packing, and LGO-CV calculations were generalized accordingly.

- Added an RTMB-based linear-model family through `lmFit()`. It replaces the latent GP with a formula-based linear predictor while retaining the GP model's likelihood template. Likelihood hyperparameters can be fixed or optimized. The resulting `lmFit` objects support prediction, coefficient summaries, standard errors, `logLik()`, `AIC()`, `extractAIC()`, `nobs()`, and model-selection workflows such as `step()`. `lmFitCV()` provides fold-based cross-validation for these models.

## Smaller new features

### Prediction and cross-validation

- Reworked prediction around the likelihood-template pipeline:
	- `type = "latent"` returns latent GP predictions.
	- `type = "terms"` adds likelihood-template terms.
	- `type = "response"` applies the likelihood link stage.
	- GP, linear-model, and LGO predictions now share the same latent-prediction post-processing through `gpPredictFromLatent()`.

- Replaced the old `CI` argument to `predict.gp()` with `conf.int` and `conf.level`. Enabling confidence intervals now automatically enables standard-error calculation.

- `gpFitCV()` now:
	- Fits and predicts a null model independently within every fold.
	- Automatically evaluates predictions through the likelihood template.
	- Accepts `pred.options`, `lmFit.options`, and `lmFit.pred.options`.
	- Stores both GP and null-model cross-validation predictions.
	- Reports likelihood-specific predictive-performance statistics.

- Added parallel processing of chunked predictions with `predict.gp(parallel = TRUE)`. Parallel prediction jobs support log files, diagnostic dumps, shortened tracebacks, and appendable per-worker logs.

- Added Bernoulli predictive-performance evaluation through `predPerfBern()`, including AUC, TSS, deviance-based and likelihood-ratio pseudo-\(R^2\), NLL, and null-model NLL.

- Added Poisson predictive-performance evaluation through `predPerfPois()`, including ordinary \(R^2\), deviance-based and likelihood-ratio pseudo-\(R^2\), NLL, and null-model NLL.

- Likelihood-template `terms()` callbacks may now return either one shared row or one row for every latent value, with explicit result validation.

### Fitting and diagnostics

- Added `convergence()` for a concise assessment of optimizer and Laplace-approximation convergence.

- Added `coef.gp()` for inspecting fitted hyperparameters together with marginal-likelihood, prior, and posterior gradients.

- Added `iter()` for examining recorded hyperparameter and Laplace iterations, including support for two-stage fits.

- Added `opt.method` to `gpFit()`, allowing the optimization method passed to `optim()` to be selected.

- Added `save.every.iter` to `gpFit()`. This saves a packed checkpoint after every hyperparameter iteration and once more when fitting finishes.

- Added `gpHyperparReadLog()` for recovering the last fitted hyperparameter vector from a fitting log.

- Hyperparameter optimization can now run without analytic gradients by setting `grad.computation = FALSE`. Gradient computation now defaults to the value of `opt.h`.

- Added likelihood-aware `f_start` handling. Likelihood templates may define an `f_start` transformation that adjusts the previous latent solution when likelihood hyperparameters change.

- Added `gpLik()`, a convenience function for running selected stages of a fitted model's likelihood template.

- Fitted models now record package/build provenance, including the package version, Git commit, dirty flag, build time, and R version when available.

- Model-run metadata now includes timing, maximum memory usage, hostname, start and end times, and BLAS/LAPACK information.

### Data handling

- `gpDataPrepare()` now uses an explicit `scale.as` argument when another dataset's scaling should be reused. This prevents accidental reuse of inappropriate training-data scaling.

- Zero-variance columns are detected during scaling, reported, and represented by zeros instead of resulting `NaN` values.

- Added `.onSubset` to `gpData()`, allowing a callback to update derived data after every `gpDataSubset()` operation.

- Added `scale_cols()` for scaling selected columns only.

### Performance and storage

- Substantially optimized calculations for non-diagonal likelihood Hessians, reducing dense intermediate allocations and avoiding several unnecessary sparse-to-dense conversions.

- Reworked construction of \(B = W^{1/2}KW^{1/2} + I\) to reduce peak memory consumption.

- Optimized common marginal-likelihood gradient calculations and reused cached matrix products more extensively.

- Posterior covariance masking is now performed before expensive dense/sparse products. In the motivating case, this reduced one operation from approximately nine seconds to around 0.1 seconds.

- `gpPack()` now removes additional reconstructible matrices, the instantiated likelihood class, and stored LGO predictions. The likelihood object is reconstructed from its saved constructor arguments by `gpUnpack()`.

- Cross-validation fold models are packed more aggressively. `lmFitCV()` no longer stores all fold models, which can otherwise become extremely large for LOO/LGO-style runs.

- Likelihood-template constructor defaults and stored constructor arguments were reorganized to avoid capturing large environments. Large user-function environments now produce a warning.

- Reused `LTinv.rW` instead of repeatedly reconstructing or solving with the Cholesky factor.

- Added `f_cov_masked`, a sparse posterior covariance matrix restricted to the structural nonzero pattern of `W`.

## Fixes

- Fixed a critical likelihood-hyperparameter gradient bug ([#1](https://github.com/telenskyt/gp/issues/1)): `d0_dhyp()` constructed its automatic-differentiation tape at `f = 0` and also evaluated it there instead of at the fitted latent values. This could produce incorrect fitted results whenever likelihood hyperparameters were optimized.

- Hardened `f_start` against several serious numerical failures:
	- A previous solution is rejected when its likelihood under the current hyperparameters is non-finite or drastically worse than the ordinary starting value.
	- Failed Cholesky decomposition while using `f_start` now resets the starting value and restarts the Laplace optimization ([#2](https://github.com/telenskyt/gp/issues/2)).
	- Non-finite Laplace objectives trigger the same fallback.
	- Missing likelihood-specific `f_start` methods are reported before optimization begins.

- Reimplemented the likelihood-hyperparameter derivative calculation in `d1_dhyp()` to prevent extreme memory consumption with non-diagonal likelihood Hessians.

- Fixed `predict.gp(cov.fit = TRUE)`, which previously failed immediately because it tested the wrong variable.

- Fixed incorrect row reindexing in `predict.gp()` for models whose GP dimension differs from the main-table dimension.

- Fixed likelihood-template prediction stages receiving unscaled rather than correctly prepared prediction data.

- Fixed `predict.lmFit(type = "terms")` failing to include fixed likelihood hyperparameters.

- Fixed prediction and sparse-matrix operations when matrices do not have dimension names.

- Fixed cross-validation training data potentially inheriting scaling from the full dataset rather than being scaled independently.

- Fixed `gpUnpack()` for maximally packed cross-validation models by refitting the final iteration when the latent solution is absent.

- Extended `gpUnpack()` to reconstruct diagonal and non-diagonal fitting matrices consistently, including Laplace-Fisher models and the numerical correction of small negative diagonal `W` values.

- Fixed an incorrect `rW` reference during unpacking.

- Fixed a typo that disabled validation that likelihood hyperparameters appear at the end of the hyperparameter table.

- Fixed `gpHyperparList()` reordering hyperparameters alphabetically within components. Original table order is now preserved.

- Improved `gpHyperparCheck()` so out-of-range hyperparameters are identified by component, parameter, index, variable, limits, and supplied value.

- Fixed accidental partial matching of `fitCV` as `fit`.

- Fixed `gpFitCV()` warning-level handling so it does not lower an explicitly stricter `options(warn)` setting.

- Fixed parallel diagnostic dumps potentially stalling indefinitely while deparsing enormous `do.call()` expressions. Dumped call labels are now safely bounded.

- Rewrote several `do.call()` sites so `match.call()` does not expand large argument objects into stored calls.

- Ensured RTMB is attached when `gp` is loaded, preventing automatic-differentiation failures in functions such as `plogis()`.

- Fixed likelihood-template function handling that could remove required environments and cause RTMB "non-numeric argument" failures. Large environments are now warned about instead.

- Fixed several sparse-matrix operations that produced unnecessary gigabyte-scale dense allocations.

## Minor changes and improvements

- Renamed the fitted quantity `mnll` to the mathematically standard `nlml`—negative log marginal likelihood.

- Added targeted warning suppression with `suppress_warnings_from()` and utilities for inspecting the warning call stack.

- Added `log.append` to `parallelJobWrapper()`.

- Error logs now use concise condition messages instead of printing enormous condition objects.

- `gpFit()` now reports the GP matrix dimension at the beginning of fitting.

- Hyperparameter priors are printed more clearly in tibbles through a custom pillar method.

- Removed an unnecessary prediction error for extra columns in `newdata`; required training columns are still checked.

- Added checks for invalid likelihood-hyperparameter indices, non-finite coefficient tables, incompatible likelihood settings, missing fitted data, and malformed template outputs.

- Stored LGO results now include the package version and the arguments used to create them, including the fold definition.

- Added instructions for accelerating matrix operations with Intel OneMKL.

- Expanded manual pages and updated the pkgdown website, including a dedicated News section.

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
