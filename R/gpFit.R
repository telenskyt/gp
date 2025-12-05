#' Fit a Gaussian process model
#'
#' @param gp GP model object.
#' @param h optional; a numeric vector of starting values of hyperparameters (untransformed, i.e. on their most natural scale). If \code{NULL},
#'   the starting values are taken from \code{gp} hyperparameter table (\code{gp$hyperpar[,"start"]}).
#' @param opt.h logical; should hyperparameters be optimised? Default TRUE; Only rarely you may want to say FALSE.
#' @param two.stage logical; if TRUE, the two-stage optimization method will be applied (see Details below).
#' @param stages numerical vector - the optimization stages to be performed. Only valid for \code{two.stage = TRUE}. Default is \code{1:2}, i.e. both stages.
#' @param stage1low the low threshold for the sigma2 parameter of the components for which the two-stage fitting applies.
#' @param use.prior logical; if TRUE, likelihood + prior (posterior) is optimized; if FALSE, only likelihood is optimized.
#' @param hessian logical; compute and return the Hessian during the optimization. Warning: may be heavy on CPU time! Default FALSE.
#' @param opt.control list of control parameters passed to the \code{optim(control = )}
#' @param verbose logical
#' @param method character; fitting method. Currently \code{"Laplace"} supported (Laplace approximation).
#' @param use_f_start logical; internal optimization - start each Laplace approximation from the optimized \code{f} vector from the previous hyperparameter iteration.
#' @param weights numeric or \code{NULL}; observation weights (not implemented).
#' @param grad.computation logical; compute gradients during the gaussian process optimization (required when
#'   \code{opt.h = TRUE}).
#' @param recursive internal option; do not use.
#'
#' @return Returns the input \code{gp} object with the model fit (stored in \code{gp$fit}) and updated hyperparameters.
#'
#' *NOTE!!* The returned object is usually large, because it contains pre-calculated matrices that are used for predictions. If you need to save the model,
#' call \code{gpPack()} first to remove these; when you need to use the model again, call \code{gpUnpack()} to restore these matrices.
#'
#' @details
#' The two-stage fitting method (\code{two.stage = TRUE}) is a special optimization procedure (Telenský et al 2026) where the optimization is split into two stages:
#'
#' 1. In the first stage, sigma2 parameters are not allowed to go too close to zero (as per the \code{stage1low} parameter). This way, the other hyper-parameters in the same component have time
#'    to optimize themselves before the sigma2 parameter could shut the whole component down to zero. (This is only applied to components
#'    that have also other hyper-parameters apart from sigma2.)
#'
#' 2. In the second stage, the restriction applied in the first stage is lifted and all the sigma2 hyper-parameters are allowed to go further down to zero.
#'
#' In some cases, where user tries to optimize the CPU time by manually splitting the procedure in multiple steps, it might be a good idea to run the stage 2 only in the last step,
#' and use only stage 1 in all the previous steps (see Vignette XXX).
#'
#' @examples
#' \dontrun{
#' gp <- gp(...) # construct gp object
#' gp <- gpFit(gp)
#' }
#'
#' @export


# prior - mean function (obsolete now)

# l - vector of hyper-parameters (on the most natural hyper-parameter scale). Called "l" because originally it were squared length-scales
#		(and some of them still are :-)). Squared length-scales = square of the l defined by Rasmussen & Williams 2006, which are length-scales at the scale of the covariates
#		note that in the paper Golding & Purse 2016 is an opposite thing, the l is defined as a sqrt of the Rasmussen & Williams 2006 length-scale (super confusion there :))
#
# theta.link - function to convert hyper-parameters to the theta - to the optimization scale used by optim() function
#
# grad.computation - computation of gradients (quite time consuming) - can be set to FALSE for FIXED models (opt.l = FALSE)
#
# opt.control - option to pass to optim(control = )

gpFit <- function (gp, h = NULL, opt.h = TRUE,
				  two.stage = TRUE,
				  stages = 1:2,
				  stage1low = 1e-2,
				  use.prior = TRUE,
				  hessian = FALSE,
#				  optim.low = -Inf, optim.up = Inf,
				  opt.control = list(factr = 1e10, maxit = 190),
				  verbose = FALSE, method = c('Laplace'),
				  use_f_start = TRUE,
				  weights = NULL,
				  grad.computation = TRUE,
				  recursive = FALSE #internal option only, don't use!
				  )
{
	method <- match.arg(method)

	if (!recursive) { # top level call, from the outside (not an internal call) - will be called just once in the beginning
		if (grad.computation == FALSE && opt.h == TRUE)
			stop("grad.computation must be TRUE when hyperparameters are optimized")

		# get all visible object as a list
		# (dela to moc slozite, stacilo by c(as.list(environment()), list(...))  viz https://stackoverflow.com/a/17244041
		# ale necham to bejt uz...)
		args <- capture.all()

		# get the expected objects
		expected_args <- names(formals(gpFit))

		# remove any unexpected arguments
		args <- args[names(args) %in% expected_args]

		cat("\n\n=====================\nABSOLUTE START\n\n")
		args2print <- args
		args2print$gp <- NULL
		args2print$h <- NULL
		args2print$recursive <- NULL # this is internal, should be hidden
		cat(sub("^list\\(", "gpFit(", deparse(args2print)))
		cat("\n(gp, h args deleted)\n\n")
		args2print <- NULL
		mstart(TRUE, id = "gpFit", mem_precise = TRUE) # start measuring time & memory
	}

	# optionally, do hyperparameter optimization (by recursively calling this function)
	if (opt.h) { # will be only set in the top level call; not in the optimise()'d ones
		stopifnot(!recursive) # so again, this will be executed just in the very beginning

		# initiate the starting value of h, if needed
		if (is.null(h)) {
			args$h <- gpHyperparExportVector(gp, "start")
		}

		if (!two.stage) {
			# pass this to optimiser
			fit <- optimise.gp(args)
			gp$hyperpar <- gpHyperparImportVector(gp, fit$h) # don't forget to import the resultant hyperparameter vector to the table
		} else {
			# two-stage fitting method!

			# determine the components to be staged
			hy <- gp$hyperpar
			staged_components <- hy %>%
				summarize(.by = component, nhyp = n_distinct(hyperpar), with_sigma2 = any(hyperpar == "sigma2")) %>%
				filter(nhyp > 1 & with_sigma2) %>%
				pull(component)

			cat("Two-stage fitting method\n")
			# fit stage 1
			if (1 %in% stages) {
				warning("testing warning, remove this one!!!")
				cat("Fitting stage 1\n")
				# set limits for STAGE 1
				old_low_limits <- gp$hyperpar$low[gp$hyperpar$hyperpar == "sigma2" & gp$hyperpar$component %in% staged_components]
				if (any(abs(old_low_limits - stage1low) < sqrt(.Machine$double.eps)))
					warning("gpFit(): some sigma2 hyperparameters seem already restricted to stage 1")
				if (any(old_low_limits > stage1low + sqrt(.Machine$double.eps)))
					stop("gpFit(): low limits for some sigma2 hyperparameters have higher value than the limit for stage 1 (which is stage1low = ", stage1low, ")")

				args$gp <- NULL # remove this duplicity now, to avoid any unnecessary copy-on-write
				gp$hyperpar$low[gp$hyperpar$hyperpar == "sigma2" & gp$hyperpar$component %in% staged_components] <- stage1low
				gpHyperparCheckAll(gp)

				# pass this to optimiser
				args$gp <- gp
				fit <- optimise.gp(args)
				args$gp <- NULL # remove this duplicity now, to avoid any unnecessary copy-on-write

				gp$hyperpar <- gpHyperparImportVector(gp, fit$h) # don't forget to import the resultant hyperparameter vector to the table

				# recover the previous limits
				gp$hyperpar$low[gp$hyperpar$hyperpar == "sigma2" & gp$hyperpar$component %in% staged_components] <- old_low_limits
			}

			if (2 %in% stages) {
				cat("Fitting stage 2\n")

				if (any(gp$hyperpar$low[gp$hyperpar$hyperpar == "sigma2" & gp$hyperpar$component %in% staged_components] >= stage1low - sqrt(.Machine$double.eps)))
					warning("gpFit(): trying to proceed with stage 2, but low limits for some sigma2 components are at the level of stage 1 (which is stage1low = ", stage1low, ") or higher!")

				# pass this to optimiser
				args$gp <- gp
				fit <- optimise.gp(args)
				args$gp <- NULL # remove this duplicity now, to avoid any unnecessary copy-on-write

				gp$hyperpar <- gpHyperparImportVector(gp, fit$h) # don't forget to import the resultant hyperparameter vector to the table
			}
		}
		# right before return, when opt.h = TRUE finishes
		gp_fit_res <- mstop(TRUE, id = "gpFit")
		fit$opt.h <- TRUE
		fit$time = gp_fit_res$timeSec
		fit$max_memory_usage_mb = gp_fit_res$maxMemMB
		fit$call <- match.call() # bacha duplicitni kod je nize!
		fit$args <- args
		fit$args$gp <- NULL # save space, at to tam neni 2x
		#fit$my_args <- c(as.list(environment())) # c(as.list(environment()), list(...)) # simple. https://stackoverflow.com/a/17244041
		# skip out of this top-level function call and return the final result
		gp$fit <- fit
		return (gp)
	} else if (!recursive) { # direct call, no optim
		# initiate the starting value of h, if needed
		if (is.null(h)) {
			h <- gpHyperparExportVector(gp, "start")
		}
		if (use_f_start) {
		#	warning("use_f_start works only for opt.h = TRUE")
			use_f_start <- FALSE
		}
	}

	if (!is.null(weights)) {
		stop("Weights are not implemented yet.")
	} else {
		weights = 1
	}

	stopifnot(!is.null(h))
	gpHyperparCheckAll(gp)
	gpHyperparCheck(gp, h)

	# fit model
	if (method == 'Laplace') {
		# by Laplace approximation
		fit <- gpFitLaplace(gp, h = h, wt = weights, e = NULL, verbose = verbose, use_f_start = use_f_start,
			grad.computation = grad.computation)
	} else {
		stop("Not implemented")
	}

	fit$hessian <- hessian
	#XXYY
	#if (is.function(model_report_AUC))
	#	fit$AUC <- model_report_AUC(fit$f, fit$y, hyperpar_decode(l))
	if (!recursive) { # case when opt.h = FALSE, called once at the top level before finishing.
		gp_fit_res <- mstop(TRUE, id = "gpFit")
		fit$opt.h <- FALSE
		fit$time = gp_fit_res$timeSec
		fit$max_memory_usage_mb = gp_fit_res$maxMemMB
		fit$call <- match.call() # bacha duplicitni kod je vyse!
		fit$args <- args
		fit$args$gp <- NULL # save space, at to tam neni 2x
		#fit$my_args <- c(as.list(environment())) # c(as.list(environment()), list(...)) # simple. https://stackoverflow.com/a/17244041
		gp$fit <- fit
		gp$hyperpar <- gpHyperparImportVector(gp, fit$h) # don't forget to import the resultant hyperparameter vector to the table
		return (gp)
	}
	return(fit)
}

capture.all <- function() {
	# capture all visible objects in the parent environment and pass to a list
	env <- parent.frame()
	object_names <- objects(env)
	objects <- lapply(object_names,
					  get,
					  envir = env)
	names(objects) <- object_names
	return (objects)
}
