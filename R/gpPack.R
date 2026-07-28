
#' Pack the model - reduce its size for memory or disk storage
#'
#' Reduce model size by deleting auxiliary data structures which speed up the computations 
#' but take a lot of space.
#'
#' @param maximum should it be taken to maximum level? Default FALSE.
#' @export
gpPack <- function(gp, maximum = FALSE)
{
	# have to use gp[["fit"]] instead of gp$fit due to the damn partial matching!! Because it matched $fitCV instead!!! Sucks as fuck!!!
	gp[["fit"]] <- gpPackFit(gp[["fit"]], maximum = maximum)
	gp$data <- NULL
	gp$lik <- NULL
		# gp$lik class is no longer permanently stored in the packed object! Only the user supplied arguments to construct it! For two reasons:
		#	- wrong design. Thanks to the unfortunate RC classes design in R, the gp$lik object instance contains class code, 
		#	  which might be subject to updates in new package versions. We don't want that. We want only user-supplied code (parts of the template)
		#	  to be stored in the object. The gp$lik class is now deleted in gpPack() and restored in gpUnpack().
		#	- the gp$lik object storing class code was too big. Not so much compressed, but still.	
	gp
}

# internal function
#
# maximum - should it be the max level of packing
# pack just the gp$fit object
gpPackFit <- function(fit, maximum = FALSE)
{
	fit$LTinv.rW <- NULL
	fit$L <- NULL
	fit$rW <- NULL
	fit$K <- NULL
	fit$W <- NULL
	if (!is.null(fit[["LGOCV"]])) {	
		fit$LGOCV$pred <- NULL
	}
	if (maximum) {
		fit$f <- NULL # nakonec i bez tohoto se da obejit!!! :-))) # ale defaultne necham tam radsi uplne to puvodni, kvuli ruzne num nestabilite....
		fit$f_cov_masked <- NULL
	} else
		rownames(fit$f) <- NULL # rownames take a lot of space! Just one column matrix with them 279872 bytes, without 31440 bytes!!!
							# And they only had 6 characters each! (nrow = 3881, ncol = 1). R is not very efficient storing them.
							# And rownames are in obsx already
	fit
}

#' Unpack the model - prepare it for computations (mostly predictions)
#'
#' Restore the auxiliary data structures, which are redundant and take significant space, but are needed for computations.
#' Note that restoring these structures can take some CPU time.
#'
#' @param compute logical - should we prepare for computations (mostly predictions)? This is the CPU-heavy part. \code{compute = FALSE} will skip this.
#'
#' @param need.K,need.L,need.LTinv.rW logical - in case of \code{compute = TRUE}, specify whether particular matrix (K, L) is needed. If only one of these matrices (K, L, need.LTinv.rW) is needed, it can 
#' significantly save memory and sometimes also CPU time. But the object depends on each other (K -> L -> LTinv.rW),
#' so if you set this to \code{FALSE} for an object on which another needed object depends on,
#' you will only save the memory size of the resultant object, but not CPU time.
#'
#' Note that the status of the \code{rW} object will be exactly the same as the \code{L} object (so it is directed by \code{need.L} argument).
#'
#' @param num.correct.W.tol numerical correction of W. Applies only to the case of diagonal W (gp$W.type == "diag"). 
#' If not NULL, values of W will be corrected for small negative numbers, within given tolerance.
#'
#' @export
gpUnpack <- function (gp, compute = TRUE, need.K = compute, need.L = FALSE, need.LTinv.rW = compute, num.correct.W.tol = 10*sqrt(.Machine$double.eps))
{
	gp$data <- gpDataPrepare(gp, gp$obsdata)
	if (!is.null(gp[["lik.constr.args"]])) {
		gp$lik <- do.call(getRefClass(gp$lik.class)$new, gp$lik.constr.args)
	} else if (is.null(gp[["lik"]])) {
		stop("Cannot reconstruct gp$lik: missing likelihood constructor arguments.")
	}

	if (is.null(gp[["fit"]]))
		return(gp)
		
	if (!compute) {
		stopifnot(!need.K)
		stopifnot(!need.L)
		stopifnot(!need.LTinv.rW)
		return(gp)
	}

	if (is.null(gp$fit[["f"]])) {
		cat("gpUnpack: $fit$f is missing: need to re-run the gpFitLaplace() for the last iteration\n")
		
		gf <- gpFit(gp, h = gpHyperparExportVector(gp, "value"), opt.h = FALSE, grad.computation = FALSE, num.correct.W.tol = num.correct.W.tol) 
			# !! there should be an option to avoid nlml computation, too, but perhaps that wouldn't help us anyway, since we need the L too...
		stopifnot(!is.null(gf[["fit"]]))
		stopifnot(!is.null(gf$fit[["f"]]))
		stopifnot(!is.null(gf$fit[["a"]]))
		stopifnot(!is.null(gf$fit[["K"]]))
		stopifnot(!is.null(gf$fit[["W"]]))
		gp$fit$f <- gf$fit$f
		gp$fit$f_cov_masked <- gf$fit$f_cov_masked
		gp$fit$a <- gf$fit$a
		gp$fit$K <- gf$fit$K
		gp$fit$W <- gf$fit$W
		if (!is.null(gf$fit[["L"]])) {
			gp$fit$L <- gf$fit$L
			gp$fit$rW <- gf$fit$rW
		}
		if (!is.null(gf$fit[["LTinv.rW"]]))
			gp$fit$LTinv.rW <- gf$fit$LTinv.rW
		return(gp)
	}
	
	mstart(id = "K")
cat("Re-creating covariance matrix... ")
	gp$fit$K <- K_matrix(gp, comp_means = TRUE)
	cat(sprintf("(%4dx%-4d): \t", nrow(gp$fit$K), ncol(gp$fit$K)))
	mstop(id = "K")
#mstop()
#mstart()
	n <- nrow(gp$fit$K)
	
	hyperpar <- gpHyperparList(gp)

	mn <- mnfun(gp, hyperpar = hyperpar)
	#gp$fit$f <- gp$fit$K %*% gp$fit$a + mn
	
	fisher <- NULL
	if (gp$fit$method == "Laplace-Fisher")
		fisher <- gp$fit$fisher.options

	W <- -d2(gp, gp$fit$f, gp$data, hyperpar, fisher)

	if (gp$W.type == "diag" && gp$fit$method != "Laplace-Fisher" && any(W < 0) && !is.null(num.correct.W.tol) && is.finite(num.correct.W.tol)) { # 2026 osetreni - jen drobne numericke chybky muzem zarovnat na 0, s ohledem na komentar nize
		W[W < 0 & W >= -num.correct.W.tol] <- 0
		warning("corrected small negative values of W (within tolerance num.correct.W.tol)")
	}
		
	if (gp$W.type == "diag" && gp$fit$method != "Laplace-Fisher" && any(W < 0)) {
		# (2026 note: going through the log files of my models, it looks like this didn't actually work. 
		# It happened in q_extend models, which didn't work out, only once it happened in normal model (O-tsc,sitesm/Picus_viridis), resulting in an error and core dump anyways.)
		stop("some W < 0, i.e. the diagonal of the hessian has negative values (range: ", paste(range(W), collapse=" to "),
			"),\n\tmeaning the hessian of the neg. log likelihood is not positive definite. Consider using method = \"Laplace-Fisher\". Increasing num.correct.W.tol can help for small negative values, but use with caution!")
		
	}
	gp$fit$W <- W

	if (!need.LTinv.rW && !need.L) {
		if (!need.K) {
			empty_K <- NA
			mostattributes(empty_K) <- attributes(gp$fit$K)
			gp$fit$K <- empty_K
		}
		return(gp)
	}

	L <- gp$fit[["L"]]
	rW <- gp$fit[["rW"]]
	# need.LTinv.rW || need.L
	if (is.null(gp$fit[["L"]]) && (need.L || is.null(gp$fit[["LTinv.rW"]]))) {
		# need to recompute L
		mstart(id = "L")
		if (gp$W.type == "diag") {
			rW <- sqrt(gp$fit$W)
			xx <- rW %*% t(rW)# * K + diag(n)
			#xx <- rW %*% t(rW) * K + diag(n)
		gc()
			xx <- xx * gp$fit$K
		gc() # ta vlozena gc() tady jsou velmi dulezita!!!
			diag(xx) <- diag(xx) + 1
		gc()
		} else {
		gc()
			rW <- chol_W(gp$fit$W)
			xx <- rW %*% gp$fit$K %*% t(rW)
		gc() # ta vlozena gc() tady jsou velmi dulezita!!! (pro nizsi konzumaci pameti)
			diag(xx) <- diag(xx) + 1
		gc()
		}
		if (!need.K) {
			empty_K <- NA
			mostattributes(empty_K) <- attributes(gp$fit$K)
			gp$fit$K <- empty_K
		}
		gc() # ta vlozena gc() tady jsou velmi dulezita!!!		
		L <- chol(xx) # see the note on L above (it's a transpose of L in R&W2006, but it's OK)
		rm(xx)
		mstop(id = "L")
		gc()
	} else {
		if (!need.K) {
			empty_K <- NA
			mostattributes(empty_K) <- attributes(gp$fit$K)
			gp$fit$K <- empty_K
		}
		gc() # ta vlozena gc() tady jsou velmi dulezita!!!		
	}
	if (need.LTinv.rW && is.null(gp$fit[["LTinv.rW"]])) {
		mstart(id = "LTinv.rW")
		stopifnot(!is.null(rW))
		stopifnot(!is.null(L))
		if (gp$W.type == "diag")
			LTinv.rW <- backsolve(L, diag(rW), transpose = TRUE) # L^T^-1 chol(W)
		else {
			rW.dense <- suppress_warnings_from(as.matrix(rW), "sparse->dense coercion: allocating vector of size", fun = "asMethod")
			LTinv.rW <- backsolve(L, rW.dense, transpose = TRUE) # L^T^-1 chol(W)
		}
		mstop(id = "LTinv.rW")
	}
	if (need.L) {
		gp$fit$L <- L
		gp$fit$rW <- rW
	}
	if (need.LTinv.rW) 
		gp$fit$LTinv.rW <- LTinv.rW

	gp
}
