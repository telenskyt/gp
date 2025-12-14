
#' Pack the model - reduce its size for memory or disk storage
#'
#' Reduce model size by deleting auxiliary data structures which speed up the computations 
#' but take a lot of space.
#'
#' @param maximum should it be taken to maximum level? Default FALSE.
#' @export
gpPack <- function(gp, maximum = FALSE)
{
	gp$fit <- gpPackFit(gp$fit, maximum = maximum)
	gp$data <- NULL
	gp
}

# internal function
#
# maximum - should it be the max level of packing
# pack just the gp$fit object
gpPackFit <- function(fit, maximum = FALSE)
{
	fit$L <- NULL
	fit$K <- NULL
	fit$W <- NULL
	if (maximum) 
		fit$f <- NULL # nakonec i bez tohoto se da obejit!!! :-))) # ale defaultne necham tam radsi uplne to puvodni, kvuli ruzne num nestabilite....
	else
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
#' @param need.K,need.L logical - in case of \code{compute = TRUE}, specify whether particular matrix (K, L) is needed. If only one of these matrices (K, L) is needed, it can 
#' significantly save memory, and in case of L matrix, also CPU time. Note that if \code{compute = TRUE}, \code{K} matrix will be calculated regardless of \code{need.K};
#' so \code{need.K = FALSE} will only save the memory size of the resultant object, but not CPU time.
#'
#' @export
gpUnpack <- function (gp, compute = TRUE, need.K = compute, need.L = compute)
{
	gp$data <- gpDataPrepare(gp, gp$obsdata)

	if (!compute) {
		stopifnot(!need.K)
		stopifnot(!need.L)
		return(x)
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

	gp$fit$W <- -(gp$fit$wt * d2(gp, gp$fit$f, gp$data, hyperpar))

	stopifnot(all(gp$fit$W >= 0))

	if (!need.L) {
		if (!need.K) {
			empty_K <- NA
			mostattributes(empty_K) <- attributes(gp$fit$K)
			gp$fit$K <- empty_K
		}
		return(gp)
	}

	# recompute L
	mstart(id = "L")
cat("Re-computing L matrix ... ")

	rW <- sqrt(gp$fit$W)
	xx <- rW %*% t(rW)# * K + diag(n)
gc()
    xx <- xx * gp$fit$K
	if (!need.K) {
		empty_K <- NA
		mostattributes(empty_K) <- attributes(gp$fit$K)
		gp$fit$K <- empty_K
	}
gc() # ta vlozena gc() tady jsou velmi dulezita!!!

	xx <- xx + diag(n)
gc()
	gp$fit$L <- tryCatch(chol(xx), error = function(x) return(NULL))
	rm(xx)
	mstop(id = "L")
gc()
#mstop()
	gp
}
