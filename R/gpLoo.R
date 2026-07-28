
# implementation based on Vehtari et al. (2016), eq 38-39; probably what the paper calls LA-LOO (with LA-L)
# folds - in case these are specified, it becomes LGO-CV (leave-group-out CV)
# G - incidence matrix (e.g. for buffered LOO-CV) - doesn't have to be groups, can be e.g. the buffer
# -> nefacha to (to folds, a G nebude taky; bez toho facha skvele)
#' @export
gpLoo <- function (gp, folds = NULL, G = NULL)
{
	W <- gp$fit$W
	if (!is.null(folds)) {
		W <- ave(W, folds, FUN = sum)
	} else if (!is.null(G)) {
		W <- G %*% W
	}
	pr <- predict(gp, type = "latent", se.fit = TRUE)
	pr.cov <- pr[,"f_SE"]^2 # diagonal of the predicted posterior distribution covariance matrix
	v.loo <- 1/(1/pr.cov - W) 
		# ?? hele, neni to v.loo
		# protoze pr.cov = (K^-1 + W)^-1 # vynecham to diag() a udelam ho az na konci
		# takze 1(1/pr.cov - W), kdyz zobecnim = (pr.cov^-1 - W)^-1 = (K^-1)^-1 = K !!! 
		# takze je to totez jako vzit diag(K)???
		#
		# jenze trik je asi v tom, ze 1/diag(A) neni totez co diag(A^-1) - ANO! Presne!
	#d1 <- d1(gp, f = gp$fit$f, hyperpar = gpHyperparList(gp)) # this is the same as gp$fit$a
	f.loo <- gp$fit$f - v.loo * gp$fit$a
	data.frame(f = f.loo, f_SE = sqrt(v.loo))
}

# generalized version for LGO-CV
# jen zoufaly pokus, asi je to blbost!! (neslo mi invertovat pr.cov, tak jsem zkusil ginv()... ale je to blbost, nefacha)
#' @importFrom MASS ginv
#' @export
gpLooGen_v0 <- function (gp, folds)
{
	pr.cov <- predict(gp, type = "latent", cov.fit = TRUE)

	mask <- outer(folds, folds, "==")*1L
	pr.cov <- pr.cov$cov * mask
	pr.cov.inv <- ginv(pr.cov)

	R <- ginv(pr.cov.inv - gp$fit$W)
	v.loo <- diag(R)

	f.loo <- gp$fit$f - v.loo * gp$fit$a
	data.frame(f = f.loo, f_SE = sqrt(v.loo))
}



#' LGO-CV prediction - an approximation of leave-group-out cross-validation
#' @param gp GP model object
#' @param fold.col either a name of a column in the main table, or a vector along factor \code{fold.fact}.
#'		The column (or the supplied vector) must be a vector of integers from \code{1} to \code{N} (\code{N} being the number of folds), 
#'		specifying the number of the cross-validation fold the given record belongs to.
#' @param fold.fact a factor along which the cross-validation folds (\code{fold.col}) are specified. Factor \code{"1"} means the folds are specified for the rows of the main table.
#'		Note that if the GP dimension is given by some real grouping factor (i.e. \code{gp$GP_factor != "1"}), then the \code{fold.fact} must be that factor.
#
#' @importFrom Matrix Matrix
#' @importFrom Matrix fac2sparse
#' @importFrom Matrix Diagonal
#' @export
#
# warning: can take longer time if only few folds (groups) chosen (then the block matrices are big and matrix operations take longer time)
#
# odvodil jsem si to sam jako analogii eq 67 ve Vehtari 2016, kdyz se to f[i] a t~[i] vezmou jako multidimenzionalni
# a pak se provede ten trik s delenim dvou gaussian PDF v eq 38, 39, ale zobecni se to na multidimenzionalni.
# a ono to konecne funguje!!!

gpLGOpred <- function (gp, fold.col, fold.fact = "1")
{
	fold.col <- gpFitCV__validate_and_get_fold_col(gp, fold.col, fold.fact)
		# validate fold.col and fold.fact arguments and get evaluated fold.col (as a vector)

	# similar code as in gpFitCV:
	if (fold.fact != gp$GP_factor) {
		# now, reindex the fold.col to the dimension of the prediction (gp$GP_factor)
		# thanks to the condition checked in gpFitCV__validate_and_get_fold_col(), the only case when reindexing might be needed is when GP_factor = "1" and fold.fact = something else (=something smaller)
		stopifnot(gp$GP_factor == "1") # consequence of the checks above
		stopifnot(fold.fact != "1") # consequence of the checks above
		# so now, we have to reindex fold.col to the main table, which is exactly corresponding to the GP dimension
		stopifnot(gp$GP_size == gpDataSize(gp$obsdata, "1"))
		
		stopifnot(gp:::gpDataHasMainTable(gp[["obsdata"]])) # has to have it in this case
		fold_idx_col <- paste0(fold.fact, "_idx")
		fold.col <- fold.col[gp$obsdata[[1]][[fold_idx_col]]]		
	}
	# fold.col is now of the dimension of the GP

	### calculate S - the posterior covariance matrix masked by the groups
	S <- predict(gp, type = "latent", cov.fit = TRUE)

	mstart(id = "rest_of_LGO")
	mask <- crossprod(fac2sparse(factor(fold.col))) # make a sparse matrix `mask` of the GP dimension, which has 1 where the elements belong to the same fold
	S <- mask(S$cov, mask) # now mask the posterior to it

	mstart(id = "rW")
	### retrieve rW - the cholesky factor of W = t(rW) %*% rW
	stopifnot(!is.null(gp[["fit"]]))
	if (gp$W.type == "diag") {
		if (!is.null(gp$fit[["rW"]]))
			rW <- Diagonal(x = gp$fit$rW)
		else {
			need(gp$fit, "W")
			rW <- Diagonal(x = sqrt(gp$fit$W))
		}
	} else {
		if (!is.null(gp$fit[["rW"]]))
			rW <- gp$fit$rW
		else {
			need(gp$fit, "W")
			warning("had to compute rW using cholesky; would be more optimal to have rW stored in the object.")
			rW <- chol_W(gp$fit$W)
		}
	}
	cat("getting the rW took ")
	mstop(id = "rW")
	
	mstart()
	M <- - rW %*% S %*% t(rW)
	diag(M) <- diag(M) + 1
	cat("M <- took ")
	mstop()

	mstart()
	LM <- Matrix::chol(M)
	cat("LM <- took ")
	mstop()
	
	mstart()
	v <- Matrix::solve(t(LM), rW %*% S, sparse = TRUE) # backsolve() for the sparse matrix
	cat("v <- Matrix::solve(t(LM), rW %*% S, sparse = TRUE) took ")
	mstop()
#	if (gp$W.type == "diag")
#		R <- diag(S) + colSums(v^2)
#	else
	mstart()
		R <- S + crossprod(v)
	cat("R <- S + crossprod(v) took ")
	mstop()
		# ??? co tady ten crossprod?? nemelo by se to maskovat na group mask? Asi z kouzelnych aritmetickych duvodu neni treba?
		# the W.type == "diag" pathway could be probably optimized a bit more - but that only in case of LOO and not LGO?
	R <- Matrix::forceSymmetric(R)
		# it won't waste time by testing the symmetry; it will just take half of it and say the rest is symmetric :)
		# note though that when you do mask(R, ...), the symmetry is gone - mask() doesn't support it now
	mstart()
	f.loo <- gp$fit$f - R %*% gp$fit$a
	cat("f.loo <- took ")
	mstop()
	cat("Computing the rest of LGO took ")
	mstop(id = "rest_of_LGO")
	list(f = as.matrix(f.loo), cov = R)
		# keep f as a column vector, same as in gpFitLaplace()
		# we will keep the masking of R by mask(R, g$fit$W) on the caller now
}

# GPStuff implementation
# vypada to ze nefunguje, ze tam maji chybu, ze to neni ekvivaletni tomu nasemu co delam podle Vehtari 2016
#' @export
gpLoo2 <- function (gp)
{
	stopifnot(is.vector(gp$fit$W))
	v <- backsolve(gp$fit$L, diag(sqrt(gp$fit$W)))
	v.loo <- 1/diag(crossprod(v)) - 1/gp$fit$W

	#d1 <- d1(gp, f = gp$fit$f, hyperpar = gpHyperparList(gp)) # this is the same as gp$fit$a
	f.loo <- gp$fit$f - v.loo * gp$fit$a
	data.frame(f = f.loo, f_SE = sqrt(v.loo))
}

#' @export
gpLooExact <- function (gp, fold.fact = g$GP_factor, ...)
{
	N <- gpDataSize(gp$data, fact = fold.fact)
	g <- gpFitCV(gp, fold.fact = fold.fact, fold.col = 1:N, log.fn = NULL, dump.fn = NULL, pred.options = list(type = "latent", se.fit = TRUE), opt.h = FALSE, grad.computation = FALSE, ...)
	g
}

# co delam, vzdyt mam fci gpFitCV!!!
if (0) {
gpLooExact <- function (gp)
{
	for (i in 1:gp$GP_size) {
		# remove i from data
		gp_i <- gp
		gp_i$data <- gpDataSubset(gp$data, fact = "1", -i) # !!! can I do -i?
		gp_i <- gpFit(gp, opt.h = FALSE, grad.computation = FALSE)
		pr_i <- predict(gp, se.fit = TRUE)
		pr.cov <- pr_i[,"f_SE"]^2 # diagonal of the predicted posterior distribution covariance matrix

		pred(gp, ...) # predict for i 
		# bundle result
	}
}
}