
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


# generalized version for LGO-CV
# opravena verze, tohle je presne to co jsem chtel!!
# stale predpoklada, ze W je diagonal! Ale dela to po skupinach. Odsud je uz jen krucek k non-diagonal verzi!
#
# warning: can take longer time if only few folds chosen (then the block matrices are big and take time to invert)
#
# odvodil jsem si to sam jako analogii eq 67 ve Vehtari 2016, kdyz se to f[i] a t~[i] vezmou jako multidimenzionalni
# a pak se provede ten trik s delenim dvou gaussian PDF v eq 38, 39, ale zobecni se to na multidimenzionalni.
# a ono to konecne funguje!!!
#
#' @importFrom Matrix Matrix
#' @export
gpLooGen <- function (gp, folds)
{
	stopifnot(is.vector(gp$fit$W))
	pr.cov <- predict(gp, type = "latent", cov.fit = TRUE)

	mask <- outer(folds, folds, "==")*1L
	mask <- Matrix(mask, sparse = TRUE)

	pr.cov <- pr.cov$cov * mask # this mask will practically make this block-diagonal w.r.t. the folds (generalization of diag()!)

	M <- solve(pr.cov - diag(1/gp$fit$W))
	R <- pr.cov - pr.cov %*% M %*% pr.cov # use matrix inversion lemma, the form from R&W 2006 page 201 top!
		# this calculates the inverse(inverse(pr.cov) - gp$fit$W)), even though the inverse(pr.cov) itself doesn't exist!!!
		# (pozn.: na toto reseni me privedla page https://math.stackexchange.com/a/1251958/15731, i kdyz s tim primo nesouvisi - ale nasel jsem ji kdyz 
		# jsem googlil: division of two multivariate normal PDFs)
		#
		# pozn: mozna by to slo spocitat nejak rychleji a lepe bez toho inversu, kdyz mam L? Ale asi ne, protoze ta pr.cov je uz po tom maskovani
	
	v.loo <- diag(R) # here we marginalize each an every single one, for this is needed the diagonality of W
					# the generalization will not do this, it will have v.loo non-diagonal, but then it will be harder to integrate the predictive density over it
	f.loo <- gp$fit$f - v.loo * gp$fit$a
	data.frame(f = f.loo, f_SE = sqrt(v.loo))
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