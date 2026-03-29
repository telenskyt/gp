
# implementation based on Vehtari et al. (2016), eq 38-39; probably what the paper calls LA-LOO (with LA-L)
#' @export
gpLoo <- function (gp)
{
	pr <- predict(gp, type = "latent", se.fit = TRUE)
	pr.cov <- pr[,"f_SE"]^2 # diagonal of the predicted posterior distribution covariance matrix
	v.loo <- 1/(1/pr.cov - gp$fit$W) 
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