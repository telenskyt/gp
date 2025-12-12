
#' zahod vse co se da "snadno" vypocitat
#' @export
gpPack <- function(gp)
{
	gp$fit$L <- NULL
	gp$fit$K <- NULL
	gp$fit$W <- NULL
	gp$data <- NULL

	#gp$fit$f <- NULL # nakonec i bez tohoto se da obejit!!! :-))) # ale necham tam radsi uplne to puvodni, kvuli ruzne num nestabilite....
	rownames(gp$fit$f) <- NULL # rownames take a lot of space! Just one column matrix with them 279872 bytes, without 31440 bytes!!!
							# And they only had 6 characters each! (nrow = 3881, ncol = 1). R is not very efficient storing them.
							# And rownames are in obsx already

	gp
}

# compute - priprava na vypocty - predikce. Je ale casove narocna viz vyse. compute = FALSE nebude delat ty casove narocne veci...
# need.K, need.L - v pripade compute = TRUE, pokud nepotrebuju nekterou z techto matic, lze to dat vedet a vyrazne to usetri pamet.
#	K se ale pocita vzdy (pro vypocet MAP) a jeji attributy (comp_means) zustanou i kdyz need.K = FALSE
#		- tak zamer s usetrenim pameti se moc nekona, protoze samotne vytvareni matice K_matrix() potrebuje (pro bigmatrix 9471x9471) whole process: 2703.5Mb
#			- musely by se udelat optimalizace navrzene v komentarich fce K_matrix v K-tsc_spat2.R, ale to by melo byt snadne
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
