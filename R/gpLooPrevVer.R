# gpLGOpred0 - previous version of computation, had problems in Lapwing nest survival model (non-diagonal W)
#' @export
gpLGOpred0 <- function (gp, fold.col, fold.fact = "1")
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

	#stopifnot(is.vector(gp$fit$W))
	pr.cov <- predict(gp, type = "latent", cov.fit = TRUE)

	mask <- crossprod(fac2sparse(factor(fold.col))) # make a sparse matrix `mask` of the GP dimension, which has 1 where the elements belong to the same fold
	pr.cov <- mask(pr.cov$cov, mask) # now mask the posterior to it

	if (gp$W.type == "diag")
		M <- solve(pr.cov - Diagonal(x = 1/gp$fit$W))
	else {
		Winv <- solve(gp$fit$W)
		M <- solve(pr.cov - Winv)
	}
	R <- pr.cov - pr.cov %*% M %*% pr.cov # use matrix inversion lemma, the form from R&W 2006 page 201 top!
		# this calculates the inverse(inverse(pr.cov) - gp$fit$W)), even though the inverse(pr.cov) itself doesn't exist!!!
		# (pozn.: na toto reseni me privedla page https://math.stackexchange.com/a/1251958/15731, i kdyz s tim primo nesouvisi - ale nasel jsem ji kdyz 
		# jsem googlil: division of two multivariate normal PDFs)
		#
		# pozn: mozna by to slo spocitat nejak rychleji a lepe bez toho inversu, kdyz mam L? Ale asi ne, protoze ta pr.cov je uz po tom maskovani
	
	if (gp$W.type == "diag")
		f.loo.cov.masked <- Diagonal(x = diag(R))
	else
		f.loo.cov.masked <- mask(R, gp$fit$W) 
	f.loo <- gp$fit$f - f.loo.cov.masked %*% gp$fit$a
	list(f = as.matrix(f.loo), f_cov_masked = f.loo.cov.masked)
		# keep f as a column vector, same as in gpFitLaplace()
}

# gpLGOpred0_nomask - previous version of computation, but without the masking to the level of W
#
#' @export
gpLGOpred0_nomask <- function (gp, fold.col, fold.fact = "1")
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

	#stopifnot(is.vector(gp$fit$W))
	pr.cov <- predict(gp, type = "latent", cov.fit = TRUE)

	mask <- crossprod(fac2sparse(factor(fold.col))) # make a sparse matrix `mask` of the GP dimension, which has 1 where the elements belong to the same fold
	pr.cov <- mask(pr.cov$cov, mask) # now mask the posterior to it

	if (gp$W.type == "diag")
		M <- solve(pr.cov - Diagonal(x = 1/gp$fit$W))
	else {
		Winv <- solve(gp$fit$W)
		M <- solve(pr.cov - Winv)
	}
	R <- pr.cov - pr.cov %*% M %*% pr.cov # use matrix inversion lemma, the form from R&W 2006 page 201 top!
		# this calculates the inverse(inverse(pr.cov) - gp$fit$W)), even though the inverse(pr.cov) itself doesn't exist!!!
		# (pozn.: na toto reseni me privedla page https://math.stackexchange.com/a/1251958/15731, i kdyz s tim primo nesouvisi - ale nasel jsem ji kdyz 
		# jsem googlil: division of two multivariate normal PDFs)
		#
		# pozn: mozna by to slo spocitat nejak rychleji a lepe bez toho inversu, kdyz mam L? Ale asi ne, protoze ta pr.cov je uz po tom maskovani
	
	f.loo <- gp$fit$f - R %*% gp$fit$a
	list(f = as.matrix(f.loo), cov = R)
		# keep f as a column vector, same as in gpFitLaplace()
}

