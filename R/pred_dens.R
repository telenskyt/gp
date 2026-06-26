

# input:
# pred - predictions from LOO/LGO-CV approximation, or real fold CV
#		- to co vylezlo z __predictor casti templaty ("s" u cejek)
#		- s covariance matrix na urovni granularity W (tj. cov matrix marginalizovana na ty komponenty likelihood)
#			-> cili musi byt uz transformovana tim predictorem
#				- u cejek tam je to "A" identity matrix, takze to posune jenom meanem
#				- u q atlasu to neni linearni transformace ???
#

# predictive density using integrate()

# !!!! slo by udelat i multi-dimensional verzi pro W.type = "bdiag" s nejakou cubature()
#	- ??? ale jak pak pres to iterovat? Neresi to nekdo? podle me mask(.., W) je k tomu blizko


pred_dens_logit__int <- function(gp, data = gp[["data"]], pred, subdivisions = 10000L, ...)
{
	#vec <- gp$negLogLik$process(data, pred)
	# vec - vector for given likelihood
	# but now, for each factor of the likelihood / block of the W, we need to integrate over it
	# since here we use integrate() which is just 1-D integration, this function will thus only work if there is 1-1 correspondence,
	# which means W is diagonal
	if (gp$W.type != "diag")
		stop("W is not diagonal - please choose a different method for calculating predictive density (e.g. sampling)")
	
	if (is.list(pred)) {
		stopifnot(!is.null(pred[["f"]]))
		stopifnot(!is.null(pred[["f_cov_masked"]]))
		stopifnot(inherits(pred[["f_cov_masked"]], "sparseMatrix"))
		stopifnot(isDiagonal(pred[["f_cov_masked"]]))
		pred <- as.matrix(data.frame(f = pred$f, f_SE = sqrt(diag(pred[["f_cov_masked"]]))))
	} else if (is.data.frame(pred)) {
		pred <- as.matrix(pred)
	} else {
		stopifnot(is.matrix(pred))
	}
	stopifnot(all(c("f", "f_SE") %in% colnames(pred)))
	# for now, we only support the case when each f[i] corresponds to one likelihood component
	#!!!! zde je jeste potreba nejaky typ checku,ze korespondence mezi in	dexy y a f je 1:1!!!!!
	#?? nejaky attribute te __process templaty?
	

	hyperpar <- gpHyperparList(gp)
	N <- gpDataSize(data, fact = gp$GP_factor)
	works_for_reindex_confirmed <- FALSE
	pd <- c()
	err <- c()		
	for (i in 1:N) {
		data_i <- gpDataSubset(data, fact = gp$GP_factor, i)
		par_i <- gp:::gpGetParForLikTemplate(gp, pred[i, "f"], data_i, hyperpar)
		stopifnot(!is.null(par_i[["f"]]))
		if (length(par_i$f) > 1) {
			stopifnot(all(par_i$f == pred[i, "f"]))
			works_for_reindex_confirmed <- TRUE
		}
		qlo <- qnorm(1e-10, pred[i, "f"], pred[i, "f_SE"])
		qhi <- qnorm(1 - 1e-10, pred[i, "f"], pred[i, "f_SE"])
		int <- integrate(function (x) { 
			par_i$f <- rep(x, length(par_i$f)) # needed for the case of reindexing
			gp$lik$stages(data_i, par_i, stages = 1:4)*dnorm(x, pred[i, "f"], pred[i, "f_SE"])
		}, qlo, qhi, subdivisions = subdivisions, ...)
		pd[i] <- int$value
		err[i] <- int$abs.error/int$value
	}
	if (works_for_reindex_confirmed)
		cat("!!! it was confirmed that it works for reindexing! CONF-123\n")
	list(value = -sum(log(pd)), err = sum(log1p(err)))
}

# predictive density using importance sampling

pred_dens_logit__IS  <- function(gp, data, pred, samp = 1000L, ...)
{
	# convert `pred` to unified form where we always have cov matrix
	if (!is.null(pred$cov))
		cov <- pred$cov
	else
		cov <- diag(pred$f_SE)
	# !!! verify that this cov is full rank!!!
	
	# generate matrix of samples from MVN(pred$f, cov); each sample is a column - "f" vector
	#F <- ...
	
	# non-vectorized version
	par3_acc <- list()
	#for (i in all columns of F) {
	for (i in 1:ncol(F)) {
		par <- pred
		par$f <- F[,i] # replace the mean estimated f with the current sample
		par2 <- template$link(data, par) 
		par3 <- template$process(data, par2)
		for (el in names(par3)) {
			par3_acc[[el]] <- cbind(par3_acc[[el]], par3[[el]]) # here I prefer cbind(), as it creates matrix instead of data.frame (like bind_cols)
																# the problem might be if one argument is a data.frame() and the other is NULL
																# but that probably will not happen here
		}
	}
	# now we have the matrix of probabilities
	template$lik(par3_acc)
	# and we can get the density(F) and calculate weights from it...
	
	
	if (0) { # vectorized - vyzaduje specialni interface: bud vectorized process & lik zvlast, nebo dohromady
	#[P; matrix(ces) for all samples] <- gp$negLogLik$process__vectorized(data, F)
		# nebo $process a volat ve for loopu a brat z toho vse krom y
		
	#[y] <- gp$negLogLik$process(data, F[,1]) # you can get the data just once
	# or like this:
	#[y] <- gp$negLogLik$process(data, pred)
	
	#gp$negLogLik$ELEMENTWISE_or_VECTORIZED_somehow__lik(data or y, P)*density(F)/weights...
	#-> prezvejkat -> return
	}
	
}


pred_dens_logit__int_naive <- function (y, f, f_SE, subdivisions = 10000L, ...)
{
	pd <- c()
	err <- c()
	for (i in 1:length(y)) {
		# integrate vraci presnost atd., to ted zahodime
		qlo <- qnorm(1e-10, f[i], f_SE[i])
		qhi <- qnorm(1 - 1e-10, f[i], f_SE[i])
		if (y[i] == 1)
			int <- integrate(function (x) plogis(x)*dnorm(x, f[i], f_SE[i]), qlo, qhi, subdivisions = subdivisions, ...)
		else if (y[i] == 0)
			int <- integrate(function (x) (1-plogis(x))*dnorm(x, f[i], f_SE[i]), qlo, qhi, subdivisions = subdivisions, ...)
		pd[i] <- int$value
		err[i] <- int$abs.error/int$value
	}
	list(value = -sum(log(pd)), err = sum(log1p(err)))
}

pred_dens_logit__samp <- function (y, f, f_SE, samp = 1000L)
{
	mstart(id = "rnorm")
	F <- matrix(rnorm(samp*length(f), f, f_SE), ncol = samp)
	mstop(id = "rnorm")
	mstart(id = "eval")
	if (1) { # this way is faster
		Y <- matrix(2*y-1, nrow = length(y), ncol = samp)
		res <- -sum(log(rowMeans(plogis(Y*F))))
	} else {
		res <- -sum(log(rowMeans(dbinom(y, size = 1, plogis(F)))))
	}
	mstop(id = "eval")
	res
}

# use laplace approximation of the tilted distribution!!!
# => not very good: 1264.397 while the integrate()'d is 1218.479
pred_dens_logit__LA <- function (gp, y, f, f_SE)
{
	Y <- 2*y-1
	peak_p <- plogis(Y*f)
	#peak_cavity <- 1/sqrt(2*pi)/f_SE
	pr <- predict(gp, type = "latent", se.fit = TRUE)
	# pr[,"f_SE"] is a square root of the diagonal of the predicted posterior distribution covariance matrix	
	#pd <- (peak_p*peak_cavity)/(1/sqrt(2*pi)/pr[,"f_SE"])
	pd <- peak_p*pr[,"f_SE"]/f_SE # simplified the above
	-sum(log(pd))
}


# verify hypothesis that f-hat (MAP f estimate) is the peak of all the tilted distributions
# yes, it is true!!!
pred_dens_logit__verify_f <- function (y, f, f_SE)
{
	f_max <- c()
	for (i in 1:length(y)) {
		# integrate vraci presnost atd., to ted zahodime
		qlo <- qnorm(1e-10, f[i], f_SE[i])
		qhi <- qnorm(1 - 1e-10, f[i], f_SE[i])
		if (y[i] == 1)
			max <- optimise(function (x) plogis(x)*dnorm(x, f[i], f_SE[i]), c(qlo, qhi), maximum = TRUE)
		else if (y[i] == 0)
			max <- optimise(function (x) (1-plogis(x))*dnorm(x, f[i], f_SE[i]), c(qlo, qhi), maximum = TRUE)
		f_max[i] <- max$maximum
	}
	f_max
}
