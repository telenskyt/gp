# newdata - ***unscaled*** data! They must be on the same scale as was the unscaled training dataset.
# when you omit newdata, it will be done on the training dataset! And it will be super fast (if CI = NULL)


# maxn - maximum `n`, i.e. the maximum dataset size to fit at once (for the unit, see the comment on x_dataset_size of the appropriate model). 
#   It doesn't affect the result; is only to speed-up the computation and use less memory by splitting the dataset -> the K (pred,pred) matrix
#	is then smaller.
#   So it only speeds it up for CI != NULL (in the backsolve function)
#	for CI = NULL its only slowdown

# CI - at the moment, for model specific stuff (model_expand_predictions()), the CIs will always be set at 0.95, or NULL.
#	if this is set to NULL, no standard errors and CIs will be calculated, so it's faster, but also no mean on the response scale will be
#	calculated - because for this, the standard error of the latent GP `f` is needed.

# pred.sims (timingy u plotu: = 1e6 => 24s; 1e5 => 12.58 ; 1e4 => 11s -> volim 1e5 jako default)
#	=> now is obsolete, we are using numerical integration. Whether or not it's used depends on CI.

# comp_missing - how to compensate for missing components, in mean and variance. Method "none" does no compensation. In this case, 
#		the interpretation of the result's variance isn't much clear to me :-) "avg" is replacing the missing components by "averaging over them" - more
#		precisely, it predicts for the average effects of these components - for these, the variance/CI corresponds to confidence interval 
#		(uncertainty around the mean predictions) and not the prediction interval - see http://www.sthda.com/english/articles/40-regression-analysis/166-predict-in-r-model-predictions-and-confidence-intervals/
#		For models which use reindexed sub-matrices (like big matrix model) be aware of how the average is weighted - see CAUTION CAU001 in the code.
#	=> tak jsem zkusil jeste jeden napad s "intercept", ale nefacha to... hazi to na diagonale zaporny cisla a odmocnina pak vrati NaN..
#		myslim ze tady by to chtelo tu mou genialni metodu na variance of a draw... kterou jsem vymyslel genialne okolo toho 8.5., ale 
#		ona na tech komponentovych maticich nekterych nefacha, protoze cholesky rekne, ze neni positive-semi definite.
#	=> so the only reasonable setting now is "avg"; unless one cares only about the prediction mean (`f`) of the latent variable, in which case
#		a simple "none" would also do the job.
#
# groupMeans - vector of the length of the predicted dataset; will be coerced to a factor

# ... - passed to model_expand_predictions()

predict.graf <-	function(object, newdata = NULL, hyperpar = NULL, type = c('response', 'latent'),
						 CI = 0.95, maxn = NULL, pred.sims = 100000, w = NULL, Kx.cache = NULL, Kxx.cache = NULL, 
						 cov = FALSE, components = NULL, comp_missing = c("avg", "intercept", "none"), groupMeans = NULL, ...)  
{
	comp_missing <- match.arg(comp_missing)
	need <- function (object, x) if (is.null(object[[x]])) stop("Model object is missing the `", x, "` element - try to call `unpack` on it")
	
	if (class(newdata) %in% c('Raster', 'RasterBrick', 'RasterStack')) {
		
		# predict to a raster
		ans <- predict.graf.raster(object = object,
								   x = newdata,
								   type = type,
								   CI = CI,
								   maxn = maxn,
								   ...)
		
		return (ans)
		
	}
	
	use_model(object)
	need(object, "x")
	need(object, "y")
	need(object, "obsx")
	
	all_model_components <- object$model.opts$components
	if (is.null(components)) 
		components <- all_model_components
	validate_components(components, within = all_model_components)	
	
	type = match.arg(type)
	# set up data
	same <- FALSE
	if (is.null(newdata)) {
		# use already set up inference data if not specified
		x_new <- object$x
		if (is.null(components) || components == object$model.opts$components)
			same <- TRUE
	} else {
		# data must not be scaled!!!
		# use some hints to detect if user by mistake supplied scaled data
		if (model_data_is_scaled(newdata))
			stop("newdata must not be scaled for predict()!! predict() will scale them. The newdata must be on the same scale as the unscaled training data.")
		# convert any ints to numerics
		obsx <- model_prepare_obsx(newdata)
		x_new <- model_prepare_x(obsx, object$x) # scale atd. a ty dalsi veci 
				# it will scale the data in exactly the way that training data set was
		
		# cunarna tyto testy ale seru na to uz @@@
		if (is.data.frame(newdata$env) &&
			!all(sapply(object$obsx$env, class) == sapply(newdata$env, class)))
			stop('newdata must be with the same elements as used for inference, or NULL')
		#if (is.data.frame(newdata$visits) &&  # toto neplati v datetime plotu...
		#	!all(sapply(object$obsx$visits, class) == sapply(newdata$visits, class)))
		#	stop('newdata must be with the same elements as used for inference, or NULL')
	}
	
	# check CI
	if (!is.null(CI)) {
		if (!(CI == 'std')) {
			if (CI >= 1 | CI <= 0) {
				stop("CI must be a number between 0 and 1, or NULL")
			}
			err <- qnorm( 1 - (1 - CI) / 2 )
		}
	}

	if (is.null(maxn)) maxn <- ceiling(x_dataset_size(x_new)  / 10) # we are splitting the newdata, not the training dataset; ceiling needed instead of round here

	if (is.null(hyperpar))
		hyperpar <- hyperpar_decode(object$ls)
		
	# first, put up together the latent prediction
	mstart(id = "whole predict", mem_precise = TRUE)
	mstart(id = "whole pred", mem_precise = TRUE)
	if (is.null(CI)) { # if CIs aren't wanted
		pred <- pred(x_new, object, hyperpar = hyperpar, std = FALSE, maxn = maxn, same = same, w = w, Kx.cache = Kx.cache, Kxx.cache = Kxx.cache, cov = cov, components = components, comp_missing = comp_missing, groupMeans = groupMeans)
	} else {
		pred <- pred(x_new, object, hyperpar = hyperpar, std = TRUE, maxn = maxn, same = same, w = w, Kx.cache = Kx.cache, Kxx.cache = Kxx.cache, cov = cov, components = components, comp_missing = comp_missing, groupMeans = groupMeans)
	}
	cat("pred() took ")
	mstop(id = "whole pred")
	if (cov) {
		pred.cov <- pred$cov
		pred <- pred$pred
	}
	
	if (type == 'latent') {
		if (is.null(CI)) {
			ans <- pred
		} else if (CI == 'std') { # if standard deviations are wanted instead
			ans <- pred
		} else { # CI wanted
			upper <- pred[, 1] + err * pred[, 2]
			lower <- pred[, 1] - err * pred[, 2]
			dCI <- cbind(lower, upper)
			colnames(dCI) <- c(paste0("f_lower_", round(100 * CI), "CI"), #paste("lower ", round(100 * CI), "% CI", sep = ""),
							   paste0("f_upper_", round(100 * CI), "CI")) #paste("upper ", round(100 * CI), "% CI", sep = "")
			ans <- cbind(pred, dCI)							   
		}	
	} else { # if not only latent wanted, predict also other stuff
		# correct for the "intercept" of the missing components (we chose not to do it in the 'latent' case)
		if (0) { # no longer done here; done in pred() now, so it can be done also for (co)variances
				# hmm, tak zpetne si rikam, ze tohle bylo skoro nejlepsi nakonec.. mohl bych to sem znova dat jako jednu z voleb
				# pro comp_missing... ale vyuzit ted uz ty predpocitany comp_means
			missing_components <- paste0(setdiff(strsplit(all_model_components, "")[[1]], strsplit(components, "")[[1]]), collapse = "")
			if (nchar(missing_components) > 0) {
				cat("calculating posterior mean for missing components ", missing_components, ", to correct for it\n")
				mstart(id = "missing", mem_precise = TRUE)
				pr <- pred(object$x, object, hyperpar = hyperpar, std = FALSE, maxn = Inf, components = missing_components) # @@@@@@@@@@@ experimental
				pred[,'f'] <- pred[,'f'] + mean(pr[,'f'])
				cat("pred() for missing components took ")
				mstop(id = "missing")			
			}
		}
		# calculate model specific response from the latent
		mstart(id = "mep", mem_precise = TRUE)
		if ("hyperpar" %in% names(formals(model_expand_predictions))) # supporting newer as well as older interface, when loading older models
			exp_pred <- model_expand_predictions(pred, hyperpar, pred.sampling = !is.null(CI), pred.sims = pred.sims, ...)
		else
			exp_pred <- model_expand_predictions(pred, pred.sampling = !is.null(CI), pred.sims = pred.sims, ...)
		ans <- cbind(pred, exp_pred)
		cat("model_expand_predictions() took ")
		mstop(id = "mep")
	}
	cat("returning memory - gc() took ")
	mstart(id = "gc")
	gc() # Return the memory. Important! :-) 
	mstop(id = "gc")
	cat("whole predict() took ")
	mstop(id = "whole predict")
	if (!cov)
		ans
	else
		list(pred = ans, cov = pred.cov)
}