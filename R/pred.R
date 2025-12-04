# toto vraci prediction na urovni toho f (resp. f*)
# tj. Algorithm 3.2 z Rasmussen & Williams 2006 (chap 3.4, pg 47), az po radek 6, radek 7 (integrace) zde neimplementovana
#
# std = FALSE, cov = FALSE - no standard errors, no variance computed  (just one column, f, is returned)
# std = TRUE,  cov = FALSE - individual variance (diagonal of the covariance matrix) is computed, and returned sqrt of it (columns "f", "f_SE" returned)
# std = TRUE,  cov = TRUE  - whole covariance matrix is computed (returns a list, a vector f and covariance matrix cov)
#								=> drive vracel pred.cov do glob. promenne
# same = TRUE in case that you predict for the training dataset (should be faster, everything will be reused)
#
# mn - mean function transformed to scale of probabilities
#
# w - calculate average of these predictions with these weights. Will do it before predicting. Which is same as after, as it is all linear,
# but probably faster.
#		if w is matrix, calculate more different averages with different weights. Each column of `w` corresponds to one average - one set of weights
#
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
#	26.6.2006 - nove zjisteno - mozna ze to comp_missing = "avg" je blbe, protoze tam muzou vyjit zaporny cisla na diagonale!
#				viz GRaF-my_modif\my_devel\pred-com_missing=avg_zaporny_cisla_na_diag\ => 27.6. jsem zjistil ze to bylo tim CAUTION CAU001,
#				vyreseno cistou implementaci
#
# groupMeans - vector of the length of the predicted dataset; will be coerced to a factor
#
# Three modes of prediction are thus possible (je jich mnohem vic, ale toto je principialni deleni):
#	1) without any SE/CI - just mean (std = FALSE and cov = FALSE). Fastest
#	2) with individual variance/SE/CI  (std = TRUE and cov = FALSE and w = NULL). Only diagonal of covariance matrix needed.
#	3) with whole cov matrix needed (cov = TRUE or (std = TRUE and w != NULL)).

pred <- function(gp, predx, std = TRUE, maxn = 250, same = FALSE, hyperpar = gpHyperparList(gp), w = NULL, Kx.cache = NULL, Kxx.cache = NULL, cov = FALSE, 
	components = NULL, recursive = FALSE, comp_missing = c("avg", "intercept", "none"), groupMeans = NULL)
{
	fit <- gp$fit
	#cat("pred(): print(mstart):\n")
	#print(mstart)
	#cat("pred(): print(mstartOptions(mem_precise)):\n")
	#print(mstartOptions("mem_precise"))
	if (identical(comp_missing, FALSE))
		comp_missing <- "none"
	else
		comp_missing <- match.arg(comp_missing)
	need <- function (obj, x) if (is.null(obj[[x]])) stop("Model object is missing the `", x, "` element - try to call `unpack` on it")
	n <- x_dataset_size(predx)
	#validate_components(components) # allowing empty components here; K_matrix will solve it if needed :)
	# and also trusting the validation from predict
	if (!recursive)
		cat("pred(..., std = ", std, ", same =", same, ", components = ", components, ")\n")
	if (same && !is.null(components) && components != model.obj$components) {
		warning("pred(): resetting same to FALSE as arg components (\"", components, "\") isn't the same as model components (\"", model.obj$components, "\")")
		same <- FALSE
	}
	if (cov)
		std <- TRUE
	if (!is.null(groupMeans) && (cov || !is.null(w)))
		stop("these options (groupMeans & w, cov) do not fit together")
	if (n > maxn && !same && is.null(Kx.cache) && is.null(Kxx.cache) && is.null(w) && !cov && is.null(groupMeans)) { # split the data set
		if (!is.null(predx$visits)) # @@@ docasne reseni; pak lze splitting opravit oro tento pripad nebo ho vyhodit uplne
			warning("Unsolved problem in this variant of splitting: subset reindexes but we are putting it together just by rbind... so it would mess up the data. Set maxn = Inf as a workaround.")
		cat("n = ", n, ", splitting dataset to chunks of size maxn = ", maxn, "\n")
		inds <- split(1:n, ceiling((1:n) / maxn))
		#print(sapply(inds, length))
		fun <- function(ind, X, fit, std, maxn, same, hyperpar, components, comp_missing) {
			pred(gp, gpDataSubset(X, ind), fit, std, maxn, same, hyperpar, components = components, recursive = TRUE, comp_missing = comp_missing)
		}    
		prediction <- lapply(inds, fun, predx, fit, std, maxn, same, hyperpar, components, comp_missing)
		prediction <- do.call('rbind', prediction)
	} else {
		mstart(id = "predmean")
		cat("\t* computing prediction (posterior mean):\n")
		if(same) {
			# if predicting back to input data re-use covariance matrix
			need(fit, "K")
			need(fit, "MAP")
			cat("\t\t* (same)\n")
			Kx <- fit$K
			prediction <- fit$MAP
			stopifnot(nrow(prediction) == ncol(Kx)) # prediction is column vector
			if (!is.null(w)) {
				stopifnot((is.matrix(w) && nrow(w) == ncol(Kx)) || (is.vector(w) && length(w) == ncol(Kx)))
				# nrow(w) == nrow(prediction)
				prediction <- crossprod(w, prediction)
			}
		} else {
			cat("\t\t* K_matrix(train, pred)... ")
			need(fit, "x")
			need(fit, "a")
			mstart(id = "Kx")
			Kx <- K_matrix(gp, hyperpar, x2 = predx, K.cache = Kx.cache, components = components) # radky: training, sloupce: predicted records
			
			# for all of the missing components, put there Covariance with average of all the predictions from the training dataset (on those components)!
			# this will create the proper compensation for the missing components in the prediction mean as well as covariance!
			all_model_components_v <- strsplit(fit$model.opts$components, "")[[1]]
			components_v <- strsplit(components, "")[[1]]
			missing_components_v <- setdiff(all_model_components_v, components_v)
			if (comp_missing != "none" && length(missing_components_v) > 0) {
				stopifnot(!is.null(attr(fit$K, 'comp_means')))
				stopifnot(all(missing_components_v %in% colnames(attr(fit$K, 'comp_means'))))
				Kt_miss <- apply(attr(fit$K, 'comp_means')[,missing_components_v, drop = FALSE], 1, sum)
				stopifnot(nrow(Kt_miss) == nrow(Kx))
				if (comp_missing == "avg")
					Kx <- Kx + matrix(rep(Kt_miss, ncol(Kx)), nrow = nrow(Kx), ncol = ncol(Kx)) # add the vector of the covariances of missing components with their mean to each column of Kx
				else if (comp_missing == "intercept") { 
					# spocitej corresponding "intercept sigma2" scale za Kt_miss a udelej z nej matici!
					intcept.sigma2_scale <- Kt_miss %*% fit$a / sum(fit$a)
					cat("pred: compensation for missing components ", missing_components_v, " is intcept.sigma2_scale = ", intcept.sigma2_scale, "\n")
					Kx <- Kx + matrix(intcept.sigma2_scale, nrow = nrow(Kx), ncol = ncol(Kx))
				}
			    #            "\t\t* K_matrix(train, pred)... "
				#cat(sprintf("\t\t* - for missing c. %6s.. ", missing_components))
				#mstart(id = "K missing")
				#Kt_miss <- K_matrix(hyperpar, fit$x, x2 = NULL
				#cat(sprintf("(%4dx%-4d): \t", nrow(Kx), ncol(Kx)))
				#mstop(id = "K missing")
			}
			cat(sprintf("(%4dx%-4d): \t", nrow(Kx), ncol(Kx)))
			mstop(id = "Kx")
			mpred <- mnfun(fit$x, fit$y, predx) # mpred - mean on the normal scale of GP
			stopifnot(length(mpred) == ncol(Kx))
			if (!is.null(w)) {
				stopifnot((is.matrix(w) && nrow(w) == ncol(Kx)) || (is.vector(w) && length(w) == ncol(Kx)))
				# nrow(w) == length(mpred)
				mpred <- crossprod(w, mpred)
			}
		}
		if (!is.null(w) && (!same || std)) {
			stopifnot((is.matrix(w) && nrow(w) == ncol(Kx)) || (is.vector(w) && length(w) == ncol(Kx)))
			cat("\t\t* K* w       : \t\t\t\t\t")
			mstart(id = "w")
			Kx <- Kx %*% w
			mstop(id = "w")
		}
		if (!same) {
			cat("\t\t* K* a + mean ... \t\t\t\t")
			mstart(id = "mul")
			prediction <- crossprod(Kx, fit$a) + mpred
				# mas-li vektor kovariancí dane predikovane veci se vsemi training (to je dany sloupec z Kx), tak skalarni soucin s 
				# vektorem fit$a (to je kouzelny vektor) ti dá predikovanou hodnotu! :-)
			
			mstop(id = "mul")
		}
		# prediction je sloupcovy vektor
		if (!is.null(groupMeans)) {
			stopifnot(length(groupMeans) == nrow(prediction))
			prediction <- t(t(tapply(prediction, groupMeans, mean))) # col. vector
		}
		colnames(prediction) <- c("f")
	cat("\t\t* returning memory - gc() took \t\t\t")
	mstart(id = "gc")
	gc() # Return the memory. Important! :-) 
	mstop(id = "gc")		
		cat("\t\t- prediction mean part took \t\t\t")
		mstop(id = "predmean")
		
		mstart(id = "std")
		if (std) {
			cat("\t* computing posterior (co)variances:\n")
			need(fit, "L")
			need(fit, "W")			
			cat("\t\t* backsolve... \t\t\t\t\t")
			mstart(id = "backsolve")
			stopifnot(all(fit$W >= 0))
			v <- backsolve(fit$L, sqrt(as.vector(fit$W)) * Kx, transpose = T)
			mstop(id = "backsolve")
			
					# vector * Matice = diag(vector) %*% Matice
			# Kxx - k(x*,x*)
			#		@@@ TODO!!! Optimization: calculate only diagonal :-) But if maxn is relatively small, this is not so bad as calculating K for the big matrix
			cat("\t\t* K_matrix(pred, pred)")
			mstart(id = "Kxx")
			if (!same) {
				cat("...  ")
				Kxx <- K_matrix(hyperpar, x1 = predx, x2 = NULL, K.cache = Kxx.cache, components = components)
				if (comp_missing != "none" && length(missing_components_v) > 0) {
					if (comp_missing == "avg")
						Kxx <- Kxx + matrix(mean(Kt_miss), nrow = nrow(Kxx), ncol = ncol(Kxx))
					else if (comp_missing == "intercept") {
						#intcept.sigma2_scale <- Kt_miss %*% fit$a / sum(fit$a)
						Kxx <- Kxx + matrix(intcept.sigma2_scale, nrow = nrow(Kxx), ncol = ncol(Kxx))
					}
				}
			} else {
				cat("same ")
				need(fit, "K")
				Kxx <- fit$K
			}
			#cat("(", paste(dim(Kxx), collapse = "x"), "): \t")
			cat(sprintf("(%4dx%-4d): \t", nrow(Kxx), ncol(Kxx)))
			mstop(id = "Kxx")
			if (!is.null(w)) {
				cat("\t\t* w^T K** w  : \t\t\t\t\t")
				mstart(id = "w") 
				Kxx <- t(w) %*% Kxx %*% w
				mstop(id = "w")
			}
			cat("\t\t* K** - v^T v: \t\t\t\t\t")
			mstart(id = "mul") 
			if (cov) { # we will need to calculate full cov matrix	
				if (!is.null(groupMeans))
					stop("these options (groupMeans & w, cov) do not fit together")
				predvar <- Kxx - crossprod(v)
				stopifnot(all(diag(predvar) >= 0, na.rm = TRUE)) # allow NA/NaN, e.g. with NaN weights
				prediction <- cbind(prediction, sqrt(diag(predvar)))
				colnames(prediction) <- c("f", "f_SE")
				if (cov) { # vrat covariancni matici celou - do glob promenne, temporary dirty hack!
					#cat(" I am here")
					prediction <- list(pred = prediction, cov = predvar)
					# drive vracel pred.cov do glob. promenne
				}
			} else {
				if (is.null(groupMeans))
					predvar <- diag(Kxx) - colSums(v^2)
				else { # groupMeans . Otestovano rucne v play_with_predictions.R - mean i variance (funguje)
					predvar <- c()
					stopifnot(length(groupMeans) == nrow(Kxx))
					grp <- as.factor(groupMeans)
					for (i in levels(grp)) {
						predvar <- c(predvar, mean(Kxx[grp == i, grp == i] - crossprod(v[, grp == i, drop = FALSE])))
					}
				}
				stopifnot(all(predvar >= 0, na.rm = TRUE))		 # allow NA/NaN, e.g. with NaN weights	
				prediction <- cbind(prediction, sqrt(predvar))
				colnames(prediction) <- c("f", "f_SE")			
			}
			mstop(id = "mul")
		}
		#stop("check it out")
	cat("\t\t* returning memory - gc() took \t\t\t")
	mstart(id = "gc")
	gc() # Return the memory. Important! :-) 
	mstop(id = "gc")		
		cat("\t\t- prediction (co)variance part took \t\t")
		mstop(id = "std")
		#cat("\n") 
		#stopifnot(all(is.finite(prediction[,"f_SE"])))
	}
	prediction
}