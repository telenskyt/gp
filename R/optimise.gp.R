# v1 - vratil jsem tam memoise(), ktery mel zakomentovany! Zrychli vypocet 2x! 
#		Nicmene je potreba moje verze sekvencniho memoise, ktery zaroven nezere tolik pameti.
#

# define the objective and gradient functions to optimise the hyperparameters
# of a GRaF model
objective <- function(theta, prior.pars, isfac, args, fun) 
{
	cat("opt: objective\n")
	cat("Hyper-parameter iteration: ", hyper_iter+1, "\n")
	# unpack theta
	h <- gpLinkInv(args$gp, theta)
	
	# set the hyperparameters
	args$h <- h
	
#	if (hyper_iter > 1) {
#		cat("\thyper_iter", hyper_iter, "head(f_start):", head(f_[,hyper_iter - 1]), "\n")
#		args$f_start <- f_[,hyper_iter - 1]
#	}
	# run the model
	m <- do.call(fun, args)
	
	# log likelihood and prior
	llik <- -m$mnll

	lpri <- sum(gpPrior(args$gp, h))
	
	# log posterior
	lpost <- llik + lpri
	
	# and objective
	if (args$use.prior)
		objective <- -lpost
	else
		objective <- m$mnll
	cat("opt: mnll:", m$mnll, ", -logPrior: ", -lpri, ", negLogPosterior:", -lpost, "\n")
	
	
	AUC <- NA
	if (!is.null(m$AUC))
		AUC <- m$AUC
	cat(" **** AUC = ", AUC, "\n")
	
	if (nrow(iter_info) < hyper_iter) {
		if (nrow(iter_info) < 1) 
			prev_iter_time <- 0
		else
			prev_iter_time <- iter_info[nrow(iter_info), 'time']
		iter_info <<- rbind(iter_info, data.frame(
			iter = hyper_iter,
			LA_iterations = m$iterations, # pojmenovano LA, ale muze to byt i EP!
			mnll = m$mnll, # marginal likelihood, ~ -log p(y|X), integral toho psi pro vsechna f
			negLogPosterior = -lpost, # posterior (marginal), ~ -log p(theta|X,y)  toto se minimalizuje v hyperparam. optimalizaci
			psi = ifelse(args$method == 'Laplace', m$psi, NA), # ~ - log p(f|X,y), vlastne posterior pro latent f, toto se minimalizuje v LA iteracich
			lastLAObjDiff = m$lastLAObjDiff, # last difference between LA iteration objective (psi, which is minimized). When positive, last iteration objective was going up, which is something to look at still.
			fStartReset = as.integer(m$f_start_was_reset), 
			AUC = AUC,
			time = (proc.time() - time_start)[3] # time_start je v optimise.graf
		))
		cat("** hyperpar iteration took ", (proc.time() - time_start)[3] - prev_iter_time, "s, total since beginning: ", (proc.time() - time_start)[3], "s\n")
		stopifnot(all(h == m$h))
		iter_h <<- rbind(iter_h, h)
		if (args$use_f_start)
			f_start <<- list(f = m$f, a = m$a, psi = m$psi)
		f_ <<- cbind(f_, m$f)
	}

	return (objective)
	
}


gradient <- function(theta, prior.pars, isfac, args, fun) {
	cat("opt: gradient\n")
	# unpack theta
	h <- gpLinkInv(args$gp, theta)
	
	# set the hyperparameters
	args$h <- h

#	if (hyper_iter > 1) {
#		cat("\thyper_iter", hyper_iter, "head(f_start):", head(f_[,hyper_iter - 1]), "\n")
#		args$f_start <- f_[,hyper_iter - 1]
#	}
	# run the model
	m <- do.call(fun, args)
	
	# gradient of llik w.r.t. h
	dLdh <- m$h_grads
	
	# gradient of h w.r.t. theta
	dhdtheta <- gpLinkInvDer(args$gp, theta)
	
	# gradient of lpri w.r.t. theta
	dpdtheta <- gpPriorGradient(args$gp, h) * dhdtheta
	
	# gradient of llik w.r.t. theta
	dLdtheta <- dLdh * dhdtheta
	
	# gradient of lpost w.r.t. theta
	dPdtheta <- dLdtheta + dpdtheta
	
	# gradient of objective w.r.t. lpost (aha tady se to obraci aby to odpovidalo MNLL!!!!
	dOdP <- -1
	
	# gradient of objective w.r.t. theta
	if (args$use.prior)
		dOdtheta <- dOdP * dPdtheta 
	else
		dOdtheta <- dOdP * dLdtheta
	
	if (nrow(iter_dmll_dh) < hyper_iter) {
		# derivace na skale log lik, nikoli negloglik
		# komentare peclive prekontrolovany 16.5.2020, muze slouzit jako dokumentace
		iter_dmll_dh <<- rbind(iter_dmll_dh, m$h_grads) # gradient of marginal log likelihood w.r.t. h (h ~ hyperparameters, SE lengthscales are squared) (mll = marginal log lik, pozor not negative)
		dpdh <- dpdtheta/dhdtheta
		iter_dP_dh <<- rbind(iter_dP_dh, m$h_grads + dpdh) # gradient of log posterior (the thing which is being maximized) w.r.t lh
															# (h ~ hyperparameters, SE lengthscales are squared)
														# toto je derivace -objective (kontroloval jsem) podle h; tato fce vraci toto * (-1) * dhdtheta (aby to byla derivace dle theta)
		iter_dp_dh <<- rbind(iter_dp_dh, dpdh) # gradient of log prior w.r.t. h
	}
			
	return (dOdtheta)
	
}

optimise.gp <- function(args) 
{
	# pass all the arguments of a call to gpFit, memoize and
	# optimise the model, and return afitted version
	
	# set optim to FALSE
	cat('opt: start\n')
	args[['opt.h']] <- FALSE
	args[['recursive']] <- TRUE

	# put all these functions into a closure of this one!
	# trochu bastlirna: vsechny tyto fce budou v closure teto fce
	# ale je to stale lepsi nez globalni promenne
	# cistsi by bylo udelat si new.env() a dat tam jen to, ale stejne 
	# ty fce musi byt robustni aby nepouzili globalni promennou...
	environment(gpFit) <- environment() # diky tomu jak je udelan memoise(), toto bude closurem i toho mem_gpFit!!! Coz potrebujem!
	environment(gpFitLaplace) <- environment()
	environment(objective) <- environment()
	environment(gradient) <- environment()
	#environment(args$gp$ll) <- globalenv()
	
	hyper_iter <- 0 # 
	iter_info <- data.frame()
#	iter_hyperpar <<- c()
#	iter_hyperpar_grads <<- c()
	iter_h <- c()
	iter_dmll_dh <- data.frame()
	iter_dP_dh <- c() 
	iter_dp_dh <- c() 
	f_ <- c() # matrix(NA, nrow = length(args$y), ncol = 1)
	f_start <- NULL

	# memoise the gpFit function 
	#print(memoise)
	#return()
	mem_gpFit <- memoise_seq_bastl(gpFit)
	#mem_gpFit <- gpFit

	
	# convert hyperparameters to optimization scale
	theta <- gpLink(args$gp, args$h)
		
	# optimisation arguments
	if (args$method == 'Laplace') {
		meth <- 'L-BFGS-B'
		grad <- gradient
	} else {
		stop("not implemented")
		meth <- 'BFGS'
		grad <- NULL
	}
		
	time_start <- proc.time()
	
	opt <- optim(theta,
				 fn = objective,
				 gr = grad,
				 prior.pars = args$theta.prior.pars,
				 isfac = isfac,
				 args = args,
				 fun = mem_gpFit,
				 hessian = args$hessian,
				 lower = gpLink(args$gp, gpHyperparExportVector(args$gp, "low")),
				 upper = gpLink(args$gp, gpHyperparExportVector(args$gp, "up")),
				 method = meth,
				 control = args$opt.control)
	
	# get the resultant hyperparameters
	h <- gpLinkInv(args$gp, opt$par)
	
	args$h <- h
	
	cat("\nAFTER optim():   fit the final model:\n")
	# fit the final model and return
	m <- do.call(mem_gpFit, args)
	
	# replace hessian with the hessian matrix or NULL
	if (args$hessian) {
		m$hessian <- opt$hessian
	} else {
		m$hessian <- NULL
	}
	
	#if (opt$convergence != 0)
	#	cat("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!!!!\n!!!!!   optim $convergence code:", opt$convergence, "\n\n")
	#cat("Before forget:\n")
	#print(gc())
	# un-memoize graf
	forget(mem_gpFit)
	#cat("After forget:\n")
	#print(gc())
	m$fit$iter_info <- iter_info
	m$fit$iter_h <- iter_h
	m$fit$iter_dmll_dh <- iter_dmll_dh
	m$fit$iter_dP_dh <- iter_dP_dh
	m$fit$iter_dp_dh <- iter_dp_dh
	m$fit$optim <- opt
	f_start <- NULL	
	return (m)
}
