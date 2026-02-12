# gpFitLaplace - find optimum of the GP for given set of hyperparameters using Laplace approximation
# ???? proc f_start neni parametr???
#
# mn - mean vector, on GP scale. If NULL, mnfun() is used.
# h - vector of hyperparameters that are being optimized (not fixed), optim.link funkci musi resit volajici, zde se uz okolo tohoto nic neresi
# wt - weights. !!! They are not implemented properly now!!! In good design they should be parameters to d0, d1, d2, etc. functions
# grad.computation - TRUE (computes all), FALSE (computes none), vector of indices - computes just those
# num.correct.W.tol - numerical correction of W - if not NULL, values of W (hessian diagonal) will be corrected for small negative numbers, within given tolerance
#
# Value:
# h_grads jsou: d marginal log likelihood / d h (!!!! pozor neni to derivace mnll, ale -mnll!!!!!!)
#
# changelog - previous modifications of graf::gpFitLaplace() by Tomas Telensky: 
# v8: 
#	- v h_grads computation jsem nechal jen jedno volani gc() hned za rm(dK_i) (ve v7 to byl overkill)
#	- 122 vars, h_grads computation time: 62.890s, total time: 320.25, mem: 1097.4Mb.
# v7:
#	- reducing memory usage: 
#		- sikovne rm() nekterych promennych, rozebrani expressions na sub-expressions a vsude prospikovano volanim gc() - to dela velky rozdil!!!
#	- 122 vars, h_grads computation time: 86.120s, total time: 343.23, mem: 1097.4Mb!!!! (druhy beh na 1330.0Mb, 86.7s, 343s) - tuto pamet to melo i jen s 8 promennymi
# v6:
#	- trying to reduce memory usage
#	- 122 vars, h_grads computation time: 60.960s, total time: 314.00s, mem: 3114.6Mb.
# v4:
#	- merged the computation of d2_i with the next equation
# 	- 122 vars, h_grads computation time: 68.280s, total time: 322.62; mem: 3314.2Mb.
# v3:
#	- optimization: totez do jednoho %*% !!!
#	- 122 vars, h_grads computation time: 87.840s, total time: 342.91, ten jeden optimalizovany radek: 0.14s !!
# v2:
#	- optimization: as.matrix(dist()) into %*% and matrix()
#	- 122 vars, h_grads computation time: 152.230s, total time: 410.67, ten jeden optimalizovany radek: 0.59s
# v1: 
#	- merged in the cov.SE.d1() code to save memory requirements for large number of covariates
#	- warning instead of print()
#	- 122 vars, h_grads computation time: 243.670s, total time: 501.23s, ten jeden optimalizovany radek: 1.24s

#
# wow! Here are derived the gradients of marginal likelihood for Likelihood hyperparameters! (https://math.stackexchange.com/q/3207960/15731)
# Groot, P., Peters, M., Heskes, T., & Ketter, W. (2014). Fast Laplace Approximation for Gaussian Processes with a Tensor Product Kernel.
#

gpFitLaplace <- function (gp, h, mn = NULL, wt = 1, e = NULL, tol = 10 ^ -6, itmax = 50,
            verbose = FALSE, use_f_start = FALSE, grad.computation = TRUE, mem_verbose = FALSE, num.correct.W.tol = 10*sqrt(.Machine$double.eps)) 
{
	print(match.call())
	wt <- 1 # currently weights are not supported
	x <- gp$data # alias for **scaled** data; but there is also y in there....
	y <- gp$data
options(scipen = 9)
#options(warn = 2) # just for debug
#cat("  ... (pro kontrolu hyper_iter = ", hyper_iter, ")\n") # tady to vychazi ze promenna hyper_iter cisluje od nuly, uz to tak necham
cat("Fitting hyperparameters: ")
print(h)	
	#hyperpar <- hyperpar_decode(h)
	gp$hyperpar <- gpHyperparImportVector(gp, h); 
	hyperpar <- gpHyperparList(gp)
	if (is.null(mn))
		mn <- mnfun(gp, hyperpar = hyperpar)	
print(gc())
	mstart(id = "graf_fit", mem_precise = TRUE)
	# create the covariance matrix
	#K <- cov.SE(x1 = x, e1 = e, e2 = NULL, l = l)
	K.cache <- K_cache(gp, hyperpar, x, cache_for = "derivative")
	K <- K_matrix(gp, hyperpar, x, K.cache = K.cache, comp_means = TRUE)

	n <- nrow(K) # size of the GP matrix, length of the whole GP "f" vector

	# an identity matrix for the calculations (oh no please)
	#eye <- diag(n) 
	if (use_f_start && exists("f_start") && !is.null(f_start) && hyper_iter >= 2) { # myslim ze dava smysl zvednout iteraci, od ktere se pouziva f_start...
		# hyper_iter podminka: ono to failovalo kdyz se fci graf() dal l.start a bylo use_f_start = TRUE, nekdy to slo do uplne numerickejch chyb
		# tak tam davam tuto podminku ze pouziju f_start az kdyz se to trochu stabilizuje
		# !! pozor, z nejakeho duvodu promenna hyper_iter tady bude cislovana od 0 (= 0 v 1. iteraci)
		cat("..using f_start\n")
		f <- f_start$f # !!!  nevim jestli to bude fungovat!	Podle me je potreba s tim se privest i `a` a obj.
		#a <- f_start$a
		obj <- Inf # proste nevime, stare `a' nelze pouzit, je nutne ho spocist znovu, stejne jako obj!
		a <- rep(0, n)
		a_doesnt_correspond_with_f <- TRUE
		#f <- mn
		#obj <- -sum(wt * d0(f, y))		
		using_f_start_now <- TRUE
	} else {
		# initialise (bacha tento kod je jeste duplikovan nize!)
		a <- rep(0, n)
		f <- mn
		a_doesnt_correspond_with_f <- FALSE
		obj <- -sum(wt * d0(gp, f, y, hyperpar))
		using_f_start_now <- FALSE
	}
	obj.old <- Inf
gc()
	# Algorithm 3.1: Mode-finding in Rasmussen & Williams 2006, pg 46
	# start newton iterations
	# optimized value obj is -psi from Rasmussen and Williams, and it is minimized here (see also psi() function for more comments on this).
	# -obj = psi(f) = log p(y|f) + log p(f|X) coz az na konstantu odpovida p(f|X,y), coz je posterior pro latent f 
	#       p(y|f) je likelihood a p(f|X) je gaussian prior.
	it <- 0
	LAiter <- c()
	f_start_was_reset <- FALSE
	while ((obj.old == Inf || obj < obj.old - tol) && it < itmax) { 
		cat("Iteration ", it+1, "...\n")
		#print(gc())
		it <- it + 1
		obj.old <- obj
		W <- -(wt * d2(gp, f, y, hyperpar))
		
		if (any(W < 0) && !is.null(num.correct.W.tol) && is.finite(num.correct.W.tol)) { # 2026 osetreni - jen drobne numericke chybky muzem zarovnat na 0, s ohledem na komentar nize
			W[W < 0 & W >= -num.correct.W.tol] <- 0
			warning("corrected small negative values of W (within tolerance num.correct.W.tol)")
		}
			
		if (any(W < 0)) {
			# (2026 note: going through the log files of my models, it looks like this didn't actually work. 
			# It happened in q_extend models, which didn't work out, only once it happened in normal model (O-tsc,sitesm/Picus_viridis), resulting in an error and core dump anyways.)
			stop("some W < 0, i.e. the diagonal of the hessian has negative values (range: ", paste(range(W), collapse=" to "),
				"),\n\tmeaning the hessian of the neg. log likelihood is not positive definite.\n\tIncreasing num.correct.W.tol can help, but use with caution!")
			
			# moje bastloidni osetreni!
			# zkusim remedy!! Posunout se po tech f tak, abych byl v te concave casti! 
			# (nemusi pak fungovat approx. marginal likelihood ale... ale to zatim neresim...)
			
			# debug using likelihood_concavity.R

				if (all(sign(d1(gp, f, y, hyperpar))[W < 0] < 0)) # @@@@!!!! this test only makes sense for TSC likelihood... remove later
					cat("according to expectation, all(sign(d1(f, y, hyperpar))[W < 0] < 0)\n")
				else 
					cat("!!! interesting! some sign(d1(f, y, hyperpar))[W < 0] is >= 0!\n")


			f0 <- as.vector(f)
			W0 <- W
			if (1) {  # lepsi remedy! Najdi hodnotu f takovou, aby W = 0!!
				for (ind in which(W < 0)) { 
					f[ind] <- uniroot(function (fi, i) { f[i] <- fi ; -d2(gp, f, y, hyperpar)[i] }, interval = c(-20,20), i = ind)$root
				}
				W <- -(wt * d2(gp, f, y, hyperpar))
				if (all(f <= f0)) { # @@@@!!!! this test only makes sense for TSC likelihood... remove later
					cat("according to expectation, all(f <= f0)\n")
				} else
					cat("!!! interesting! some f > f0 !\n")
				
			} else { # starsi remedy
				my_step <- 0
				while (any(W < 0) && my_step < 10) {
					# v mistech kde je W zaporne (druha der. kladna), se posunu o 1 smerem podle prvni derivace
					# je-li kladna, prictu 1; je-li zaporna, odectu 1
					f <- ifelse(W >= 0, f, f + 1 * sign(d1(gp, f, y, hyperpar)))
					if (all(sign(d1(gp, f, y, hyperpar))[W < 0] < 0)) # @@@@!!!! this test only makes sense for TSC likelihood... remove later
						cat("according to expectation, all(sign(d1(f, y, hyperpar))[W < 0] < 0)\n")
					else 
						cat("!!! interesting! some sign(d1(f, y, hyperpar))[W < 0] is >= 0!\n")
					
					# a ted prepocitej W:
					W <- -(wt * d2(gp, f, y, hyperpar))
					
					my_step <- my_step + 1
				}
			}
			if (any(W < 0))
				stop("Remedy for W < 0 didn't work")
			warning("Hyp iter ", if (exists("hyper_iter")) hyper_iter+1 else NA, ", LA iter ", it, ": log lik non-concave, bastl-remedy (in ", my_step, " steps)")
			if (0) {
			cat("f0:\n")
			print(f0)
			cat("f:\n")
			print(f)
			cat("The changes to `f` made:\n")
			print(f - f0)
			}
			chgs <- data.frame(f0 = f0, f = f, f_diff = f - f0, W0 = W0, W = W)
			cat("The changes to `f` made (max 12 shown):\n")
			print(head(chgs[chgs$f_diff != 0,], 12))
			a_doesnt_correspond_with_f <- TRUE
		}
		rW <- sqrt(W)
		cf <- f - mn
		#L <- tryCatch(chol(rW %*% t(rW) * K + diag(n)),
		#			error = function(x) return(NULL)) # nechapu proc tryCatch() proboha!! Ja o te chybe potrebuju vedet!
		L <- chol(rW %*% t(rW) * K + diag(n))
			# caution! This L is a transpose of the L in R&W2006. Here, the L <- chol(B) means t(L) %*% L == B, in R&W2006 it is L %*% t(L) == B
			# but I've checked this code of finding mode and also the predictions and it's used correctly (it takes this into account). -- Tomas 02/2020
		b <- W * cf + wt * d1(gp, f, y, hyperpar)
		a_new_proposed <- b - rW * backsolve(L, forwardsolve(t(L), rW * (K %*% b)))
		rm(L) # jeste vice usetri pameti!
gc()
		if (a_doesnt_correspond_with_f) {
			a <- as.vector(a_new_proposed)
			cat("a_doesnt_correspond_with_f, doing just simple method\n")

		} else {
			adiff <- a_new_proposed - a 
			dim(adiff) <- NULL
			# find optimum step size using Brent's method
			res <- optimise(psiline, c(0, 2), gp, adiff, a, as.matrix(K), y, mn, wt, hyperpar)
			a <- a + res$minimum * adiff # toto `a` je takovy kouzelny vektor... K %*% a dá f a Kx (pro predikce) dá f predikcí!
		}
		f <- K %*% a + mn
		a_doesnt_correspond_with_f <- FALSE
		obj <- psi(gp, a, f, mn, y, wt, hyperpar)
		LAiter <- rbind(LAiter, data.frame(LAiter = it, obj = obj, using_f_start_now = using_f_start_now))
		cat("\tminimized obj = -psi = ", obj, "\n")
		if (!is.finite(obj)) {	
			# tato vec se v jednom pripade stala (obj = Inf) kvuli tomu, ze f melo v jednom miste extremni hodnotu (cca 200) 
			# a d0(f,y) ve fci psi() vratil -Inf
			# kdyz se vypne f_start, uz se to nedeje; a to i v jinych pripadech; Pustime tedy cely cyklus znova, bez f_start
			warning("Hyperpar iter", hyper_iter+1, ", LA iter:", it, ": obj = ", obj, " !!!!")
			message("LA iterations up to now:\n")
			write(capture.output(print(LAiter)), stderr()) # super clumsy but that's R sometimes...			
			if (using_f_start_now) { 
				# vypni f_start a inicializuj jako by nebyl
				warning(" f_start was TRUE; resetting it to FALSE and restarting the whole LA optimisation") # this is OK, no problem, klidne casem vyhodit z warningu tuto vetev
				# initialise (bacha duplicitni kod s kodem vyse!)
				a <- rep(0, n)
				f <- mn
				a_doesnt_correspond_with_f <- FALSE
				obj <- -sum(wt * d0(gp, f, y, hyperpar))
				using_f_start_now <- FALSE

				obj.old <- Inf
				f_start_was_reset <- TRUE
			} else { # f_start nebyl nastaven
				# nevime, zda to vubec nastavalo; tento pripad zatim vubec nebyl nalezen
				# zde je treba vyhodit chybu, protoze by to stejne spadlo ve funkci optim
				# a pokud toto bude problem, bude asi potreba to resit vice principialneji a delat nejaky numericky osetreni ve fci psi()
				options(error = recover) # tady chci debug rovnou :-)
				stop("gpFitLaplace: -PSI = obj is not finite:", obj, ", and f_start not used")
				# bohuzel spousta promennych co pouzivam jako globalnich jsou jen v environmentu v optimise.graf()
			}
		} else if (obj > obj.old + tol) {
			cat("\n!!!!!!!!!!!!!!\nobj > obj.old + tol !! obj.old =", obj.old, ", obj =", obj, "\n!!!!!!!!!!!!!!!\n\n")
			# obj se nezmensil ale naopak zvetsil (zhorsil)
			# v tom pripade toto bude posledni iterace, dle soucasne definice while; nebylo by ale lepsi vzit si teda tu predchozi, 
			# kdyz se to zhorsilo?
			# da se toto oznacit za uspesnou konvergenci?
			
			# toto se stavalo kdyz obj vysel Inf, ale i s konecnymi hodnotami a tam je to nejspis v pohode (koukal jsem na priklad 
			# kdy to jen malicko stouplo a to muze byt znak toho ze jsme dokonvergovali a jsme na minimu a min uz to nejde)
			# ale je potreba se podivat na to do logu, jestli to nekde nebude delat problem...
			warning("Hyperpar iter", hyper_iter+1, ", LA iter:", it, ": obj - obj.old = ", obj - obj.old, " > tol !!!")
			message("LA iterations up to now:\n")
			write(capture.output(print(LAiter)), stderr()) # super clumsy but that's R sometimes...
			# v iteracich se to pozna tak ze lastLAObjDiff bude > 0
		}
		
		#print(gc())
	}
  
	# recompute key components
	cf <- f - mn
	W <- -(wt * d2(gp, f, y, hyperpar))
	if (any(W < 0) && !is.null(num.correct.W.tol) && is.finite(num.correct.W.tol)) { # 2026 osetreni - jen drobne numericke chybky muzem zarovnat na 0, s ohledem na komentar nize
		W[W < 0 & W >= -num.correct.W.tol] <- 0
		warning("corrected small negative values of W (within tolerance num.correct.W.tol)")
	}	
	if (!all(W >= 0)) {
		# set options(error = recover) and debug it using likelihood_concavity.R
		stop("even in the optimum f, some W < 0, i.e. the diagonal of the hessian has negative values (range: ", paste(range(W), collapse=" to "),
			"),\n\tmeaning the hessian of the neg. log likelihood is not positive definite, meaning log likelihood isn't concave function in the optimum f.",
			"\n\tIncreasing num.correct.W.tol can help, but use with caution!")
	}
	mstart(id = "L")	
	rW <- sqrt(W)
	xx <- rW %*% t(rW)# * K + diag(n)
	#xx <- rW %*% t(rW) * K + diag(n)
gc()
    xx <- xx * K
gc() # ta vlozena gc() tady jsou velmi dulezita!!!
	xx <- xx + diag(n)
gc()
    L <- tryCatch(chol(xx),
                  error = function(x) return(NULL)) # see the note on L above (it's a transpose of L in R&W2006, but it's OK)
	rm(xx)
	cat("Computing L matrix took ")	
	mstop(id = "L")
gc()

	# `a` bylo spocitano relativne slozite, ale v tuto chvili doiterovalo do hodnoty d1(f, y, hyperpar) ! Viz (3.17), (3.21) a alg. 3.2 v R&W2006 :)

    # return marginal negative log-likelihood, (3.32) in R&W 2006 (pg 48), see alg. 3.1 (pg 46)
    mnll <- (a %*% cf)[1, 1] / 2 + sum(log(diag(L))) - sum((wt * d0(gp, f, y, hyperpar))) 
		# (a %*% cf)[1, 1] = t(a) %*% cf ... otestoval jsem ze to tak je
gc()
	cat("gpFitLaplace(): fit took ")
	mstop(id = "graf_fit")

    # vector to store gradients
    h_grads <- rep(NA, length(h))

	grad.computation.idx <- 1:length(h)
	if (is.numeric(grad.computation)) {
		grad.computation.idx <- grad.computation
		grad.computation <- TRUE
	}
if (grad.computation) {
    # Get 1) partial gradients of the marginal LL wrt kernel hyperparameters h (Algorithm 5.1 in Rasmussen & Williams 2006, pg 126)
	# And newly also 2) gradients of marginal likelihood for likelihood hyperparameters! 
	# Groot, P., Peters, M., Heskes, T., & Ketter, W. (2014). Fast Laplace Approximation for Gaussian Processes with a Tensor Product Kernel.)
	# (O tom paperu jsem se dozvedel odsud: https://math.stackexchange.com/q/3207960/15731)
	# tam taky pak zpropagovat muj package, pokud ho udelam!!! Rict ze thanks, ze jsem to diky tomu implementoval a tady to je :-))

	cat("Computing gradients...\n")
	# pro analyzu vyuziti pameti:
	gc_rv <- gc()
	if (mem_verbose) {
		print(gc_rv)
		vardim <- t(sapply(ls(), function (x) if (is.null(dim(get(x)))) c(length(get(x)), 0, object.size(get(x))) else c(dim(get(x)), object.size(get(x)))))
		str(K.cache)
		print(memobj(K.cache))
		
		colnames(vardim) <- c("rows", "cols", "object.size")
		print(vardim)
	}
	
	mstart(id = "graf_grad", mem_precise = TRUE)
    
	mstart(id = "graf_grad_common", mem_precise = TRUE)
    # gradient components
    W12 <- matrix(rep(rW, n), n)
gc()
	fw <- forwardsolve(t(L), diag(rW)) 
gc()
    R <- W12 * backsolve(L, fw)
rm(fw)
gc()
    W12 <- W12 * K
gc()
    colSums_C2 <- colSums(forwardsolve(t(L), W12)^2)
    rm(W12)
gc()
	cat("gpFitLaplace(): common gradient computations took ")
	mstop(id = "graf_grad_common")

#	print(gc())
    # partial gradients of the kernel
    #dK <- cov.SE.d1(x, e, h)
	#gc1 <- gc()
	#t1 <- proc.time()
	#Rprof(line.profiling = TRUE, memory.profiling = TRUE)
	#K <- cov.SE(x1 = x, e1 = e, h = h) # ?? neni to totez co nahore? !!! ANO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    # rate of change of likelihood w.r.t. the mode
	diag_posterior_cov <- (diag(K) - colSums_C2) # diagonal of (K^-1 + W)^-1, the posterior covariance matrix
    s2 <- diag_posterior_cov / 2 * d3(gp, f, y, hyperpar)

	#KR <- K %*% R
	
    for (i in grad.computation.idx) {
		#cat("  computing marginal likelihood gradient for hyperparameter i = ", i, "\n")
		if (mem_verbose)
			mstart(id = "graf_grad_i", mem_precise = TRUE)

		if (gp$hyperpar[gpHyperparIdx(gp, i), "component"] != ".lik") { # kernel hyperparameter
			# gradient for each parameter h[i]: dK/dl[i] (this is what cov.SE.d1() computes, but optimized memory- and cpu- wise)
			dK_i <- dK_dhi(gp, hyperpar, x, i, K.cache)
		  
			s1 <- (t(a) %*% dK_i %*% a) / 2 - sum(R * dK_i) / 2 # s1 (eq 5.22) from Alg. 5.1
			b <- dK_i %*% d1(gp, f, y, hyperpar)
			rm(dK_i) 
gc()
			#grad <- s1 + t(s2) %*% (cbind(KR, b) %*% c(-b, 1)) # h_grads computation time: 160.560s; mem: 4198.2Mb. # 2x slower and more memory
			h_grads[i] <- s1 + t(s2) %*% (b - K %*% (R %*% b))
		} else { # likelihood hyperparameter! Finally possible!! Celebration!! This I implemented using formulas in Groot et al 2014! Otestovano a proslo!
			s1 <- d0_dhyp(gp, f, y, hyperpar, i) + t(diag_posterior_cov) %*% d2_dhyp(gp, f, y, hyperpar, i) / 2
			b <- K %*% d1_dhyp(gp, f, y, hyperpar, i)
			h_grads[i] <- s1 + t(s2) %*% (b - K %*% (R %*% b))
		}
		der_wrt <- gp$hyperpar[gpHyperparIdx(gp, i),] # which var we derivate w.r.t to :-)
		if (mem_verbose) {
			cat("gradient computation for hyperpar ", der_wrt$hyperpar, "(i = ", i, ") took ")
			mstop(id = "graf_grad_i")
		}
    }

	#Rprof(NULL)
	#t2 <- proc.time()
	#gc2 <- gc()
	#cat(sprintf("h_grads computation time: %.3fs\n", (t2 - t1)[3]))
	cat("gpFitLaplace(): whole gradients computation took ")
	mstop(id = "graf_grad")
	#print(gc())
	#cat(sprintf("h_grads computation mem: %.1fMb.\n", sum(gc2[,6] - gc1[,6])))
}
   
	# get local variable dimensions, to be able to read code better and maybe optimize it
	vardim <- t(sapply(ls(), function (x) if (is.null(dim(get(x)))) c(length(get(x)), 0) else dim(get(x))))
	colnames(vardim) <- c("rows", "cols")
	
    if(verbose ) cat(paste("  ", it, "Laplace iterations\n"))
	if (it == itmax) {
		warning("timed out, don't trust the inference!") # !!!! make this a warning, at least!
	}
	# vector h is already imported to gp$hyperpar
    fit <- list(h = h, f = f, a = a, W = W, L = L, K = K,
                e = e, mnll = mnll, wt = wt, psi = obj, tol = tol,
                h_grads = h_grads, vardim = vardim, iterations = it, itmax = itmax, 
				LAiter = LAiter, lastLAObjDiff = obj - obj.old, f_start_was_reset = f_start_was_reset)
	fit
}
  
# res <- optimise(psiline, c(0, 2), gp, adiff, a, as.matrix(K), y, d0, mn, wt, hyperpar)
psiline <- function(s, gp, adiff, a, K, y, mn = 0, wt, hyperpar) 
{
	a <- a + s * as.vector(adiff)
	f <- K %*% a + mn
	psi(gp, a, f, mn, y, wt, hyperpar)
}

# tato fce psi() dava -psi (psi s obracenym znamenkem) toho psi from Rasmussen & Williams 2006, 
# psi(f) = log p(y|f) + log p(f|X) coz odpovida p(f|X,y), az na integral p(y|X) coz je konstanta vuci f
# je to tedy vlastne posterior pro latent f! p(y|f) je likelihood a p(f|X) je gaussian prior!
# jen s obracenym znaminkem
# toto se minimalizuje v LA iteracich

psi <- function(gp, a, f, mn, y, wt, hyperpar) 
{
	0.5 * t(a) %*% (f - mn) - sum(wt * d0(gp, f, y, hyperpar))
}
