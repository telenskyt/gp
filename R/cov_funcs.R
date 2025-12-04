
# todo:
# - zautomatizovat len, 
cov_funcs <- list(
	"1" = list( # !!! resit toto pomoci specialni entry? Mit na to spec fci? Ono vlastne na to nejaka funkcionalita treba je.
				# toto je matice J
		hyperpar = list(),
		scaling = NULL
	),
	cov.SE = list(
		hyperpar = list(
			ls = list(
				#len = "ncol",
				start = 200,
				low = 0.1,
				up = 200,
				optim.link = "cubicrootlog",
				prior = quote(lasso_ls(x, magn = 3 * 2.077589)) # env.ls: 10, date/time ls: 3
				#prior.lambda = 
			)
		),
		scaling = TRUE # scaling should be performed on the matrix
	),
	cov.Matern = list(
		hyperpar = list(
			ds = list( # distance-scale - like a length-scale (not squared), but only one dimension 
				#len = 1,
				start = 10, # in km
				low = 2.5, # in km; the point at which covariance is 0.44
				up = 1000, # in km
				optim.link = "log",
				prior = quote(ls_pcmatern_lp(x, lambda = 90, d = 2))
				#prior.lambda = 
			)
		),
		scaling = FALSE # scaling doesn't have to be performed on the matrix (= will not be by default)
	),
	cov.I = list(
		hyperpar = list(),
		scaling = NA # there is no matrix here
	),
	cov.I.factor = list(
		hyperpar = list(),
		scaling = NA # there is no matrix here
	),
	cov.NN = list(
		hyperpar = list(
			sigma2_diag = list(
				#len = "ncol",
				start = 1, # ??? tady nevim kde zacit? Davam kompromis mezi slope a intercept u cov.NN.add
				low = 1e-7,
				up = 100,
				optim.link = "log",
				prior = quote(sigma2_exp_lp(x, lambda = 1))
				#prior.lambda = 
			)
		),
		scaling = TRUE
	),
	cov.NN.add = list(
		hyperpar = list(
			sigma2_int = list(
				#len = "ncol",
				start = 0.1,
				low = 1e-7,
				up = 100,
				optim.link = "log",
				prior = quote(sigma2_exp_lp(x, lambda = 0.5))
				#prior.lambda = 
			),
			sigma2_slope = list(
				#len = "ncol",
				start = 10,
				low = 1e-7,
				up = 100,
				optim.link = "log",
				prior = quote(sigma2_exp_lp(x, lambda = 1.5))
				#prior = 
				#prior.lambda = 
			)			
		),
		scaling = TRUE
	)
)
