# cov.SE(x1, x2, e1, e2, ls)
# x1 - matrix of covariates for the training dataset (rows - data records, cols - covariates)
# x2 - matrix of covariates for the predicted dataset (rows - data records, cols - covariates)
# e1, e2: (no longer supported) An optional matrix of standard deviations associated with x1, x2. If this is missing, covariates are assumed to be measured without error.
# ls - length-scales as defined by Rasmussen & Williams 2006, which are length-scales at the scale of the covariates
#
# original memory requirements: n1*n2*n3
# new memory requirements: n1*n2 (since v1)
#
# changelog:
# v3 - single matrix multiplication with added rows and cols (the trick I applied in gpFitLaplace.R as well); further speedups and memory savings for 1,3,4,5, see cov.SE.test_LOG.txt
#
# argument variants:
# 1) no x2, no e1/e2
# 2) no x2, with e1 (and e2 gets automatically initialized to e1) 
# 3) with x2, no e1/e2
# 4) with x2, with e1
# 5) with x2, with e1 and e2
#
# output: 
#	- covariance matrix of dimensions nrow(x1) x nrow(x1) (variants 1,2)
#	- covariance matrix of dimensions nrow(x1) x nrow(x2) (variants 3,4,5)
#
# at distance ls the covariance is 0.6065307 (at diagonal = 1)
cov.SE <- function(x1, x2 = NULL, e1 = NULL, e2 = NULL, ls) 
{
  n1 <- nrow(x1)
  n2 <- ifelse(is.null(x2), n1, nrow(x2))
  n3 <- ncol(x1)
  
#  print(dim(x1))
 # print(n1)
  #print(n2)
  #print(n3)
  
  stopifnot(all(ls > 0)) # not sure about this... test it!
  
  # with error matrices
  if (!is.null(e1)) { 
    # run through each covariate
	stop("no longer supported")
    
	if (is.null(x2)) {
		e2 <- e1
	}	
    sumdiffs <- 0
    denom <- 1
	ones1 <- t(rep(1, n2))
	ones2 <- t(rep(1, n1))
    for (i in 1:n3) {

	  E1 <- e1[, i] %*% ones1 # matrix n1 x n2: e1[row, i]  (all columns the same)
	  if (!is.null(e2)) {
	    E2 <- t(e2[, i] %*% ones2) # matrix n1 x n2: e2[col, i] (all rows the same)
	  } else {
		E2 <- 0
	  }

      err <- E1 + E2 # matrix n1 x n2: e1[row, i] + e2[col, i]
	  if(is.null(x2)) { # variant 2)
		dists_i <- dist(x1[, i]) ^ 2
		if (!exists("lower", inherits = FALSE))
			lower <- lower.tri(E1)		
        err <- err[lower] # save only lower portion for speed up	
		#dists_i <- cbind(x1[, i], x1[, i]^2, 1)  %*% rbind(-2*x1[, i], 1, x1[, i]^2) # !!! nechapu proc test tohoto selze, kdyz toto da stejny vysledek jako dist()^2
	  } else { # variants 4), 5)
		#dists_i <-   x1[, i] ^ 2 %*% t(rep(1, n2)) +
		#	rep(1, n1) %*% t(x2[, i] ^ 2) - 2 * x1[, i] %*% t(x2[, i]) # matrix n1 x n2: (x1[row, i] - x2[col, i])^2	  
		dists_i <- cbind(x1[, i], x1[, i]^2, 1)  %*% rbind(-2*x2[, i], 1, x2[, i]^2)
	  }
      sumdiffs <- sumdiffs + dists_i / (err + ls[i]^2)
      denom <- denom * (1 + err / ls[i]^2) # I don't know how to vectorize this; if we could vectorize this, we wouldn't need to compute the whole err matrix in each iteration
    }
    # inverse kronecker delta
    ikds <- as.numeric(sumdiffs > 0) # !!! here, what if you have two records with the same covariates, then you don't get the IKD!
    diag(ikds <- 1)
    denom <- sqrt(denom) * ikds
    K <- exp(-0.5 * sumdiffs) / denom
    
  } else {
	# without error matrices
	sumdiffs <- my_squared_dist(x1, x2, ls) 
    K <- exp(-0.5 * sumdiffs)  # to matrix?
  }
  
  if('dist' %in% class(sumdiffs)) {
    K <- as.matrix(K)
    diag(K) <- 1
  }
  stopifnot("matrix" %in% class(K))
  K
}
