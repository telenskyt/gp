# cov.NN(x1, x2, l)
# Neural Network Covariance Function - best optimized version with just one matrix multiplication!
# as defined in Rasmussen & Williams, 2006, eq. 4.29, with diagonal Sigma matrix
#
# x1 - matrix of covariates for the training dataset (rows - data records, cols - covariates)
# x2 - matrix of covariates for the predicted dataset (rows - data records, cols - covariates)
# sigma2_diag - diagonal of the Sigma matrix (contains "small sigmas"^2 on the diagonal), 
#		a vector of length ncol(x1)+1
#
# !!!2025 thoughts: meditating on whether it does need a global sigma2 multiplication (which is not part of this function)
#		- it is not necessary for zeroing the matrix, since zeros sigma2_diag will accomplish that, but it is needed 
#		for scaling it up, since maximum value in the resultant cov.NN() matrix is 1.
# 
#
# argument variants:
# 1) no x2
# 3) with x2

# output: 
#	- covariance matrix of dimensions nrow(x1) x nrow(x1) (variants 1)
#	- covariance matrix of dimensions nrow(x1) x nrow(x2) (variants 3)
#
# This is internal version cov.NN.b, the fastest, even faster than cov.NN.2b! See test_cov_NN.R
#
# cov.NN.b, 4851 sites, length(sigma2) = 96: 1.5s, lepsi!
# cov.NN.b, 4851 sites, length(sigma2) = 2: 0.59s, lepsi!

cov.NN <- function(x1, x2 = NULL, sigma2_diag)  # verze kterou pojmenovavam cov.NN.b
{
	stopifnot(ncol(x1) + 1 == length(sigma2_diag))
	xt1 <- as.matrix(cbind(1, x1)) # x tilde, augmented	
	if (is.null(x2)) {
		xt2 <- xt1
	} else {
		stopifnot(ncol(x2) + 1 == length(sigma2_diag))
		xt2 <- as.matrix(cbind(1, x2))
	}	
	#Sigma <- diag(sigma2_diag)
	
	Sigma.cols <- matrix(sigma2_diag, nrow = nrow(xt1), ncol = length(sigma2_diag), byrow = TRUE)

	# I decomposed it so that only one matrix multiplication is done. Main idea: prepare everything you need from x1, 
	# then from x2, and then perform one final matrix multiplication	
	clen2 <- 2 * xt1 * Sigma.cols / sqrt(1 + 2 * xt1^2 %*% sigma2_diag)[,1]
	clen3 <-     xt2              / sqrt(1 + 2 * xt2^2 %*% sigma2_diag)[,1]
	K <- 2/pi * asin(clen2 %*% t(clen3))	
	K	
}
