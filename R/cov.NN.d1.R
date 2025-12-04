# cov.NN(x1, x2, l)
# Neural Network Covariance Function - best optimized version with just two matrix multiplications!
#
# x1 - matrix of covariates for the training dataset (rows - data records, cols - covariates)
# x2 - matrix of covariates for the predicted dataset (rows - data records, cols - covariates)
#
# argument variants:
# 1) no x2
# 3) with x2

# output: 
#	- covariance matrix of dimensions nrow(x1) x nrow(x1) (variants 1)
#	- covariance matrix of dimensions nrow(x1) x nrow(x2) (variants 3)
#
# internal version cov.NN.b.d1.R. Based on cov.NN.b(), computed in a similar manner.

cov.NN.d1 <- function(x1, x2 = NULL, sigma2_diag, der_wrt, der_wrt_i) # verze kterou interne pojmenovavam cov.NN.b.d1
{
	stopifnot(der_wrt == "sigma2_diag")
	stopifnot(der_wrt_i <= length(sigma2_diag) && der_wrt_i >= 1)
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

	JA <- sqrt(1 + 2 * xt1^2 %*% sigma2_diag)[,1]
	JB <- sqrt(1 + 2 * xt2^2 %*% sigma2_diag)[,1]
	
	# I decomposed it so that only one matrix multiplication is done. Main idea: prepare everything you need from x1, 
	# then from x2, and then perform one final matrix multiplication
	clen2 <- 2 * xt1 * Sigma.cols / JA
	clen3 <-     xt2              / JB

	# calculate the derivative of what is inside the asin() function:
	# here I also decomposed it in the same way as above, so that only one matrix multiplication is done
	d_clen2 <- cbind(2 * xt1 * Sigma.cols * (-xt1[,der_wrt_i]^2/JA^3), 2*xt1[,der_wrt_i]/JA, clen2)
	d_clen3 <- cbind(clen3, clen3[,der_wrt_i], xt2*(-xt2[,der_wrt_i]^2/JB^3)) #-clen3*clen3[,der_wrt_i]/JB
	K <- 2/pi * 1/sqrt(1 - (clen2 %*% t(clen3))^2)*(d_clen2 %*% t(d_clen3))
	K	
}

