
# from original legacy code, from model_expand_predictions(): d:\tomas\gp\acr\vypocet\GRaF-my_modif\R\model\K-tsc_spat2.R 

############ Predictions
# simplest possible version!
#' Expand the predictions by quantities derived using the inverse link function
#'
#' Gets predictions of `f` from the Gaussian Process on the input, and from this computes also prediction for derived quantities,
#' in this case, using an inverse of given link function.
#'
#' @param pred data.frame or matrix of basic predictions to be expanded (`f`, etc.), as would be returned by predict.gp(type = "latent")
#' @param link character - name of link function to take inverse of - same as used by make.link()
#' @param parname character - name of the new parameter
#' @returns Expanded predictions.
#' @export
predict_expand_link <- function (pred, link = "logit", parname)
{
	stopifnot(is.character(link))
	stopifnot(is.character(parname) && nchar(parname) > 0)

	link_str <- make.link2(link)
	link <- link_str$linkfun
	inv.link <- link_str$linkinv
		
	pred_e <- pred[,colnames(pred) != "f_SE", drop = FALSE]
	colnames(pred_e) <- sub("^f", parname, colnames(pred_e)) 
	pred_e <- inv.link(pred_e)
	cbind(pred, pred_e)
}


# zatim bez pridani CI, SE!
predict_expand_fun <- function(gp, newdata, pred, pred_fun, parname)
{
	hyperpar <- gpHyperparList(gp)
	par <- c(hyperpar[[".lik"]], list(f = pred[,'f']))
	xx <- pred_fun(newdata, par)
	if (is.vector(xx)) {
		xx <- t(t(xx))
		colnames(xx) <- parname
	}
	cbind(pred, xx)
}
