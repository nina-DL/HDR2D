# Main function: HDR.2d --------------------

#' This function computes an HDR for an observed sample, given a certain coverage probability
#'
#' @param sn data sample of dimension n x d, with d = 2
#' @param measure the measure to be used for computing the HDR. Should be one of the following: "KDE", "DE.NPCop", "DE.PCop", "KNN", "CDF.KNN", "CDF.e", "CDF.NPCop", "CDF.PCop".
#' @param coverage_prob the coverage probability in (0,1). By default, this is 0.95
#' @param build_plot logical, stating whether to plot the HDR classification or not. By default, this is T
#' @param k an integer constant with value 1 to n representing the number of neighbors to consider. By default it is = NULL and the rule of thumb is used.
#' @param H the bandwidth hyperparameter, in the form of a 2 x 2 matrix. By default it is = NULL and the optimal bandwidth of Chacón, J.E. et al. (2011) is used.
#' @param eps a bidimensional hyperparameter with positive elements, representing the half-length of the rectangular neighborhood. By default it is = NULL and the rule of thumb is used.
#' @param method method used for estimating the copula with kde; default is "TLL2nn" (see kdecopula package for more details).
#' @param margin_family a two-dimensional vector specifying the two marginal families in the form 'norm', 't', or 'mixnorm'. No other specifications are allowed for the marginals. Default is c('norm', 'norm').
#' @param cdf_x1 an estimate of the cdf of the first marginal. By default it is set to NULL and the empirical cdf is taken
#' @param cdf_x2 an estimate of the cdf of the first marginal. By default it is set to NULL and the empirical cdf is taken
#' @param fit_cop a fitted copula model. By default it is set to NULL and it is estimated with kde.
#' @param xlab name of the x axis
#' @param ylab name of the y axis
#' @param cex points dimension
#' @param col points color
#' @param pch points type
#' @param ...
#'
#' @return   in_out: a vector of 0/1 indicating which data point is inside (1) or outside (0) the HDR; HDR_in: the set of data points from sn classified as HDR data points; HDR_out: the set of data points from sn classified as non-HDR data points
#' @export
#'
#'
#' @examples
#' # Generate some bivariate data
#' R = matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2, byrow = TRUE)
#' draws_2d = MASS::mvrnorm(n = 2000, mu = c(0, 0), Sigma = R)
#' # Estimate the HDR with KDE
#' HDR.2d(sn = draws_2d, measure = "KDE", coverage_prob = 0.95, build_plot = T)

HDR.2d = function(sn, measure = "KDE", coverage_prob = 0.95, build_plot = T,
                  k = NULL, H = NULL, eps = NULL, method = "TLL2nn",
                  margin_family = c('norm', 'norm'),
                  cdf_x1 = NULL, cdf_x2 = NULL, fit_cop = NULL,
                  xlab = expression(X[1]), ylab = expression(X[2]),
                  cex=0.5, col = "darkmagenta", pch = 4, ...){

  n = nrow(sn)

  hdr_mv = data.frame(Order = 1:n, X = sn[,1], Y = sn[,2])
  hdr_mv = hdr_mv[order(hdr_mv$Order),]

  if(measure == "KDE"){
    est_density = KDE.2d(sn, H = H)
  } else if(measure == "DE.NPCop"){
    est_density = DE.NPCop.2d(sn, method = method)
  } else if(measure == "DE.PCop"){
    est_density = DE.PCop.2d(sn, margin_family = margin_family)
  } else if(measure == "KNN"){
    est_density = KNN.2d(sn, k = k)
  } else if(measure == "CDF.KNN"){
    est_density = CDF.KNN.2d(sn, k = k, cdf_x1 = cdf_x1, cdf_x2 = cdf_x2)
  } else if(measure == "CDF.e"){
    est_density = CDF.e.2d(sn, eps = eps)
  } else if(measure == "CDF.NPCop"){
    est_density = CDF.NPCop.2d(sn, eps = eps, cdf_x1 = cdf_x1, cdf_x2 = cdf_x2, fit_cop = fit_cop)
  } else if(measure == "CDF.PCop"){
    est_density = CDF.PCop.2d(sn, eps = eps, margin_family = margin_family)
  } else {
    stop("The measure was not specified or has an incorrect specification.\n Please specify one of the following: KDE, DE.NPCop, DE.PCop, KNN, CDF.KNN, CDF.e, CDF.NPCop, CDF.PCop")
  }

  hdr_mv$fxy = est_density

  # Neighborhood-quantile method
  cutoff_level = quantile(hdr_mv$fxy, probs = (1-coverage_prob), na.rm = T)

  hdr_mv$ok = as.integer(hdr_mv$fxy>cutoff_level)

  # plot parameters
  hdr_mv$col = sapply(hdr_mv$ok, FUN = function(x, col) ifelse(x==0, col, "cadetblue"), col = col)
  hdr_mv$pch = sapply(hdr_mv$ok, FUN = function(x, pch) ifelse(x==0, pch, 20), pch = pch)

  if(build_plot == T){
    HDR_title = bquote("Estimated "*.(coverage_prob*100)*"% HDR with "*.(measure))
    HDR_plot = plot(hdr_mv$X, hdr_mv$Y,
                    xlim=c(min(hdr_mv$X), max(hdr_mv$X)), ylim=c(min(hdr_mv$Y), max(hdr_mv$Y)),
                    cex=cex, pch = hdr_mv$pch, col=hdr_mv$col,
                    main = HDR_title, bty="n", xlab = xlab, ylab = ylab, ...)
  } else {HDR_plot = NA}

  return(list(in_out = hdr_mv$ok, HDR_in = sn[hdr_mv$ok==1,], HDR_out = sn[hdr_mv$ok==0,], HDR_plot = HDR_plot))
}
