# |- KDE.2d (M0:KDE) ---------------

#' This measure represents M0:KDE (concentration measure)
#'
#' @param sn data sample of dimension n x d, with d = 2
#' @param H the bandwidth hyperparameter, in the form of a 2 x 2 matrix
#'
#' @return the density estimate with M0:KDE distance for each point in the sample
#' @export
#'
#' @examples
#' # Generate some bivariate data
#' R = matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2, byrow = TRUE)
#' draws_2d = MASS::mvrnorm(n = 2000, mu = c(0, 0), Sigma = R)
#' KDE.2d(sn = draws_2d)
KDE.2d = function(sn, H = NULL){

  sn = as.matrix(sn)
  n = dim(sn)[1]
  if(dim(sn)[2] != 2) stop("Sample data sn must be bidimensional, with the two dimensions on the columns")

  if(is.null(H)){
    # OPTIMAL bandwidth: Chacón, J.E., Duong, T., Wand, M.P. (2011), Asymptotics for General Multivariate Kernel Density Derivative Estimators, Statistica Sinica, 21, 807–840.
    H = (4/(dim(sn)[2]+4))^(1/(dim(sn)[2]+6))*dim(sn)[1]^(-1/(dim(sn)[2]+6))*abs(cov(sn))^(1/2)
  }

  fit_kde = ks::kde(sn, eval.points = sn, H = H) # from ks package
  return(dist = fit_kde$estimate)
}


# |- DE.NPCop.2d (M0-NPCop:DE) ---------------

#' This function computes M0-NPCop:DE (concentration measure)
#'
#' @param sn data sample of dimension n x d, with d = 2
#' @param method method used for estimating the copula with kde; default is "TLL2nn". See kdecopula package for more details: https://cran.r-project.org/web/packages/kdecopula/vignettes/kdecopula.pdf
#'
#' @return the density estimate with M0-NPCop:DE distance for each point in the sample
#' @export
#'
#' @examples
#' # Generate some bivariate data
#' R = matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2, byrow = TRUE)
#' draws_2d = MASS::mvrnorm(n = 2000, mu = c(0, 0), Sigma = R)
#' DE.NPCop.2d(sn = draws_2d)
DE.NPCop.2d = function(sn, method = "TLL2nn"){

  sn = as.matrix(sn)
  n = dim(sn)[1]
  if(dim(sn)[2] != 2) stop("Sample data sn must be bidimensional, with the two dimensions on the columns")

  U = copula::pobs(sn)
  fit_cop = kdecopula::kdecop(U, method = method)

  fit_x <- kde1d::kde1d(sn[,1]) # estimate marginal 1
  fit_y <- kde1d::kde1d(sn[,2]) # estimate marginal 2

  # combine copula and marginals
  npC_est = kdecopula::dkdecop(U, obj = fit_cop)*kde1d::dkde1d(sn[,1], fit_x)*kde1d::dkde1d(sn[,2], fit_y)

  return(dist = npC_est)
}


#' This function fits a parametric model for each of the two marginals and for the copula model
#'
#' @param sn data sample of dimension n x d, with d = 2
#' @param margin_family a two-dimensional vector specifying the two marginal families in the form 'norm', 't', or 'mixnorm'. No other specifications are allowed for the marginals.
#'
#' @return margin_family: a two-dimensional vector specifying the two marginal families; margin_params: a two-dimensional list specifying the estimated parameters of each of the two marginal families; copula_model: the estimated copula model
#' @export
#'
#' @examples
#' # Generate some bivariate data
#' R = matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2, byrow = TRUE)
#' draws_2d = MASS::mvrnorm(n = 2000, mu = c(0, 0), Sigma = R)
#' Pfit_margins_cop(sn = draws_2d)
Pfit_margins_cop = function(sn, margin_family = c('norm', 'norm')){

  sn = as.matrix(sn)
  n = dim(sn)[1]
  if(dim(sn)[2] != 2) stop("Sample data sn must be bidimensional, with the two dimensions on the columns")

  margin_params = list()
  U = matrix(NA, nrow = dim(sn)[1], ncol = dim(sn)[2])
  for(m in 1:length(margin_family)){
    if(margin_family[m]=="norm"){
      fit_m = MASS::fitdistr(sn[,m],"normal")
      margin_params[[m]] = list(mean = fit_m$estimate[1], sd = fit_m$estimate[2])
      U[,m] = pnorm(sn[,m], mean = margin_params[[m]]$mean, sd = margin_params[[m]]$sd)
    } else if (margin_family[m]=="t"){
      fit_m = MASS::fitdistr(sn[,m],"t", start = list(m=mean(sn[,m]),s=sd(sn[,m]), df=3), lower=c(-1, 0.001,1))
      margin_params[[m]] = list(df = fit_m$estimate[3])
      U[,m] = pt(sn[,m], df = margin_params[[m]]$df)
    } else if (margin_family[m]=="mixnorm"){
      fit_m = mclust::Mclust(sn[,m], modelNames = mclustModelNames("V"), G = 2)$param
      margin_params[[m]] = list(mean = fit_m$mean, sd = sqrt(fit_m$variance$sigmasq), pro = fit_m$pro)
      U[,m] = KScorrect::pmixnorm(sn[,m], mean = margin_params[[m]]$mean, sd = margin_params[[m]]$sd, pro = margin_params[[m]]$pro)
    } else {
      stop("Execution stopped. Only 'norm', 't', and 'mixnorm' specifications are allowed for the marginals")
    }
  }

  estCopula = VineCopula::BiCopSelect(U[,1], U[,2])
  estCopula2 = VC2copula::BiCop2copula(family = estCopula$family, par = estCopula$par, par2 = estCopula$par2)

  return(list(margin_family = margin_family, margin_params = margin_params, copula_model = estCopula, copula_model2 = estCopula2))
}

#' This function computes M0-PCop:DE (concentration measure)
#'
#' @param sn data sample of dimension n x d, with d = 2
#' @param margin_family a two-dimensional vector specifying the two marginal families in the form 'norm', 't', 'beta', or 'mixnorm'. No other specifications are allowed for the marginals.
#'
#' @return the density estimate with M0-PCop:DE distance for each point in the sample
#' @export
#'
#' @examples
#' # Generate some bivariate data
#' R = matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2, byrow = TRUE)
#' draws_2d = MASS::mvrnorm(n = 2000, mu = c(0, 0), Sigma = R)
#' DE.PCop.2d(sn = draws_2d)
DE.PCop.2d = function(sn, margin_family = c('norm', 'norm')){

  sn = as.matrix(sn)
  n = dim(sn)[1]
  if(dim(sn)[2] != 2) stop("Sample data sn must be bidimensional, with the two dimensions on the columns")

  parametric_fit = Pfit_margins_cop(sn = sn, margin_family = margin_family)

  mvd_est = suppressMessages(copula::mvdc(copula = parametric_fit$copula_model2,
                                          margins = parametric_fit$margin_family,
                                          paramMargins = parametric_fit$margin_params))

  pC_est = copula::dMvdc(sn, mvdc = mvd_est)

  return(dist = pC_est)
}

# |- KNN.2d (M1:kNN-Eucl) ---------------

#' This function computes M1:kNN-Eucl in Definition 3 (sparsity measure)
#'
#' @param sn data sample of dimension n x d, with d = 2
#' @param k an integer constant with value 1 to n representing the number of neighbors to consider in the summation; by default the rule of thumb is used
#'
#' @return the M1:kNN-Eucl distance for each point in the sample
#' @export
#'
#' @examples
#' # Generate some bivariate data
#' R = matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2, byrow = TRUE)
#' draws_2d = MASS::mvrnorm(n = 2000, mu = c(0, 0), Sigma = R)
#' KNN.2d(sn = draws_2d)
KNN.2d = function(sn, k = NULL){

  sn = as.matrix(sn)
  n = dim(sn)[1]

  if(dim(sn)[2] != 2) stop("Sample data sn must be bidimensional, with the two dimensions on the columns")

  if(is.null(k)==T){
    k = round(sqrt(n*.5))
  } else {
    if(k < 1) stop("Parameter k must be a positive integer")
    if(k >= n ) warning("Parameter k should be smaller than the sample size n. \n The maximum value corresponding to k = n-1 neighbors is taken.")
    if(as.integer(k) != k) warning("Parameter k is not an integer. The floor value is considered")
  }

  dist = c()
  for(i in 1:n){
    x = sn[i,]
    if(length(x) != 2) stop("Data point x must be one single bidimensional point")

    all_dist = mapply(function(xk1, xk2) sqrt((x[1]-xk1)^2+(x[2]-xk2)^2), sn[,1], sn[,2])
    dist[i] = 1/sum(sort(all_dist)[2:(k+1)], na.rm = T)
  }

  return(dist)
}

# |- CDF.KNN.2d (M2:kNN-CDF) ---------------

#' This function computes M2:kNN-CDF distance in Definition 4 (sparsity measure)
#'
#' @param data sample of dimension n x d, with d = 2
#' @param k an integer constant with value 1 to n representing the number of neighbors to consider in the summation
#' @param cdf_x1 an estimate of the cdf of the first marginal. By default it is set to NULL and the empirical cdf is taken
#' @param cdf_x2 an estimate of the cdf of the first marginal. By default it is set to NULL and the empirical cdf is taken
#'
#' @return the M2:kNN-CDF distance for each point in the sample
#' @export
#'
#' @examples
#' # Generate some bivariate data
#' R = matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2, byrow = TRUE)
#' draws_2d = MASS::mvrnorm(n = 2000, mu = c(0, 0), Sigma = R)
#' CDF.KNN.2d(sn = draws_2d)
CDF.KNN.2d = function(sn, k = NULL, cdf_x1 = NULL, cdf_x2 = NULL){

  sn = as.matrix(sn)
  n = dim(sn)[1]
  if(dim(sn)[2] != 2) stop("Sample data sn must be bidimensional, with the two dimensions on the columns")

  if(is.null(k)==T){
    k = 30
  } else {
    if(k < 1) stop("Parameter k must be a positive integer")
    if(k >= n ) warning("Parameter k should be smaller than the sample size n. \n The maximum value corresponding to k = n-1 neighbors is taken.")
    if(as.integer(k) != k) warning("Parameter k is not an integer. The floor value is considered")
  }

  if(is.null(cdf_x1)){
    cdf_x1 = ecdf(sn[,1])
  }

  if(is.null(cdf_x2)){
    cdf_x2 = ecdf(sn[,2])
  }

  dist = c()
  for(i in 1:n){
    x = sn[i,]
    if(length(x) != 2) stop("Data point x must be one single bidimensional point")

    all_dist = mapply(function(xk1, xk2) sqrt((x[1]-xk1)^2+(x[2]-xk2)^2), sn[,1], sn[,2])
    idx = order(all_dist)[1:(k+1)]

    k_CDFdist = mapply(function(xk1, xk2) sqrt((cdf_x1(x[1])-cdf_x1(xk1))^2+(cdf_x2(x[2])-cdf_x2(xk2))^2), sn[idx,1], sn[idx,2])

    num = sort(all_dist)[2:(k+1)]
    den = k_CDFdist[2:(k+1)]

    dist[i] = 1/sum(num/den, na.rm = T)
  }

  return(dist)
}

# |- CDF.e.2d (M3:e-CDF) ---------------

#' This function computes M3:e-CDF distance in Definition 5 (concentration measure)
#'
#' @param sn data sample of dimension n x d, with d = 2
#' @param eps a bidimensional hyperparameter with positive elements, representing the half-length of the rectangular neighborhood
#'
#' @return the M3:e-CDF distance for each point in the sample
#' @export
#'
#' @examples
#' # Generate some bivariate data
#' R = matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2, byrow = TRUE)
#' draws_2d = MASS::mvrnorm(n = 2000, mu = c(0, 0), Sigma = R)
#' CDF.e.2d(sn = draws_2d)
CDF.e.2d = function(sn, eps = NULL){

  sn = as.matrix(sn)
  n = dim(sn)[1]
  if(dim(sn)[2] != 2) stop("Sample data sn must be bidimensional, with the two dimensions on the columns")

  if(is.null(eps)==T){
    eps = rep(exp(2.13-0.3*log(n)), 2)
  } else {
    if(any(eps <= 0)) stop("Parameter eps must be positive")
    if(length(eps) != 2 ) warning("Parameter eps should be bidimensional. \n The first element is used as the half-length of the square neighborhood.")
  }

  dist = c()
  for(i in 1:n){
    x = sn[i,]
    if(length(x) != 2) stop("Data point x must be one single bidimensional point")

    colnames(sn) = c("V1", "V2")
    a_CDF_mv = mltools::empirical_cdf(x=data.table::data.table(sn), ubounds = data.table::data.table(V1 = x[1]+eps[1], V2 = x[2]+eps[2]))$CDF
    b_CDF_mv = mltools::empirical_cdf(x=data.table::data.table(sn), ubounds = data.table::data.table(V1 = x[1]+eps[1], V2 = x[2]-eps[2]))$CDF
    c_CDF_mv = mltools::empirical_cdf(x=data.table::data.table(sn), ubounds = data.table::data.table(V1 = x[1]-eps[1], V2 = x[2]+eps[2]))$CDF
    d_CDF_mv = mltools::empirical_cdf(x=data.table::data.table(sn), ubounds = data.table::data.table(V1 = x[1]-eps[1], V2 = x[2]-eps[2]))$CDF

    dist[i] = a_CDF_mv - b_CDF_mv - c_CDF_mv + d_CDF_mv
  }

  return(dist)
}

# |- CDF.NPCop.2d (M3-NPCop:e-CDF) ---------------

#' This function computes M3-NPCop:e-CDF (concentration measure)
#'
#' @param sn data sample of dimension n x d, with d = 2
#' @param eps a bidimensional hyperparameter with positive elements, representing the half-length of the rectangular neighborhood
#' @param cdf_x1 an estimate of the cdf of the first marginal. By default it is set to NULL and the empirical cdf is taken
#' @param cdf_x2 an estimate of the cdf of the first marginal. By default it is set to NULL and the empirical cdf is taken
#' @param fit_cop a fitted copula model. By default it is set to NULL and is estimated with kde
#' @param method method used for estimating the copula with kde; default is "TLL2nn". See kdecopula package for more details: https://cran.r-project.org/web/packages/kdecopula/vignettes/kdecopula.pdf
#'
#' @return the M3-NPCop:e-CDF distance for each point in the sample
#' @export
#'
#' @examples
#' # Generate some bivariate data
#' R = matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2, byrow = TRUE)
#' draws_2d = MASS::mvrnorm(n = 2000, mu = c(0, 0), Sigma = R)
#' CDF.NPCop.2d(sn = draws_2d)
CDF.NPCop.2d = function(sn, eps = NULL, cdf_x1 = NULL, cdf_x2 = NULL, fit_cop = NULL, method = "TLL2nn"){

  sn = as.matrix(sn)
  n = dim(sn)[1]
  if(dim(sn)[2] != 2) stop("Sample data sn must be bidimensional, with the two dimensions on the columns")

  if(is.null(cdf_x1)){
    cdf_x1 = ecdf(sn[,1])
  }

  if(is.null(cdf_x2)){
    cdf_x2 = ecdf(sn[,2])
  }

  if(is.null(eps)==T){
    eps = rep(exp(1.74-0.26*log(n)), 2)
  } else {
    if(any(eps <= 0)) stop("Parameter eps must be positive")
    if(length(eps) != 2 ) warning("Parameter eps should be bidimensional. \n The first element is used as the half-length of the square neighborhood.")
  }

  if(is.null(fit_cop)==T){
    U = copula::pobs(sn)
    fit_cop = kdecopula::kdecop(U, method = method)
  }

  dist = c()
  for(i in 1:n){
    x = sn[i,]
    if(length(x) != 2) stop("Data point x must be one single bidimensional point")

    Ux_r = cdf_x1(x[1]+eps[1])
    Ux_l = cdf_x1(x[1]-eps[1])
    Uy_r = cdf_x2(x[2]+eps[2])
    Uy_l = cdf_x2(x[2]-eps[2])

    a_CDF_mv = kdecopula::pkdecop(c(Ux_r, Uy_r), obj = fit_cop)
    b_CDF_mv = kdecopula::pkdecop(c(Ux_r, Uy_l), obj = fit_cop)
    c_CDF_mv = kdecopula::pkdecop(c(Ux_l, Uy_r), obj = fit_cop)
    d_CDF_mv = kdecopula::pkdecop(c(Ux_l, Uy_l), obj = fit_cop)

    dist[i] = a_CDF_mv - b_CDF_mv - c_CDF_mv + d_CDF_mv
  }

  return(dist)
}

# |- CDF.PCop.2d (M3-PCop:e-CDF) ---------------

#' This function computes M3-PCop:e-CDF (concentration measure)
#'
#' @param sn data sample of dimension n x d, with d = 2
#' @param margin_family a two-dimensional vector specifying the two marginal families in the form 'norm', 't', or 'mixnorm'. No other specifications are allowed for the marginals.
#' @param eps a bidimensional hyperparameter with positive elements, representing the half-length of the rectangular neighborhood
#'
#' @return the density estimate with M0-PCop:DE distance for each point in the sample
#' @export
#'
#' @examples
#' # Generate some bivariate data
#' R = matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2, byrow = TRUE)
#' draws_2d = MASS::mvrnorm(n = 2000, mu = c(0, 0), Sigma = R)
#' CDF.PCop.2d(sn = draws_2d)
CDF.PCop.2d = function(sn, margin_family = c('norm', 'norm'), eps = NULL){

  sn = as.matrix(sn)
  n = dim(sn)[1]
  if(dim(sn)[2] != 2) stop("Sample data sn must be bidimensional, with the two dimensions on the columns")

  if(is.null(eps)==T){
    eps = rep(exp(1.60 - 0.41*log(n)), 2)
  } else {
    if(any(eps <= 0)) stop("Parameter eps must be positive")
    if(length(eps) != 2 ) warning("Parameter eps should be bidimensional. \n The first element is used as the half-length of the square neighborhood.")
  }

  dist = c()
  for(i in 1:n){
    x = sn[i,]
    if(length(x) != 2) stop("Data point x must be one single bidimensional point")

    parametric_fit = Pfit_margins_cop(sn = sn, margin_family = margin_family)

    U_r = U_l = c()
    for(m in 1:length(margin_family)){
      if(margin_family[m]=="norm"){
        U_r[m] = pnorm(x[m]+eps[m], mean = parametric_fit$margin_params[[m]]$mean, sd = parametric_fit$margin_params[[m]]$sd)
        U_l[m] = pnorm(x[m]-eps[m], mean = parametric_fit$margin_params[[m]]$mean, sd = parametric_fit$margin_params[[m]]$sd)
      } else if (margin_family[m]=="t"){
        U_r[m] = pt(x[m]+eps[m], df = parametric_fit$margin_params[[m]]$df)
        U_l[m] = pt(x[m]-eps[m], df = parametric_fit$margin_params[[m]]$df)
      } else if (margin_family[m]=="mixnorm"){
        U_r[m] = KScorrect::pmixnorm(x[m]+eps[m], mean = parametric_fit$margin_params[[m]]$mean, sd = parametric_fit$margin_params[[m]]$sd, pro = parametric_fit$margin_params[[m]]$pro)
        U_l[m] = KScorrect::pmixnorm(x[m]-eps[m], mean = parametric_fit$margin_params[[m]]$mean, sd = parametric_fit$margin_params[[m]]$sd, pro = parametric_fit$margin_params[[m]]$pro)
      } else {
        stop("Execution stopped. Only 'norm', 't', and 'mixnorm' specifications are allowed for the marginals")
      }
    }

    a_CDF_mv = VineCopula::BiCopCDF(U_r[1], U_r[2], parametric_fit$copula_model)
    b_CDF_mv = VineCopula::BiCopCDF(U_r[1], U_l[2], parametric_fit$copula_model)
    c_CDF_mv = VineCopula::BiCopCDF(U_l[1], U_r[2], parametric_fit$copula_model)
    d_CDF_mv = VineCopula::BiCopCDF(U_l[1], U_l[2], parametric_fit$copula_model)

    dist[i] = a_CDF_mv - b_CDF_mv - c_CDF_mv + d_CDF_mv
  }

  return(dist)
}
