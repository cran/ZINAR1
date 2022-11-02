#' @title Parameter Estimation for ZINAR(1) Models
#'
#' @description This function uses the EM algorithm to find the maximum likelihood estimates of a ZINAR(1) model.
#'
#' @usage
#' EST_ZINAR(y,init = NULL,tol = 1e-05,iter = 1000,model,innovation,desc = FALSE)
#'
#' @param y A vector containing a discrete non-negative time series dataset.
#' @param init A vector containing the initial parameters estimates to maximize the likelihood function. If not informed, uses Yule-Walker method to calculate.
#' @param tol Tolerance for the convergence of the algorithm. Defaults to 1e-5.
#' @param iter Maximum number of iterations of the algorithm. Defaults to 1000.
#' @param model Must be "zinar", if the innovation have Zero-Inflated distribution, and "inar", otherwise.
#' @param innovation Must be "Po" if Poisson, "NB" if Negative binomial or "GI" if Gaussian inverse.
#' @param desc TRUE to plot the exploratory graphs. Defaults to FALSE.
#'
#'
#' @examples
#'
#' # Estimates the parameters of an INAR(1) and a ZINAR(1) models with Poisson innovations
#' # for the monthly number of drug offenses recorded from January 1990 to December 2001
#' # in Pittsburgh census tract 2206.
#'
#' data(PghTracts)
#'
#' y=ts(PghTracts$DRUGS,start=c(1990,1),end=c(2001,12),frequency=12)
#'
#' Inar1 = EST_ZINAR(y, init = c(0.3,0.5,2), model = "inar", innovation = "Po",desc = TRUE)
#'
#' ZIPInar1 = EST_ZINAR(y, init = c(0.3,0.5,2), model = "zinar", innovation = "Po",desc = TRUE)
#'
#' @return Returns a list containing the parameters estimates and the number of interactions.
#'
#' @references  Aldo M.; Medina, Francyelle L.; Jales, Isaac C.; Bertail, Patrice. First-order integer valued AR processes with zero-inflated innovations. Cyclostationarity: Theory and Methods, Springer Verlag - 2021, v. 1, p. 19-40.
#'
#' @export

EST_ZINAR = function(y,init=NULL,tol= 1e-5,iter = 1000,model,innovation,desc = FALSE){

  if(!(is.null(dim(y)))) stop("y must be a numeric vector or ts")
  if(length(y)<2) stop("length of y must be greather than 2")
  if(sum(y%%1!=0)>0 | sum(y<0)>0) stop("y must be a positive integer")
  if(tol < 0) stop("tolerance for the convergence of the algorithm must be positive")
  if(!(is.numeric(tol))) stop("the tolerance must be numeric")
  if(!(is.numeric(iter))) stop("the maximum number of iterations must be numeric")
  if(!(model %in% c("inar","zinar"))) stop("model must be set to 'inar' or 'zinar'")
  if(!(innovation %in% c("Po","NB","GI"))) stop("family must be set to Poisson ('Po'), Negative Binomial ('NB') or Gaussian Inverse ('GI')")
  if(iter%%1!=0|iter<1) stop("maximum number of iterations must be a positive integer")
  if(!(desc == TRUE | desc == FALSE)) stop("desc must be 'TRUE' to plot the exploratory graphs or 'FALSE' otherwise")


  if(is.null(init)){

    results = EM(y,initial = YW(y = y,model = model,family = innovation),tol = tol,int = iter,model = model,family = innovation)

    print(results$parameters)

  } else {

    if(innovation %in% c("NB","GI") & !(is.numeric(init))) stop("the initial values must be a numeric vector containing the initial parameters estimates of alpha, rho, mu and phi")
    if(innovation == "Po" & !(is.numeric(init))) stop("init must be a numeric vector containing the initial parameters estimates of alpha, rho and lambda")
    if(init[1]<0 | init[1]>=1) stop("initial value of the parameter alpha must be in [0,1)")
    if(init[2]<0 | init[2]>=1) stop("initial value of the parameter rho must be in [0,1)")
    if(innovation %in% c("NB","GI") & !(length(init)==4)) stop("init must be a vector containing the initial parameters estimates of alpha, rho, mu and phi")
    if(innovation == "Po" & !(length(init) == 3)) stop("init must be a vector containing the initial parameters estimates of alpha, rho and lambda")
    if(innovation == "Po" & init[3] < 0) stop("initial value of the parameter lambda must be positive")
    if(innovation %in% c("NB","GI") & (init[3] < 0 | init[4] < 0)) stop("initial values of the parameters mu and phi must be positive")

    results = EM(y,initial = init,tol = tol,int = iter,model = model,family = innovation)

    print(results$parameters)

  }

  if(desc == TRUE){

    explore_zinar1(y)

  }

  return(results)

}
