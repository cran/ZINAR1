#' @title Simulate values for ZINAR(1)
#'
#' @description This function generates realizations of a ZINAR(1) process.
#'
#' @usage
#' SIM_ZINAR(n, alpha, rho, th, innovation)
#'
#' @param n Number of realizations of the ZINAR(1) process.
#' @param alpha The probability of an element remaining in the process. The parameter alpha must be in [0,1].
#' @param rho The probability of the innovation be from the state zero. The parameter rho must be in [0,1].
#' @param th Is equal the value of the parameter lambda, if the innovations follow a Zero-Inflated Poisson (ZIP) distribution, and is a vector containing the values of the parameters (mu,phi), if the innovations follow a Zero-Inflated Negative Binomial (ZINB) or Zero-Inflated Inverse Gaussian (ZIPIG) distribution.
#' @param innovation Must be "Po" if Poisson, "NB" if Negative binomial or "GI" if Gaussian inverse.
#'
#' @examples
#'
#' # Simulates values for ZIPInar1 model and estimate its parameters.
#'
#' set.seed(5)
#'
#' model = "zinar"
#' innv = "Po"
#' y = SIM_ZINAR(n = 500,alpha = 0.3,rho = 0.5,th = 3,innovation = innv)
#' ZIPInar1 = EST_ZINAR(y,model=model,innovation=innv,desc = TRUE)
#'
#'
#' @return Returns a numeric vector representing a realization of a ZINAR(1) process.
#'
#' @references  Aldo M.; Medina, Francyelle L.; Jales, Isaac C.; Bertail, Patrice. First-order integer valued AR processes with zero-inflated innovations. Cyclostationarity: Theory and Methods, Springer Verlag - 2021, v. 1, p. 19-40.
#'
#' @export

SIM_ZINAR = function(n,alpha,rho,th,innovation){

  if(n%%1!=0|n<1) stop("length of the time serie must be a positive integer")
  if(!(is.numeric(alpha))) stop("parameter alpha must be numeric")
  if(alpha<0 | alpha>=1) stop("the value of the parameter alpha must be in [0,1)")
  if(!(is.numeric(rho))) stop("parameter rho must be numeric")
  if(rho<0 | rho>=1) stop("the value of the parameter rho must be in [0,1)")
  if(innovation %in% c("NB","GI") & !(is.numeric(th))) stop("the parameter th must be a vector containing the positive values of the parameters mu and phi")
  if(innovation %in% c("NB","GI") & (th[1] < 0 | th[2] < 0)) stop("the parameter th must be a vector containing the positive values of the parameters mu and phi")
  if(innovation == "Po" & !(is.numeric(th))) stop("the parameter th must be a positive number containing the value of the parameter lambda")
  if(innovation == "Po" & (th[1] < 0)) stop("the parameter th must be a positive number containing the value of the parameter lambda")
  if(!(innovation %in% c("Po","NB","GI"))) stop("innovation must be set to Poisson ('Po'), Negative binomial ('NB') or Gaussian inverse ('GI')")
  if(innovation %in% c("NB","GI") & !(length(th) == 2)) stop("the parameter th must be a vector containing the values of the parameters mu and phi")
  if(innovation == "Po" & !(length(th) == 1)) stop("the parameter th must be a scalar containing the value of the parameter lambda")

  result = rzinar1(n=n,a=alpha,r=rho,th=th,family=innovation)

  return(result)
}


