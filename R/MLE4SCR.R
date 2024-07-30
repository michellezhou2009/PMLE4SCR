#' One-stage maximum likelihood estimation (MLE) for the copula-based analysis of semi-competing risks (SCR) data
#'
#' @description
#'
#' @param data a data frame containing variables names in the formulas: \code{T.fmla}, \code{D.fmla}, and \code{formula} in \code{copula.control}.
#' @param time a character string specifying the variable name in the \code{data} for the time to a non-terminal event.
#' @param death a character string specifying the variable name in the \code{data} for the time to death.
#' @param status_time a character string specifying the variable name in the \code{data} for the censoring indicator of the non-terminal event: \code{1} indicates the non-terminal time is observed, and \code{0} indicates censored.
#' @param status_death a character string specifying the variable name in the \code{data} for the censoring indicator of death: \code{1} indicates the death time is observed, and \code{0} indicates censored.
#' @param T.fmla an object of class \code{\link[stats]{formula}}: a symbolic description of the regression model to be fitted for the marginal distribution of the non-terminal event time.
#' @param D.fmla an object of class \code{\link[stats]{formula}}: a symbolic description of the regression model to be fitted for the marginal distribution of the death time.
#' @param Gfun a list of two components \code{T} and \code{D}, both character strings specifying the link function in the SPT model for the non-terminal event and death, respectively: "PH" (proportional hazards as the default value) or "PO" (proportional odds).
#' @param copula.family a character string: "Clayton", "Gumbel", or "Frank", specifying the copula family to be used for the dependence between the bivariate event times.
#' @param copula.control a list of two components: \code{link} (a character string: "identity", "log", or "log-1", specifying the link function for the regression model of the copula parameter; if \code{link = NULL} (default), the link function will be the default function for the specified copula family: "log" for Clayton, "identity" for Frank, and "log-1" for Gumbel), and \code{formula} (an object of class \code{\link[stats]{formula}}: a symbolic description of the regression model to be fitted for the copula parameter under the specified link function; if \code{formula = ~ 1} (default), the copula parameter is a constant)
#' @param initial a numerical value or a vector of numerical values for the initial values of the copula parameter or the regression coefficients for the copula parameter.
#'
#' @import trust rlang tidyverse VineCopula
#'
#' @export
#'
#' @return a list of the following components:
#' \itemize{
#'   \item{\code{gamma}: }{a data frame containing the estimate and robust standard error (SE) of the copula parameter or regression coefficients for the copula parameter.}
#'  \item{\code{gamma.cov}: }{a matrix containing the variance or variance-covariance matrix of the copula parameter or regression coefficients for the copula parameter.}
#'  \item{\code{betaT}: }{a data frame containing the estimate and robust SE of the regression coefficients for the marginal distribution of the non-terminal event time.}
#'  \item{\code{dLambdaT}: }{a data frame containing the estimate and robust SE of the jump size of the baseline function for the marginal distribution of the non-terminal event time.}
#'  \item{\code{thetaT.cov}: }{a matrix containing the variance-covariance matrix of the parameters for the marginal distribution of the non-terminal event time.}
#'  \item{\code{betaD}: }{a data frame containing the estimate and robust SE of the regression coefficients for the marginal distribution of the death time.}
#'  \item{\code{dLambdaD}: }{a data frame containing the estimate and robust SE of the jump size of the baseline function for the marginal distribution of the death time.}
#'  \item{\code{thetaD.cov}: }{a matrix containing the variance-covariance matrix of the parameters for the marginal distribution of the death time.}
#'  \item{\code{naive}}{a list containing naive estimates of \code{betaT}, \code{dLambdaT}, and \code{gamma} by estimating the marginal distribution of the non-terminal event time without the information on the between-event dependence.}
#'  \item{\code{call}: }{a list containing the specified values of input arguments \code{time}, \code{death}, \code{status_time}, \code{status_death}, \code{T.fmla}, \code{D.fmla}, \code{copula.family}, and the following two components:}
#'    \itemize{
#'      \item{\code{copula.link}: }{a list containing three R functions: "h.fun" (the link function used for the copula parameter), "dot.h.fun" (the first-order derivative of "h.fun"), and "ddot.h.fun" (the second-order derivative of "h.fun").}
#'      \item{\code{copula.fmla}: }{the specified value of \code{fomula} of the input argument \code{copula.control}.}
#'    }
#'  \item{\code{Par2Tau}: }{a list containing two R functions: "tau.alpha" (transformation from the copula parameter to Kendall's tau), and "Dtau.alpha" (the first-order derivative of "tau.alpha" function)}
#' }
#' @examples
#' \code{data(BMT, package = "SemiCompRisks")}
#' \code{data = BMT %>%
#'         mutate(g = factor(g, levels = c(2, 3, 1),
#'                     labels = c("AML-low", "AML-high", "ALL")))}
#' \code{myfit = PMLE4SCR(data, time = "T2", death = "T1",
#'                         status_time = "delta2", status_death = "delta1",
#'                         T.fmla = ~ g, D.fmla = ~ g,
#'                         copula.family = "Clayton",
#'                         copula.control = list(link = "identity", formula = ~ g),
#'                         initial = c(2, 0, 0))}
#' \code{myfit$gamma}
#' \code{myfit$betaT}
#'
MLE4SCR = function(data, time, death, status_time, status_death,
                    T.fmla = ~ 1, D.fmla = ~ 1,
                    Gfun = list(T = "PH", D = "PH"),
                    copula.family,
                    copula.control = list(link = NULL, formula = ~ 1),
                    initial){
  
  if (!rlang::is_formula(T.fmla))
    stop("Argument \"T.fmla\" needs to be a formula.")
  if (!rlang::is_formula(D.fmla))
    stop("Argument \"D.fmla\" needs to be a formula.")
  link.fun = copula.control$link; copula.fmla = copula.control$formula
  if (!rlang::is_formula(copula.fmla))
    stop("Argument \"formula\" of \"copula.control\" needs to be a formula.")
  
  allvars = list(time, death, status_time, status_death)
  allnm = c("time", "death", "status_time", "status_death")
  tmp.index = sapply(allvars, function(x) is.character(x))
  if (any(tmp.index == FALSE))
    stop(paste0("Values of arguments: \"",
                paste0(allnm[which(tmp.index == FALSE)], collapse = ", "),
                "\" needs to be characters."))
  var.nm = unique(c(all.vars(T.fmla), all.vars(D.fmla),
                    all.vars(copula.fmla)))
  allvars = c(unlist(allvars), var.nm);
  allnm = c(allnm, var.nm)
  exists.index = allvars %in% colnames(data)
  if (any(exists.index == FALSE))
    stop(paste0("Variables: ",
                paste0(allvars[exists.index == FALSE], collapse = ", ")),
         " do not exists in the data.")
  copula.initial = initial
  dat = data[, allvars]; colnames(dat) = allnm
  N = nrow(dat); col.nm = colnames(dat)
  Xi = dat[, "time"]; Ci = dat[, "death"];
  deltaT = dat[, "status_time"]; deltaD = dat[, "status_death"]
  Zmat.T = model.matrix(T.fmla, dat)[ , -1, drop = F]; n.bT = ncol(Zmat.T)
  Zmat.D = model.matrix(D.fmla, dat)[ , -1, drop = F]; n.bD = ncol(Zmat.D)
  Wmat = model.matrix(copula.fmla, dat); n.gamma = ncol(Wmat)
  tk = sort(unique(Xi[deltaT == 1])); n.tk = length(tk)
  dk = sort(unique(Ci[deltaD == 1])); n.dk = length(dk)
  
  switch(copula.family,
         "Clayton" = {
           copula.index = 3; copula.lwr = 0; copula.upr = 28
           tau.alpha = function(alpha) alpha / (alpha + 2)
           Dtau.alpha = function(alpha){}
           body(Dtau.alpha) = D(expression(alpha/(alpha + 2)), "alpha")
           link.default = "log"
         },
         "Frank" = {
           copula.index = 5;
           copula.lwr = 0; copula.upr = 50
           tau.alpha = function(alpha){
             fun0 = function(t) t / (exp(t) - 1)
             sapply(alpha, function(a){
               1 - 4 / a + (4 / a ^ 2) *
                 integrate(fun0, lower = 0, upper = a)$value
             })
           }
           Dtau.alpha = function(alpha){
             fun0 = function(t) t / (exp(t) - 1)
             sapply(alpha, function(a){
               4 / a ^ 2 + 4 / (a * (exp(a) - 1)) -
                 (8 / a ^ 3) * integrate(fun0, lower = 0, upper = a)$value
             })
           }
           link.default = "identity"
         },
         "Gumbel" = {
           copula.index = 4;
           copula.lwr = 1; copula.upr = 17
           tau.alpha = function(alpha) 1 - 1 / alpha
           Dtau.alpha = function(alpha) 1 / (alpha ^ 2)
           link.default = "log-1"
         }
  )
  if (is.null(link.fun)) link.fun = link.default
  switch(link.fun,
         "identity" = {
           copula.link = list(h.fun = function(x) {x},
                              dot.h.fun = function(x) {rep(1, length(x))},
                              ddot.h.fun = function(x) {rep(0, length(x))})
         },
         "log" = {
           copula.link = list(h.fun = function(x) {exp(x)},
                              dot.h.fun = function(x) {exp(x)},
                              ddot.h.fun = function(x) {exp(x)})
         },
         "log-1" = {
           copula.link = list(h.fun = function(x) {exp(x) + 1},
                              dot.h.fun = function(x) {exp(x)},
                              ddot.h.fun = function(x) {exp(x)})
         }
  )
  control = list(copula.lwr = copula.lwr, copula.upr = copula.upr)
  
  # Initial value for time to death ----
  fitD = fitSPT(dat, time = "death", status = "status_death",
                formula = D.fmla, Gfun = "PH")
  betaD = as.vector(fitD$beta$est)
  dLambdaD = as.vector(fitD$dLambda$est)
  Psi.thetaD = do.call(cbind, fitD$Psi.theta)
  Psi.uD = predict.fitSPT(fitD, dat)$Psi.surv
  thetaD.est = c(betaD, dLambdaD)
  thetaD.cov = fitD$varcov$robust
  thetaD.se = sqrt(diag(thetaD.cov))
  
  # Initial value for time to the non-terminal event ----
  fitT = fitSPT(dat, time = "time", status = "status_time",
                formula = T.fmla, Gfun = "PH")
  betaT.naive = as.vector(fitT$beta$est)
  dLambdaT.naive  = as.vector(fitT$dLambda$est)
  thetaT.naive = c(betaT.naive, dLambdaT.naive)
  
  # Initial value for copula
  objfun<- function(x){
    f = lln.fun(theta = x, thetaT = thetaT.naive, thetaD = thetaD.est,
                Xi, Ci, deltaT, deltaD,
                Zmat.T, Zmat.D, copula.index, Gfun,
                copula.link, Wmat, control)
    g = dlln.fun(theta = x, thetaT = thetaT.naive, thetaD = thetaD.est,
                 Xi, Ci, deltaT, deltaD,
                 Zmat.T, Zmat.D, copula.index, Gfun,
                 copula.link, Wmat, control)
    B = ddlln.fun(theta = x, thetaT = thetaT.naive, thetaD = thetaD.est,
                  Xi, Ci, deltaT, deltaD,
                  Zmat.T, Zmat.D, copula.index, Gfun,
                  copula.link, Wmat, control)
    list(value = f, gradient = g, hessian = B)
  }
  est.res <- trust(objfun, copula.initial, 5, 100, iterlim = 300,
                   minimize= FALSE, blather = T)
  if (!est.res$converged)
    stop("Error: Maximizing the log-likelihood function for a naive estimate
         did not converge.")
  gamma.naive = est.res$argument
  para.naive = data.frame(type = "copula", para = colnames(Wmat), time = NA,
                          est = gamma.naive)
  
  para.naive = rbind(
    para.naive,
    data.frame(type = "betaT", para = colnames(Zmat.T), time = NA,
               est = thetaT.naive[c(1 : n.bT)])
  )
  para.naive = rbind(
    para.naive,
    data.frame(type = "dLambdaT", para = "dLambdaT", time = tk,
               est = thetaT.naive[n.bT + c(1 : n.tk)])
  )
  
  # MLE ----
  MLE.ini = c(gamma.naive, thetaT.naive, thetaD.est)
  objfun<- function(x){
    f = lln.fun(theta = x, thetaT = NULL, thetaD = NULL, 
                Xi, Ci, deltaT, deltaD, 
                Zmat.T, Zmat.D, copula.index, Gfun, 
                copula.link, Wmat, control)
    g = dlln.fun(theta = x, thetaT = NULL, thetaD = NULL, 
                 Xi, Ci, deltaT, deltaD, 
                 Zmat.T, Zmat.D, copula.index, Gfun, 
                 copula.link, Wmat, control)
    B = ddlln.fun(theta = x, thetaT = NULL, thetaD = NULL, 
                  Xi, Ci, deltaT, deltaD, 
                  Zmat.T, Zmat.D, copula.index, Gfun, 
                  copula.link, Wmat, control) 
    list(value = f, gradient = g, hessian = B)
  }
  
  est.res <- trust(
    objfun, MLE.ini, 5, 100, iterlim = 300, 
    minimize= FALSE, blather = T)
  if (!est.res$converged) stop("Error: PMLE did not converge.")
  theta.est = est.res$argument
  Imat = - ddlln.fun(theta = theta.est, thetaT = NULL, thetaD = NULL, 
                     Xi, Ci, deltaT, deltaD, 
                     Zmat.T, Zmat.D, copula.index, Gfun, 
                     copula.link, Wmat, control) 
  Psi =  dll.fun(theta = theta.est, thetaT = NULL, thetaD = NULL, 
                 Xi, Ci, deltaT, deltaD, 
                 Zmat.T, Zmat.D, copula.index, Gfun, 
                 copula.link, Wmat, control) 
  Vmat = crossprod(Psi, Psi) / N
  theta.cov = solve(Imat) %*% Vmat %*% solve(Imat) / N
  theta.cov.model = solve(Imat) / N
  

  gamma.summary = data.frame(
    para = colnames(Wmat),
    est = theta.est[1 : n.gamma],
    se = sqrt(diag(theta.cov))[1 : n.gamma])
  
  betaT.summary = data.frame(
    para = colnames(Zmat.T),
    est = theta.est[n.gamma + c(1 : n.bT)],
    se = sqrt(diag(theta.cov))[n.gamma + c(1 : n.bT)])
  
  dLambdaT.summary = data.frame(
    time = tk,
    est = theta.est[n.gamma + n.bT + c(1 : n.tk)],
    se = sqrt(diag(theta.cov))[n.gamma + n.bT + c(1 : n.tk)])
  
  betaD.summary = data.frame(
    para = colnames(Zmat.D),
    est = thetaD.est[c(1 : n.bD)],
    se = sqrt(diag(thetaD.cov))[c(1 : n.bD)])
  
  dLambdaD.summary = data.frame(
    time = dk,
    est = thetaD.est[n.bD + c(1 : n.dk)],
    se = sqrt(diag(thetaD.cov))[n.bD + c(1 : n.dk)])
  
  
  gamma.cov = theta.cov[1 : n.gamma, 1 : n.gamma, drop = F]
  thetaT.cov = theta.cov[n.gamma + c(1 : (n.bT + n.tk)),
                         n.gamma + c(1 : (n.bT + n.tk))]
  
  list(gamma = gamma.summary, gamma.cov = gamma.cov,
       betaT = betaT.summary, dLambdaT = dLambdaT.summary,
       thetaT.cov = thetaT.cov,
       betaD = betaD.summary, dLambdaD = dLambdaD.summary,
       thetaD.cov = thetaD.cov,
       naive = list(betaT = betaT.naive, dLambdaT = dLambdaT.naive,
                    gamma = gamma.naive),
       call = list(time = time, death = death,
                   status_time = status_time, status_death = status_death,
                   T.fmla = T.fmla, D.fmla = D.fmla,
                   copula.family = copula.family, copula.link = copula.link,
                   copula.fmla = copula.fmla),
       Par2Tau = list(tau.alpha = tau.alpha, Dtau.alpha = Dtau.alpha))
}
