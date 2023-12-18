PMLE4SCR = function(data, time, death, status_time, status_death,
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

  # Stage I of PMLE: estimate the marginal for death ----
  fitD = fitSPT(dat, time = "death", status = "status_death",
                formula = D.fmla, Gfun = "PH")
  betaD = as.vector(fitD$beta$est)
  dLambdaD = as.vector(fitD$dLambda$est)
  Psi.thetaD = do.call(cbind, fitD$Psi.theta)
  Psi.uD = predict.fitSPT(fitD, dat)$Psi.surv
  thetaD.est = c(betaD, dLambdaD)
  thetaD.cov = fitD$varcov$robust
  thetaD.se = sqrt(diag(thetaD.cov))

  # Initial value for the marginal of relapse for stage II of PMLE ----
  fitT = fitSPT(dat, time = "time", status = "status_time",
                formula = T.fmla, Gfun = "PH")
  betaT = as.vector(fitT$beta$est)
  dLambdaT  = as.vector(fitT$dLambda$est)
  thetaT.naive = c(betaT, dLambdaT)

  # Naive estimate
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

  # PMLE ----
  PMLE.ini = c(gamma.naive, thetaT.naive)
  objfun<- function(x){
    f = lln.fun(theta = x, thetaT = NULL, thetaD = thetaD.est,
                Xi, Ci, deltaT, deltaD,
                Zmat.T, Zmat.D, copula.index, Gfun,
                copula.link, Wmat, control)
    g = dlln.fun(theta = x, thetaT = NULL, thetaD = thetaD.est,
                 Xi, Ci, deltaT, deltaD,
                 Zmat.T, Zmat.D, copula.index, Gfun,
                 copula.link, Wmat, control)
    B = ddlln.fun(theta = x, thetaT = NULL, thetaD = thetaD.est,
                  Xi, Ci, deltaT, deltaD,
                  Zmat.T, Zmat.D, copula.index, Gfun,
                  copula.link, Wmat, control)
    list(value = f, gradient = g, hessian = B)
  }

  est.res <- trust(objfun, PMLE.ini, 5, 100, iterlim = 300,
                   minimize= FALSE, blather = T)
  if (!est.res$converged) stop("Error: PMLE did not converge.")
  theta.est = est.res$argument
  Imat = - ddlln.fun(theta = theta.est, thetaT = NULL, thetaD = thetaD.est,
                     Xi, Ci, deltaT, deltaD,
                     Zmat.T, Zmat.D, copula.index, Gfun,
                     copula.link, Wmat, control)
  dll.k = ddll.fun2(theta = theta.est, thetaD = thetaD.est,
                    Xi, Ci, deltaT, deltaD,
                    Zmat.T, Zmat.D, copula.index, Gfun,
                    copula.link, Wmat, control)
  dll.k[is.na(dll.k)] = 0
  Psi =  dll.fun(theta = theta.est, thetaT = NULL, thetaD = thetaD.est,
                 Xi, Ci, deltaT, deltaD,
                 Zmat.T, Zmat.D, copula.index, Gfun,
                 copula.link, Wmat, control)
  Psi[is.na(Psi)] = 0
  Psi.theta = Psi + Psi.uD %*% dll.k / N
  Vmat = crossprod(Psi.theta, Psi.theta) / N
  theta.cov = solve(Imat) %*% Vmat %*% solve(Imat) / N


  para.est = data.frame(type = "copula", para = colnames(Wmat), time = NA,
                        est = theta.est[1 : n.gamma],
                        se = sqrt(diag(theta.cov))[1 : n.gamma])

  para.est = rbind(
    para.est,
    data.frame(type = "betaT", para = colnames(Zmat.T), time = NA,
               est = theta.est[n.gamma + c(1 : n.bT)],
               se = sqrt(diag(theta.cov))[n.gamma + c(1 : n.bT)])
  )
  para.est = rbind(
    para.est,
    data.frame(type = "dLambdaT", para = "dLambdaT", time = tk,
               est = theta.est[n.gamma + n.bT + c(1 : n.tk)],
               se = sqrt(diag(theta.cov))[n.gamma + n.bT + c(1 : n.tk)])
  )
  para.est = rbind(
    para.est,
    data.frame(type = "betaD", para = colnames(Zmat.D), time = NA,
               est = thetaD.est[c(1 : n.bD)],
               se = sqrt(diag(thetaD.cov))[c(1 : n.bD)])
  )
  para.est = rbind(
    para.est,
    data.frame(type = "dLambdaD", para = "dLambdaD", time = dk,
               est = thetaD.est[n.bD + c(1 : n.dk)],
               se = sqrt(diag(thetaD.cov))[n.bD + c(1 : n.dk)])
  )

  theta.cov = list(
    copula = theta.cov[1 : n.gamma, 1 : n.gamma, drop = F],
    thetaT = theta.cov[n.gamma + c(1 : (n.bT + n.tk)),
                       n.gamma + c(1 : (n.bT + n.tk))],
    thetaD = thetaD.cov)

  out.PMLE = list(est = para.est, cov = theta.cov)

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

  est.res <- trust(objfun, MLE.ini, 5, 100, iterlim = 300,
    minimize= FALSE, blather = T)
  if (!est.res$converged) stop("Error: MLE did not converge.")
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

  para.est = data.frame(
    type = "copula", para = colnames(Wmat), time = NA,
    est = theta.est[1 : n.gamma],
    se = sqrt(diag(theta.cov))[1 : n.gamma],
    se.model = sqrt(diag(theta.cov.model))[1 : n.gamma])

  para.est = rbind(
    para.est,
    data.frame(
      type = "betaT", para = colnames(Zmat.T), time = NA,
      est = theta.est[n.gamma + c(1 : n.bT)],
      se = sqrt(diag(theta.cov))[n.gamma + c(1 : n.bT)],
      se.model = sqrt(diag(theta.cov.model))[n.gamma + c(1 : n.bT)])
  )
  para.est = rbind(
    para.est,
    data.frame(
      type = "dLambdaT", para = "dLambdaT", time = tk,
      est = theta.est[n.gamma + n.bT + c(1 : n.tk)],
      se = sqrt(diag(theta.cov))[n.gamma + n.bT + c(1 : n.tk)],
      se.model = sqrt(diag(theta.cov.model))[n.gamma + n.bT + c(1 : n.tk)])
  )
  para.est = rbind(
    para.est,
    data.frame(
      type = "betaD", para = colnames(Zmat.D), time = NA,
      est = theta.est[n.gamma + n.bT + n.tk + c(1 : n.bD)],
      se = sqrt(diag(theta.cov))[n.gamma + n.bT + n.tk + c(1 : n.bD)],
      se.model =
        sqrt(diag(theta.cov.model))[n.gamma + n.bT + n.tk + c(1 : n.bD)])
  )
  para.est = rbind(
    para.est,
    data.frame(
      type = "dLambdaD", para = "dLambdaD", time = dk,
      est = theta.est[n.gamma + n.bT + n.tk + n.bD + c(1 : n.dk)],
      se = sqrt(diag(theta.cov))[n.gamma + n.bT + n.tk + n.bD + c(1 : n.dk)],
      se.model =
        sqrt(diag(theta.cov.model))[n.gamma + n.bT + n.tk + n.bD + c(1 : n.dk)])
  )
  theta.cov = list(
    copula = theta.cov[1 : n.gamma, 1 : n.gamma, drop = F],
    thetaT = theta.cov[n.gamma + c(1 : (n.bT + n.tk)),
                       n.gamma + c(1 : (n.bT + n.tk))],
    thetaD = theta.cov[n.gamma + n.bT + n.tk + c(1 : (n.bD + n.dk)),
                       n.gamma + n.bT + n.tk + c(1 : (n.bD + n.dk))])

  out.MLE = list(est = para.est, cov = theta.cov)

  list(PMLE = out.PMLE, MLE = out.MLE, naive = para.naive,
       copula.family = copula.family, copula.link = copula.link,
       Tau2Par = list(tau.alpha = tau.alpha, Dtau.alpha = Dtau.alpha))
}
