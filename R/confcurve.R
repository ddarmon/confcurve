#' @export
bootcurve = function(data, statistic, B = 2000, formula = NULL){
  if (length(B) > 1){
    if (is.null(formula)){
      t0 = statistic(data, 1:nrow(data))
    }else{
      t0 = statistic(formula, data, 1:nrow(data))
    }

    boot.out = list(t0 = t0, t = B)

    B = nrow(B)
  }else{
    if (is.null(formula)){
      boot.out = boot(data = data, statistic = statistic, R = B)
    }else{
      boot.out = boot(data = data, statistic = statistic, R = B, formula = formula)
    }
  }

  p = length(boot.out$t0) # Dimension of the parameter vector

  # Compute bias adjustment:
  z0 = rep(0, p)

  for (i in 1:p){
    # Handle possible NAs:

    comp = boot.out$t[, i] <= boot.out$t0[i]
    B.wo.na = sum(!is.na(comp))

    z0[i] = qnorm(sum(comp, na.rm = TRUE)/(B.wo.na))
  }

  n = nrow(data)

  if (is.null(n)){
    n = length(data)
    # Reshape data into a matrix

    data = matrix(data, nrow = n)
  }

  # Compute acceleration adjustment:
  u = matrix(rep(0, n*p), nrow = n)

  n1 <- sqrt(n * (n - 1))

  for (i in seq_len(n)) {
    if (is.null(formula)){
      u[i, ] <- statistic(data[-i, ], seq_len(n-1))
    }else{
      u[i, ] <- statistic(formula, data[-i, ], seq_len(n-1))
    }

  }
  t. <- sweep(-u, 2, colMeans(u), "+") * (n - 1)
  a <- (1 / 6) * colSums(t.^3) / (colSums(t.^2))^1.5

  Gn = list()

  for (i in 1:p){
    Gn[[i]] = stats::ecdf(boot.out$t[, i])
  }

  return(list(t0 = boot.out$t0, t = boot.out$t, Gn = Gn, z0 = z0, a = a))
}

#' @export
confdist = function(bc, theta, param){
  Gn = bc$Gn[[param]]
  Phi.invs = qnorm(Gn(theta))

  # The BCa confidence distribution
  Hn = pnorm((Phi.invs - bc$z0[param])/(1 + bc$a[param]*(Phi.invs - bc$z0[param])) - bc$z0[param])

  return(Hn)
}

#' @export
confdens = function(bc, param){
  # density.out = density(bc$t[, param], bw = "SJ") # Seems to undersmooth
  density.out = density(bc$t[, param], bw = "bcv", n = 1024) # Seems to oversmooth, which in this case is good.

  gn.percentile = density.out$y
  thetas = density.out$x

  z0 = bc$z0[param]; a = bc$a[param]

  Gn = bc$Gn[[param]]

  w = function(theta, Gn, z0, a){
    ztheta = qnorm(Gn(theta)) - z0

    bca.fac = dnorm(ztheta/(1+a*ztheta) - z0)/((1 + a*ztheta)^2*dnorm(ztheta + z0))

    return(bca.fac)
  }

  gn.bca = gn.percentile*w(thetas, Gn, z0, a)

  # Reweight so sums to 1, at the given discretization of theta.
  gn.bca = gn.bca/(sum(gn.bca, na.rm = TRUE)*diff(thetas[1:2]))

  return(list(theta = thetas, gn.perc = gn.percentile, gn.bca = gn.bca))
}

#' @export
confpvalue = function(object, theta, param = 1){
  if(class(object) == 'lm'){
    lm.info = summary(object)

    df = lm.info$df[2]

    b = lm.info$coefficients[param, 1]
    se.b = lm.info$coefficients[param, 2]

    t.obs = (b - theta)/se.b

    P = 2*pt(-abs(t.obs), df)

    return(P)
  } else{
    # A modification that prevents NAN P-values:

    n.leq = sum(object$t[, param] <= theta)

    if (n.leq == 0){
      cat('Warning: True bootstrap P-value likely smaller than reported P-value, since the null parameter value is smaller than any of the bootstrap parameter estimates. Reporting Percentile Bootstrap-based P-value.\n')

      Gn = 1/(length(object$t[, param]) + 1)

      return(2*min(Gn, 1 - Gn))
    }else if (n.leq == length(object$t[, param])){
      cat('Warning: True bootstrap P-value likely smaller than reported P-value, since the null parameter value is larger than any of the bootstrap parameter estimates. Reporting Percentile Bootstrap-based P-value.\n')

      Gn = 1/(length(object$t[, param]) + 1)

      return(2*min(Gn, 1 - Gn))
    }else{
      # The standard definition:
      Gn = object$Gn[[param]]
      Phi.invs = qnorm(Gn(theta))

      # The BCa confidence distribution
      Hn = pnorm((Phi.invs - object$z0[param])/(1 + object$a[param]*(Phi.invs - object$z0[param])) - object$z0[param])

      return(2*min(Hn, 1 - Hn))
    }
  }
}

#' @export
confcurve = function(bc, conf.level, param, warn.na = TRUE){
  alpha = (1 - conf.level)/2
  zalpha.l = qnorm(alpha)
  zalpha.r = -zalpha.l

  z0 = bc$z0[param]
  a = bc$a[param]

  zadj.l = z0 + zalpha.l
  zadj.r = z0 + zalpha.r

  Ps.l = pnorm(z0 + zadj.l/(1 - a*zalpha.l))
  Ps.r = pnorm(z0 + zadj.r/(1 - a*zalpha.r))

  num.na = sum(is.na(bc$t[, param]))

  if(num.na > 0){
    if (warn.na) cat(sprintf("\n\nWarning: %g NAs present in the bootstrapped estimates. If %g is large relative to B = %g, this may indicate that bootstrapping will not work for this data set.\n\n", num.na, num.na, nrow(bc$t)))
  }

  cc.l = quantile(x = bc$t[, param], probs = Ps.l, na.rm = TRUE)
  cc.u = quantile(x = bc$t[, param], probs = Ps.r, na.rm = TRUE)

  return(list(cc.l = cc.l, cc.u = cc.u, conf.level = conf.level))
}

#' @export
plot.confcurve = function(object, cs = seq(0.001, 0.999, by = 0.001), conf.level = 0.95, param = 1, col = 'black', xlim = NULL, xlab = NULL, add = FALSE){
  if(class(object) == 'lm'){
    confcurve.out = confcurve.lm(object = object, conf.level = cs, param = param)

    if (length(conf.level) == 1){
      ci = confcurve.lm(object = object, conf.level = conf.level, param = param)
    }else{
      ci = list()

      for (cl.ind in 1:length(conf.level)){
        cl = conf.level[cl.ind]
        ci[[cl.ind]] = confcurve.lm(object = object, conf.level = cl, param = param)
      }
    }

    if(is.null(xlab)) xlab = names(object$coefficients)[param]
  }else{
    confcurve.out = confcurve(bc = object, conf.level = cs, param = param)

    if (length(conf.level) == 1){
      ci = confcurve(bc = object, conf.level = conf.level, param = param, warn.na = FALSE)
    }else{
      ci = list()

      for (cl.ind in 1:length(conf.level)){
        cl = conf.level[cl.ind]

        ci[[cl.ind]] = confcurve(bc = object, conf.level = cl, param = param, warn.na = FALSE)
      }
    }

    if(is.null(xlab)) xlab = names(object$t0)[param]
  }

  if(is.null(xlab)){
    xlab = 'Parameter'
  }

  if (is.null(xlim)){
    xlim = range(confcurve.out$cc.l, confcurve.out$cc.u)
  }

  if (add){
    lines(confcurve.out$cc.l, confcurve.out$conf.level, type = 'l', xlim = xlim, col = col)
  }else{
    plot(confcurve.out$cc.l, confcurve.out$conf.level, type = 'l', xlim = xlim, ylab = 'Confidence Level', xlab = xlab, col = col)
  }

  lines(confcurve.out$cc.u, confcurve.out$conf.level, type = 'l', col = col)

  if (length(conf.level) == 1){
    segments(x0 = ci$cc.l, x1 = ci$cc.u, y0 = conf.level, lwd = 2, col = col)

    if(class(object) == 'lm'){
      points(x = object$coefficients[param], y = conf.level, pch = 16, cex = 2, col = col)
    }else{
      points(x = object$t0[param], y = conf.level, pch = 16, cex = 2, col = col)
    }
  }else{
    for (cl.ind in 1:length(conf.level)){
      cl = conf.level[cl.ind]
      ci.cur = ci[[cl.ind]]

      segments(x0 = ci.cur$cc.l, x1 = ci.cur$cc.u, y0 = cl, lwd = 2, col = col)

      if(class(object) == 'lm'){
        points(x = object$coefficients[param], y = cl, pch = 16, cex = 2, col = col)
      }else{
        points(x = object$t0[param], y = cl, pch = 16, cex = 2, col = col)
      }
    }
  }

  # return(list(confcurve = confcurve.out, ci = ci))
}

#' @export
confcurve.lm = function(object, conf.level, param){
  lm.info = summary(object)

  df = lm.info$df[2]

  alpha = (1 - conf.level)/2

  pm = qt(1-alpha, df)*lm.info$coefficients[param, 2]

  cc.l = lm.info$coefficients[param, 1] - pm
  cc.u = lm.info$coefficients[param, 1] + pm

  return(list(cc.l = cc.l, cc.u = cc.u, conf.level = conf.level))
}

#' @export
plot.lm.coef = function(object, conf.level = 0.95, cex = 1){
lm.summary = summary(object)
alpha = 1 - conf.level

b = lm.summary$coefficients[, 1]

pm = qt(alpha/2, lm.summary$df[2], lower.tail = FALSE)*lm.summary$coefficients[, 2]

lcb = b - pm
ucb = b + pm

effects.df = data.frame(coef.name = rownames(lm.summary$coefficients), point.estimate = lm.summary$coefficients[, 1], lcb = lcb, ucb = ucb)

gf_pointrangeh(coef.name ~ point.estimate + lcb + ucb, data = effects.df, cex = cex) %>%
  gf_vline(xintercept = ~ 0, lty = 2) %>%
  gf_labs(x = "Coefficient Value", y = "Regressor Name")
}

#' @export
confcurve.or = function(ys, ns, conf.level = 0.95, or.upper = NA, plot = FALSE){
  y1 = ys[1]; y2 = ys[2]
  n1 = ns[1]; n2 = ns[2]

  m1 = y1 + y2

  if (is.na(or.upper)){
    or.upper = 30
  }


  psis = seq(0.01, or.upper, by = 0.01)

  H = rep(0, length(psis))

  for (psi.ind in 1:length(psis)){
    psi = psis[psi.ind]
    p.i = dnoncenhypergeom(x = NA, n1 = n1, n2 = n2, m1 = m1, psi = psi)

    sum.inds = which(p.i[, 1] > y2)

    p.i = p.i[, 2]

    H[psi.ind] = sum(p.i[sum.inds]) + 0.5*p.i[sum.inds[1] - 1]
  }

  H.fun = approxfun(psis, H)

  first.pos.ind = which(H - 0.5 > 0)[1]

  or.median.est = uniroot(function(x) H.fun(x) - 0.5, interval = c(psis[first.pos.ind - 1], psis[first.pos.ind]))$root

  cc = abs(1 - 2*H)

  # Find confidence interval using confidence curve:
  # neg.inds = which(cc - conf.level < 0)
  # first.neg.ind = neg.inds[1]
  # last.neg.ind = tail(neg.inds, 1)
  #
  #   cc.fun = approxfun(psis, cc)
  #
  #   lr = uniroot(function(x) cc.fun(x) - conf.level, interval = c(psis[first.neg.ind - 1], psis[first.neg.ind]))$root
  #   rr = uniroot(function(x) cc.fun(x) - conf.level, interval = c(psis[last.neg.ind], psis[last.neg.ind + 1]))$root

  # Find confidence interval using confidence distribution

  alpha = 1 - conf.level
  ad2 = alpha/2

  first.pos.ind = which(H - ad2 > 0)[1]

  lr = uniroot(function(x) H.fun(x) - ad2, interval = c(psis[first.pos.ind - 1], psis[first.pos.ind]))$root

  first.pos.ind = which(H - (1-ad2) > 0)[1]

  rr = uniroot(function(x) H.fun(x) - (1-ad2), interval = c(psis[first.pos.ind - 1], psis[first.pos.ind]))$root

  if (plot){
    plot(psis, cc, type = 'l'); segments(x0 = c(lr), x1 = c(rr), y0 = c(0.95), y1 = c(0.95), lwd = 2, col = 'blue')
    abline(v = or.median.est, lty = 3, col = 'black')
    abline(v = 1)

#     plot(psis, H, type = 'l')
#     abline(v = lr); abline(h = ad2)
#     abline(v = rr); abline(h = 1-ad2)
  }

  return(list(ci = cbind(lr, rr), or.median.est = or.median.est))
}
