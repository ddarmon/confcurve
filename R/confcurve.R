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
bootcurve.lm = function(formula, data, B = 2000){
  statistic = function(formula, data, indices){
    d = data[indices, ]

    return(coefficients(lm(formula, data = d)))
  }

  return(bootcurve(data, statistic, B = B, formula = formula))
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
plot.confcurve = function(object, param = 1, cs = seq(0.001, 0.999, by = 0.001), col = 'black', conf.level = NULL, xlim = NULL, xlab = NULL, add = FALSE){
  if(class(object) == 'lm'){
    confcurve.out = confcurve.lm(object = object, conf.level = cs, param = param)

    if (!is.null(conf.level)){
      if (length(conf.level) == 1){
        ci = confcurve.lm(object = object, conf.level = conf.level, param = param)
      }else{
        ci = list()

        for (cl.ind in 1:length(conf.level)){
          cl = conf.level[cl.ind]
          ci[[cl.ind]] = confcurve.lm(object = object, conf.level = cl, param = param)
        }
      }
    }

    if(is.null(xlab)) xlab = names(object$coefficients)[param]
  }else{
    confcurve.out = confcurve(bc = object, conf.level = cs, param = param)

    if (!is.null(conf.level)){
      if (length(conf.level) == 1){
        ci = confcurve(bc = object, conf.level = conf.level, param = param, warn.na = FALSE)
      }else{
        ci = list()

        for (cl.ind in 1:length(conf.level)){
          cl = conf.level[cl.ind]

          ci[[cl.ind]] = confcurve(bc = object, conf.level = cl, param = param, warn.na = FALSE)
        }
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

  if (!is.null(conf.level)){
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
  }
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
  if(class(object) == 'lm'){
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
  }else{
    point.est = object$t0

    lcb = rep(0, length(point.est))
    ucb = rep(0, length(point.est))

    for (param in 1:length(point.est)){
      ci.obj = confcurve(object, conf.level = conf.level, param = param)

      lcb[param] = ci.obj$cc.l
      ucb[param] = ci.obj$cc.u
    }

    effects.df = data.frame(coef.name = names(point.est), point.estimate = point.est, lcb = lcb, ucb = ucb)

    gf_pointrangeh(coef.name ~ point.estimate + lcb + ucb, data = effects.df, cex = cex) %>%
      gf_vline(xintercept = ~ 0, lty = 2) %>%
      gf_labs(x = "Coefficient Value", y = "Regressor Name")
  }
}

#' @export
confcurve.or = function(ys, ns, conf.level = 0.95, n = 10000, xlim = NULL, plot.ci = FALSE, show.normal.approximation = FALSE){
  ## See page 236 of *Confidence, Likelihood, Probability*
  ## for the form of the confidence distribution
  ## for the odds ratio.

  y1 = ys[1]; y2 = ys[2]
  n1 = ns[1]; n2 = ns[2]

  p1 = y1/n1
  p2 = y2/n2

  or = (p2*(1-p1))/(p1*(1-p2))

  # Let rho be the odds ratio, and
  #     psi be the log-odds ratio:
  #
  # psi = log(rho)

  log.or = log(or)

  kappa.s = 1/y1 + 1/(n1 - y1) + 1/y2 + 1/(n2 - y2)

  m1 = y1 + y2

  # Use normal approximation to determine a
  # reasonable upper-bound for the odds ratio
  # rho:

  log.or.upper = qnorm(0.999, mean = log.or, sd = sqrt(kappa.s))

  or.upper = exp(log.or.upper)

  rhos = seq(0.01, or.upper, length.out = n)

  H = rep(0, length(rhos))

  for (rho.ind in 1:length(rhos)){
    rho = rhos[rho.ind]
    p.i = dnoncenhypergeom(x = NA, n1 = n1, n2 = n2, m1 = m1, psi = rho)

    sum.inds = which(p.i[, 1] > y2)

    p.i = p.i[, 2]

    H[rho.ind] = sum(p.i[sum.inds]) + 0.5*p.i[sum.inds[1] - 1]
  }

  cc = abs(1 - 2*H)

  H.norm = pnorm(log(rhos), mean = log.or, sd = sqrt(kappa.s))

  cc.norm = abs(1 - 2*H.norm)

  ci.out = ci.or(ys, ns, conf.level)

  lr = ci.out$ci[1]
  rr = ci.out$ci[2]

  or.median.est = ci.out$or.median.est

    if(is.null(xlim)){
      xlim = c(0, or.upper)
    }

    plot(rhos, cc, type = 'l', xlim = xlim, xlab = expression('Odds Ratio' ~ rho), ylab = expression(cc(rho)))
    if(plot.ci) segments(x0 = c(lr), x1 = c(rr), y0 = c(conf.level), y1 = c(conf.level), lwd = 2, col = 'blue', xlab = expression("Odds Ratio" ~ rho), ylab = expression(cc(rho)))

    if (show.normal.approximation){
      lines(rhos, cc.norm, col = 'red')
    }
    abline(v = or.median.est, lty = 3, col = 'black')
    abline(v = 1)

    if (show.normal.approximation){
      legend('bottomright', legend = c('Exact CC', 'Normal Approx. CC'), col = c('black', 'red'))
    }

  return(list(ci = cbind(lr, rr), or.median.est = or.median.est))
}

#' @export
ci.or = function(ys, ns, conf.level = 0.95){
  ## See page 236 of *Confidence, Likelihood, Probability*
  ## for the form of the confidence distribution
  ## for the odds ratio.

  y1 = ys[1]; y2 = ys[2]
  n1 = ns[1]; n2 = ns[2]

  p1 = y1/n1
  p2 = y2/n2

  or = (p2*(1-p1))/(p1*(1-p2))

  # Let rho be the odds ratio, and
  #     psi be the log-odds ratio:
  #
  # psi = log(rho)

  log.or = log(or)

  kappa.s = 1/y1 + 1/(n1 - y1) + 1/y2 + 1/(n2 - y2)

  m1 = y1 + y2

  # The confidence distribution:

  Hfun = function(rho){
    p.i = dnoncenhypergeom(x = NA, n1 = n1, n2 = n2, m1 = m1, psi = rho)

    sum.inds = which(p.i[, 1] > y2)

    p.i = p.i[, 2]

    H = sum(p.i[sum.inds]) + 0.5*p.i[sum.inds[1] - 1]

    return(H)
  }

  # Values for finding lower confidence bound,
  # median, and upper confidence bound using
  # Brent's root finding method:

  alpha = 1-conf.level
  ad2 = alpha/2

  lower.probs = c(0.001, 0.25, 0.75)
  upper.probs = c(0.5, 0.75, 0.999)

  # The desired values of H(rho) from
  # the confidence distribution.

  zero.vals = c(ad2, 0.5, 1-ad2)

  # The values of rho that solve H(rho) = c

  interval.vals = rep(0, 3)

  for (i in 1:3){
    log.or.lb = qnorm(lower.probs[i], mean = log.or, sd = sqrt(kappa.s))
    log.or.ub = qnorm(upper.probs[i], mean = log.or, sd = sqrt(kappa.s))

    or.lb = exp(log.or.lb)
    or.ub = exp(log.or.ub)

    root.out = brentDekker(fun = function(x) Hfun(x) - zero.vals[i], a = or.lb, b = or.ub)$root

    interval.vals[i] = root.out
  }

  return(list(ci = c(interval.vals[1], interval.vals[3]), or.median.est = interval.vals[2]))
}

#' @export
confcurve.TukeyHSD = function(object, ordered = FALSE, conf.level = 0.95, which.term = 1, xlim = NULL, dc = 0.01, ncol = 3){
  tukey.out = TukeyHSD(object, ordered = ordered, conf.level = 0.95)

  which.diff = 1
  ndiffs = nrow(tukey.out[[which.term]])

  rnames = rownames(tukey.out[[which.term]])

  cs = seq(0, 1-dc, by = dc)

  conf.level = 0.95

  cc = array(NA, dim = c(ndiffs, length(cs), 3))

  cind.at.conf.level = conf.level / dc + 1

  if((cind.at.conf.level - trunc(cind.at.conf.level)) != 0){
    cat(sprintf('WARNING: Asking for a confidence interval at a finer resolution than the confidence curve. Please choose dc so that conf.level / dc is an integer value.'))
    cind.at.conf.level = NULL
  }

  for (c.ind in 1:length(cs)){
    c = cs[c.ind]

    tukey.out = TukeyHSD(object, ordered = ordered, conf.level = c)

    for (which.diff in 1:ndiffs){
      cc[which.diff, c.ind, ] = tukey.out[[which.term]][which.diff, 1:3]
    }
  }

  cc.range = range(cc)

  cex.use = 1

  if (is.null(xlim)){
    xlim = cc.range
  }

  par(mfrow = c(ceiling(ndiffs/ncol), ncol), mar=c(5,5,2,1), cex.lab = cex.use, cex.axis = cex.use)
  for (which.diff in 1:ndiffs){
    plot(cc[which.diff, , 2], cs, type = 'l', xlim = xlim, xlab = rnames[which.diff], ylab = 'Confidence Curve')
    lines(cc[which.diff, , 3], cs)
    abline(v = cc[which.diff, 1, 1])
    abline(v = 0, lty = 2)
    if(!is.null(cind.at.conf.level))
      segments(x0 = cc[which.diff, cind.at.conf.level, 2], x1 = cc[which.diff, cind.at.conf.level, 3], y0 = conf.level, lwd = 2)
  }
}

#' @export
confcurve.ScheffeTest = function(object, conf.level = 0.95, which.term = 1, xlim = NULL, dc = 0.01, ncol = 3){
  scheffe.out = ScheffeTest(object, conf.level = 0.95)

  which.diff = 1
  ndiffs = nrow(scheffe.out[[which.term]])

  rnames = rownames(scheffe.out[[which.term]])

  cs = seq(0, 1-dc, by = dc)

  conf.level = 0.95

  cc = array(NA, dim = c(ndiffs, length(cs), 3))

  cind.at.conf.level = conf.level / dc + 1

  if((cind.at.conf.level - trunc(cind.at.conf.level)) != 0){
    cat(sprintf('WARNING: Asking for a confidence interval at a finer resolution than the confidence curve. Please choose dc so that conf.level / dc is an integer value.'))
    cind.at.conf.level = NULL
  }

  for (c.ind in 1:length(cs)){
    c = cs[c.ind]

    scheffe.out = ScheffeTest(object, conf.level = c)

    for (which.diff in 1:ndiffs){
      cc[which.diff, c.ind, ] = scheffe.out[[which.term]][which.diff, 1:3]
    }
  }

  cc.range = range(cc)

  cex.use = 1

  if (is.null(xlim)){
    xlim = cc.range
  }

  par(mfrow = c(ceiling(ndiffs/ncol), ncol), mar=c(5,5,2,1), cex.lab = cex.use, cex.axis = cex.use)
  for (which.diff in 1:ndiffs){
    plot(cc[which.diff, , 2], cs, type = 'l', xlim = xlim, xlab = rnames[which.diff], ylab = 'Confidence Curve')
    lines(cc[which.diff, , 3], cs)
    abline(v = cc[which.diff, 1, 1])
    abline(v = 0, lty = 2)
    if(!is.null(cind.at.conf.level))
      segments(x0 = cc[which.diff, cind.at.conf.level, 2], x1 = cc[which.diff, cind.at.conf.level, 3], y0 = conf.level, lwd = 2)
  }
}
