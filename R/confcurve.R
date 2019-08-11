#' @export
bootcurve = function(data, statistic, B = 2000){
  if (length(B) > 1){
    t0 = statistic(data, 1:nrow(data))

    boot.out = list(t0 = t0, t = B)

    B = nrow(B)
  }else{
    boot.out = boot(data = data, statistic = statistic, R = B)
  }

  p = length(boot.out$t0) # Dimension of the parameter vector

  # Compute bias adjustment:
  z0 = rep(0, p)

  for (i in 1:p){
    z0[i] = qnorm(sum(boot.out$t[, i] <= boot.out$t0[i])/(B))
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
    u[i, ] <- statistic(data[-i, ], seq_len(n-1))
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
confpvalue = function(bc, theta, param){

  # A modification that prevents NAN P-values:

  n.leq = sum(bc$t[, param] <= theta)

  if (n.leq == 0){
    cat('Warning: True bootstrap P-value likely smaller than reported P-value, since the null parameter value is smaller than any of the bootstrap parameter estimates. Reporting Percentile Bootstrap-based P-value.\n')

    Gn = 1/(length(bc$t[, param]) + 1)

    return(2*min(Gn, 1 - Gn))
  }else if (n.leq == length(bc$t[, param])){
    cat('Warning: True bootstrap P-value likely smaller than reported P-value, since the null parameter value is larger than any of the bootstrap parameter estimates. Reporting Percentile Bootstrap-based P-value.\n')

    Gn = 1/(length(bc$t[, param]) + 1)

    return(2*min(Gn, 1 - Gn))
  }else{
    # The standard definition:
    Gn = bc$Gn[[param]]
    Phi.invs = qnorm(Gn(theta))

    # The BCa confidence distribution
    Hn = pnorm((Phi.invs - bc$z0[param])/(1 + bc$a[param]*(Phi.invs - bc$z0[param])) - bc$z0[param])

    return(Hn)
  }

}

#' @export
confcurve = function(bc, conf.level, param){
  alpha = (1 - conf.level)/2
  zalpha.l = qnorm(alpha)
  zalpha.r = -zalpha.l

  z0 = bc$z0[param]
  a = bc$a[param]

  zadj.l = z0 + zalpha.l
  zadj.r = z0 + zalpha.r

  Ps.l = pnorm(z0 + zadj.l/(1 - a*zalpha.l))
  Ps.r = pnorm(z0 + zadj.r/(1 - a*zalpha.r))

  cc.l = quantile(x = bc$t[, param], probs = Ps.l)
  cc.u = quantile(x = bc$t[, param], probs = Ps.r)

  return(list(cc.l = cc.l, cc.u = cc.u, conf.level = conf.level))
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
