library(confcurve)

n = 100

B = 2000

xlim.use = c(-0.5, 0.5)
# xlim.use = c(-2, 2)

theta0 = 0
thetaa = 0.5

thetas = seq(-1, 1, by = 0.01)
cs = seq(0.001, 0.999, by = 0.001)

p = 5

x = matrix(rnorm(n*p), nrow = n, ncol = p)
x = cbind(rep(1, n), x)

eps = rnorm(n)
# eps = rt(n, 2)
# eps = rexp(n) - 1

beta.true = c(0, 1, 2, rep(0, p-2))

y = beta.true[1] + beta.true[2]*x[, 2] + beta.true[3]*x[, 3] + eps

X = cbind(x, y)

statistic = function(x, d){
  p = ncol(x) - 1
  lm.out = lm.fit(x[d, 1:p], x[d, p+1])

  return(coefficients(lm.out))
}

# boot.out = boot(data = x, statistic = statistic, R = B)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Stand-alone code for generic confidence
# distribution and confidence curve construction.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bootcurve.out = bootcurve(X, statistic, B = B)

par(mfrow = c(2, 3))
for(param in 1:(p+1)){
  confdist.out = confdist(bc = bootcurve.out, thetas, param = param)
  confcurve.out = confcurve(bc = bootcurve.out, conf.level = cs, param = param)

  plot(confcurve.out$cc.l, confcurve.out$conf.level, type = 'l', xlim = xlim.use)
  lines(confcurve.out$cc.u, confcurve.out$conf.level, type = 'l')

  # lines(thetas, abs(2*confdist.out-1), type = 'l')

  conf.level = 0.95
  ci = confcurve(bc = bootcurve.out, conf.level = conf.level, param = param)
  segments(x0 = ci$cc.l, x1 = ci$cc.u, y0 = conf.level, lwd = 2)
  points(x = bootcurve.out$t0[param], y = conf.level, pch = 16, cex = 2)
  abline(v = beta.true[param], lty = 2)

  show(confpvalue(bc = bootcurve.out, theta = theta0, param = param))
}

lm.out = lm(X[, p+2] ~ X[, 1:(p+1)])

par(mfrow = c(2, 3))
for(param in 1:(p+1)){
confcurve.lm.out = confcurve.lm(object = lm.out, conf.level = cs, param = param)

plot(confcurve.lm.out$cc.l, confcurve.lm.out$conf.level, type = 'l', xlim = xlim.use)
lines(confcurve.lm.out$cc.u, confcurve.lm.out$conf.level, type = 'l')
conf.level = 0.95
ci = confcurve.lm(object = lm.out, conf.level = conf.level, param = param)
segments(x0 = ci$cc.l, x1 = ci$cc.u, y0 = conf.level, lwd = 2)
points(x = bootcurve.out$t0[param], y = conf.level, pch = 16, cex = 2)
abline(v = beta.true[param], lty = 2)
}

par(mfrow = c(2, 3))
for(param in 1:(p+1)){
  confdist.out = confdist(bc = bootcurve.out, thetas, param = param)
  confcurve.out = confcurve(bc = bootcurve.out, conf.level = cs, param = param)

  plot(confcurve.out$cc.l, confcurve.out$conf.level, type = 'l', col = 'blue', xlim = range(confcurve.out$cc.l, confcurve.out$cc.u))
  lines(confcurve.out$cc.u, confcurve.out$conf.level, type = 'l', col = 'blue')

  confcurve.lm.out = confcurve.lm(object = lm.out, conf.level = cs, param = param)

  lines(confcurve.lm.out$cc.l, confcurve.lm.out$conf.level, type = 'l', xlim = xlim.use, col = 'red')
  lines(confcurve.lm.out$cc.u, confcurve.lm.out$conf.level, type = 'l', col = 'red')

  abline(v = beta.true[param], lty = 2)
  abline(h = 0.95)
}

cars.statistic = function(df, d){
  return(lm(dist ~ speed, df[d, ])$coefficients)
}

bootcurve.out = bootcurve(cars, cars.statistic, B = B)

param = 1
confcurve.out = confcurve(bc = bootcurve.out, conf.level = cs, param = param)

par(mfrow = c(1, 1), mar=c(5,5,2,1), cex.lab = 2, cex.axis = 2)
plot(confcurve.out$cc.l, confcurve.out$conf.level, type = 'l', xlim = range(confcurve.out$cc.l, confcurve.out$cc.u), col = 'blue')
lines(confcurve.out$cc.u, confcurve.out$conf.level, type = 'l', col = 'blue')

lm.out = lm(dist ~ speed, cars)

confcurve.lm.out = confcurve.lm(object = lm.out, conf.level = cs, param = param)

lines(confcurve.lm.out$cc.l, confcurve.lm.out$conf.level, type = 'l', col = 'red')
lines(confcurve.lm.out$cc.u, confcurve.lm.out$conf.level, type = 'l', col = 'red')

chickwts.statistic = function(df, d){
  return(lm(weight ~ feed, df[d, ])$coefficients)
}

bootcurve.out = bootcurve(chickwts, chickwts.statistic, B = B)

plot(weight ~ feed, chickwts)

par(mfrow = c(2, 3))

cat('\n\n')

lm.out = lm(weight ~ feed, chickwts)

show(cbind(summary(lm.out)$coefficients[, 1]))

for(param in 1:6){
  confcurve.out = confcurve(bc = bootcurve.out, conf.level = cs, param = param)

  plot(confcurve.out$cc.l, confcurve.out$conf.level, type = 'l', xlim = range(confcurve.out$cc.l, confcurve.out$cc.u), col = 'blue')
  lines(confcurve.out$cc.u, confcurve.out$conf.level, type = 'l', col = 'blue')

  confcurve.lm.out = confcurve.lm(object = lm.out, conf.level = cs, param = param)

  lines(confcurve.lm.out$cc.l, confcurve.lm.out$conf.level, type = 'l', col = 'red')
  lines(confcurve.lm.out$cc.u, confcurve.lm.out$conf.level, type = 'l', col = 'red')
  abline(v = 0, lty = 2)
  abline(h = 0.95)

  show(confpvalue(bc = bootcurve.out, theta = 0, param = param))
}

show(cbind(summary(lm.out)$coefficients[, 4]))

library(alr4)

plot(Late ~ Early, ftcollinssnow)

ftcollinssnow.statistic = function(df, d){
  return(lm(Late ~ Early, df[d, ])$coefficients)
}

bootcurve.out = bootcurve(ftcollinssnow, ftcollinssnow.statistic, B = B)

plot(Late ~ Early, ftcollinssnow)

par(mfrow = c(2, 1))

cat('\n\n')

lm.out = lm(Late ~ Early, ftcollinssnow)

show(cbind(summary(lm.out)$coefficients[, 1]))

for(param in 1:2){
  confcurve.out = confcurve(bc = bootcurve.out, conf.level = cs, param = param)

  plot(confcurve.out$cc.l, confcurve.out$conf.level, type = 'l', xlim = range(confcurve.out$cc.l, confcurve.out$cc.u), col = 'blue')
  lines(confcurve.out$cc.u, confcurve.out$conf.level, type = 'l', col = 'blue')

  confcurve.lm.out = confcurve.lm(object = lm.out, conf.level = cs, param = param)

  lines(confcurve.lm.out$cc.l, confcurve.lm.out$conf.level, type = 'l', col = 'red')
  lines(confcurve.lm.out$cc.u, confcurve.lm.out$conf.level, type = 'l', col = 'red')
  abline(v = 0, lty = 2)
  abline(h = 0.95)

  confdens.out = confdens(bootcurve.out, param)

  plot(confdens.out$theta, confdens.out$gn.perc, type = 'l', col = 'red', lwd = 2)
  lines(confdens.out$theta, confdens.out$gn.bca, col = 'blue', lwd = 2)

  show(confpvalue(bc = bootcurve.out, theta = 0, param = param))
}

show(summary(lm.out))
