library(confcurve)

par(mar=c(5,5,2,1), cex.lab = 2, cex.axis = 2)

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

eps = rnorm(n)
# eps = rt(n, 2)
# eps = rexp(n) - 1

beta.true = c(0, 1, 2, rep(0, p-2))

y = beta.true[1] + beta.true[2]*x[, 2] + beta.true[3]*x[, 3] + eps

dafr = data.frame(x = x, y = y)

statistic = function(formula, data, indices){
  d = data[indices, ]
  lm.out = lm(formula, data = d)

  return(coefficients(lm.out))
}

# boot.out = boot(data = x, statistic = statistic, R = B)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Stand-alone code for generic confidence
# distribution and confidence curve construction.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bootcurve.out = bootcurve(dafr, statistic, B = B, formula = y ~ .)

par(mfrow = c(2, 3))
for(param in 1:(p+1)){
  confdist.out = confdist(bc = bootcurve.out, thetas, param = param)
  confcurve.out = confcurve(bc = bootcurve.out, conf.level = cs, param = param)

  plot.confcurve(bootcurve.out, param = param)

  show(confpvalue(bc = bootcurve.out, theta = theta0, param = param))
}

lm.out = lm(y ~ ., data = dafr)

par(mfrow = c(2, 3))
for(param in 1:(p+1)){
  plot.confcurve(lm.out, param = param)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Cars data set:
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bootcurve.out = bootcurve(cars, statistic, B = B, formula = dist ~ speed)

par(mfrow = c(2, 1))
for (param in 1:2){
  confcurve.out = confcurve(bc = bootcurve.out, conf.level = cs, param = param)

  plot.confcurve(bootcurve.out, param = param)

  lm.out = lm(dist ~ speed, cars)

  confcurve.lm.out = confcurve.lm(object = lm.out, conf.level = cs, param = param)

  plot.confcurve(lm.out, param = param, add = TRUE)
}

par(mfrow = c(2, 2))
plot(lm.out)
par(mfrow = c(1, 1))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Chicken Weights Data Set
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bootcurve.out = bootcurve(chickwts, statistic, B = B, formula = weight ~ feed)

plot(weight ~ feed, chickwts)

par(mfrow = c(2, 3))

cat('\n\n')

lm.out = lm(weight ~ feed, chickwts)

show(cbind(summary(lm.out)$coefficients[, 1]))

for(param in 1:6){
  confcurve.out = confcurve(bc = bootcurve.out, conf.level = cs, param = param)

  plot.confcurve(bootcurve.out, param = param)
  plot.confcurve(lm.out, param = param, add = TRUE)

  confcurve.lm.out = confcurve.lm(object = lm.out, conf.level = cs, param = param)

  show(confpvalue(bc = bootcurve.out, theta = 0, param = param))
}

par(mfrow = c(2, 2))
plot(lm.out)
par(mfrow = c(1, 1))

show(cbind(summary(lm.out)$coefficients[, 4]))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Fort Collins Snow Data Set
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(alr4)

plot(Late ~ Early, ftcollinssnow)

bootcurve.out = bootcurve(ftcollinssnow, statistic, B = B, formula = Late ~ Early)

par(mfrow = c(2, 1))

cat('\n\n')

lm.out = lm(Late ~ Early, ftcollinssnow)

show(cbind(summary(lm.out)$coefficients[, 1]))

for(param in 1:2){
  confcurve.out = confcurve(bc = bootcurve.out, conf.level = cs, param = param)

  plot.confcurve(bootcurve.out, param = param)
  plot.confcurve(lm.out, param = param, add = TRUE)

  confcurve.lm.out = confcurve.lm(object = lm.out, conf.level = cs, param = param)

  confdens.out = confdens(bootcurve.out, param)

  plot(confdens.out$theta, confdens.out$gn.perc, type = 'l', col = 'red', lwd = 2)
  lines(confdens.out$theta, confdens.out$gn.bca, col = 'blue', lwd = 2)

  show(confpvalue(bc = bootcurve.out, theta = 0, param = param))
}

show(summary(lm.out))

par(mfrow = c(2, 2))
plot(lm.out)
par(mfrow = c(1, 1))
