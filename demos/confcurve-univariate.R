library(confcurve)

n = 100

B = 2000

statistic = function(x, d){
  return(median(x[d]))
}

# statistic = function(x, d){
#   return(mean(x[d]))
# }

# statistic = function(x, d){
#   return(sd(x[d]))
# }

# statistic = function(x, d){
#   return(mean(abs(x[d] - mean(x[d]))))
# }

# statistic = function(x, d){
#   return(sd(x[d])/mean(x[d]))
# }

# For mean
theta0 = 0
thetaa = -0.1

thetas = seq(-1, 1, by = 0.01)
cs = seq(0.001, 0.999, by = 0.001)

x = 2*(runif(n)-0.5) + thetaa

# For sd:
# theta0 = 1
# thetaa = 1.25
#
# x = thetaa*rnorm(n) + thetaa

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Stand-alone code for generic confidence
# distribution and confidence curve construction.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bootcurve.out = bootcurve(x, statistic, B = B)

confdist.out = confdist(bc = bootcurve.out, thetas, param = 1)
confcurve.out = confcurve(bc = bootcurve.out, conf.level = cs, param = 1)

plot.confcurve(bootcurve.out, xlim = c(-2, 2))
abline(v = theta0, lty = 2)
abline(v = thetaa, lty = 1)

show(confpvalue(bc = bootcurve.out, theta = theta0))
