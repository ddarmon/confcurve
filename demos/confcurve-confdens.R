library(confcurve)

n = 100

B = 5000

# statistic = function(x, d){
#   return(median(x[d]))
# }

# statistic = function(x, d){
#   return(mean(x[d]))
# }

statistic = function(x, d){
  return(sd(x[d]))
}

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
# Stand-alone code for generic confidence density
# reconstruction:
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bootcurve.out = bootcurve(x, statistic, B = B)

cdens.out = confdens(bootcurve.out, param = 1)

plot(cdens.out$theta, cdens.out$gn.perc, type = 'l', col = 'red', xlim = c(0, 1))
lines(cdens.out$theta, cdens.out$gn.bca, type = 'l', col = 'blue')
