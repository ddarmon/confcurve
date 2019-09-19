library(confcurve)

# n1 = 15
# n2 = 15

n1 = 100
n2 = 100

theta1 = 0.5
theta2 = 0.5

or.true = theta1*(1-theta2)/(theta1*(1-theta2))

y1 = rbinom(1, n1, theta1)
y2 = rbinom(1, n1, theta1)

p1 = y1/n1; p2 = y2/n2

odds1 = p1*(1 - p2)/(p2*(1 - p1))
odds2 = p2*(1 - p1)/(p1*(1 - p2))

# (00*11)/(01*10)

odds3 = (n1-y1)*(y2)/(y1*(n2-y2))

x = c(y1, y2, n1-y1, n2-y2)

fisher.test(x = matrix(x, nrow = 2, byrow = TRUE))

show(odds1)
show(odds2)
show(odds3)

ns = c(n1, n2)
ys = c(y1, y2)

cc.out = confcurve.or(ys, ns, or.upper = 20, plot = TRUE)

cc.out$ci
cc.out$or.median.est
