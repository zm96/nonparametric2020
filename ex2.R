set.seed(2020)
x <- (seq(from=1, to=201, by=1) - 100) / 100
y <- 3 * x + cos(4 * pi * x ) + rnorm(201, 0, 0.2) 
truey <- 3 * x + cos(4 * pi * x)

library(KernSmooth)
# 窗宽
h <- dpill(x, y)
# NW核估计
nw <- ksmooth(x, y, kernel="normal", bandwidth = h)
# 局部线性估计
llr <- locpoly(x, y, bandwidth = h)
# 局部二次估计
lpr <- locpoly(x, y, degree = 2, bandwidth = h)

plot(x, y)
lines(x, truey, col = "black")
lines(nw, col = 'blue',lty = 2)
lines(llr, col = 'red', lty = 3)
lines(lpr, col = 'green', lty = 4)
legend("topleft", c("true", "nw", "llr", "lpr"), col=c("black", "blue","red","green"),lty=c(1, 2, 3, 4),ncol=2)



