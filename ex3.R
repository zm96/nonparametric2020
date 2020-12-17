### q1
##Generate data
set.seed(2020)
U <- (seq(from=1, to=200, by=1) - 100) / 100
beta <- data.frame(beta0 = U^2, beta1 = sin(pi*U), beta2 = cos(pi*U))
X1 = rep(1,times=200)
X2 = rnorm(200)
X3 = rnorm(200)
W = rep(1,times=200)
error <- rnorm(200, 0, 0.2)
X <- data.frame(X1, X2, X3)
y <- apply(X*beta, 1, sum) + error

## https://github.com/Harrindy/VCMforPB
library(Rcpp)
library(RcppArmadillo)
sourceCpp("VCM.cpp")
source("Functions.R")
u_grid=seq(-1,1,length=200)
lmfit=VCM(u_grid,Y,U,X,W,c(0.01,1))
lmfit$h
par(mfrow=c(ceiling(nrow(lmfit$fit)),1))
par(mar=c(2,6,.1,.1))
for(k in 1:nrow(lmfit$fit))
{
  plot(u_grid,lmfit$fit[k,],type="l",xlab="u",ylab=expression(beta(u)))
}

### q2
### https://webcache.googleusercontent.com/search?q=cache:cLMwGCLwHTcJ:https://s0cran0r-project0org.icopy.site/web/packages/tvReg/vignettes/tvReg.html+&cd=1&hl=zh-CN&ct=clnk&gl=us
library(catdata)
library(tvReg)
data(aids)
attach(aids)

data <- data.frame(y = cd4, X1 = drugs, X2 = partners, X3 = packs, X4 = cesd, X5 = age, z = time)
model.tvLC <- tvLM(y ~ 0 + X1 + X2 + X3 + X4 + X5, z = time, data = data, est = "lc", tkernel = "Gaussian")
model.tvLM <- tvLM(y ~ 0 + X1 + X2 + X3 + X4 + X5, z = time, data = data, est = "ll", tkernel = "Gaussian")

sort.index <- sort.int(time, index.return = TRUE)$ix
## beta1
plot(time[sort.index], model.tvLM$coefficients[sort.index, 1], type = "l", main = "", ylab = expression(beta[1]), 
     xlab = expression(time[t]), col = 1)
lines(time[sort.index], model.tvLC$coefficients[sort.index, 1], col = 2)
legend("topright", c("tvlm", "tvlc"), col = c(1, 2), bty = "n", lty = 1, cex = 0.5)
## beta2
plot(time[sort.index], model.tvLM$coefficients[sort.index, 2], type = "l", main = "", ylab = expression(beta[2]), 
     xlab = expression(time[t]), col = 1)
lines(time[sort.index], model.tvLC$coefficients[sort.index, 2], col = 2)
legend("topright", c("tvlm", "tvlc"), col = c(1, 2), bty = "n", lty = 1, cex = 0.5)
## beta3
plot(time[sort.index], model.tvLM$coefficients[sort.index, 3], type = "l", main = "", ylab = expression(beta[3]), 
     xlab = expression(time[t]), col = 1)
lines(time[sort.index], model.tvLC$coefficients[sort.index, 3], col = 2)
legend("topright", c("tvlm", "tvlc"), col = c(1, 2), bty = "n", lty = 1, cex = 0.5)
## beta4
plot(time[sort.index], model.tvLM$coefficients[sort.index, 4], type = "l", main = "", ylab = expression(beta[4]), 
     xlab = expression(time[t]), col = 1)
lines(time[sort.index], model.tvLC$coefficients[sort.index, 4], col = 2)
legend("topright", c("tvlm", "tvlc"), col = c(1, 2), bty = "n", lty = 1, cex = 0.5)
## beta5
plot(time[sort.index], model.tvLM$coefficients[sort.index, 5], type = "l", main = "", ylab = expression(beta[5]), 
     xlab = expression(time[t]), col = 1)
lines(time[sort.index], model.tvLC$coefficients[sort.index, 5], col = 2)
legend("topright", c("tvlm", "tvlc"), col = c(1, 2), bty = "n", lty = 1, cex = 0.5)
