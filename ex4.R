library(MASS)

Ker <- function(u){
  v <- numeric(length(u))
  if(abs(u)<=1){
    .75 * (1 - u^2)
  }else{
    0
  }
}

J <- 4
n <- 100

set.seed(2020)
# example2
x <- matrix(runif(n*J/2, -1, 1) ,nrow=n)
X <- cbind(x[,c(2,1)],x)

sigma <- matrix(nrow=J,ncol=J)
d <- c(0.2 ,0.3, 0.1, 0.4)
for (i in 1:J){
  for (j in 1:J){
    if (i == j){
      sigma[i, j] <- d[i]^2 
    }else{
      sigma[i, j] <- 0.6 * d[i] * d[j]
    }
  }
}

e <- matrix(nrow=n,ncol=J)
for (i in 1:n){
  e[i,] <- mvrnorm(1, rep(0, J), sigma)
}

Y <- 1 - 60 * X * exp(-20 * X^2) + e
truey <-  1 - 60 * X * exp(-20 * X^2)
y_flat <- as.vector(Y)
x_flat <- as.vector(X)

span <- 0.075
llc1 <- loess(y_flat ~ x_flat, degree = 1, span = span)
sort.index <- sort.int(x_flat, index.return = TRUE)$ix
plot(x_flat,y_flat)
lines(llc1$x[sort.index],llc1$fitted[sort.index],type='l')

res <- llc1$residuals
dim(res) <- c(n,J)

D<- diag(chol(sigma))^2
sortindex <- sort.int(1/D, index.return = TRUE)$ix

X <- X[,sortindex]
Y <- Y[,sortindex]
res <- res[,sortindex]
D <- c(0.16, 0.0576, 0.022, 0.0051)
X1 <- as.vector(t(X[,2:J]))

Y1 <- as.vector(t(Y[,2:J]))
Fa <- matrix(nrow=n*(J-1), ncol=J*(J-1)/2) 
count <- 1
for (i in 1:n){
  for (j in 2:J){
    Fa[count,] <- c(rep(0, (j-2)*(j-1)/2) ,res[i, 1:j-1], rep(0, (J-1)*J/2-(j-1)*j/2))
    count <- count + 1
    }
}

invG <- diag(rep(D[2:J],n))
h <- 0.03
S <- matrix(nrow=n*(J-1), ncol=n*(J-1)) 
count <- 1
for (i in 1:n){
  for (j in 2:J){
    A <- cbind(rep(1, n*(J-1)), X1-X[i,j])
    w <- (X1 - rep(X[i,j], n*(J-1))) / h
    dim(w) <- c(n*(J-1), 1)
    W <- apply(w, 1, Ker)/h
    W <- diag(W) %*% invG
    S[count,] <- matrix(c(1,0),nrow=1) %*% ginv(t(A) %*% W %*% A) %*% t(A) %*% W
    count <- count + 1
  }
}

In <- diag(n*(J-1))
fi <- ginv(t(Fa)%*%t(In-S)%*%invG%*%(In-S)%*%Fa)%*%t(Fa)%*%t(In-S)%*%invG%*%(In-S)%*%Y1
Y_star <- Y1 - Fa %*% fi
dim(Y_star) <- c(J-1,n)
Y2 <- cbind(Y[,1], t(Y_star))

invG0 <- diag(rep(D,n))
beta <- matrix(nrow=n*J, ncol=1)
count <- 1
for (i in 1:n){
  for (j in 1:J){
    A <- cbind(rep(1, n*J), as.vector(t(X))-X[i,j])
    w <- (as.vector(t(X)) - rep(X[i,j], n*J)) / h
    dim(w) <- c(n*J, 1)
    W <- apply(w, 1, Ker)/h
    W <- diag(W) %*% invG0
    B <- ginv(t(A) %*% W %*% A) %*% t(A) %*% W %*% as.vector(t(Y2))
    beta[count,] <- B[1,]
    count <- count + 1
  }
}

mean((as.vector(t(truey))-beta)^2) # 0.024