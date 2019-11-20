# read dataset
setwd("C:/MY_DATAS/OneDrive - Data ScienceTech Institute/Documents/RProjects/20191031-ContinuousOpt")
datatmp <- read.table("toy_data.txt",sep="\t",header=T,blank.lines.skip=T)
print(str(datatmp))

# convert to numbers
data <- matrix(as.numeric(unlist(datatmp)),nrow=nrow(datatmp))

# dimension of the problem
n <- dim(data)[2]-1

# number of observations
p <- dim(data)[1]

# initial guess
a0 <- rep(0,n+1)

# computes the cost function
# a = array with a0, a1..an
residual <- function(a) {
 residual <- 0
 for (j in 1:p) {
  residual <- residual + .5*(data[j,1:n]%*%a[1:n]+a[n+1]-data[j,n+1])^2
 }
 return(residual)
}

# computes the gradient
gradient <- function(a) {
 gradient <- rep(0,n+1)
 for (j in 1:p) {
  for (i in 1:n) {
   gradient[i] = gradient[i] + data[j,i]*(data[j,1:n]%*%a[1:n]+a[n+1]-data[j,n+1])
  }
  gradient[n+1] =  gradient[n+1] + data[j,1:n]%*%a[1:n]+a[n+1]-data[j,n+1]
 }
 return(gradient)
}

# check gradient with finite increments
aref <- rep(0.2,n+1)
eps <- 1e-5
direction <- runif(n=n+1,min=-1,max=1)
gradient(aref)%*%direction
(residual(aref+eps*direction)-residual(aref))/eps

# fixed step gradient algorithm
rho <- 0.01
niter <- 500
solution <- matrix(,nrow=niter+1,ncol=n+1)
values <- rep(0,niter+1)
solution[1,] <- a0
values[1] <- residual(solution[1,])
for (k in 1:niter) {
 solution[k+1,] <- solution[k,]-rho*gradient(solution[k,])
 values[k+1] <- residual(solution[k+1,])
}

# plot solutions
dev.new()
plot(values)
dev.new()
matplot(solution)
dev.new()
index <- which.max(abs(solution[niter+1,]))
plot(data[,index],data[,n+1])
lines(data[,index],solution[niter+1,index]*data[,index]+solution[niter+1,n+1])
dev.new()
index2 <- which.min(abs(solution[niter+1,]))
plot(data[,index2],data[,n+1])
lines(data[,index2],solution[niter+1,index2]*data[,index2]+solution[niter+1,n+1])

# built-in linear regression
L <- lm(y~.,data=datatmp)
summary(L)

# comparison
print('built-in:')
c(as.numeric(L$coefficients)[2:11],as.numeric(L$coefficients)[1])
residual(c(as.numeric(L$coefficients)[2:11],as.numeric(L$coefficients)[1]))
print('gradient descent:')
solution[niter+1,]
values[niter+1]

# gradient with projection (we want positive coefficients)
rho <- 0.01
niter <- 500
solutionwp <- matrix(,nrow=niter+1,ncol=n+1)
valueswp <- rep(0,niter+1)
solutionwp[1,] <- a0
valueswp[1] <- residual(solution[1,])
for (k in 1:niter) {
 solutionwp[k+1,] <- solutionwp[k,]-rho*gradient(solutionwp[k,])
 solutionwp[k+1,] <- solutionwp[k+1,]*(solutionwp[k+1,]>0)
 valueswp[k+1] <- residual(solutionwp[k+1,])
}
solutionwp[niter+1,]
valueswp[niter+1]

# conjugate gradient
Atmp <- data
Atmp[,n+1] <- 1
A <- t(Atmp)%*%Atmp
b <- t(Atmp)%*%data[,n+1]
niter <- 20
solutioncg <- matrix(,nrow=niter+1,ncol=n+1)
directioncg <- matrix(,nrow=niter+1,ncol=n+1)
valuescg <- rep(0,niter+1)
solutioncg[1,] <- a0
directioncg[1,] <- gradient(solutioncg[1,])
valuescg[1] <- residual(solutioncg[1,])
for (k in 1:niter) {
 step <- (t(gradient(solutioncg[k,]))%*%directioncg[k,])/(t(A%*%directioncg[k,])%*%directioncg[k,])
 solutioncg[k+1,] <- solutioncg[k,]-step*directioncg[k,]
 beta <- (t(gradient(solutioncg[k+1,]))%*%(A%*%directioncg[k,]))/(t(A%*%directioncg[k,])%*%directioncg[k,])
 directioncg[k+1,] <- gradient(solutioncg[k+1,])-beta*directioncg[k,]
 valuescg[k+1] <- residual(solutioncg[k+1,])
}
solutioncg[niter+1,]
valuescg[niter+1]

# check convergence after n iterations
# compare with ".Machine$double.eps" (epsilon machine)
