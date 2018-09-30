library(mcmc)
library(MASS)
library(KernSmooth)
library(fields)
library(pscl)
library(mvtnorm)
library(SparseM)


#################################
##############function#############
################################

## Read and manipulate basic data

n              <- 1000                                                                    ## number of observations
n.knot         <- 40                                                                      ## number of knots
s              <- cbind(runif(n),runif(n))                                                ## observations
knot           <- cbind(runif(n.knot),runif(n.knot))                                      ## knots
dist.knot      <- rdist(knot)                                                             ## distance matrix with knots
dist.loc       <- rdist(s)                                                                ## distance matrix with locations
dist.loc.knot  <- rdist(s,knot)                                                           ## distance matrix with knots and locations

nu             <- 3                                                                       
tap1           <- 0.30                                                                    ## tapering range in the first dimension
tap2           <- 0.30                                                                    ## tapering range in the second dimension
tap12          <- 0.60                                                                    ## tapering range in the cross-covariance

T1             <- matrix(((1+(nu+1)*as.matrix(dist.loc/tap1))*
                  ifelse(((1-as.matrix(dist.loc/tap1))>0),1-as.matrix(dist.loc/tap1),0)^
                  (nu+1)),n,n)                                                            ## tapering matrix in dimension 1

T2             <- matrix(((1+(nu+1)*as.matrix(dist.loc/tap2))*
                  ifelse(((1-as.matrix(dist.loc/tap2))>0),1-as.matrix(dist.loc/tap2),0)^
                  (nu+1)),n,n)                                                            ## tapering matrix in dimension 2

T12            <- 0.8*matrix(((1+(nu+1)*as.matrix(dist.loc/tap12))*
                  ifelse(((1-as.matrix(dist.loc/tap12))>0),1-as.matrix(dist.loc/tap12),0)^
                  (nu+1)),n,n)                                                            ## tapering matrix in cross-covariance

tap.mat                            <- matrix(0,2*n,2*n)
tap.mat[seq(1,2*n,2),seq(1,2*n,2)] <- T1
tap.mat[seq(2,2*n,2),seq(2,2*n,2)] <- T2
tap.mat[seq(1,2*n,2),seq(2,2*n,2)] <- T12
tap.mat[seq(2,2*n,2),seq(1,2*n,2)] <- T12


X         <- kronecker(rep(1,n),diag(2))                                                  ## X matrix
sm1       <- 0.5                                                                          ## smoothness parameter in dimension 1
sm2       <- 1.5                                                                          ## smoothness parameter in dimension 2
sm12      <- (sm1+sm2)/2                                                                  ## smoothness parameter in cross-covariance
phi.tr    <- 1.5                                                                          ## range parameter
sigma.tr1 <- 2                                                                            ## spatial variance in dimension 1
sigma.tr2 <- 1.5                                                                          ## spatial variance in dimension 2
rho.tr    <- 0.3                                                                          
tau.tr1   <- 0.2                                                                          ## error variance in dimension 1 
tau.tr2   <- 0.1                                                                          ## error variance in dimension 2
beta.tr   <- c(4,-1)                                                                      ## fixed effects

D1        <- sigma.tr1*Matern(dist.loc, range=1/phi.tr, smoothness=sm1)
D2        <- sigma.tr2*Matern(dist.loc, range=1/phi.tr, smoothness=sm2)
D12       <- sqrt(sigma.tr1*sigma.tr2)*rho.tr*
                  Matern(dist.loc, range=1/phi.tr, smoothness=sm12)

var.mat                            <- matrix(0,2*n,2*n)
var.mat[seq(1,2*n,2),seq(1,2*n,2)] <- D1
var.mat[seq(2,2*n,2),seq(2,2*n,2)] <- D2
var.mat[seq(1,2*n,2),seq(2,2*n,2)] <- D12
var.mat[seq(2,2*n,2),seq(1,2*n,2)] <- D12

W <- rmvnorm(1,rep(0,2*n),var.mat)                                                        ## multivariate spatial random effect
Y <- rnorm(2*n,X%*%beta.tr+c(W),kronecker(rep(1,n),sqrt(c(tau.tr1,tau.tr2))))             ## multivariate response


## Permutation matrix

P                           <- matrix(0,2*n,2*n)
sequence1                   <- seq(1,2*n-1,2)
sequence2                   <- seq(2,2*n,2)
seq.in                      <- numeric()
seq.in[sequence1]           <- c(1:n)
seq.in[sequence2]           <- c(1:n)+n
P[cbind(c(1:(2*n)),seq.in)] <-1

## Initialize

a       <- 2
b       <- 1
al      <- 2
bet     <- 1
u.phi1  <- 8
l.phi1  <- 1
f       <- 5000
beta    <- matrix(0,f,2)
phi1    <- vector()
sigma1  <- numeric()
sigma2  <- numeric()
rho12   <- numeric()
count   <- vector()
count1  <- vector()
count2  <- vector()
psi     <- matrix(NA,f,2)



## initial values

beta[1,]   <- c(mean(Y[seq(1,2*n,2)]),mean(Y[seq(2,2*n,2)]))
phi1[1]    <- 2
sigma1[1]  <- 2
sigma2[1]  <- 2
rho12[1]   <- 0.3
psi[1,]    <- c(0.3, 0.3)
count[1]   <- 0
count1[1]  <- 0
count2[1]  <- 0
start      <- 4901
end        <- f

## MCMC iteration

for(i in 2:f){


    ## beta update
    DD1                                          <- sigma1[i-1]*Matern(dist.loc, range=1/phi1[i-1], smoothness=sm1)
    DD2                                          <- sigma2[i-1]*Matern(dist.loc, range=1/phi1[i-1], smoothness=sm2)
    DD12                                         <- sqrt(sigma1[i-1]*sigma2[i-1])*rho12[i-1]*
                                                    Matern(dist.loc, range=1/phi1[i-1], smoothness=sm12)

    R                                            <- matrix(0,2*n,2*n)
    R[seq(1,2*n,2),seq(1,2*n,2)]                 <- DD1
    R[seq(2,2*n,2),seq(2,2*n,2)]                 <- DD2
    R[seq(1,2*n,2),seq(2,2*n,2)]                 <- DD12
    R[seq(2,2*n,2),seq(1,2*n,2)]                 <- DD12

     R_star1                                     <- sigma1[i-1]*Matern(dist.knot, range=1/phi1[i-1], smoothness=sm1)
     R_star2                                     <- sigma2[i-1]*Matern(dist.knot, range=1/phi1[i-1], smoothness=sm2)
     R_star12                                    <- rho12[i-1]*sqrt(sigma1[i-1]*sigma2[i-1])*
                                                    Matern(dist.knot, range=1/phi1[i-1], smoothness=sm12)
     R_star                                      <- matrix(0,2*n.knot,2*n.knot)
     R_star[seq(1,2*n.knot,2),seq(1,2*n.knot,2)] <- R_star1
     R_star[seq(2,2*n.knot,2),seq(2,2*n.knot,2)] <- R_star2
     R_star[seq(1,2*n.knot,2),seq(2,2*n.knot,2)] <- R_star12
     R_star[seq(2,2*n.knot,2),seq(1,2*n.knot,2)] <- R_star12
     R_star.inv                                  <- chol2inv(chol(R_star))


     R_cv1                                       <- sigma1[i-1]*Matern(dist.loc.knot, range=1/phi1[i-1], smoothness=sm1)
     R_cv2                                       <- sigma2[i-1]*Matern(dist.loc.knot, range=1/phi1[i-1], smoothness=sm2)
     R_cv12                                      <- rho12[i-1]*sqrt(sigma1[i-1]*sigma2[i-1])*
                                                    Matern(dist.loc.knot, range=1/phi1[i-1], smoothness=sm12)

     R_cv                                        <- matrix(0,2*n,2*n.knot)
     R_cv[seq(1,2*n,2),seq(1,2*n.knot,2)]        <- R_cv1
     R_cv[seq(2,2*n,2),seq(2,2*n.knot,2)]        <- R_cv2
     R_cv[seq(1,2*n,2),seq(2,2*n.knot,2)]        <- R_cv12
     R_cv[seq(2,2*n,2),seq(1,2*n.knot,2)]        <- R_cv12
     diagon.tpp                                  <- (R-R_cv%*%R_star.inv%*%t(R_cv))*tap.mat
     dg.tp                                       <- as.matrix.csr(diagon.tpp+diag(rep(psi[i-1,],n)))
     dgg                                         <- chol2inv(chol(as.matrix(dg.tp)))
     dgg.ps                                      <- dgg%*%R_cv
     Sigv.inv                                    <- dgg-dgg.ps%*%chol2inv(chol(R_star+t(R_cv)%*%dgg.ps))%*%t(dgg.ps)
     Sigma_bet                                   <- chol2inv(chol(t(X)%*%Sigv.inv%*%X))
     mu_bet                                      <- Sigma_bet%*%(t(X)%*%Sigv.inv%*%Y)
     beta[i,]                                    <- mvrnorm(1,mu_bet,Sigma_bet)                  ## fixed effects updated

## update variance

     target<-function(theta){
          trans.psi1 <- exp(theta[1])
          trans.psi2 <- exp(theta[2])
          dg.tp      <- as.matrix.csr(diagon.tpp+diag(rep(c(trans.psi1,trans.psi2),n)))
          dgg        <- chol2inv(chol(as.matrix(dg.tp)))
          expr1      <- dgg%*%R_cv
          Sigv       <- chol2inv(chol(R_star+t(R_cv)%*%expr1))
          diff       <- as.numeric(Y)-as.matrix(kronecker(rep(1,n),beta[i,]))
          diff.ps    <- t(diff)%*%expr1
          logdet     <- sum(log(eigen(Sigv, symmetric = TRUE, only.values = TRUE)$values))+
                        sum(log(eigen(dgg, symmetric = TRUE, only.values = TRUE)$values))

          out        <- logdet/2-(t(diff)%*%dgg%*%diff-diff.ps%*%Sigv%*%t(diff.ps))/2+
                        theta[1]+theta[2]+log(densigamma(trans.psi1,al,bet))+log(densigamma(trans.psi2,a,b))
          out
     }

## metrop stuff
     NITER       <- 1
     BatchLength <- 1

## tuning
     tuning      <- diag(rep(.01,2))

## starting values
     inits      <- log(psi[i-1,])
     metrop.out <- metrop(target, initial=inits, nbatch=NITER, blen=BatchLength, scale=tuning)
     psi[i,]    <- exp(metrop.out$batch[1,])                                                     ## error variances updated
     count[i]   <- ((i-1)*count[i-1]+metrop.out$accept)/i                                        ## acceptance rate

## update sigma1,sigma2,rho,phi1

     u.rho <- .8
     l.rho <- .05

     target1<-function(theta){

         trans.sig1                                  <- exp(theta[1])
         trans.sig2                                  <- exp(theta[2])
         trans.phi1                                  <- (u.phi1-l.phi1)*(exp(theta[3])/(1+exp(theta[3])))+l.phi1
         trans.rho                                   <- exp(theta[4])/(1+exp(theta[4]))*(u.rho-l.rho)+l.rho


         DD1                                         <- trans.sig1*Matern(dist.loc, range=1/trans.phi1, smoothness=sm1)
         DD2                                         <- trans.sig2*Matern(dist.loc, range=1/trans.phi1, smoothness=sm2)
         DD12                                        <- sqrt(trans.sig1*trans.sig2)*trans.rho*
                                                        Matern(dist.loc, range=1/trans.phi1, smoothness=sm12)

         R                                           <- matrix(0,2*n,2*n)
         R[seq(1,2*n,2),seq(1,2*n,2)]                <- DD1
         R[seq(2,2*n,2),seq(2,2*n,2)]                <- DD2
         R[seq(1,2*n,2),seq(2,2*n,2)]                <- DD12
         R[seq(2,2*n,2),seq(1,2*n,2)]                <- DD12

         R_star1                                     <- trans.sig1*Matern(dist.knot, range=1/trans.phi1, smoothness=sm1)
         R_star2                                     <- trans.sig2*Matern(dist.knot, range=1/trans.phi1, smoothness=sm2)
         R_star12                                    <- trans.rho*sqrt(trans.sig1*trans.sig2)*
                                                        Matern(dist.knot, range=1/trans.phi1, smoothness=sm12)
         R_star                                      <- matrix(0,2*n.knot,2*n.knot)
         R_star[seq(1,2*n.knot,2),seq(1,2*n.knot,2)] <- R_star1
         R_star[seq(2,2*n.knot,2),seq(2,2*n.knot,2)] <- R_star2
         R_star[seq(1,2*n.knot,2),seq(2,2*n.knot,2)] <- R_star12
         R_star[seq(2,2*n.knot,2),seq(1,2*n.knot,2)] <- R_star12
         R_star.inv <- chol2inv(chol(R_star))

         R_cv1                                       <- trans.sig1*Matern(dist.loc.knot, range=1/trans.phi1, smoothness=sm1)
         R_cv2                                       <- trans.sig2*Matern(dist.loc.knot, range=1/trans.phi1, smoothness=sm2)
         R_cv12                                      <- trans.rho*sqrt(trans.sig1*trans.sig2)*
                                                        Matern(dist.loc.knot, range=1/trans.phi1, smoothness=sm12)

         R_cv                                        <- matrix(0,2*n,2*n.knot)
         R_cv[seq(1,2*n,2),seq(1,2*n.knot,2)]        <- R_cv1
         R_cv[seq(2,2*n,2),seq(2,2*n.knot,2)]        <- R_cv2
         R_cv[seq(1,2*n,2),seq(2,2*n.knot,2)]        <- R_cv12
         R_cv[seq(2,2*n,2),seq(1,2*n.knot,2)]        <- R_cv12

         diagon.tpp                                  <- (R-R_cv%*%R_star.inv%*%t(R_cv))*tap.mat
         dg.tp                                       <- as.matrix.csr(diagon.tpp+diag(rep(psi[i,],n)))
         dgg                                         <- chol2inv(chol(as.matrix(dg.tp)))

         expr1                                       <- dgg%*%R_cv
         Sigv                                        <- chol2inv(chol(R_star+t(R_cv)%*%expr1))
         diff                                        <- as.numeric(Y)-as.matrix(kronecker(rep(1,n),beta[i,]))
         diff.ps                                     <- t(diff)%*%expr1

         logdet                                      <- sum(log(eigen(Sigv, symmetric = TRUE, only.values = TRUE)$values))+
                                                        sum(log(eigen(dgg, symmetric = TRUE, only.values = TRUE)$values))-
                                                        sum(log(eigen(R_star, symmetric = TRUE, only.values = TRUE)$values))

          out                                        <- logdet/2-(t(diff)%*%dgg%*%diff-diff.ps%*%Sigv%*%t(diff.ps))/2+theta[3]-2*log(1+exp(theta[3]))
                                                        theta[1]+theta[2]+log(densigamma(trans.sig1,al,bet))+log(densigamma(trans.sig2,al,bet))+
                                                        theta[4]-2*log(1+exp(theta[4]))
          out
    }



## tuning
    tuning     <- diag(rep(.02,4))


##starting values
    inits      <- c(log(sigma1[i-1]),log(sigma2[i-1]),log((phi1[i-1]-l.phi1)/(u.phi1-phi1[i-1])),
                   log((rho12[i-1]-l.rho)/(u.rho-rho12[i-1])))
    metrop.out <- metrop(target1, initial=inits, nbatch=NITER, blen=BatchLength, scale=tuning)
    sigma1[i]  <- exp(metrop.out$batch[1,1])                                                           ## error variance dimension 1 updated
    sigma2[i]  <- exp(metrop.out$batch[1,2])                                                           ## error variance dimension 2 updated
    phi1[i]    <- exp(metrop.out$batch[1,3])/(1+exp(metrop.out$batch[1,3]))*(u.phi1-l.phi1)+l.phi1     ## range parameter updated  
    rho12[i]   <- exp(metrop.out$batch[1,4])/(1+exp(metrop.out$batch[1,4]))*(u.rho-l.rho)+l.rho        ## correlation parameter updated
    count2[i]  <- ((i-1)*count2[i-1]+metrop.out$accept)/i                                              ## acceptance rate

    if((i %% 10 ==0)&&(i>=10)){ cat("betta, psi , sigma,rho,phi",c(beta[i,],psi[i,],sigma1[i],sigma2[i],rho12[i],phi1[i]),"\n")}

}

