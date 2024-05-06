get.SRDOR_IPW.naive <- function(X1, X2, delta1, delta2, tau){
  ##############################################################
  #            Preparations                                    #
  ##############################################################
  # sample size
  n <- length(X1)

  # refined the data to incorporate the
  # restriction time tau
  X1[X1>X2] <- X2[X1>X2]
  X1 <- pmin(X1, tau)
  delta1[X1==tau] <- 1
  X2 <- pmin(X2, tau)
  delta2[X2==tau] <- 1

  ## define delta1 to indicate response or event without response
  delta1[delta2==1] <- 1

  # estimate censoring survival function
  # G_C(t)
  deltaC <- 1-delta2 # indicator for censoring
  # fit KM curves for G_C(t)
  fitC <- survival::survfit(survival::Surv(X2, deltaC) ~ 1)
  tC <- c(0, fitC$time)
  survC <- c(1, fitC$surv)
  cdfC <- ecdf(tC) # empr dist of timeC
  # can be used to get G_C(t) from survC for any t
  nC <- length(tC)
  if (survC[nC]==0){survC[nC]=survC[nC-1]}


  Dt <- pmax(X2 - X1, 0) # should require x>=X1 to begin with
  # ordered unique values of Dt
  t1 <- sort(unique(c(0,Dt)))
  m1 <- length(t1)

  ##############################################################
  #            Naive IPCW                                      #
  ##############################################################

  ## evaluate G_C(X2) ##
  X2r <- round(cdfC(X2)*nC) # ranks of X2 in tC
  GCX2 <- survC[X2r] # n-vector of G_C(X2)

  DtMAT <- outer(Dt,t1,">") # n by m matrix of I(Dt > t1)
  # dim(DtMAT)
  # ep <- 1e-12

  # naive IPW
  surv1MAT <- matrix(delta2/GCX2, n, m1)*DtMAT
  surv1 <- colMeans(surv1MAT)

  ##############################################################
  ###    Variance estimation for naive IPCW         ###
  ##############################################################
  ## censoring martingale M_C
  Ut <- sort(unique(X2[deltaC==1]))
  l <- length(Ut)

  eeq <- function(y1,y2){
    abs(y1-y2)<1e-8
  }

  dNC <- outer(X2,Ut,eeq)*matrix(rep(deltaC,l),n,l)

  R2MAT <- outer(X2,Ut,">=")
  pi2 <- colMeans(R2MAT)

  dLambdaC <- colMeans(dNC)/pi2

  dMC <- dNC - R2MAT*matrix(rep(dLambdaC,each=n),n,l)

  ## matrices hold the IF's of surv1 and surv2
  IF1CMAT <- matrix(0,n,m1)

  for (j in 1:l){
    ## naive estimator
    zetatu <-  colMeans(surv1MAT*matrix(rep(R2MAT[,j],m1),n,m1))
    int1tu <- zetatu/pi2[j]
    IF1CMAT <- IF1CMAT + outer(dMC[,j],int1tu)
    }

  IF1MAT <- surv1MAT - matrix(rep(surv1,each=n),n,m1) + IF1CMAT

  V1 <- colMeans(IF1MAT^2)/n
  se1 <- sqrt(V1)

  ## estimate the mean DOR and its SE

  gap1 <- c(t1[-1], tau)-t1
  mu1 <- sum(surv1*gap1)

  mu.v1 <- mean((apply(t(t(IF1MAT)*gap1), 1, sum))^2)/n

  mu.se1 <- sqrt(mu.v1)

  return(list(td=t1,est.curve=surv1,est.mu=mu1,se.curve=se1,se.mu=mu.se1))
}
