get.SRDOR_IPW.improved <- function(X1, X2, delta1, delta2, tau, mono=TRUE, M=100){
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
  #            Improved IPCW                                   #
  ##############################################################

  # Vector of t2 at which S_IPW changes value
  tdn <- unique(pmax(as.vector(outer(tC,X1[delta1==1],"-")),0))
  ## t's at which G_C(X_1+t) changes value

  ## combine with t1
  t2 <- unique(sort(c(t1,tdn)))
  m2 <- length(t2)

  if (m2>M){
    t2 <- seq(0,tau,length=M)
    m2 <- M
  }


  X1tMAT <- outer(X1,t2,"+")

  X1tMATr <- round(cdfC(X1tMAT)*nC)

  GCX1t <- matrix(survC[X1tMATr],n,m2)


  DtMAT2 <- outer(Dt,t2,">") # n by m matrix of I(Dt > t1)
  surv2MAT <- DtMAT2/GCX1t

  surv2 <- colMeans(surv2MAT)

  if (mono){
    t2m <- c(-0.1,t2)

    Q <- c(0,cumsum((1-surv2)*diff(t2m)))
    surv2m <-  1- CM(t2m,Q)

    surv2m <- pmin(1,pmax(surv2m,0))
  }else{
    surv2m <- NULL
  }
  ##############################################################
  ###    Variance estimation for improved IPCW         ###
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

  IF2CMAT <- matrix(0,n,m2)


  for (j in 1:l){
    ## improved IPW
    kappatu <- colMeans(surv2MAT*(X1tMAT>=Ut[j]))
    int2tu <- kappatu/pi2[j]
    IF2CMAT <- IF2CMAT + outer(dMC[,j],int2tu)
  }

  IF2MAT <- surv2MAT - matrix(rep(surv2,each=n),n,m2) + IF2CMAT


  V2 <- colMeans(IF2MAT^2)/n
  se2 <- sqrt(V2)

  ## estimate the mean DOR and its SE

  gap2=c(t2[-1], tau)-t2
  mu2=sum(surv2*gap2)

  mu.v2=mean((apply(t(t(IF2MAT)*gap2), 1, sum))^2)/n

  mu.se2=sqrt(mu.v2)

  return(list(td=t2,est.curve=surv2,est.mu=mu2,se.curve=se2,se.mu=mu.se2))
}
