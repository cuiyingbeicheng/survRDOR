get.SRDOR_twostep <- function(td,X1,X2,delta1,delta2,h = 0.075,tau){
  # process data for tau
  Y1 <- pmin(X1, tau)
  Y2 <- pmin(X2, tau)
  d1 <- delta1
  d2 <- delta2
  d1[Y1==tau] <- 1
  d2[Y2==tau] <- 1

  # estimation
  res.est <- get.S_NP_tD(td,Y1,Y2,d1,d2,lower.bound = NULL,h,tau)
  res.var <- get.varS_NP_tD(td,Y1,Y2,d1,d2,h,tau)
  return(list(est.curve=res.est$S.NP_tD,est.mu=res.est$mu,se.curve=sqrt(res.var$var.curve),se.mu = sqrt(res.var$var.mu)))
}
