get.survC <- function(s,X2,delta2){
  surv_C <- survival::Surv(time = X2, event = as.numeric(delta2 == 0))
  csurv_C.fit <- survival::survfit(surv_C~1,type='kaplan-meier')
  ind.trun_C <- (csurv_C.fit$surv >= 0)
  if(min(csurv_C.fit$time[ind.trun_C]) == 0){
    csurv_C.fit$time <- csurv_C.fit$time[ind.trun_C];csurv_C.fit$time[1] <- 0
    csurv_C.fit$surv <- csurv_C.fit$surv[ind.trun_C];csurv_C.fit$surv[1] <- 1
  }else{
    csurv_C.fit$time <- c(0,csurv_C.fit$time[ind.trun_C])
    csurv_C.fit$surv <- c(1,csurv_C.fit$surv[ind.trun_C])
  }
  L.tpt_C = length(csurv_C.fit$time)


  num1_C <- length(s)
  csurv_C.indx <- apply( s >= t(array(rep(csurv_C.fit$time, num1_C), c(L.tpt_C, num1_C))), 1, sum)
  csurv_C.x <- csurv_C.fit$surv[csurv_C.indx]
  return(csurv_C.x)
}
