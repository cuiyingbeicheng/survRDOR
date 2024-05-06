get.S_R_t1 <- function(t1,X1,delta1){
  surv_R <- survival::Surv(time = X1, event = as.numeric(delta1 == 1))
  csurv_R.fit <- survival::survfit(surv_R~1,type='kaplan-meier')
  ind.trun_R <- (csurv_R.fit$surv >= 0)
  csurv_R.fit$time <- c(0,csurv_R.fit$time[ind.trun_R])
  csurv_R.fit$surv <- c(1,csurv_R.fit$surv[ind.trun_R])
  L.tpt_R = length(csurv_R.fit$time)

  num1_R <- length(t1)
  csurv_R.indx <- apply( t1 >= t(array(rep(csurv_R.fit$time, num1_R), c(L.tpt_R, num1_R))), 1, sum)
  csurv_R.x <- csurv_R.fit$surv[csurv_R.indx]
  return(csurv_R.x)
}
