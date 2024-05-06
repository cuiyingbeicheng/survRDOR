get.S_td_t1_vec <- function(t1,td,X1,X2,delta1,delta2,h){
  num1 <- length(t1)

  csurv.x <- NULL
  for(i in 1:num1){
    kern.weight <- stats::dnorm((X1-t1[i])/h) * delta1
    kern.weight <- kern.weight/sum(kern.weight)

    surv <- survival::Surv(time = X2 - X1, event = as.numeric(delta2 == 1))
    csurv.fit <- survival::survfit(surv~1,type='kaplan-meier',weights = kern.weight)
    ind.trun <- (csurv.fit$surv >= 0)
    csurv.fit$time <- c(0,csurv.fit$time[ind.trun])
    csurv.fit$surv <- c(1,csurv.fit$surv[ind.trun])
    L.tpt = length(csurv.fit$time)

    csurv.indx <- sum(td[i] >= csurv.fit$time)
    csurv.x <- rbind(csurv.x,csurv.fit$surv[csurv.indx])
  }

  return(csurv.x)
}
