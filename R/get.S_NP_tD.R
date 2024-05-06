get.S_NP_tD <- function(td,X1,X2,delta1,delta2, lower.bound = NULL,h,tau){
  t1.grid <- summary(survival::survfit(survival::Surv(X1, delta1)~1))$time
  S.D_R <- get.S_td_t1(t1 = t1.grid,td,X1,X2,delta1,delta2,h = h)
  F.R <- 1 - get.S_R_t1(t1 = t1.grid,X1,delta1)

  if(is.null(lower.bound)){
    S.NP_tD <- apply(rbind(S.D_R,td), 2, func<-function(x,F.R){
      nx <- length(x)
      temp_tD <- x[nx]; temp_S.D_R <- x[-nx]
      temp_res <- sum((t1.grid <= tau - temp_tD) * (t1.grid >= 0) *
                        temp_S.D_R * (F.R - c(0,F.R[-(nx-1)])))
      return(temp_res)
    },F.R = F.R)
  }else{
    S.NP_tD <- NULL
    nlb <- length(lower.bound)
    for(lbi in 1:nlb){
      S.NP_tD <- rbind(S.NP_tD,apply(rbind(S.D_R,td), 2, func<-function(x,F.R){
        nx <- length(x)
        temp_tD <- x[nx]; temp_S.D_R <- x[-nx]
        temp_res <- sum((t1.grid <= tau - temp_tD) * (t1.grid >= lower.bound[lbi]) *
                          temp_S.D_R * (F.R - c(0,F.R[-(nx-1)])))
        return(temp_res)
      },F.R = F.R))
    }
  }

  mu <- sum(S.NP_tD * (c(td[-1],tau)-td))

  return(list(S.NP_tD=S.NP_tD,mu=mu))
}
