get.varS_NP_tD <- function(td,X1,X2,delta1,delta2,h = 0.075,tau){
  ns <- length(X1)
  ntd <- length(td)
  t1.grid <- summary(survival::survfit(survival::Surv(X1, delta1)~1))$time

  S.DR.td <- get.S_td_t1(t1 = X1,td = td,X1,X2,delta1,delta2,h = h)
  S.T2mT1.td <- array(rep(get.S_td_t1_vec(t1 = X1,td = X2 - X1,X1,X2,delta1,delta2,h = h),ntd),c(ns,ntd))
  G.C.T2 <- array(rep(get.survC(s = X2,X2,delta2),ntd),c(ns,ntd))
  G.C.T1 <- array(rep(get.survC(s = X1,X2,delta2),ntd),c(ns,ntd))
  S.R.T1 <- array(rep(get.S_R_t1(t1 = X1,X1,delta1),ntd),c(ns,ntd))

  ind1.matrix <- 1 * (array(rep(X1,ntd),c(ns,ntd)) < (tau - t(array(rep(td,ns),c(ntd,ns)))))
  ind2.matrix <- array(rep(delta2,ntd),c(ns,ntd))
  ind3.matrix <- array(rep(delta1,ntd),c(ns,ntd))
  ind4.matrix <- 1 * (array(rep(X2-X1,ntd),c(ns,ntd)) <= (t(array(rep(td,ns),c(ntd,ns)))))
  ind.part1 <- ind1.matrix*ind2.matrix*ind4.matrix
  S.T2mT1.td <- ifelse(ind.part1 == 0, 1, S.T2mT1.td)
  G.C.T2 <- ifelse(ind.part1 == 0, 1, G.C.T2)
  if(sum(S.T2mT1.td==0)>0){
    #warning(paste0(sum(S.T2mT1.td==0),"/",length(S.T2mT1.td)," points removed due to zero denominator!"))
    S.T2mT1.td <- ifelse(S.T2mT1.td == 0,NA,S.T2mT1.td)
  }else{
    S.T2mT1.td <- S.T2mT1.td
  }
  res.part1 <- (1/(ns^2)) * apply(ind.part1*(S.DR.td^2)/((S.T2mT1.td^2)*(G.C.T2^2)),2,sum,na.rm = TRUE)

  intgral.T1 <- get.S_NP_tD(td,X1,X2,delta1,delta2,lower.bound = X1,h = h,tau = tau)$S.NP_tD

  ind.part2 <- ind1.matrix*ind3.matrix
  S.R.T1 <- 1 * (ind.part2 == 0) + S.R.T1 * (ind.part2 != 0)#ifelse(ind.part2 == 0, 1, S.R.T1)
  G.C.T1 <- 1 * (ind.part2 == 0) + G.C.T1 * (ind.part2 != 0)#ifelse(ind.part2 == 0, 1, G.C.T1)
  if(sum(S.R.T1==0)>0){
    #warning(paste0(sum(S.R.T1==0),"/",length(S.R.T1)," points removed due to zero denominator!"))
    S.R.T1.d <- ifelse(S.R.T1 == 0,NA,S.R.T1)
  }else{
    S.R.T1.d <- S.R.T1
  }
  res.part2 <- (1/(ns^2)) * apply(ind.part2*((S.DR.td*S.R.T1-intgral.T1)^2)/((S.R.T1.d^2)*(G.C.T1^2)),2,sum,na.rm = TRUE)

  res <- res.part1 + res.part2

  I1.mu.mat <- matrix(0, ns, ntd)
  for(m in 1:ns){
    G.C.I1mu.temp <- get.survC(s = X1[m]+td,X2,delta2)
    for(b in 1:ntd){
      if(delta2[m]==1 & td[b]<= tau-X1[m] & S.DR.td[m, b]>0){
        subgroup1 <- (td>=td[b] & td<=tau-X1[m])
        if(sum(subgroup1)>0){
          probc=G.C.I1mu.temp[b]
          td1=td[subgroup1]
          gap1=c(td1[-1], tau-X1[m])-td1
          I1.mu.mat[m,b]=sum(gap1*S.DR.td[m,subgroup1])^2/(S.DR.td[m, b]*probc)^2
        }
      }
    }
  }


  time1.jump <- t1.grid
  F.R <- 1 - get.S_R_t1(t1 = time1.jump,X1,delta1)
  prob0 <- F.R - c(0,F.R[-length(F.R)])
  prob0.X1 <- rep(0, ns)
  for(m in 1:ns)
  {if(delta1[m]==1) prob0.X1[m]=prob0[which.min(abs(time1.jump-X1[m]))]}


  I2.mu.mat=rep(0, ns)
  for(m in 1:ns)
  {if(delta1[m]==1 && S.R.T1[m,1]>0 && G.C.T1[m,1]>0){
    lower= G.C.T1[m,1]^2*S.R.T1[m,1]^2

    subgroup1=(td < tau-X1[m])
    td1=td[subgroup1]
    gap1=c(td1[-1], tau-X1[m])-td1
    upper1=S.R.T1[m,1]*sum(S.DR.td[m, subgroup1]*gap1)

    upper2=sum(intgral.T1[m,subgroup1]*gap1)

    I2.mu.mat[m]=(upper1-upper2)^2/lower}
  }

  res.mu=mean(I2.mu.mat)/ns+mean(I1.mu.mat)/ns

  return(list(var.curve=res,var.mu = res.mu))
}
