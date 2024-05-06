CM<-function(G,Q){
  m<-length(G)
  y<-rep(NA,m-1)
  j<-1
  while(j<m){
    slopes<-rep(NA,m-j)
    for (k in (j+1):m){
      slopes[k-j]<-(Q[k]-Q[j])/(G[k]-G[j])
    }
    j1<-which.min(slopes)+j
    y[j:(j1-1)]<-min(slopes)
    j<-j1
  }
  return(y)
}
