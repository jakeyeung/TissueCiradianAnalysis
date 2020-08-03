# Felix Naef
# 2015-10-28
# Optimal partitioning of data on an interval
# coded as per http://ieeexplore.ieee.org/xpls/icp.jsp?arnumber=1381461

make.block=function(a,b) outer(a,b,"+")

S=function(M)
{
  M=as.matrix(M)
  n=nrow(M)
  m=ncol(M)
  
  b=apply(M,2,mean)
  a=apply(M-matrix(b,n,m, byrow=T), 1, mean)
  
  # list(a,b)
  sum((M-make.block(a,b))^2)
}

seg=function(M, min.size=5, max.size=15)
  # max size may need tuning
{
  n=nrow(M)
  m=ncol(M)
  
  if(m<2*min.size) res=list(s=S(M), br=m)
  else
  {
    s=S(M)
    br=m
    for(l in min.size:min(max.size,m-min.size))
    {
      # cat(l, "\n")
      ind=1:l
      f1=seg(M[, ind])
      f2=seg(M[,-ind])
      ff=f1$s+f2$s
      # cat(ff, "\n")
      if(ff<s){br=c(f1$br, f2$br); s=ff}
    }
    res = list(s=s, br=br)
  }
  res
}

S2=function(j, M, min.size=5, max.size=20, BIC=T, s2)
  # s2 really needs to be estimated well.
  # too large and you are too stringent with the number of allowed params (not enough blocks)
  # too small and you are allow too many params (too many blocks)
{
  n=nrow(M)
  m=ncol(M)
  M=as.matrix(M[,j:m])
  m=ncol(M)
  if(m<min.size | m>max.size){res=Inf}
  else{
    b=apply(M,2,mean)
    a=apply(M-matrix(b,n,m, byrow=T), 1, mean)
    
    # list(a,b)
    res=sum((M-make.block(a,b))^2)
    if(BIC) res = res + s2 * (n+m-1)*log(n*m)
  }
  res
}

seg.dp=function(M, min.size=5, max.size=20, BIC=T, sigma2=0)
{
  # if BIC is TRUE: returns BIC value
  # if BIC is FALSE: returns RSS (an estimate of variance you can think)
  if(is.null(sigma2) & BIC) sigma2=var(as.vector(M))
  print(paste("Sigma2:", sigma2))
  n=nrow(M)
  m=ncol(M)
  
  scores	= rep(0,m)
  change	= rep(0,m)
  
  for (k in 1:m)
  {
    if(k==1) s = c(0)
    if(k>1)  s = c(0, scores[1:(k-1)])
    s.tmp = s + sapply(1:k, S2, M=as.matrix(M[,1:k]), min.size=min.size, max.size=max.size, BIC=BIC, s2=sigma2)
    jmin = which.min(s.tmp)
    scores[k] = s.tmp[jmin]
    change[k] = jmin
  }
  
  #backtrack
  bk=change[m]
  while(bk[1]>1)
  {
    bk=c(change[bk[1]-1],bk)
  }
  sizes=diff(c(bk, ncol(M)+1))
  scores=scores/(n*m)
  list(sizes=sizes, start.pos=bk, changes=change, scores=scores)
}

# #main
# 
# n=100
# M1=make.block(rnorm(n), rnorm(12))
# M2=make.block(rnorm(n), rnorm(5))
# M3=make.block(rnorm(n), rnorm(5))
# M4=make.block(rnorm(n), rnorm(12))
# M5=make.block(rnorm(n), rnorm(5))
# 
# # M=cbind(M1,M2)
# M=cbind(M1, M2, M3, M4, M5)
# M=M + rnorm(length(M))/10
# 
# # res=seg(M)
# # print(res)
# 
# res=seg.dp(M, sigma2=NULL, BIC=T)
# print(res)
