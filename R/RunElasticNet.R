# RunElasticNet.R
# from Jingkui

RunElasticNet <- function(y, x, plot.out, standardize = TRUE, alpha){
  intercept=0
  standardize=standardize ### standardize matrix of motif occurrence makes more sense because the absolute number of motif occurrence is not precise.
  standardize.response=FALSE
  alpha = alpha
  binary.matrix = FALSE
  if(binary.matrix){x = x >0; standardize=FALSE}
  
  ### use Cross-validation to select tuning paprameter 
  cv.fit=cv.glmnet(x,y,family='mgaussian',grouped=TRUE, alpha=alpha, nlambda=500, standardize=standardize, standardize.response=standardize.response, intercept=intercept)
  plot(cv.fit)
  #cv.fit$lambda
  
  optimal = which(cv.fit$lambda==cv.fit$lambda.min)
  #optimal = which(cv.fit$lambda==cv.fit$lambda.1se)
  fit=glmnet(x,y,alpha=alpha, lambda=cv.fit$lambda,family='mgaussian', type.multinomial=c("grouped"), standardize=standardize, standardize.response=standardize.response, intercept=intercept)
  #fit=glmnet(x,y,alpha=alpha, family='mgaussian', type.multinomial=c("grouped"), standardize=standardize, standardize.response=standardize.response, intercept=intercept)
  #plot(fit)
  colnames(x)[which(fit$beta[[1]][,optimal]!=0)]
  
  
  # Plot results ------------------------------------------------------------
  
  fit.pol2 = fit
  optimal.pol2 = optimal
  ### processed the results of Lasso
  ###
  kk = c(which(fit.pol2$beta[[1]][,optimal.pol2]!=0),which(fit.pol2$beta[[2]][,optimal.pol2]!=0))
  kk = unique(kk)
  names.pol2 = colnames(x)[kk]
  a.pol2 = fit.pol2$beta[[1]][,optimal.pol2][kk]
  b.pol2 = fit.pol2$beta[[2]][,optimal.pol2][kk]
  aa = a.pol2 
  bb = b.pol2 
  period = 24
  keep.pol2 = c()
  for(n in 1:length(names.pol2))
  {
    phase=period/(2*pi)*atan2(bb[n], aa[n])
    if(phase<0) phase=phase+period
    if(phase>period) phase=phase-period
    amp = sqrt(aa[n]^2+bb[n]^2)
    keep.pol2 = rbind(keep.pol2, c(phase,amp))
  }
  rownames(keep.pol2) = names.pol2
  colnames(keep.pol2) = c('phase','ampl')
  
  infer = keep.pol2 # infered motifs by elastic-net
  
  ##### Plot phases of inferred Motifs
  motif.amp = 7.5;
  phase.m = infer[,1]
  
  o1 = order(phase.m)
  phase.m = phase.m[o1]
  motif.names = rownames(infer)[o1]
  
  amp = motif.amp
  motif.a = amp*cos(2*pi/24*phase.m)
  motif.b = amp*sin(2*pi/24*phase.m)
  CC = (motif.a -1i*motif.b) * exp(1i*pi/2)
  motif.aa = Re(CC)
  motif.bb = Im(CC)
  
  amp = motif.amp + 0.25
  motif.a = amp*cos(2*pi/24*phase.m)
  motif.b = amp*sin(2*pi/24*phase.m)
  CC = (motif.a -1i*motif.b) * exp(1i*pi/2)
  motif.txt1.aa = Re(CC)
  motif.txt1.bb = Im(CC)
  
  amp = motif.amp-0.25
  motif.a = amp*cos(2*pi/24*phase.m)
  motif.b = amp*sin(2*pi/24*phase.m)
  CC = (motif.a -1i*motif.b) * exp(1i*pi/2)
  motif.txt2.aa = Re(CC)
  motif.txt2.bb = Im(CC)
  
  rmm=max(abs(c(motif.aa,motif.bb)))+2.0
  rr=c(-rmm,rmm)
  xlim = rr
  ylim = rr
  
  pdf(paste(plot.out, alpha, '.pdf', sep=''), width=5, height=5)
  par(cex = 0.5, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,5,0.8)+0.1, tcl = -0.3)
  plot(motif.aa, motif.bb, main=paste('alpha = ', alpha, sep=''), type='n', xlim=xlim, ylim=ylim, axes=F, xlab='', ylab='', lwd=2, pch=23, col='black',bg='green',cex=2.0)
  
  abline(v=0,h=0, col='darkgray',lwd=2.0)
  phi=seq(0,2*pi,len=1000)
  #lines(rm*cos(phi), rm*sin(phi), col='darkgray', lwd=2)
  lines((motif.amp-0.5)*cos(phi), (motif.amp-0.5)*sin(phi), col='darkgray', lwd=2)
  rainbow = rainbow(length(motif.names),s = 0.85, v = 0.85)
  for(n in 1:length(motif.names))
  {
    points(motif.aa[n], motif.bb[n], pch=21, cex=1.5, col='black',bg=rainbow[n])
    
    if(phase.m[n]<=6) srt = 90-phase.m[n]/24*360;
    if(phase.m[n]>6 & phase.m[n]<=12) srt = 90-phase.m[n]/24*360;
    if(phase.m[n]>12 & phase.m[n]<=18) srt = 270-phase.m[n]/24*360;
    if(phase.m[n]>=18) srt = 270-phase.m[n]/24*360;
    
    if(phase.m[n]<=12) {
      if(n%%2==1){
        adj = 0; 
        pos=1;
        text(motif.txt1.aa[n], motif.txt1.bb[n], motif.names[n], cex=0.7, col='black', srt=srt, adj=adj)
      }else{
        adj = 1;
        text(motif.txt2.aa[n], motif.txt2.bb[n], motif.names[n], cex=0.7, col='black', srt=srt, adj=adj)
      }
      
    }else{
      if(n%%2==1)
      {
        adj = 1; 
        pos=3.0
        text(motif.txt1.aa[n], motif.txt1.bb[n], motif.names[n], cex=0.7, col='black', srt=srt, adj=adj)
      }else{
        adj = 0; 
        pos=1;
        text(motif.txt2.aa[n], motif.txt2.bb[n], motif.names[n], cex=0.7, col='black', srt=srt, adj=adj)
      }
    }
  }
  dev.off()
}