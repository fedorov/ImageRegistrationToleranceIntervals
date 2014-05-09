
                                        #
# Mark Vangel, Massachusetts General Hospital
# vangel@nmr.mgh.harvard.edu
# 2013
#

library(locfit)
dyn.load('../mproot-Mac/mproot.so')
source('../mproot-Mac/mproot.R')

# load the data
dat<-read.table('all.sims.2.6.12.csv',header=T,sep=',')
dat<-dat[dat$Type=="BSpline" & dat$Landmark==1,]


logit<-function(p){return(log(p/(1-p)))}
ilogit<-function(x){return(exp(x)/(exp(x)+1))}
samples.all<-1:10
for(skip in 1:10){
  var.within<-NULL; p<-NULL; fit.tbl.1<-NULL; s2e<-NULL; s2b<-NULL
  samples<-samples.all[-skip]
  for(i in samples){
    idx<-dat$Sample==i & dat$Landmark=="1" & dat$Type=="BSpline"
    y<-dat$Success[idx]
    x<-dat$Misalign[idx]
    o<-order(x)
    x<-x[o]
    y<-y[o]
    y<-!y

    fit<-locfit(y~lp(x),family='binomial')
    # predicted probability of failure for each of the misalignment values
    yhat<-predict(fit,x)
    fit.tbl.1<-cbind(fit.tbl.1,yhat)
    # p is the logit transformed probability
    p<-cbind(p,logit(yhat))
    # standard error
    s2e<-cbind(s2e,predict(fit,x,se.fit=TRUE)$se.fit^2)
  }

nsamp<-length(samples)
# s2b is across sample standard error
  for(i in 1:500){
     s2b<-c(s2b,mproot(p[i,], s2e[i,], nsamp-1)$root)
}

   w<-NULL
   for(i in 1:500){
     if(s2b[i]!=0){
        w<-rbind(w,s2b[i]/(s2b[i]+s2e[i,]))}
     if(s2b[i]==0){
        w<-rbind(w,1/s2e[i,])}
}
   w.sum<-matrix(apply(w,1,'sum'),500,length(samples))
   w.norm<-w/w.sum
   p.hat<-apply(w.norm*p,1,'sum')
   s2b<-matrix(s2b,500,length(samples))
   se.mean<-sqrt(apply((p-p.hat)^2/(s2b+s2e)^2,1,'sum'))/
              apply(1/(s2b+s2e),1,'sum')
   k<-qnorm(0.95)
   ucb<-ilogit(p.hat+k*se.mean)
   lcb<-ilogit(p.hat-k*se.mean)
   p.hat<-ilogit(p.hat)
   par(mfrow=c(2,1))

#sd<-sqrt(apply(s2e,1,'mean')+s2b[,1])
   sd<-se.mean *sqrt(length(samples))
   k<-qnorm(0.95)
   ucb<-ilogit(logit(p.hat)+k*sd)
   lcb<-ilogit(logit(p.hat)-k*sd)
   tol.lim<-approx(ucb[x>1],x[x>1],xout=0.1)$y
   fail<-!dat$Success[dat$Sample==skip]
   par<-dat$Misalign[dat$Sample==skip]
   p.val<-sum(fail[par<=tol.lim])/sum(par<=tol.lim)
   print(c(skip,tol.lim,p.val))
   print(paste('Tolerance interval:',tol.lim))
   }   


