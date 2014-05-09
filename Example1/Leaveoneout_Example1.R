
dyn.load('../mproot-Mac/mproot.so')
source('../mproot-Mac/mproot.R')

logit<-function(p){return(log(p/(1-p)))}
ilogit<-function(x){return(exp(x)/(exp(x)+1))}
dat.1<-read.table('summary_bs10_3000iter_0.2inc.csv',header=TRUE,sep=',',
                  na.strings='-1')
dat.2<-read.table('summary_bs10_3000iter.csv',header=TRUE,sep=',',
                  na.string='-1')

dat<-rbind(dat.1,dat.2)
#  Init.Err proportional to the amount of the initial deformation
#  Res.Err  mesaured misalignment after registration
# 10 subjects, 16 initial misalignments, 10 experiments/misalignment
#  Experiment involves regist. with a different seed  ** 160 seeds **
#
dat$CaseID<-paste("C",dat$CaseID,sep='')
dat$CaseID<-as.factor(dat$CaseID)
#dat$Init.Err<-round(dat$Init.Err) <-- Mark: I am not sure why this was here?
dat<-dat[dat$Init.Err>0 & !is.na(dat$Res.Err),]
cases<-levels(dat$CaseID)
dat$Failure<-rep(FALSE,dim(dat)[1])
h<-3.0

for(case in cases){
    idx<-dat$CaseID==case
    dat$Failure[idx]<-dat$Res.Err[idx]>0.5 
}

cases.all<-cases
for(skip in 1:length(cases.all)){
  y<-NULL; s2e<-NULL
  cases<-cases.all[-skip]
  for(case in cases){
    idx<-dat$CaseID==case
    fit<-glm(Failure~Init.Err,family='binomial',data=dat,
             subset=idx)
    x<-seq(1,15,0.25); x<-data.frame(x); names(x)<-"Init.Err"
    pred<-predict(fit,newdata=x,type='response',se.fit=TRUE)
    y<-cbind(y,pred$fit)
    s2e<-cbind(s2e,pred$se.fit^2)
  }

  n.x<-dim(y)[1]; n.case<-dim(y)[2]; s2b<-NULL
  y<-log(y/(1-y))
  for(i in 1:n.x){
    z<-y[i,]
    var<-s2e[i,]
    root<-mproot(z,var,n.case-1)
    s2b<-c(s2b,root$root)
  }


  s2b<-matrix(s2b,n.x,n.case)
  w<-s2b/(s2b+s2e)
  w[s2b==0]<-1/s2e[s2b==0]
  w.sum<-matrix(apply(w,1,'sum'),n.x,n.case)
  w<-w/w.sum
  y.hat<-matrix(apply(y*w,1,'sum'),n.x,n.case)
  se.mean<-sqrt(apply((y-y.hat)^2/(s2b+s2e)^2,1,'sum'))/
              apply(1/(s2b+s2e),1,'sum')
  x<-x[,1]; s2b<-s2b[,1]; y.hat<-y.hat[,1]
  sd<-sqrt(s2b+se.mean^2)
  k<-qnorm(0.95)
  ucb<-ilogit(y.hat+k*sd); lcb<-ilogit(y.hat-k*sd); y.hat<-ilogit(y.hat)

  tol.lim<-approx(ucb,x,xout=0.1)$y

  print(paste('Tolerance limit: ',tol.lim))

  ii<-dat$CaseID==cases.all[skip]
  fail<-dat$Failure[ii]
  par<-dat$Init.Err[ii]
  par.vals<-unique(par[par<=tol.lim])
  par.vals<-unique(par)
  p.fail<-NULL
  for(jj in 1:length(par.vals)){
    p.fail<-c(p.fail,sum(fail[par==par.vals[jj]])/length(fail[par==par.vals[jj]]))
  }
  print(paste(" Proportion failed for case ",cases.all[skip]))
  o<-order(par.vals)
  print(cbind(par.vals[o],p.fail[o]))
}
