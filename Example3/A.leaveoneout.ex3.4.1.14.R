#
# Mark Vangel, Massachusetts General Hospital
# vangel@nmr.mgh.harvard.edu
# 2013
#
library(locfit)
dyn.load('~/Dropbox/Registration.Prostate/2013-TMIRevisionWork/IntraopCasesSimulations-2013-10-22/mproot.so')
source('~/Dropbox/Registration.Prostate/2013-TMIRevisionWork/IntraopCasesSimulations-2013-10-22/mproot-Mac/mproot.R')
# load the data
dat1<-read.table(
    'Dir14mD0MD10mrX-15MrX15mrY-15MrY15mrZ-15MrZ15dD2dr3-Case1.csv',
    header=TRUE,sep=',',na.strings="NA")
dat2<-read.table(
    'Dir14mD0MD10mrX-15MrX15mrY-15MrY15mrZ-15MrZ15dD2dr3-Case2.csv',
    header=TRUE,sep=',',na.strings="NA")
dat3<-read.table(
    'Dir14mD0MD10mrX-15MrX15mrY-15MrY15mrZ-15MrZ15dD2dr3-Case3.csv',
    header=TRUE,sep=',',na.strings="NA")
dat4<-read.table(
    'Dir14mD0MD10mrX-15MrX15mrY-15MrY15mrZ-15MrZ15dD2dr3-Case4.csv',
    header=TRUE,sep=',',na.strings="NA")
dat5<-read.table(
    'Dir14mD0MD10mrX-15MrX15mrY-15MrY15mrZ-15MrZ15dD2dr3-Case5.csv',
    header=TRUE,sep=',',na.strings="NA")
dat6<-read.table(
    'Dir14mD0MD10mrX-15MrX15mrY-15MrY15mrZ-15MrZ15dD2dr3-Case6.csv',
    header=TRUE,sep=',',na.strings="NA")
dat7<-read.table(
    'Dir14mD0MD10mrX-15MrX15mrY-15MrY15mrZ-15MrZ15dD2dr3-Case7.csv',
    header=TRUE,sep=',',na.strings="NA")
dat8<-read.table(
    'Dir14mD0MD10mrX-15MrX15mrY-15MrY15mrZ-15MrZ15dD2dr3-Case8.csv',
    header=TRUE,sep=',',na.strings="NA")
dat9<-read.table(
    'Dir14mD0MD10mrX-15MrX15mrY-15MrY15mrZ-15MrZ15dD2dr3-Case9.csv',
    header=TRUE,sep=',',na.strings="NA")
dat10<-read.table(
    'Dir14mD0MD10mrX-15MrX15mrY-15MrY15mrZ-15MrZ15dD2dr3-Case10.csv',
    header=TRUE,sep=',',na.strings="NA")

n<-dim(dat1)[1]
dat1$Sample<-rep(1,n)
dat2$Sample<-rep(2,n)
dat3$Sample<-rep(3,n)
dat4$Sample<-rep(4,n)
dat5$Sample<-rep(5,n)
dat6$Sample<-rep(6,n)
dat7$Sample<-rep(7,n)
dat8$Sample<-rep(8,n)
dat9$Sample<-rep(9,n)
dat10$Sample<-rep(10,n)

nm<-names(dat1); dat1<-dat1[,nm!="Error.4"]
nm<-names(dat4); dat4<-dat4[,nm!="Error.4"]
nm<-names(dat5); dat5<-dat5[,nm!="Error.4"]
nm<-names(dat6); dat6<-dat6[,nm!="Error.4"]


dat<-rbind(dat1,dat2,dat3,dat4,dat5,dat6,dat7,dat8,dat9,dat10)
sample.ids<-c(1,2,3,4,5,6,7,8,9,10)
dat<-rbind(dat1,dat2,dat3,dat4,dat5,dat6,dat7,dat8,dat10)
sample.ids<-c(1,2,3,4,5,6,7,8,10)
#
t<-3
dat$Error.1<- is.na(dat$Error.1) | ((!is.na(dat$Error.1) & (dat$Error.1 >t)))
dat$Error.2<- is.na(dat$Error.2) | ((!is.na(dat$Error.2) & (dat$Error.1 >t)))
dat$Error.3<- is.na(dat$Error.3) | ((!is.na(dat$Error.3) & (dat$Error.1 >t)))
max.rot<-7
dat<-subset(dat, abs(dat$Rot.X)<max.rot & abs(dat$Rot.Y)<max.rot & abs(dat$Rot.Z)<max.rot)
dat$Failure<-dat$Error.1 | dat$Error.2 | dat$Error.3 #| dat$Error.4
dat$Success<-!dat$Failure
dat$Rotation<-paste(dat$Rot.X,dat$Rot.Y,dat$Rot.Z,sep='.')
dat$Direction<-as.factor(dat$Direction)
dat$Rotation<-as.factor(dat$Rotation)
levels(dat$Direction)<-c(
    "+x","-x", "+y","-y", "+z","-z",  
    "+x+y+z", "+x+y-z", "+x-y+z", "-x+y+z",
     "-x-y-z", "+x-y-z","-x+y-z","-x-y+z")
#
n.samp<-500; n.x<-100; n.case<-length(sample.ids)
disp<-c(0,2,4,6,8,10)
tbl<-NULL
for(sample in sample.ids){
   for(d in disp){
       y<-sample(dat$Success[dat$Displacement==d &
                             dat$Sample==sample],
                 n.samp,replace=TRUE)
       tbl<-rbind(tbl,cbind(sample,d,y))
   }}
tbl<-data.frame(tbl)
names(tbl)<-c("Sample","Misalign","Success")
dat<-tbl


logit<-function(p){return(log(p/(1-p)))}
ilogit<-function(x){return(exp(x)/(exp(x)+1))}

pdf(paste('Illustration3-',n.samp,'samp-',length(sample.ids),'cases.pdf',sep=''))

#plot(0,0,type='n',xlim=c(1,10),ylim=c(0,1),
#     xlab='misalignment magnitude, mm',
#     ylab='probability of failure')


cases.all<-sample.ids
for(skip in cases.all){
  var.within<-NULL; p<-NULL; fit.tbl.1<-NULL; s2e<-NULL; s2b<-NULL
  n.case<-length(sample.ids)-1
  x0<-seq(0,n.case,length=n.x)
  cases<-cases.all[cases.all!=skip]
  for(i in 1:n.case){
    idx<-dat$Sample==sample.ids[sample.ids!=skip][i]
    y<-dat$Success[idx]
    x<-dat$Misalign[idx]
    o<-order(x)
    x<-x[o]
    y<-y[o]
    y<-!y

    fit<-locfit(y~lp(x),family='binomial')
    # predicted probability of failure for each of the misalignment values
    yhat<-predict(fit,x0)
    fit.tbl.1<-cbind(fit.tbl.1,yhat)
    # p is the logit transformed probability
    p<-cbind(p,logit(yhat))
    # standard error
    s2e<-cbind(s2e,predict(fit,x0,se.fit=TRUE)$se.fit^2)

    # include lines for individual cases
    pred<-predict(fit,newdata=x0,se.fit=TRUE,type='response')
#    lines(x0,pred$fit,type='l')
}
# s2b is across sample standard error
   for(i in 1:n.x){
     s2b<-c(s2b,mproot(p[i,], s2e[i,], n.case-1)$root)
   }

   w<-NULL
   for(i in 1:n.x){
     if(s2b[i]!=0){
        w<-rbind(w,s2b[i]/(s2b[i]+s2e[i,]))}
     if(s2b[i]==0){
        w<-rbind(w,1/s2e[i,])}
}
   w.sum<-matrix(apply(w,1,'sum'),n.x,n.case)
   w.norm<-w/w.sum
   p.hat<-apply(w.norm*p,1,'sum')
   p.hat<-matrix(p.hat,n.x,n.case); s2b<-matrix(s2b,n.x,n.case)
   se.mean<-sqrt(apply((p-p.hat)^2/(s2b+s2e)^2,1,'sum'))/
              apply(1/(s2b+s2e),1,'sum')
   s2b<-s2b[,1]; p.hat<-p.hat[,1]
   sd<-sqrt(s2b+se.mean^2)
   k<-qnorm(0.95)
   ucb<-ilogit(p.hat+k*sd)
   lcb<-ilogit(p.hat-k*sd)
   tol.lim<-approx(ucb[x0>1],x0[x0>1],xout=0.1)$y
   fail<-!tbl$Success[tbl$Sample==skip]
   mis<-tbl$Misalign[tbl$Sample==skip]
   p.test<-sum(fail[mis<=tol.lim])/sum(mis<=tol.lim)
   print(c(skip,tol.lim,p.test))
   p.hat<-ilogit(p.hat)
   plot(x0,p.hat,type='l',lwd=2,
        ylim=range(ucb,lcb),
        xlab='misalignment norm, mm',
        ylab='probability of failed registration')
#     main=paste('Average Probability of Failed Registration with\n95% Upper Prediction Limit,\nand (0.90, 0.95) Lower Tolerance Limit = ',round(tol.lim,1)))
   lines(x0,ucb,lty=2,lwd=2)
   lines(x0,lcb,lty=2,lwd=2)
   abline(h=0.1,col='blue',lwd=2)
   abline(v=tol.lim,col='blue',lwd=2)
dev.off()
}
      



    

    
