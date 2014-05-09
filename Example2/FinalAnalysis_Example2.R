#
# Mark Vangel, Massachusetts General Hospital
# vangel@nmr.mgh.harvard.edu
# 2013
#
#
#     Local likelihood regression library 
library(locfit)
#
#    -- Object file corresponding to Fotran root-finding subroutine
dyn.load('../mproot-Linux/mproot.so')
#
#    -- R wrapper which calls external root-finder
source('../mproot-Mac/mproot.R')
#
#     -- Definition of logit and inverse-logit functions
logit<-function(p){return(log(p/(1-p)))}
ilogit<-function(x){return(exp(x)/(exp(x)+1))}

#    -- Load the data
dat<-read.table('all.sims.2.6.12.csv',header=T,sep=',')
#
#    -- BSpline only. Do not use multiple landmarks.
dat<-dat[dat$Type=="BSpline" & dat$Landmark==1,]

#
#       p:   ith column will be logit-scale prob. failure for ith sample
#       s2e: ith column will be with-sample SE^2 of logit for ith sample
#       s2b: between-sample variance. Identical column for each sample.
p<-NULL; s2e<-NULL; s2b<-NULL
#
#      -- Initialize Plot #1
pdf('Figure1.pdf')
plot(0,0,type='l',lwd=2,
     ylim=range(0,1),
     xlim=range(0,10),
     xlab='misalignment norm, mm',
     ylab='probability of failed fegistration')
#
#     y: probability of failure for ith case
#     x: misalignment norm for ith case
for(i in 1:10){
    idx<-dat$Sample==i 
    y<-dat$Success[idx]
    x<-dat$Misalign[idx]
#
#   -- Sort x, carry along y, negate y to turn successes into failures.    
    o<-order(x); x<-x[o]; y<-y[o];   y<-!y
#
#      -- Nonparametric logistic regression fit
    fit<-locfit(y~lp(x),family='binomial')
#
#    -- Predicted probability of failure for each of the misalignment values
    yhat<-predict(fit,x)
#
#   -- Logit transform fit, build matrix colummnwise
    p<-cbind(p,logit(yhat))
#
#    -- Similarly, build matrix of SE^2 of fitted logit probabilities
    s2e<-cbind(s2e,predict(fit,x,se.fit=TRUE)$se.fit^2)
#
#      -- Add fitted line to plot
    lines(x,yhat,lwd=2)
  }
dev.off()

nsamp<-10
#
#      -- Obtain between-case variance (s2b) for x_i, accross cases (row of y)
#           10 samples, 500 values of misalignment norm.
for(i in 1:500){
#
#      -- mproot is wrapper function which calls FORTRAN root-finder.
#           [i,] indicates the nsamp values in the ith row.    
  s2b<-c(s2b,mproot(p[i,], s2e[i,], nsamp-1)$root)
}
#
#      -- Weight matrix
w<-NULL
for(i in 1:500){
  if(s2b[i]!=0){
     w<-rbind(w,s2b[i]/(s2b[i]+s2e[i,]))}
  if(s2b[i]==0){
     w<-rbind(w,1/s2e[i,])}
}
w.sum<-matrix(apply(w,1,'sum'),500,10)
w.norm<-w/w.sum
#
#      -- Mean probability of failure, averaging over samples, is
#            a weighted mean, which is constructed here.
p.hat<-apply(w.norm*p,1,'sum')
#
#      -- Reshape the between-variance vector into a matrix by columns
s2b<-matrix(s2b,500,10)
#
#      -- Standard error of the mean, Rukhin and Vangel, JASA, 1999.
se.mean<-sqrt(apply((p-p.hat)^2/(s2b+s2e)^2,1,'sum'))/
              apply(1/(s2b+s2e),1,'sum')
#
#      -- Upper and lower confidence limits on inverse logit of mean
#           prob. of failure, suing normal approximation.
k<-qnorm(0.95)
ucb<-ilogit(p.hat+k*se.mean)
lcb<-ilogit(p.hat-k*se.mean)
p.hat<-ilogit(p.hat)
#
#     -- Initialize plot #2
pdf('Figure2.pdf')
#
#     -- Create plot with mean prob. of failure curve
plot(x,p.hat,type='l',lwd=2,
     ylim=range(ucb,lcb),
     xlab='misalignment norm, mm',
     ylab='probability of failed fegistration',
     main='Average Probability of Failed Registration with\n95% Upper Confidence Limit on Mean')
#
#      -- Add curves for upper and lower confidence limits
lines(x,ucb,lty=2,lwd=2)
lines(x,lcb,lty=2,lwd=2)
dev.off()
#
#       -- Upper and lower prediction intervals 
sd<-sqrt(se.mean^2+s2b[,1])
#sd<-se.mean *sqrt(10)
k<-qnorm(0.95)
ucb<-ilogit(logit(p.hat)+k*sd)
lcb<-ilogit(logit(p.hat)-k*sd)
#
#     -- Get tolerance limit by linearly interpolating intersection of
#          horizontal line at 0.1 with upper prediction interval. We
#          restrict to x>1 so that we find the largest intersecting x
#          value.
tol.lim<-approx(ucb[x>1],x[x>1],xout=0.1)$y
#
#    -- Initialize plot #3
pdf('Figure3.pdf')
#
#     -- Plot mean probability of failure curve
plot(x,p.hat,type='l',lwd=2,
     ylim=range(ucb,lcb),
     xlab='misalignment norm, mm',
     ylab='probability of failed registration',
     main=paste('Average Probability of Failed Registration with\n95% Upper Prediction Limit,\nand (0.90, 0.95) Lower Tolerance Limit = ',round(tol.lim,1)))
#
#      -- Upper and lower prediction limits
lines(x,ucb,lty=2,lwd=2)
lines(x,lcb,lty=2,lwd=2)
print(paste('Tolerance interval:',tol.lim))
#
#       -- Horizontal line at 0.1; vertical line at tolerance limit
abline(h=0.1,col='blue',lwd=2,lty=3)
abline(v=tol.lim,col='blue',lwd=2,lty=3)
dev.off()

