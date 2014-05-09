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
#
#     -- Fine-resolution data
dat.1<-read.table('summary_bs10_3000iter_0.2inc.csv',header=TRUE,sep=',',
                  na.strings='-1')
#
#     -- Coarse-resoution data
dat.2<-read.table('summary_bs10_3000iter.csv',header=TRUE,sep=',',
                  na.string='-1')

#dat<-read.table('radial_bs10.csv',header=TRUE,sep=',',
#                na.strings='-1')
dat<-rbind(dat.1,dat.2)
#
#     -- Some notes on the data:
#
#        Init.Err proportional to the amount of the initial deformation
#        Res.Err  mesaured misalignment after registration
#        10 subjects, 16 initial misalignments, 10 experiments/misalignment
#      (Registration simulation with different seeds, 160 seeds in all.)
#
#
#     -- Categorical misalignment factor
dat$Misalign<-paste("M",round(dat$Init.Err),sep='.')
dat$Misalign<-factor(dat$Misalign,ordered=TRUE,
         levels=c(paste("M",0:15,sep='.')))
#
#     -- Categorical case-ID factor
dat$CaseID<-paste("C",dat$CaseID,sep='')
dat$CaseID<-as.factor(dat$CaseID)
cases<-levels(dat$CaseID)
#
#     -- Initial errors are very nearly integers, so we round to integers
dat$Init.Err<-round(dat$Init.Err)
#
#     -- Omit records with missing residual error, or zero initial error
dat<-dat[dat$Init.Err>0 & !is.na(dat$Res.Err),]
#
#     -- Recreate factor after deleting records
dat$Misalign<-as.factor(as.character(dat$Misalign))
#
#     -- Define "failure" in terms of a threshold using the parameter "h"
dat$Failure<-rep(FALSE,dim(dat)[1])
#
#      -- Initialize plot
pdf('Illustration1.1.pdf')
plot(0,0,type='n',xlim=c(1,15),ylim=c(0,1),lwd=2,
     xlab='misalignment magnitude, mm',
     ylab='probability of failure')

#
#      -- For each case, residual error > 0.5 is a failure   
for(case in cases){
    idx<-dat$CaseID==case
    dat$Failure[idx]<-dat$Res.Err[idx]>0.5 
}

#
#      -- y: Matrix of logistic regression fits, one column for each case
#         s2e: correspoding matrix of within-case squared stand. errors (SE)
y<-NULL; s2e<-NULL
for(case in cases){
    idx<-dat$CaseID==case
#
#      -- Simple logistic regression with Init.Err as covariate, for "case"
    fit<-glm(Failure~Init.Err,family='binomial',data=dat,
             subset=idx)
#
#      -- Evaluate fit at 1,1.25,...,15
    x<-seq(1,15,0.25); x<-data.frame(x); names(x)<-"Init.Err"
    pred<-predict(fit,newdata=x,type='response',se.fit=TRUE)
#
#      -- Retrieve fit and SE of fit from prediction object
    y<-cbind(y,pred$fit)
    s2e<-cbind(s2e,pred$se.fit^2)
#
#      -- Add fitted line to plot
    lines(x[,1],pred$fit,type='l',lwd=2)
}

#
#      -- Obtain between-case variance (s2b) for x_i, accross cases (row of y)
n.x<-dim(y)[1]; n.case<-dim(y)[2]; s2b<-NULL; y<-log(y/(1-y))
for(i in 1:n.x){
    z<-y[i,]
    var<-s2e[i,]
#
#       -- Wrapper for FORTRAN root-finder
    root<-mproot(z,var,n.case-1)
    s2b<-c(s2b,root$root)
}
#
#      -- Close plot
dev.off()
#
#      -- Reshape s2b into a n.x by n.case matrix by columns. Values are
#         repeated to fill the matrix, so each column is identical
s2b<-matrix(s2b,n.x,n.case)
#
#      -- Weight matrix: s2b same for each column (case), but s2e is not
w<-s2b/(s2b+s2e)
#
#      -- Limiting case of s2b=0
w[s2b==0]<-1/s2e[s2b==0]
#
#      -- Weight sum for each x, reshaped into matrix of equal columns
w.sum<-matrix(apply(w,1,'sum'),n.x,n.case)
#
#      -- Standardize weights to sum to 1
w<-w/w.sum
#
#      -- Weighted mean prob. of failure for each x, averaging over
#          cases, again reshaped into matrix of equal columns.
y.hat<-matrix(apply(y*w,1,'sum'),n.x,n.case)
#
#      -- Standard error of the mean, Rukhin and Vangel, JASA, 1999.
se.mean<-sqrt(apply((y-y.hat)^2/(s2b+s2e)^2,1,'sum'))/
              apply(1/(s2b+s2e),1,'sum')
#
#      -- Since columns of x, s2b, y.hat are the same, use only first column
x<-x[,1]; s2b<-s2b[,1]; y.hat<-y.hat[,1]
#
#      -- Standard deviation for prediction interval
sd<-sqrt(s2b+se.mean^2)
#
#      -- Upper and lower prediction limits on inverse logit of mean
#           prob. of failure, suing normal approximation.
k<-qnorm(0.95)
ucb<-ilogit(y.hat+k*sd); lcb<-ilogit(y.hat-k*sd); y.hat<-ilogit(y.hat)
#
#      -- Initialize second plot
pdf('Illustration1.2.pdf')
#
#      -- Fitted mean prob of failure
plot(x,y.hat,ylim=range(lcb,ucb),type='l',lwd=2,
     xlab='displacement magnitude, mm',ylab='probability of failure',
     main='Tolerance Limit Plot\nk=1.5, 50% Failures, (.90,.95) LTL = 7.7')
#
#      -- Upper and lower prediction limits onf mean prob of failure
lines(x,lcb,lty=2,lwd=2); lines(x,ucb,lty=2,lwd=2)
#
#      -- Tolerance limit is the interpolated value of the intersection
#           of "ucb" with horizontal line 0.1.
tol.lim<-approx(ucb,x,xout=0.1)$y
print(paste('Tolerance limit: ',tol.lim))
#
#     -- Horizontal line at 0.1
abline(h=0.1,lty=3,col='blue',lwd=2)
#
#      -- Vertical line at tolerance limit
lines(c(tol.lim,tol.lim),c(-0.05,1),col='blue',lty=3,lwd=2)
dev.off()
