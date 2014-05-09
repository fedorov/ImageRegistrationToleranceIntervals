#
#     Local likelihood regression library 
library(locfit)
#
#    -- Object file corresponding to Fotran root-finding subroutine
dyn.load('../mproot-Linux/mproot.so')
#
#    -- R wrapper which calls external root-finder
source('../mproot-Linux/mproot.R')
#
#    -- Load simulation results: one for each of 10 samples
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
#
#     -- Before combining these 10 data frames into one, we create a
#         coulumn for sample number in each data frame.
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
#
#      -- Delete a column which appears in data frames for only
#            some of the samples.
nm<-names(dat1); dat1<-dat1[,nm!="Error.4"]
nm<-names(dat4); dat4<-dat4[,nm!="Error.4"]
nm<-names(dat5); dat5<-dat5[,nm!="Error.4"]
nm<-names(dat6); dat6<-dat6[,nm!="Error.4"]
#
#      -- Join all 10 data frames into one, by rows.
dat<-rbind(dat1,dat2,dat3,dat4,dat5,dat6,dat7,dat8,dat9,dat10)
sample.ids<-c(1,2,3,4,5,6,7,8,9,10)  # all samples
# uncomment to get the plot excluding the outlier sample
dat<-rbind(dat1,dat2,dat3,dat4,dat5,dat6,dat7,dat8,dat10)
sample.ids<-c(1,2,3,4,5,6,7,8,10)    # omit outlier (sample 9)
#
#     -- Failures are defined to be either missing values (nonconvergence),
#         or else errors > t, where threshold t was chosen by examining
#         within-sample errors.
t<-3
dat$Error.1<- is.na(dat$Error.1) | ((!is.na(dat$Error.1) & (dat$Error.1 >t)))
dat$Error.2<- is.na(dat$Error.2) | ((!is.na(dat$Error.2) & (dat$Error.1 >t)))
dat$Error.3<- is.na(dat$Error.3) | ((!is.na(dat$Error.3) & (dat$Error.1 >t)))
#
#      -- Omit cases with rotations > 7 
max.rot<-7
dat<-subset(dat, abs(dat$Rot.X)<max.rot & abs(dat$Rot.Y)<max.rot & abs(dat$Rot.Z)<max.rot)
#
#      -- Define failure, and success as the conjugate
dat$Failure<-dat$Error.1 | dat$Error.2 | dat$Error.3 #| dat$Error.4
dat$Success<-!dat$Failure
#
#      -- Create rotation and direction factors.  For the rotation
#            factor, the levels are denoted "x.y.z", where x,y,z
#            are rotation angles around x,y, and z, respectively.
dat$Rotation<-paste(dat$Rot.X,dat$Rot.Y,dat$Rot.Z,sep='.')
dat$Direction<-as.factor(dat$Direction)
dat$Rotation<-as.factor(dat$Rotation)
#
#     -- Names for levels of rotation factor; in the data file
#          they are just set to 1,...,14.
levels(dat$Direction)<-c(
    "+x","-x", "+y","-y", "+z","-z",  
    "+x+y+z", "+x+y-z", "+x-y+z", "-x+y+z",
     "-x-y-z", "+x-y-z","-x+y-z","-x-y+z")
#
#      -- "n.samp" combinations of computer experiment factors will
#             be chosen at random for each combination of sample (values
#             in "sample.ids") with initial displacement (0,2,4,6,8 or
#             10.
n.samp<-500; n.x<-100; n.sample<-length(sample.ids)
disp<-c(0,2,4,6,8,10)
tbl<-NULL
for(sample in sample.ids){
  for(d in disp){
#
#      -- The "sample" function below selcts "n.samp" elements
#            with replacement from the first argument.
#       
#        NOTE: Each time this function is called it will use
#           a randomly changed seed, so the plots will be
#           slightly different each time this script is run.
#           To avoid this behavior, call the function "set.seed"
#           before the first call to "sample".
#       
       y<-sample(dat$Success[dat$Displacement==d &
                             dat$Sample==sample],
                 n.samp,replace=TRUE)
       tbl<-rbind(tbl,cbind(sample,d,y))
   }}
#
#      -- By random sampling from the large simulation experiment,
#           we create a new dataset with 500 values for each combination
#           of sample and initial displacement.  We then overwrite the
#           original data frame "dat".
tbl<-data.frame(tbl)
names(tbl)<-c("Sample","Misalign","Success")
dat<-tbl
#
#
#      -- p: Matrix of logistic regression fits, one column for each sample
#         s2e: correspoding matrix of within-sample squared stand. errors (SE).
#         s2b: matrix of between-sample variances (of rank 1: all columns the same)
p<-NULL; s2e<-NULL; s2b<-NULL
#
#     -- Definition of logit and inverse-logit functions
logit<-function(p){return(log(p/(1-p)))}
ilogit<-function(x){return(exp(x)/(exp(x)+1))}
#
#    -- Initialize plot
pdf(paste('Illustration3-',n.samp,'samp-',length(sample.ids),'cases.pdf',sep=''))
#
#     -- Evaluate fitted curves along "n.x" equally spaced
#          initial displacements.
x0<-seq(0,max(dat$Misalign),length=n.x)
#
#     -- Get prob. of failure curves for each sample.
for(i in 1:n.sample){
#
#     -- Extract data for sample i
    idx<-dat$Sample==sample.ids[i]
    y<-dat$Success[idx];  x<-dat$Misalign[idx]
#
#   -- Sort x, carry along y, negate y to get failure from success.
    o<-order(x);  x<-x[o];  y<-y[o];  y<-!y
#
#   -- Local likelihood weighted logistic regression
    fit<-locfit(y~lp(x),family='binomial')
#    
#   -- Predicted probability of failure for each of the misalignment values
    yhat<-predict(fit,x0)
#
#    -- Matrix "p" is logit prob of failure; one column for each sample
    p<-cbind(p,logit(yhat))
#
#    -- Within-sample SE^2 for logit prob of failure    
    s2e<-cbind(s2e,predict(fit,x0,se.fit=TRUE)$se.fit^2)
}
#
#    -- Calculate estimate of between-samplevariance for
#          each misalighment value using Fortran routine
#          called by R wrapper function "mproot".
for(i in 1:n.x){
  s2b<-c(s2b,mproot(p[i,], s2e[i,], n.sample-1)$root)
}
#
#      -- Weight matrix. Each row corresponds to a different
#           misalignment norm; each column corresponds to a
#           difference sample.
w<-NULL
for(i in 1:n.x){
  if(s2b[i]!=0){
     w<-rbind(w,s2b[i]/(s2b[i]+s2e[i,]))}
  if(s2b[i]==0){
     w<-rbind(w,1/s2e[i,])}
}
#
#      -- Weight sum for each x, reshaped into matrix of equal columns
w.sum<-matrix(apply(w,1,'sum'),n.x,n.sample)
#
#      -- Standardize weights to sum to 1
w.norm<-w/w.sum
#
#      -- Weighted mean prob. of failure for each x, averaging over
#          cases, reshaped into matrix of equal columns. The
#          vector s2b is also reshaped into a matrix in the
#          same way.
p.hat<-apply(w.norm*p,1,'sum')
p.hat<-matrix(p.hat,n.x,n.sample); s2b<-matrix(s2b,n.x,n.sample)
#
#      -- Standard error of the mean, Rukhin and Vangel, JASA, 1999.
se.mean<-sqrt(apply((p-p.hat)^2/(s2b+s2e)^2,1,'sum'))/
              apply(1/(s2b+s2e),1,'sum')
#
#      -- Since columns of x, s2b, p.hat are the same, use only first column
s2b<-s2b[,1]; p.hat<-p.hat[,1]
#
#      -- Standard deviation for prediction interval
sd<-sqrt(s2b+se.mean^2)
#
#      -- Upper and lower prediciton limits, calculated on
#           logit scale and then tranfrmed into probabilities.
k<-qnorm(0.95)
ucb<-ilogit(p.hat+k*sd)
lcb<-ilogit(p.hat-k*sd)
#
#      -- Tolerance limit is the interpolated value of the intersection
#           of "ucb" with horizontal line 0.1.  We find the root in
#           the interval x0>1 in order to attempt to find the
#           largest intesection x value, in the event that there are
#           more than one such values.
#
tol.lim<-approx(ucb[x0>1],x0[x0>1],xout=0.1)$y
#
#      -- The mean probability curve was determined on the logit scale;
#            here we transform back.
p.hat<-ilogit(p.hat)
#
#      -- Plot mean prbability of failure curve
plot(x0,p.hat,type='l',lwd=2,
     ylim=range(ucb,lcb),
     xlab='misalignment norm, mm',
     ylab='probability of failed registration',
     main=paste('Average Probability of Failed Registration with\n95% Upper Prediction Limit,\nand (0.90, 0.95) Lower Tolerance Limit = ',round(tol.lim,1)))
#
#      -- Upper and lower prediction limits
lines(x0,ucb,lty=2,lwd=2)
lines(x0,lcb,lty=2,lwd=2)
#
#      -- Horizontal line at 0.1; vertical line at tolerance limit
abline(h=0.1,col='blue',lwd=2)
abline(v=tol.lim,col='blue',lwd=2)
dev.off()
      



    

    
