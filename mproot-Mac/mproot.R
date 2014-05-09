mproot<-function(x,t2,rhs){
nlab <-length(x)
maxit<-100
tol  <-1.e-7
root <-0
niter<-0
ier  <-0
z<-.Fortran("mproot",as.double(x),
                    as.double(t2),
                    as.integer(nlab),
                    as.double(tol),
                    as.double(rhs),
                    root=as.double(root),
                    as.integer(maxit),
                    niter=as.integer(niter),
                    ier=as.integer(ier))
mproot<-list(root=z$root, niter=z$niter, ier=z$ier)
return(mproot)}


