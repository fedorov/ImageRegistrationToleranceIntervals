      subroutine mproot
     $ (x, t, nlab, tol, rhs, root, maxit, niter, ier)
      implicit double precision (a-h, o-z)
      dimension x(nlab), t(nlab)
c
c     -- If root is negative or nonexistent, set it to 
c        zero.
      y  = 0
      f0 = f(x, t, nlab, y)
      if (f0 .lt. rhs) then
         root = 0
         ier  = 1
         return
      end if
c
c    -- Loop until convergence
      niter = 0
 30   continue
         niter = niter +1
         fd = fderiv (x, t, nlab, y, wsum)
         ynew = y +(f(x, t, nlab, y)-rhs)/fd
         if (abs (ynew-y) .le. tol .and.
     $            niter   .le. maxit) then
            root = ynew
            ier  = 0
            return
         else if (niter .gt. maxit) then
            root = ynew
            write (*,*) ' Failure in mproot : ',ynew
            write (*,*) ' ts : ',(t(k),k=1,nlab)
            ier  = 2
            return
         else
           y = ynew
           go to 30
         end if
      end
      double precision function f(x, t, nlab, y)
      implicit double precision (a-h, o-z)
      dimension x(nlab), t(nlab)
c
c     -- Calculate the mean estimate
      wsum = 0
      xhat = 0
      do 10 i=1, nlab
         wght = 1/(t(i)+y)
         wsum = wsum +wght
         xhat = xhat +x(i)*wght
 10   continue
      xhat = xhat /wsum
c
c    -- Evaluate the function
      f = 0
      do 20 i=1, nlab
         f = f +1/(t(i) +y) *(x(i) -xhat)**2
 20   continue
      return
      end
      double precision function fderiv(x, t, nlab, y, wsum)
      implicit double precision (a-h, o-z)
      dimension x(nlab), t(nlab)
c
c     -- Calculate the mean estimate
      wsum = 0
      xhat = 0
      do 10 i=1, nlab
         wght = 1/(t(i)+y)
         wsum = wsum +wght
         xhat = xhat +x(i)*wght
 10   continue
      xhat = xhat /wsum
c
c    -- Evaluate the function
      fderiv = 0
      do 20 i=1, nlab
         fderiv = fderiv +1/(t(i) +y)**2 *(x(i) -xhat)**2
 20   continue
      return
      end
