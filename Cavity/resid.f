      subroutine resid( n, np2, y, fy, visc, c, d, f, g, w )

      implicit none

c --------------------------------------------------------------------
c This subroutine is called by the interface subroutine fcav. It does 
c the work of evaluating the residual for the time independent 
c Navier-Stokes problem.
c --------------------------------------------------------------------

      integer n
      integer np2

      double precision visc

      double precision c(*)
      double precision d(n,4)
      double precision f(n,*)
      double precision fy(n,*)
      double precision g(np2,*)
      double precision w(np2,*)
      double precision y(n,*)

      integer i
      integer im1
      integer ip1
      integer j
      integer jm1
      integer jp1
      integer np1

      double precision gimjp, gijp, gipjp
      double precision gimj,  gij,  gipj
      double precision gimjm, gijm, gipjm
      double precision h
      double precision h2
      double precision t
      double precision wijp, wipjp
      double precision wimj,  wij,  wipj
      double precision wimjm, wijm
      double precision s, s1, s2, s3, s4, s5, s6

      double precision cvmgz
      external cvmgz

      save

      np1 = n + 1
      h = 1.d0/dfloat(np1)
      h2 = h*h
c
c     load grid array with interior points.
c
      do 20 j = 2, np1
         jm1 = j - 1
         do 10 i = 2, np1
            g(i,j) = y(i-1,jm1)
10       continue
20    continue
      call wload( n, np1, np2, h, g, d, c, w )
c
c     residual at interior points.
c
      do 40 j = 2, np1
         jp1 = j + 1
         jm1 = j - 1
         do 30 i = 2, np1
            ip1 = i + 1
            im1 = i - 1
            gij = g(i,j)
            gijp = g(i,jp1)
            gijm = g(i,jm1)
            gimj = g(im1,j)
            gimjp = g(im1,jp1)
            gimjm = g(im1,jm1)
            gipj = g(ip1,j)
            gipjp = g(ip1,jp1)
            gipjm = g(ip1,jm1)
c
c           linear term
c
            t = 20.d0*gij - 8.d0*( gimj + gipj + gijm + gijp )
     *          + 2.d0*( gimjp + gipjp + gimjm + gipjm )
            if ( i .eq. 2 ) then
               t = t + gij + 2.d0*h*d(jm1,4)
            else
               t = t + g(i-2,j)
            end if
            if ( i .eq. np1 ) then
               t = t + gij + 2.d0*h*d(jm1,2)
            else 
               t = t + g(i+2,j)
            end if
            if ( j .eq. 2 ) then
               t = t + gij + 2.d0*h*d(im1,1)
            else
               t = t + g(i,j-2)
            end if
            if ( j .eq. np1 ) then
               t = t + gij + 2.d0*h*d(im1,3)
            else 
               t = t + g(i,j+2)
            end if
            t = t*visc/h2
c
c           nonlinear term.
c
            wij = w(i,j)
            wipj = w(ip1,j)
            wimj = w(im1,j)
            wijp = w(i,jp1)
            wipjp = w(ip1,jp1)
            wijm = w(i,jm1)
            wimjm = w(im1,jm1)
            s1 = ( wipj - wij )*( gipjp - gipj )
     *           -  ( wipjp - wipj )*( gipj - gij )
            s2 = ( wipjp - wijp )*( gijp - gij )
     *           -  ( wijp - wij )*( gipjp - gijp )
            s3 = ( wij - wimj )*( gijp - gij )
     *           -  ( wijp - wij )*( gij - gimj )
            s4 = ( wij - wimj )*( gimj - gimjm )
     *           -  ( wimj - wimjm )*( gij - gimj )
            s5 = ( wijm - wimjm )*( gij - gijm )
     *           -  ( wij - wijm )*( gijm - gimjm )
            s6 = ( wipj - wij )*( gij - gijm )
     *           -  ( wij - wijm )*( gipj - gij )
            s = ( s1 + s2 + s3 + s4 + s5 + s6 )/(6.d0)
            fy(im1,jm1) = -t - s + h2*f(im1,jm1)
30       continue
40    continue
c
c     include influence of boundary conditions for the poisson
c     problem to be solved next.
c
      do 50 i = 1, n
         ip1 = i + 1
         fy(i,1) = fy(i,1) + g(ip1,1)
         fy(1,i) = fy(1,i) + g(1,ip1)
         fy(n,i) = fy(n,i) + g(np2,ip1)
         fy(i,n) = fy(i,n) + g(ip1,np2)
 50   continue
      return
c --------------------------------------------------------------------
c End of subroutine resid. 
c --------------------------------------------------------------------
      end

