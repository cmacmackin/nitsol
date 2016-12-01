      subroutine bihsl( n, np2, lwpd, visc, c, d, g, wp, v, ier )

      implicit none

c----------------------------------------------------------------------
c This subroutine provides an interface to Bjorstad's dbihar.
c On output, z contains the solution. 
c
c Note: iflag is set to 3 below, resulting in the Cholesky option in 
c the fast biharmonic solver. 
c
c----------------------------------------------------------------------

      integer ier
      integer lwpd
      integer n
      integer np2

      double precision visc

      double precision c(*)
      double precision d(n,4)
      double precision g(np2,*)
      double precision v(n,*)
      double precision wp(lwpd)

c NOTE:  Need lwpd .ge 
c     11*n + 20*(n+3) + 19     if iflag = 2 (dbihar CG option)
c or
c     10*n + (n+1)**2/2 + 19   if iflag = 3. (dbihar Cholesky option) 

      integer i
      integer iflag
      integer itcg
      integer j
      integer jm1
      integer lwp
      integer np1

      double precision a0
      double precision alpha
      double precision b0
      double precision beta 
      double precision c0
      double precision d0
      double precision tol
      double precision visci

      save

      iflag = 3
      lwp = 10*n + ((n+1)**2)/2 + 19
      ier = 0
      if (lwp .gt. lwpd) then 
         write (6,*) 'Need larger lwpd in bihsl.'
         ier = -1
         return
      endif
      np1 = n + 1
      a0 = 0.0d0
      b0 = 1.0d0
      c0 = 0.0d0
      d0 = 1.0d0
      visci = 1.0d0/visc
      alpha = 0.0d0
      beta = 0.0d0
c
c     load grid array with interior points.
c
      do 20 j = 2, np1
         jm1 = j - 1
         do 10 i = 2, np1
            g(i,j) = -v(i-1,jm1)*visci
10       continue
20    continue
      tol = 1.d-6
      call dbihar ( a0, b0, n, d(1,4), d(1,2), d(1,1), d(1,1), c0, d0,
     *              n, g, np2, alpha, beta, iflag, tol, itcg, wp, lwp )
      if (iflag .lt. 0) then
         ier = iflag
         return
      endif
c
c if a high number of cg iterations occurs, reset iflag.
c
      if (itcg .ge. 15) then
         write(6,9000) itcg
         iflag = 2
      endif
c
c     load answers back into array v.
c
      do 40 j = 2, np1
         jm1 = j - 1
         do 30 i = 2, np1
            v(i-1,jm1) = g(i,j)
30       continue
40    continue

      return
c --------------------------------------------------------------------
c End of subroutine bihsl. 
c --------------------------------------------------------------------

9000    format(' **** itcg = ',i2,' iflag reset ****')

      end
