      subroutine init2(n, np2, y, re, c, d, f, g, visc, wk)

      implicit none

c --------------------------------------------------------------------
c This subroutine is called by the interface subroutine init and does 
c the work of loading the tail of y with the data for the Navier-Stokes 
c problem. 
c --------------------------------------------------------------------

      integer n
      integer np2

      double precision re
      double precision visc

      double precision c(4)
      double precision d(n,4)
      double precision f(n,n)
      double precision g(np2,np2)
      double precision wk(np2,np2)
      double precision y(n,*)

      integer i
      integer j
      integer np1

      double precision zero,       half,        one
      parameter      ( zero=0.0d0, half=0.50d0, one=1.0d0 )

      save
c
c     fill viscosity
c
      visc = 1/re
c
c     fill normal derivatives at the corners.
c     (corners are ordered: sw,se,ne,nw)
c
      c(1) = zero
      c(2) = zero
      c(3) = half
      c(4) = half
c
c     fill normal derivatives on the sides.
c     (sides are ordered s,e,n,w)
c
      do 10 j = 1, n
         d(j,1) = zero
         d(j,2) = zero
         d(j,3) = one
         d(j,4) = zero
10    continue
c
c     fill forcing term.
c
      do 30 j = 1, n
         do 20 i = 1, n
            f(i,j) = zero
20       continue
30    continue
c
c     fill dirichlet conditions on auxiliary grid.
c
      np1 = n + 1
c
      do 40 i = 1, np2
         g(i,1) = zero
         g(i,np2) = zero
40    continue
c
      do 50 j = 2, np1
         g(1,j) = zero
         g(np2,j) = zero
50    continue
c
      do 70 j = 1, n
         do 60 i = 1, n
            y(i,j) = zero
 60      continue
 70   continue
      return
c --------------------------------------------------------------------
c End of subroutine init2. 
c --------------------------------------------------------------------
      end
