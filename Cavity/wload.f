      subroutine wload( n, np1, np2, h, g, d, c, w )

      implicit none

c --------------------------------------------------------------------
c This subroutine is called by resid to calculate the 5-spot difference 
c scheme at the i,j grid point, including the boundary.
c --------------------------------------------------------------------

      integer n
      integer np1
      integer np2

      double precision h

      double precision c(*)
      double precision d(n,4)
      double precision g(np2,*)
      double precision w(np2,*)

      double precision two,       three,       four,       six 
      parameter      ( two=2.0d0, three=3.0d0, four=4.0d0, six=6.0d0 )

      integer i
      integer j

      double precision h2

      save
c
c     calculate w = g   + g  .
c                    xx    yy
c
      h2 = h*h
      i = 1
        j = 1
          w(i,j) = (three/two)*(two*g(i,j)-g(i+1,j)-g(i,j+1))/h2
     *             - (three/h)*c(1)
        do 10 j = 2, np1
          w(i,j) = (four*g(i,j)-two*g(i+1,j)-g(i,j-1)-g(i,j+1) )/h2
     *             - (two/h)*d(j-1,4)
 10       continue
        j = np2
          w(i,j) = three*(two*g(i,j)-g(i+1,j)-g(i,j-1))/h2
     *             - (three/h)*c(4)
      do 30 i = 2, np1
        do 20 j = 2, np1
          w(i,j) = (four*g(i,j)-g(i-1,j)-g(i+1,j)-g(i,j-1)-g(i,j+1))/h2
 20       continue
 30     continue
      i = np2
        j = 1
          w(i,j) = three*(two*g(i,j)-g(i-1,j)-g(i,j+1) )/h2
     *             - (six/h)*c(2)
        do 40 j = 2, np1
          w(i,j) = (four*g(i,j)-two*g(i-1,j)-g(i,j-1)-g(i,j+1))/h2
     *             - (two/h)*d(j-1,2)
 40       continue
        j = np2
          w(i,j) = (three/two)*(two*g(i,j)-g(i,j-1)-g(i-1,j))/h2
     *             - (three/h)*c(3)
      j = 1
        do 50 i = 2, np1
          w(i,j) = (four*g(i,j)-two*g(i,j+1)-g(i+1,j)-g(i-1,j))/h2
     *             - (two/h)*d(i-1,1)
 50       continue
      j = np2
        do 60 i = 2, np1
          w(i,j) = (four*g(i,j)-two*g(i,j-1)-g(i-1,j)-g(i+1,j))/h2
     *             - (two/h)*d(i-1,3)
 60       continue
      return
c --------------------------------------------------------------------
c End of subroutine wload. 
c --------------------------------------------------------------------
      end
