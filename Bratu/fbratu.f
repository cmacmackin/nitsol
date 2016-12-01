      subroutine fbratu(n, xcur, fcur, rpar, ipar, itrmf)
c --------------------------------------------------------------------
c This subroutine evaluates the nonlinear residual for the Bratu test 
c problem. 
c --------------------------------------------------------------------
      implicit none
      integer i, itrmf, j, j1, j2, n, ipar(*), nx, ny 
      double precision cl, cr, h2l, xcur(n), fcur(n), rpar(*)

c  Set local variables from user-supplied parameters.

      nx = ipar(1)
      ny = ipar(2)
      cl = rpar(1) 
      cr = rpar(2)
      h2l = rpar(3)

c  Evaluate nonlinear residual.

      do 100 j = 1, ny 
         j1 = (j - 1)*nx + 2 
         j2 = j*nx - 1
         do 110 i = j1, j2
            fcur(i) = cr*xcur(i+1) + cl*xcur(i-1) - 4.d0*xcur(i) 
     $           + h2l*dexp(xcur(i)) 
 110     continue
         j1 = j1 - 1
         fcur(j1) = cr*xcur(j1+1) - 4.d0*xcur(j1) + h2l*dexp(xcur(j1)) 
         j2 = j2 + 1 
         fcur(j2) = cl*xcur(j2-1) - 4.d0*xcur(j2) + h2l*dexp(xcur(j2)) 
         if (j .ne. 1) then 
            do 120 i = j1, j2 
               fcur(i) = fcur(i) + xcur(i-nx)
 120        continue
         endif
         if (j .ne. ny) then 
            do 130 i = j1, j2 
               fcur(i) = fcur(i) + xcur(i+nx)
 130        continue
         endif
 100  continue 

c  Set  termination flag for success.

      itrmf = 0

      return
c --------------------------------------------------------------------
c End of subroutine fbratu. 
c --------------------------------------------------------------------
      end
