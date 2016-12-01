      subroutine jacvcav(n, xcur, fcur, ijob, v, z, rpar, ipar, itrmjv)

      implicit none

c --------------------------------------------------------------------
c This subroutine evaluates Jacobian-vector products and preconditioner 
c solves for the time independent Navier-Stokes problem. This is an 
c interface routine. 
c --------------------------------------------------------------------

      integer ijob
      integer itrmjv
      integer n

      integer ipar(*)

      double precision fcur(n)
      double precision rpar(*)
      double precision v(n)
      double precision xcur(n)
      double precision z(n)

      integer ier

      if (ijob .eq. 0) then 
         write (6,*) ' Analytic Jacobian-vectpr product not implemented' 
         itrmjv = 2
         go to 900
      endif
      call dcopy(n, v, 1, z, 1)
      call psol(n, z, rpar, ipar, ier)
      if (ier .lt. 0) then 
         itrmjv = 2
         go to 900
      endif
      itrmjv = 0
c
c All returns made here.
c
 900  continue
      return
c --------------------------------------------------------------------
c End of subroutine jacvcav. 
c --------------------------------------------------------------------
      end
