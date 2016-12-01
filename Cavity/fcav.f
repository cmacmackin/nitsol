      subroutine fcav(neq, y, fvec, rpar, ipar, itrmf)

      implicit none

c --------------------------------------------------------------------
c This subroutine evaluates the nonlinear residual for the time independent
c Navier-Stokes problem. This is an interface routine. 
c --------------------------------------------------------------------

      integer itrmf
      integer neq

      integer ipar(*)

      double precision fvec(neq)
      double precision rpar(*)
      double precision y(neq)

      integer lc
      integer ld
      integer lf
      integer lg
      integer lvisc
      integer lwk
      integer n
      integer np2

      save

      n = ipar(1)
      np2 = ipar(2)
      lc = 1
      ld = lc + 4
      lf = ld + 4*n
      lg = lf + n*n
      lvisc = lg + np2*np2
      lwk = lvisc + 1
      call resid( n, np2, y, fvec, rpar(lvisc), rpar(lc), rpar(ld), 
     $     rpar(lf), rpar(lg), rpar(lwk))

      itrmf = 0

      return
c --------------------------------------------------------------------
c End of subroutine fcav. 
c --------------------------------------------------------------------
      end
