      subroutine init(neq, y, re, rpar, ipar)

      implicit none

c --------------------------------------------------------------------
c This subroutine loads the tail of y with the data for the 
c Navier-Stokes problem. This is an interface routine. 
c --------------------------------------------------------------------

      integer neq

      integer ipar(*)

      double precision re

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
      call init2(n, np2, y, re, rpar(lc), rpar(ld), rpar(lf), 
     $     rpar(lg), rpar(lvisc), rpar(lwk))
      return
c --------------------------------------------------------------------
c End of subroutine init. 
c --------------------------------------------------------------------
      end

