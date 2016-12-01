      subroutine psol (neq, v, rpar, ipar, ier)

      implicit none

c-----------------------------------------------------------------------
c This subroutine is called by the interface subroutine jacvcav. It calls 
c another interface routine which in turn calls the biharmonic solver dbihar 
c to compute v <-- (nu*a)-inverse*v. 
c-----------------------------------------------------------------------

      integer ier
      integer neq

      integer ipar(*)

      double precision rpar(*)
      double precision v(neq)

      integer lc
      integer ld
      integer lf
      integer lg
      integer lvisc
      integer lwk
      integer lwpd
      integer n
      integer np1
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
      np1 = n + 1
      lwpd = max(11*n + 20*(n+3),10*n + ((n+1)**2)/2) + 20
      call bihsl ( n, np2, lwpd, rpar(lvisc), rpar(lc), rpar(ld), 
     $     rpar(lg), rpar(lwk), v, ier )
      call dscal ( neq, dfloat(np1**2), v, 1 )
      return
c --------------------------------------------------------------------
c End of subroutine psol. 
c --------------------------------------------------------------------
      end
