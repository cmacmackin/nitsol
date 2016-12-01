      subroutine fpm(n, xcur, fcur, rpar, ipar, itrmf)
c --------------------------------------------------------------------
c This subroutine evaluates the nonlinear residual for the flow in
c porous media test problem. 
c --------------------------------------------------------------------
      implicit none
      integer i, itrmf, j, k, n, ipar(*), nx, ny 
      double precision bll, bur, src, coef, xcur(n), fcur(n), rpar(*) 

      double precision four
      parameter      ( four=4.0d0 )

      double precision phi, psi, v

c  Statement functions for calculating the nonlinear residual.

      phi(v) = v*v
      psi(v) = v*v*v

c  Set local variables from user-supplied parameters.

      bll = rpar(1)
      bur = rpar(2)
      src = rpar(3)
      coef = rpar(4)

      nx = ipar(1)
      ny = ipar(2)

c  Evaluate nonlinear residual.

      k = 1

c  j=1:
c     i=1:

      fcur(k) = phi(xcur(k+nx))
     &              + phi(xcur(k+1)) + coef*psi(xcur(k+1))
     &              - four*phi(xcur(k)) 
     &              + phi(bll) - coef*psi(bll) 
     &              + phi(bll)
      k = k + 1

      do i = 2, nx-1
         fcur(k) = phi(xcur(k+nx))
     &           + phi(xcur(k+1)) + coef*psi(xcur(k+1))
     &           - four*phi(xcur(k)) 
     &           + phi(xcur(k-1)) - coef*psi(xcur(k-1))
     &           + phi(bll)
         k = k + 1
      end do

c     i=nx:

      fcur(k) = phi(xcur(k+nx))
     &        + phi(bur) + coef*psi(bur)
     &        - four*phi(xcur(k)) 
     &        + phi(xcur(k-1)) - coef*psi(xcur(k-1))
     &        + phi(bll)
      k = k + 1

      do j = 2, ny-1
         fcur(k) = phi(xcur(k+nx))
     &           + phi(xcur(k+1)) + coef*psi(xcur(k+1))
     &           - four*phi(xcur(k)) 
     &           + phi(bll) - coef*psi(bll)
     &           + phi(xcur(k-nx))
         k = k + 1

         do i = 2, nx-1
            fcur(k) = phi(xcur(k+nx))
     &              + phi(xcur(k+1)) + coef*psi(xcur(k+1))
     &              - four*phi(xcur(k)) 
     &              + phi(xcur(k-1)) - coef*psi(xcur(k-1))
     &              + phi(xcur(k-nx))
            k = k + 1
         end do

         fcur(k) = phi(xcur(k+nx))
     &           + phi(bur) + coef*psi(bur)
     &           - four*phi(xcur(k)) 
     &           + phi(xcur(k-1)) - coef*psi(xcur(k-1))
     &           + phi(xcur(k-nx))
         k = k + 1
      end do

c  j=ny:
c     i=1:

      fcur(k) = phi(bur)
     &        + phi(xcur(k+1)) + coef*psi(xcur(k+1))
     &        - four*phi(xcur(k)) 
     &        + phi(bll) - coef*psi(bll)
     &        + phi(xcur(k-nx))
      k = k + 1
      
      do i = 2, nx-1
         fcur(k) = phi(bur)
     &           + phi(xcur(k+1)) + coef*psi(xcur(k+1))
     &           - four*phi(xcur(k)) 
     &           + phi(xcur(k-1)) - coef*psi(xcur(k-1))
     &           + phi(xcur(k-nx))
         k = k + 1
      end do

c     i=nx:

      fcur(k) = phi(bur)
     &        + phi(bur) + coef*psi(bur)
     &        - four*phi(xcur(k)) 
     &        + phi(xcur(k-1)) - coef*psi(xcur(k-1))
     &        + phi(xcur(k-nx))
      k = k + 1

c  Add source term.

      fcur(1) = fcur(1) + src

c  Set  termination flag for success.

      itrmf = 0 

      return
c --------------------------------------------------------------------
c End of subroutine fpm 
c --------------------------------------------------------------------
      end
