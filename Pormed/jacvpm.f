      subroutine jacvpm(n,xcur,fcur,ijob,v,z,rpar,ipar,itrmjv)
c --------------------------------------------------------------------
c This subroutine evaluates Jacobian-vector products and preconditioner 
c solves for the flow in porous media test problem. 
c --------------------------------------------------------------------
      implicit none
      integer n, ijob, ipar(*), itrmjv
      double precision xcur(n), fcur(n), v(n), z(n), rpar(*)
      double precision bll, bur, src, tol
      integer i, j, k, fill
      integer a, l, u, c, ja, il, jl, iu, ju, jnz, col, nnzlu

      double precision four
      parameter      ( four=4.0d0 )

c  For obtaining information about the state of the nonlinear iteration:
c
      include '../Nitsol/nitinfo.h'

c  Common blocks for sharing information with ILUT factorization routine.

      integer nx, ny
      common /grid/ nx, ny

c  Parameterized information for the incomplete factors.

      include 'ilut_pormed.h'

c  Storage needed to provide information on current
c  approximate solution to the incomplete factorization.

      include 'scratch.h'

c  Statement functions for evaluating Jacobian.

      double precision phipr, psipr, w

      phipr(w) = 2.d0*w
      psipr(w) = 3.d0*w*w

c  Set local variables from user-supplied parameters.

      nx = ipar(1)      
      ny = ipar(2)
      fill = ipar(3)

      bll = rpar(1) 
      bur = rpar(2)
      src = rpar(3)
      coef = rpar(4)
      tol = rpar(5)

c  Case:  evaluate J*v:

      if (ijob .eq. 0) then

c  This approach incurs a higher cost in the Jacobian-vector product
c  but saves storage.  Finite-difference Jacobian may be cheaper!
 
         k = 1

c  j=1:
c     i=1:

         z(k) = phipr(xcur(k+nx))*v(k+nx)
     &        + (phipr(xcur(k+1)) + coef*psipr(xcur(k+1)))*v(k+1)
     &        - four*phipr(xcur(k))*v(k)
         k = k + 1

         do i = 2, nx-1
            z(k) = phipr(xcur(k+nx))*v(k+nx)
     &           + (phipr(xcur(k+1)) + coef*psipr(xcur(k+1)))*v(k+1)
     &           - four*phipr(xcur(k))*v(k)
     &           + (phipr(xcur(k-1)) - coef*psipr(xcur(k-1)))*v(k-1)
            k = k + 1
         end do

c     i=nx:

         z(k) = phipr(xcur(k+nx))*v(k+nx)
     &        - four*phipr(xcur(k))*v(k)
     &        + (phipr(xcur(k-1)) - coef*psipr(xcur(k-1)))*v(k-1)
         k = k + 1

         do j = 2, ny-1
            z(k) = phipr(xcur(k+nx))*v(k+nx)
     &           + (phipr(xcur(k+1)) + coef*psipr(xcur(k+1)))*v(k+1)
     &           - four*phipr(xcur(k))*v(k)
     &           + phipr(xcur(k-nx))*v(k-nx)
            k = k + 1
         
            do i = 2, nx-1
               z(k) = phipr(xcur(k+nx))*v(k+nx)
     &              + (phipr(xcur(k+1)) + coef*psipr(xcur(k+1)))*v(k+1)
     &              - four*phipr(xcur(k))*v(k) 
     &              + (phipr(xcur(k-1)) - coef*psipr(xcur(k-1)))*v(k-1)
     &              + phipr(xcur(k-nx))*v(k-nx)
               k = k + 1
            end do
         
            z(k) = phipr(xcur(k+nx))*v(k+nx)
     &           - four*phipr(xcur(k))*v(k)
     &           + (phipr(xcur(k-1)) - coef*psipr(xcur(k-1)))*v(k-1)
     &           + phipr(xcur(k-nx))*v(k-nx)
            k = k + 1
         end do

c  j=ny:
c     i=1:

         z(k) = (phipr(xcur(k+1)) + coef*psipr(xcur(k+1)))*v(k+1)
     &        - four*phipr(xcur(k))*v(k)
     &        + phipr(xcur(k-nx))*v(k-nx)
         k = k + 1
         
         do i = 2, nx-1
            z(k) = (phipr(xcur(k+1)) + coef*psipr(xcur(k+1)))*v(k+1)
     &           - four*phipr(xcur(k))*v(k)
     &           + (phipr(xcur(k-1)) - coef*psipr(xcur(k-1)))*v(k-1)
     &           + phipr(xcur(k-nx))*v(k-nx)
            k = k + 1
         end do
         
c     i=nx:
         
         z(k) = - four*phipr(xcur(k))*v(k)
     &        + (phipr(xcur(k-1)) - coef*psipr(xcur(k-1)))*v(k-1)
     &        + phipr(xcur(k-nx))*v(k-nx)
         k = k + 1
      endif

c  Case:  evaluate P(inverse)v:

      if (ijob .eq. 1) then 

c  Partition ipar, rpar to access ILUT information.  Parameters
c  are all set in the include file ilut_pormed.h.

         a = 6
         l = a + STENCIL_SZ
         u = l + MAX_FACTR
         c = u + MAX_FACTR
             
         ja = 4
         il = ja + STENCIL_SZ
         jl = il + MAX_SZP1
         iu = jl + MAX_FACTR
         ju = iu + MAX_SZP1
         jnz = ju + MAX_FACTR
         col = jnz + MAX_SZP1

c  Check to see if P needs to be updated.

         if ( newstep .eq. 0 ) then

c  Workspace for the factorization is obtained from ipar, rpar.
c  It includes storage for one row of the Jacobian and compressed
c  sparse row representations of the incomplete factors.  Note that
c  no Jacobian is actually stored.  This version of ILUT only requires
c  one row of the coefficient matrix at a time, which is done in an
c  auxiliary routine called getrow.  The following copy is needed
c  to provide getrow with information about the current approximate
c  solution.

            do i = 1, n
               myxcur(i) = xcur(i)
            end do
            nnzlu = MAX_FACTR
            call ilutf( rpar(a), n, STENCIL_SZ, ipar(ja), tol, fill,
     &                  rpar(l), ipar(il), ipar(jl), rpar(u), ipar(iu),
     &                  ipar(ju), ipar(jnz), ipar(col), rpar(c), nnzlu )

c  Reset newstep so that P is only updated at
c  the start of a linear solution cycle.

            newstep = 1

         endif

c  Now, do the solve.
         
         call m1i( v, n, rpar(l), ipar(il), ipar(jl), MAX_FACTR, z )
         
         call m2i( z, n, rpar(u), ipar(iu), ipar(ju), MAX_FACTR )
      
      endif

c  Of course, everything always completes successfully.  (This optimism
c  is more of a reflection of the lack of error reporting capabilities
c  in the current versions of the incomplete factorization routines.)

      itrmjv = 0

      return
c --------------------------------------------------------------------
c End of subroutine jacvpm 
c --------------------------------------------------------------------
      end
