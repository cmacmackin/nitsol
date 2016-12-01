      subroutine jacvbratu(n,xcur,fcur,ijob,v,z,rpar,ipar,itrmjv)
c --------------------------------------------------------------------
c This subroutine evaluates Jacobian-vector products and preconditioner 
c solves for the Bratu test problem. 
c --------------------------------------------------------------------
      implicit none
      integer n, ijob, ipar(*), itrmjv
      double precision xcur(n), fcur(n), v(n), z(n), rpar(*)
      double precision cl, cr, dexp, h2l
      integer i, j, j1, j2, nx, ny 

c  Set local variables from user-supplied parameters.

      nx = ipar(1)
      ny = ipar(2)
      cl = rpar(1) 
      cr = rpar(2)
      h2l = rpar(3)

c  Perform requested operations.

      if (ijob .eq. 0) then

c  Analytic Jacobian-vector product.

         do 100 j = 1, ny 
            j1 = (j - 1)*nx + 2 
            j2 = j*nx - 1
            do 110 i = j1, j2
               z(i) = cr*v(i+1) + cl*v(i-1) - 
     $              (4.d0 - h2l*dexp(xcur(i)))*v(i)  
 110        continue
            j1 = j1 - 1
            z(j1) = cr*v(j1+1) - (4.d0 - h2l*dexp(xcur(j1)))*v(j1)  
            j2 = j2 + 1 
            z(j2) = cl*v(j2-1) - (4.d0 - h2l*dexp(xcur(j2)))*v(j2) 
            if (j .ne. 1) then 
               do 120 i = j1, j2 
                  z(i) = z(i) + v(i-nx)
 120           continue
            endif
            if (j .ne. ny) then 
               do 130 i = j1, j2 
                  z(i) = z(i) + v(i+nx)
 130           continue
            endif
 100     continue 

c  Set termination flag.

         itrmjv = 0

      endif

      if (ijob .eq. 1) then 

c  Precondition.  First copy input into output.

         do 200 i = 1, n
            z(i) = v(i)
 200     continue

c  Call interface routine that interfaces with the
c  fast Poisson preconditioner.  Note use of rpar
c  to provide storage for the fast Poisson solver.

         call psolb(n, z, nx, ny, rpar(4), rpar((nx+2)*(ny+2)+4))

c  Set termination flag.

         itrmjv = 0

      endif

      return
c --------------------------------------------------------------------
c End of subroutine jacvbratu. 
c --------------------------------------------------------------------
      end
