      subroutine psolb(n, x, nx, ny, rhs, wk, itrmjv)
c --------------------------------------------------------------------
c This subroutine is an interface between jacvbratu and the fast Poisson 
c solver hwscrt from Fishpack, which is used as a preconditioner for 
c the Bratu test problem. 
c --------------------------------------------------------------------
      implicit none
      integer i, idimf, ier, itrmjv, j, mx, my, n, nx, ny
      double precision x(n), rhs(nx+2,ny+2), wk(n)
      double precision a, b, c, d, dum(1), elmbda, h, pertrb 

c  Set parameters needed by hwscrt.

      mx = nx + 1
      my = ny + 1 
      ier = 0 
      h = 1.d0/dfloat(mx) 
      a = 0.d0 
      b = dfloat(mx)*h 
      c = 0.d0 
      d = dfloat(my)*h 

c  Initialize rhs for hwscrt.

      do  j = 1, my + 1
         rhs(1,j) = 0.d0
         rhs(mx+1,j) = 0.d0
      end do
      do i = 2, mx 
         rhs(i,1) = 0.d0
         rhs(i,my+1) = 0.d0
         do j = 2, my
            rhs(i,j) = x((j-2)*nx + i - 1)
         end do
      end do

c  Do the solve.

      idimf = nx + 2
      elmbda = 0.d0
      ier = 0
      call hwscrt (a, b, mx, 1, dum, dum, c, d, my, 1, dum, dum, 
     +             elmbda, rhs, idimf, pertrb, ier, wk) 
      if (ier .ne. 0) then 
         write (6,*) ' In psolb after hwscrt call, ier =', ier
         itrmjv = 2
         return
      endif

c  Copy output from hwscrt.

      do i = 2, mx 
         do j = 2, my 
            x((j-2)*nx + i - 1) = rhs(i,j)
         end do
      end do

      return
c --------------------------------------------------------------------
c End of subroutine psolb.  
c --------------------------------------------------------------------
      end
