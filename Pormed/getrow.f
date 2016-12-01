      subroutine getrow( i, max_nz, a, ja, nnz )

      implicit none

c  Purpose-  Provide a row by row description of the Jacobian of the
c            flow in porous medium model problem, in CSR format.

c  Arguments:
c    Input:   i      - integer
c                      row for which elements are desired.
c             max_nz - integer
c                      maximum number of nonzeros expected in row.
c    Output:  a      - double precision array of length max_nz.
c                      nonzero values in row.
c             ja     - integer array of length max_nz
c                      corresponding column indices.
c             nnz    - integer
c                      number of nonzeros provided in a and ja.

c  Dummy arguments-

      integer i
      integer max_nz
      integer nnz

      integer ja(max_nz)

      double precision a(max_nz)

c  Common blocks.
c  (See documentation of jacvpm for an explanation.)

      integer nx, ny
      common /grid/ nx, ny

      include 'ilut_pormed.h'

      include 'scratch.h'

      double precision four
      parameter      ( four=4.0d0 )

c  Local variables-

      integer grid_i
      integer grid_j

c  Statement functions-

      double precision phipr, psipr, u

      phipr(u) = 2.d0*u
      psipr(u) = 3.d0*u*u

c  Start of executable code-

c  Given a row index i, find the corresponding coordinates on an nx by ny 
c  grid, and store the 5 point stencil centered at that grid point.

      grid_j = (i-1)/nx + 1
      grid_i = i - ((i-1)/nx)*nx 

      nnz = 0

c  South coupling.

      if ( grid_j .gt. 1 ) then
         nnz = nnz + 1
         a(nnz) = phipr(myxcur(i-nx))
         ja(nnz) = i - nx
      endif

c  West coupling.

      if ( grid_i .gt. 1 ) then
         nnz = nnz + 1
         a(nnz) = phipr(myxcur(i-1)) - coef*psipr(myxcur(i-1))
         ja(nnz) = i - 1
      endif

c  Center of stencil.

      nnz = nnz + 1
      a(nnz) = -four*phipr(myxcur(i))
      ja(nnz) = i

c  East coupling.

      if ( grid_i .lt. nx ) then
         nnz = nnz + 1
         a(nnz) = phipr(myxcur(i+1)) + coef*psipr(myxcur(i+1))
         ja(nnz) = i + 1
      endif

c  North coupling.

      if ( grid_j .lt. ny ) then
         nnz = nnz + 1
         a(nnz) = phipr(myxcur(i+nx))
         ja(nnz) = i + nx
      endif

      return

c  End of getrow.

      end
