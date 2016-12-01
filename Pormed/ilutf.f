      subroutine ilutf( a, n, nnzrow, ja, tol, fill, l, il, jl,
     :                                   u, iu, ju, jnz, ind, cr, nnzf )

c  Purpose-  Incomplete LU factorization of a sparse matrix A
c            with fill and drop tolerance.  This is a variant
c            that only reads one row of A at a time.

c  Author-
c    Michael Pernice
c    Utah Supercomputing Institute
c    May, 1994

c  Revisions-
c    February 1995
c      -essentially a complete rewrite to mimic Saad's strategy.
c       Added features are reduced workspace requirements and
c       elimination of the need to provide all of A in CSR format.

c  Arguments-
c    Input:   a    - double precision array of length nnzrow
c                    Array for 1 row of nonzero elements of coefficient
c                     matrix A.
c             n    - integer
c                    Order of linear system.
c             nnzrow - integer
c                    Maximum number of nonzero elements in each row of A.
c             ja   - integer array of length nnzrow
c                    Array for column indices of elements of a.
c             tol  - double precision
c                    Drop tolerance.  Multipliers L_ik are dropped
c                    when | L_ik | < tol*|| A(i,:) ||.  Elements from the
c                    reduced rows of A are dropped by the same criterion.
c             fill - integer
c                    Fill level.  Each row of the triangular factor has
c                    at most (fill + number of nonzeros in A(i,1:i))
c                    nonzero elements.  (Similarly for the elements in
c                    the reduced rows of A.)
c             cr   - double precision array of length n
c                    Workspace for dense copies of sparse rows of A.
c             nnzf - integer
c                    Estimate of number of nonzero elements in the
c                    sparse incomplete triangular factor.  If A has
c                    at most k nonzero elements per row, then
c                    nnza <= k*n, and nnzf > (k+fill)*n should be
c                    sufficient.  Overwritten on output.

c    Output:   l   - double precision array of length nnzf
c                    Nonzero elements of sparse lower triangular factor,
c                    in compressed sparse row format.
c             il   - integer array of length n+1
c                    Array of row pointers in the sparse triangular
c                    factor.
c             jl   - integer array of length nnzf
c                    Array of column indices in the sparse triangular
c                    factor.
c             u    - double precision array of length nnzf
c                    Nonzero elements of sparse upper triangular factor,
c                    in compressed sparse row format.
c             iu   - integer array of length n+1
c                    Array of row pointers into u:  u(iu(i)) is the
c                    first element of row i in u.
c             ju   - double precision array of length nnzf
c                    Array of column indices of retained elements of
c                    the reduced rows of A.
c             nnzf - integer
c                    On output, actual number of nonzeros in both factors.

c  To use this routine, you must provide a routine to provide the nonzero
c  values of each row in A.  The syntax and function must be
c
c               call getrow( i, max_nz, a, ja, nnz )
c
c  where i is the row index and max_nz is the maximum number of nonzeros
c  expected in row i of A.  On return, a contains the nonzero elements of
c  row i of A, and ja the corresponding column indices, and nnz nonzero
c  values were actually provided.  This eliminates the need to store the
c  A in CSR format to begin with.

c  Todo:
c  Mostly flag error conditions.  Things to flag:
c     - lenl and/or lenu may grow too large during either the
c       unpacking or reduction phase (when fill-in occurs)
c     - a zero or small pivot may be encountered
c       A fixup is used, but a better one is probably available;
c       see comments below.
c     - inadequate space for the lower or upper triangle factor
c       may be provided
c     -exit if rownorm = 0

c  Dummy arguments-

      integer fill
      integer n
      integer nnzrow
      integer nnzf

      integer il(n+1)
      integer ind(n)
      integer iu(n+1)
      integer ja(nnzrow)
      integer jl(nnzf)
      integer jnz(n)
      integer ju(nnzf)

      double precision tol

      double precision a(nnzrow)
      double precision l(nnzf)
      double precision u(nnzf)
      double precision cr(n)

c  Parameters-

      double precision zero,       one
      parameter      ( zero=0.0d0, one=1.0d0 )

c  Local variables-

      integer i
      integer j
      integer jj
      integer jak
      integer juj
      integer k
      integer keepl
      integer keepu
      integer lenl
      integer lenu
      integer lptr
      integer nl
      integer nnz
      integer pos
      integer row
      integer uptr

      double precision fact
      double precision rownorm
      double precision rowtol
      double precision tmp

c  Externals-

      external qsplit
      external dsortja

c  Intrinsics-

      intrinsic abs
      intrinsic dble
      intrinsic min
      intrinsic sqrt

c  Start of executable code-

c  Initialize.

      do 10 i = 1, n
         cr(i) = zero
         jnz(i) = 0
         ind(i) = 0
 10   continue

c  Initialize pointers to next open location
c  in triangular factor and reduced rows of A.

      lptr = 1
      uptr = 1

c  Reduce matrix by rows (ikj nesting of symmetric Gaussian elimination).

      do 110 i = 1, n

c  Store pointers to next open locations in il and iu.

         il(i) = lptr
         iu(i) = uptr

c  Make a dense copy of row i of A and calculate the average magnitude of
c  nonzero values.  The nonzero elements are mapped to contiguous locations
c  in cr.  Arrays jnz and ind give a likewise compressed picture of the 
c  sparsity pattern of the row, and are updated as fill-in occurs.

         rownorm = zero
         lenl = 0
         lenu = 1
         jnz(i) = i
         call getrow( i, nnzrow, a, ja, nnz )
         do 20 k = 1, nnz
            tmp = a(k)
            jak = ja(k)
            if ( jak .lt. i ) then
               lenl = lenl + 1
               jnz(jak) = lenl
               ind(lenl) = jak
               cr(lenl) = tmp
            else if ( jak .eq. i ) then
               cr(i) = tmp
            else if ( jak .gt. i ) then
               lenu = lenu + 1
               jnz(jak) =  i + lenu - 1
               ind(i+lenu-1) = jak
               cr(i+lenu-1) = tmp
            endif
            rownorm = rownorm + abs(tmp)
 20      continue
         rownorm = rownorm/dble(lenu+lenl)

c  Set the tolerance for dropping elements and
c  the number of elements to keep in each factor.

         keepl = lenl+fill
         keepu = lenu+fill-1
         rowtol = rownorm*tol

c  Main reduction loop:  calculate row i of L and U.  Use only
c  those previous rows that correspond to nonzero elements in
c  the first part of the row.  This loop assumes that ind(k) is
c  in increasing order.  nl counts the number of elements going
c  into this row of the lower triangular factor and k would be
c  k = 1, ..., i-1 in standard dense ikj LU factorization.

          k = 0
          nl = 0
 30       k = k + 1

c  When there are no more elements left to eliminate jump out of loop.

          if ( k .gt. lenl ) goto 60

c  Find multiplier and column position it operates on.

          row = ind(k)

c  Must be sure that row is smallest value of ind(k:lenl).

          jj = k
          do 40 j = k+1, lenl
             if ( ind(j) .lt. row ) then
                jj = j
                row = ind(j)
             endif
 40       continue

c  If it wasn't, swap elements in description of active row.

          if ( jj .ne. k ) then

c     Elements in ind and cr...

             j = ind(jj)
             ind(jj) = ind(k)
             ind(k) = j

             tmp = cr(jj)
             cr(jj) = cr(k)
             cr(k) = tmp

c     ... and elements in jnz.

             jnz(row) = j
             jnz(ind(jj)) = jj

          endif

          fact = cr(k)*u(iu(row))
          jnz(row) = 0
          if ( abs(fact) .lt. rowtol ) goto 60

c  Reduce current row, operating on nonzeros only.

          do 50 j = iu(row)+1, iu(row+1)-1
             tmp = fact*u(j)
             juj = ju(j)
             pos = jnz(juj)
             if ( juj .ge. i ) then

c  Put result in upper triangular part of work row.

                if ( pos .eq. 0 ) then 

c  Fill element.  (Flag lenu > n-i+1 as error, and exit.)

                   lenu = lenu + 1
                   ind(i+lenu-1) = juj
                   jnz(juj) = i + lenu - 1
                   cr(i+lenu-1) = -tmp
                else
                   cr(pos) = cr(pos) - tmp
                endif
             else 

c  Put result in lower triangular part of work row.

                if ( pos .eq. 0 ) then

c  Fill element.  (Flag lenl > i-1 as error, and exit.)

                   lenl = lenl + 1
                   ind(lenl) = juj
                   jnz(juj) = lenl
                   cr(lenl) = -tmp
                else
                   cr(pos) = cr(pos) - tmp
                endif
              endif
 50        continue

c  Record factor and active row; lower triangle is built from
c  left in cr.  pos is the corresponding column pointer.
c  (Careful!  I may *not* be able to re-use cr this way.

           nl = nl + 1
           cr(nl) = fact
           ind(nl) = row
           goto 30

 60     continue

c  Sort the nl elements of row i and keep keepl largest of them.

         keepl = min(nl, keepl)
         call qsplit( cr, ind, nl, keepl )

c  Sort the kept elements so the column indices are increasing.

         call dsortja( keepl, cr, ind )

c  Store the result and zero out workspace.

         do 70 j = 1, keepl
            l(lptr) = cr(j)
            cr(j) = zero
            jl(lptr) = ind(j)
            jnz(ind(j)) = 0
            ind(j) = 0
            lptr = lptr + 1
 70      continue

c  Zero out the rest of the workspace.

         do 80 j = keepl+1, lenl
            cr(j) = zero
            jnz(ind(j)) = 0
            ind(j) = 0
 80      continue

c  Pick out the diagonal element, store its reciprocal.
c  Note:  A seemingly better fixup strategy is described in
c  Kershaw, J. Comp. Phys., 38 (1980) p. 114, but it is not
c  yet clear how to adapt it to these data structures.  

          if ( cr(i) .eq. zero ) then
             print*, ' ILUTF:  zero pivot encountered!'
             cr(i) = rownorm
          end if
          u(uptr) = one/cr(i)
          ju(uptr) = i
          uptr = uptr + 1
          cr(i) = zero
          jnz(i) = 0

c  Sort the lenu elements of row and keep keepu of them.

         keepu = min(keepu, lenu-1)
         if ( i .lt. n ) then
            call qsplit( cr(i+1), ind(i+1), lenu-1, keepu )

c  Sort the kept elements so the column indices are increasing.

            call dsortja( keepu, cr(i+1), ind(i+1) )
         endif

         do 90 j = 1, keepu
            u(uptr) = cr(i+j)
            cr(i+j) = zero
            ju(uptr) = ind(i+j)
            jnz(ind(i+j)) = 0
            ind(i+j) = 0
            uptr = uptr + 1
 90      continue

         do 100 j = keepu+1, lenu-1
            cr(i+j) = zero
            jnz(ind(i+j)) = 0
            ind(i+j) = 0
 100     continue

 110  continue

      il(n+1) = lptr
      iu(n+1) = uptr
      nnzf = lptr+uptr
      
      return

c  End of ilutf.

      end
