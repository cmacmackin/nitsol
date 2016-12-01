
c  This file contains auxiliary routines used by ILUTF:
c    qsplit  - quick-sort split of a real array
c              (from SPARSPAK)
c    dsortja - sort column indexes of a row stored in CSR format
c              (from SPARSPAK)
c    m1i     - lower triangular solve, with matrix in CSR format
c              (adapted from QMRPACK)
c    m2i     - upper triangular solve, with matrix in CSR format
c              (adapted from QMRPACK)

        subroutine qsplit  (a, ind, n, ncut)
        real*8 a(n)
        integer ind(n), n, ncut
c-----------------------------------------------------------------------
c     does a quick-sort split of a real array.
c     on input a(1:n). is a real array
c     on output a(1:n) is permuted such that its elements satisfy:
c
c     abs(a(i)) .ge. abs(a(ncut)) for i .lt. ncut and
c     abs(a(i)) .le. abs(a(ncut)) for i .gt. ncut
c
c     ind(1:n) is an integer array which permuted in the same way as a(*).
c-----------------------------------------------------------------------
        real*8 tmp, abskey
        integer itmp, first, last
        integer j, mid
c-----
        first = 1
        last = n
        if (ncut .lt. first .or. ncut .gt. last) return
c
c     outer loop -- while mid .ne. ncut do
c
 1      mid = first
        abskey = abs(a(mid))
        do 2 j=first+1, last
           if (abs(a(j)) .gt. abskey) then
              mid = mid+1
c     interchange
              tmp = a(mid)
              itmp = ind(mid)
              a(mid) = a(j)
              ind(mid) = ind(j)
              a(j)  = tmp
              ind(j) = itmp
           endif
 2      continue
c
c     interchange
c
        tmp = a(mid)
        a(mid) = a(first)
        a(first)  = tmp
c
        itmp = ind(mid)
        ind(mid) = ind(first)
        ind(first) = itmp
c
c     test for while loop
c
        if (mid .eq. ncut) return
        if (mid .gt. ncut) then
           last = mid-1
        else
           first = mid+1
        endif
        goto 1
c----------------end-of-qsplit------------------------------------------
c-----------------------------------------------------------------------
        end

      subroutine dsortja( n, arow, jrow )
c
c     purpose:
c     this is a sparskit subroutine to sort the column indices of a row
c     of a matrix stored in compressed sparse row (csr) format in order
c     (increasing).  the routine is given the  array of elements of the
c     row in arow, with the corresponding column indices in jrow.  both
c     arrays are of length n;  the elements are given  so that they can
c     be kept in correspondence with the elements in arow.
c     the routine uses heapsort  to carry out the sorting;  the code is
c     copied verbatim from numerical recipes (minus the bugs).
c
c     parameters:
c     n    = the length of the row (input).
c     arow = the row to be sorted (input/output).
c     jrow = the array of matching column indices (input/output).
c
c     noel m. nachtigal
c     october 28, 1990

      intrinsic abs
c
      integer n, jrow(n)
      double precision arow(n)
c
c     local variables.
c
      integer i, j, jtmp, k, l
      double precision dtmp
c
      if (n.le.1) return
c
      l = n / 2 + 1
      k = n
 10   if (l.gt.1) then
         l    = l - 1
         dtmp = arow(l)
         jtmp = jrow(l)
      else
         dtmp    = arow(k)
         jtmp    = jrow(k)
         arow(k) = arow(1)
         jrow(k) = jrow(1)
         k       = k - 1
         if (k.le.1) then
            arow(1) = dtmp
            jrow(1) = jtmp
            return
         end if
      end if
      i = l
      j = l + l
 20   if (j.le.k) then
         if (j.lt.k) then
            if (jrow(j).lt.jrow(j+1)) j = j + 1
         end if
         if (jtmp.lt.jrow(j)) then
            arow(i) = arow(j)
            jrow(i) = jrow(j)
            i       = j
            j       = j + j
         else
            j = k + 1
         end if
         go to 20
      end if
      arow(i) = dtmp
      jrow(i) = jtmp
      go to 10
c
      end

      subroutine m1i( x, n, l, il, jl, lenl, y )

c  Dummy Arguments-

      integer lenl
      integer n

      integer il(n+1)
      integer jl(lenl)

      double precision l(lenl)
      double precision x(n)
      double precision y(n)

c  Parameters-

      double precision zero
      parameter      ( zero=0.0d0 )

c  Local Variables- 

      integer i
      integer j
      integer k

      double precision dtmp

c  Start of executable code-

c  Loop over the rows.

      do 20 i = 1, n

c  Loop over the columns, up to the diagonal.

         dtmp = zero
         do 10 k = il(i), il(i+1)-1
            j    = jl(k)
            dtmp = dtmp + l(k) * y(j)
 10      continue

c     compute x(i). l(i,i) = 1.0.

         y(i) = x(i) - dtmp
 20   continue

      return
      end

      subroutine m2i( x, n, u, iu, ju, lenu )

      integer lenu
      integer n

      integer iu(n+1)
      integer ju(lenu)

      double precision u(lenu)
      double precision x(n)

c  Parameters-

      double precision zero
      parameter      ( zero=0.0d0 )

c  Local variables-

      integer i
      integer j
      integer k

      double precision dtmp

c  Start of executable code-

c  Loop over the rows.

      do 20 i = n, 1, -1
         dtmp = zero
         do 10 k = iu(i)+1, iu(i+1)-1
            j    = ju(k)
            dtmp = dtmp + u(k) * x(j)
 10      continue
         x(i) = (x(i) - dtmp)*u(iu(i))
 20   continue

      return
      end


