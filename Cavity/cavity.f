      program cavity

c This is a sample program for use with the NITSOL package.  It can be
c used to test an installation of NITSOL and to demonstrate use of the
c package.  For information on building the program, see the ReadMe
c file in the parent directory or consult the NITSOL User Guide.

c This program sets up and solves the streamfunction formulation of the
c driven cavity problem

c    1/Re( Delta^2 u ) - (u_y Delta u_x + u_x Delta u_y = 0 in Omega

c where Delta^2 is the biharmonic operator and Omega = [0,1]x[0,1].  The
c boundary conditions are given by

c                           u = 0, u_y = 1
c                      +----------------------+
c                      |                      |
c                      |                      |
c                      |                      | 
c                u = 0 |                      | u = 0
c                      |                      | 
c                      |                      |
c              u_x = 0 |                      | u_x = 0
c                      |                      |
c                      |                      |
c                      |                      |
c                      +----------------------+
c                             u = u_y = 0

c The problem is discretized using piecewise linear finite elements
c resulting in a system of nonlinear equations.  Code implementing this
c was provided by P. Brown, LLNL.

c The program first prompts for problem parameters:
c   nx (= ny) = the number of mesh points in (0,1) on each axis
c   Re = the Reynolds number 

c The program next prompts for some package options:
c   ikrysol - choice of Krylov iterative method
c   irpre   - optionally precondition the linearized systems

c Note:  If GMRES is chosen, the program next prompts for a maximum
c        Krylov subspace dimension (the restart value).  This being
c        Fortran, a maximum value for this is coded in a parameter
c        statement below.  If the program warns that too large a value
c        has been selected, and a larger value is desired, edit this
c        driver to change the MAXKD parameter below and recompile.
c Note:  Preconditioning is done using a fast biharmonic solver due to
c        Petter Bjorstadt.  This solver requires n to be odd and .ge. 3; 
c        the comments in the dbihar code note that "the method is 
c        somewhat faster if n+1 is a product of small primes".  Source
c        code is included in this distribution, and can also be obtained
c        from netlib.
c Note:  Analytic Jacobian-vector products have not been implemented for
c        this example.

c The program next prompts for the order of the finite-difference formula
c to be used in Jacobian-vector products.  NITSOL checks to make sure the
c input value is legal, and silently restores the default values if an
c illegal value is encountered.

c The program finally prompts for output options:
c   iplvl  - the amount of detail requested
c   ipunit - unit number to which output is sent

c Note:  The information in the prompt is self-explanatory.  NITSOL
c        checks to make sure these values are legal, and silently
c        restores the default values if an illegal value is encountered.
c
c Note:  The ipar array is used to pass nx = ipar(1) and nx+2 to the user-
c supplied routines.  The rpar array is used to pass the Reynolds number
c and provide workspace for the user-supplied routines.
c --------------------------------------------------------------------

      implicit none

      integer     MAXKD
      parameter ( MAXKD=200 )

      integer     MAXNX,    MAXNP2,         MAXN
      parameter ( MAXNX=64, MAXNP2=MAXNX+2, MAXN=MAXNX*MAXNX )

      integer     LRWORK
c>>> Alternative parameter statements for different circumstances -- HFW. 
c The following is always safe but may require a little unnecessary storage.
      parameter ( LRWORK=MAXN*(MAXKD+10)+MAXKD*(MAXKD+3))
c The following can be used if the compiler allows the "max". 
c      parameter ( LRWORK=max(11*MAXN,MAXN*(MAXKD+5)+MAXKD*(MAXKD+3)) )
c The following can be used if MAXKD > 5.
c      parameter ( LRWORK=MAXN*(MAXKD+5)+MAXKD*(MAXKD+3))

      integer     LMAX
c>>> Commented out and replaced -- HFW. 
c      parameter ( LMAX=max(MAXNP2*MAXNP2,
c     &                     max(11*MAXN+20*(MAXN+3),
c     &                         10*MAXN+((MAXN+1)**2)/2)+20) )
      parameter ( LMAX=11*MAXN+20*(MAXN+3)+20 )

      integer     LRPAR
      parameter ( LRPAR=LMAX+MAXNP2*MAXNP2+MAXN+4*MAXNX+5 )

      integer i
      integer ind
      integer iterm
      integer itrmf
      integer n
      integer neq

      integer info(6)
      integer input(10)
      integer ipar(2)

      double precision fnrm
      double precision ftol
      double precision Re
      double precision rlftol
      double precision stptol

      double precision rpar(LRPAR)
      double precision rwork(LRWORK)
      double precision y(MAXN)

      double precision dnrm2
      external dnrm2

      double precision ddot
      external ddot

      external  fcav, jacvcav

      character*8 method(0:2)
      data method /'GMRES','BiCGSTAB','TFQMR'/

c --------------------------------------------------------------------
c For printing:
      include '../Nitsol/nitprint.h'

c For internal NITSOL parameters:
      include '../Nitsol/nitparam.h'
c --------------------------------------------------------------------

c  Start of executable code-

      rlftol = 1.d-7
      stptol = 1.d-8
 10   continue

c --------------------------------------------------------------------
c Initialize. 
c --------------------------------------------------------------------
c Initialize all inputs to zero (=> default options). 

      do i = 1, 10
         input(i) = 0
      end do

c Reset particular inputs as desired. 

      write(6,*) ' Type problem parameters nx, Re (Reynold''s number): '
 20   read(5,*) n, Re
      if ( n .gt. MAXNX ) then
         write(6,*) ' Problem size too large.'
         write(6,*) ' nx must be less than ', MAXNX
         goto 20
      endif
      neq = n*n
      ipar(1) = n
      ipar(2) = n+2

      write(6,800)
      read(5,*) input(3), input(5) 

      if (input(3) .eq. 0) then 
         write(6,*) ' Type kdmax = maximum Krylov subspace dimension:'
 30      read(5,*) input(4) 
         if ( input(4) .gt. MAXKD ) then
            write(6,*) ' Maximum Krylov subspace dimension too large'
            write(6,*) ' Input value must be less than ', MAXKD
            goto 30
         endif
      endif

      write(6,*)' Type ifdord = order of finite-difference formula:'
      read(5,*) input(8) 

      write(6,810)
      read (5,*) iplvl, ipunit 
      if ( ipunit .lt. 6 ) ipunit = 6

c Complete setup. 

      call init(neq, y, Re, rpar, ipar)
      call fcav(neq, y, rwork, rpar, ipar, itrmf)
      fnrm = dnrm2(neq, rwork, 1)
      ftol = rlftol*fnrm

c Write out setup. 

      write (6,*)
      write (6,*) 
      write (6,820) method(input(3))
      if ( input(5) .eq. 0 ) then
         write (6,830)
      else
         write (6,840)
      endif
      write (6,850) input(8)
      write (6,860) n, neq, Re, stptol, ftol
      write (6,870) fnrm

c --------------------------------------------------------------------
c Call nitsol.
c --------------------------------------------------------------------

      call nitsol(neq, y, fcav, jacvcav, ftol, stptol, 
     $     input, info, rwork, rpar, ipar, iterm, ddot, dnrm2)

c --------------------------------------------------------------------
c Write results.
c --------------------------------------------------------------------

      fnrm = dnrm2(neq, rwork, 1)

      write(6,*) 
      write(6,880) iterm
      write(6,890) fnrm
      write(6,900) info(1)
      write(6,910) info(2)
      write(6,920) info(3)
      write(6,930) info(4)
      write(6,940) info(5)
      write(6,950) info(6)

c  Prompt for another case.

      write(6,*) ' To go, type 0; to stop, type 1:'
      read (5,*) ind 
      if(ind .eq. 0) go to 10

      stop

 800  format(' Type ikrysl, irpre (0-1), and ijacv (0-1):',
     $/'-------------------------------------------------------------',
     $/'ikrysl = 0 => GMRES',
     $/'         1 => BiCGSTAB', 
     $/'         2 => TFQMR', 
     $/'irpre  = 0 => no right preconditioning',
     $/'         1 => right preconditioning',
     $/'-------------------------------------------------------------')
 810  format(' Type iplvl (informational printout level), ', 
     $     'and ipunit (printout unit):'
     $/ ' iplvl = 0 => no printout, '
     $/ '       = 1 => iteration number and F-norm, '
     $/ '       = 2 => ... + stats, step-norm, lin model norm,'
     $/ '       = 3 => ... + some Krylov method and backtrack info.'
     $/ '       = 4 => ... + more Krylov method and backtrack info.')
 820  format(' Solve driven cavity problem using Newton-', a8)
 830  format(' No preconditioner is used')
 840  format(' Preconditioner is biharmonic solver with Cholesky ',
     &                                                 'decomposition.')
 850  format(' Use ',i2,'-th order finite differences for Jacobian')
 860  format(' Problem parameters:',/,
     $'   nx:', t9, i3, t21, 'n:', t24, i6, t36, 'Reynolds number: ',
     $f6.1, /, '   ftol:', t10, 1pe9.3,  t21, 'stptol:', t30, 1pe9.3 )
 870  format(' Initial f-norm:               ', 1pe9.3)
 880  format(' Termination flag iterm:       ', i9)
 890  format(' Final f-norm:                 ', t36, 1pe9.3)
 900  format(' No. function evaluations:     ', i9)
 910  format(' No. J*v evaluations:          ', i9) 
 920  format(' No. P(inverse)*v evaluations: ', i9)
 930  format(' No. linear iterations:        ', i9)
 940  format(' No. nonlinear iterations:     ', i9)
 950  format(' No. backtracks:               ', i9)
c --------------------------------------------------------------------
c End of main program. 
c --------------------------------------------------------------------
      end
