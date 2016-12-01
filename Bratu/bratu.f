      program bratu

c This is a sample program for use with the NITSOL package.  It can be
c used to test an installation of NITSOL and to demonstrate use of the
c package.  For information on building the program, see the ReadMe
c file in the parent directory or consult the NITSOL User Guide.

c This program sets up and solves a modified Bratu problem 

c     Delta u + d*u_x + lambda*exp(u) = 0 in Omega, 
c                                   u = 0 on (bdry) Omega,

c where Delta is the Laplacian and Omega = [0,1]x[0,1].

c The problem is discretized using second-order centered differences
c on a uniform grid, resulting in a system of nonlinear equations.

c The program first prompts for problem parameters:
c   nx (= ny) = the number of mesh points in (0,1) on each axis
c   d, lambda = problem parameters

c Note:  The size of the problem is n=nx*ny.  This being Fortran,
c        a maximum value for nx is coded in a parameter statement
c        below.  If the program warns that too large a value has
c        been selected, and a larger value is desired, edit this
c        driver to change the MAXNX parameter below and recompile.
c        If a value of ny that is different from nx is desired,
c        substantial modifications to the driver and other 
c        user-supplied routines are required.
c Note:  For d=0 no solutions exist for lambda > 6.95... and two
c        solutions exist if 0 < lambda < 6.95.  For d > 0, there
c        are no known bounds on lambda that guarantee that a 
c        solution exists (though solutions for lambda > 6.95...
c        can be found).  If d > 0 and lambda is too large, there
c        may be no solution.  This is manifested by failure of
c        the nonlinear iterations.

c The program next prompts for some package options:
c   ikrysol - choice of Krylov iterative method
c   irpre   - optionally precondition the linearized systems
c   ijacv   - finite-difference or analytic Jacobian-vector products

c Note:  If GMRES is chosen, the program next prompts for a maximum
c        Krylov subspace dimension (the restart value).  This being
c        Fortran, a maximum value for this is coded in a parameter
c        statement below.  If the program warns that too large a value
c        has been selected, and a larger value is desired, edit this
c        driver to change the MAXKD parameter below and recompile.
c Note:  Preconditioning is done using a fast Poisson solver from
c        FISHPAK.
c Note:  If finite-difference Jacobian-vector products are requested,
c        the program prompts for the desired order of finite-difference
c        formula.  NITSOL checks to make sure the input value is legal.

c The program next prompts for choice of forcing term.  Legal input
c values are:
c    0 - choice 1
c    1 - choice 2 with default coefficient and exponent
c    2 - choice 2 with user-specified coefficient and exponent
c    3 - user-specified forcing term
c See the User Guide for an explanation of the role of the forcing term.

c Note:  If choice 2 with user-specified coefficient and exponent is
c        selected, the program next prompts for the desired values.
c        NITSOL checks to make sure these values are legal, and
c        silently restores the default value if an illegal value is
c        encountered.

c The program finally prompts for output options:
c   iplvl  - the amount of detail requested
c   ipunit - unit number to which output is sent

c Note:  The information in the prompt is self-explanatory.  NITSOL
c        checks to make sure these values are legal, and silently
c        restores the default values if an illegal value is encountered.

c Note: The ipar array is used to pass nx = ipar(1) and ny = ipar(2) 
c       to the user-supplied routines.  The rpar array is used to pass 
c       three scalar PDE/discretization parameters plus work arrays of
c       dimension (nx+2)*(ny+2) and nx*ny to the user-supplied routines. 
c --------------------------------------------------------------------

      implicit none

      integer     MAXKD
      parameter ( MAXKD=50 )

      integer     MAXNX,     MAXNY,     MAXN
      parameter ( MAXNX=128, MAXNY=128, MAXN=MAXNX*MAXNY ) 

      integer     LRPAR
      parameter ( LRPAR=(MAXNX+2)*(MAXNY+2)+MAXNX*MAXNY+5 )

      integer     LRWORK 
c>>> Alternative parameter statements for different circumstances -- HFW. 
c The following is always safe but may require a little unnecessary storage.
      parameter ( LRWORK=MAXN*(MAXKD+10)+MAXKD*(MAXKD+3))
c The following can be used if the compiler allows the "max". 
c      parameter ( LRWORK=max(11*MAXN,MAXN*(MAXKD+5)+MAXKD*(MAXKD+3)) )
c The following can be used if MAXKD > 5.
c      parameter ( LRWORK=MAXN*(MAXKD+5)+MAXKD*(MAXKD+3))

      integer i
      integer ind
      integer iterm
      integer itrmf
      integer nx
      integer ny 
      integer n

      integer info(6)
      integer input(10)
      integer ipar(2)

      double precision cl
      double precision cr
      double precision d
      double precision fnrm
      double precision ftol
      double precision h
      double precision h2l
      double precision lambda
      double precision rlftol 
      double precision stptol

      double precision rpar(LRPAR)
      double precision rwork(LRWORK) 
      double precision x(MAXN)

      double precision ddot, dnrm2
      external ddot, dnrm2

      external fbratu, jacvbratu

      character*8 method(0:2)
      data method /'GMRES','BiCGSTAB','TFQMR'/

      character*8 forcing(0:3)
      data forcing /'Choice 1','Choice 2','Choice 2','Constant'/

c --------------------------------------------------------------------
c For printing:
      include '../Nitsol/nitprint.h'

c For internal NITSOL parameters:
      include '../Nitsol/nitparam.h'
c --------------------------------------------------------------------

c  Start of executable code-

      rlftol = 1.d-6
      stptol = 1.d-6
 10   continue

c --------------------------------------------------------------------
c Initialize. 
c --------------------------------------------------------------------
c Initialize all inputs to zero (=> default options). 

      do 20 i = 1, 10
         input(i) = 0
 20   continue

c Reset particular inputs as desired.

 30   write(6,*) ' Type problem parameters nx, d, lambda:'
      read (5,*) nx, d, lambda
      if ( nx .gt. MAXNX ) then
         write(6,*) ' Problem size too large.'
         write(6,*) ' nx must be less than ', MAXNX
         write(6,*) '   (nx = 0 to exit.)'
         goto 30
      else if ( nx .eq. 0 ) then
         stop
      endif

      write(6,800)
      read(5,*) input(3), input(5), input(2) 
      if (input(3) .eq. 0) then 
 40      write(6,*) ' Type maximum Krylov subspace dimension:'
         read(5,*) input(4)
         if ( input(4) .gt. MAXKD ) then
            write(6,*) ' Maximum Krylov subspace dimension too large'
            write(6,*) ' Input value must be less than ', MAXKD
            goto 40
         endif
      endif

      if (input(2) .eq. 0) then 
         write(6,*)' Type ifdord = order of finite-difference formula:'
         read(5,*) input(8) 
      endif

      write(6,*) ' Choice of forcing term:'
      read (5,*) input(10) 
      if ( input(10) .lt. 0 .or. input(10) .gt. 3 ) input(10) = 0
      if (input(10) .eq. 2) then 
         write (6,*) ' Type alpha, gamma:'
         read (5,*) choice2_exp, choice2_coef
      endif
      if (input(10) .eq. 3) then 
         write (6,*) ' Type fixed eta:'
         read (5,*) etafixed
      endif

      write(6,810) 
      read (5,*) iplvl, ipunit 
      if ( ipunit .lt. 6 ) ipunit = 6

c Complete setup. 

      ny = nx
      ipar(1) = nx
      ipar(2) = ny
      n = nx*ny
      h = 1.d0/dfloat(nx + 1) 
      cl = 1.d0 - h*d/2.d0 
      cr = 1.d0 + h*d/2.d0 
      h2l = h*h*lambda
      rpar(1) = cl
      rpar(2) = cr
      rpar(3) = h2l

c  Initial approximation.

      do i = 1, n
         x(i) = 0.0d0
      end do

      call fbratu(n, x, rwork(1), rpar, ipar, itrmf)
      fnrm = dnrm2(n, rwork(1), 1)
      ftol = rlftol*fnrm

c Write out setup. 

      write(6,*)
      write(6,*)
      write(6,820) method(input(3))
      if ( input(3) .eq. 0 ) then
         write(6,821) input(4)
      endif
      write(6,830) forcing(input(10))
      if ( input(10). eq. 1 .or. input(10) .eq. 2 ) then
         write(6,831) choice2_exp, choice2_coef
      else if ( input(10) .eq. 3 ) then
         write(6,832) etafixed
      endif
      if ( input(5) .eq. 0 ) then
         write(6,840)
      else
         write(6,850)
      endif
      if ( input(2) .eq. 0 ) then
         write(6,860) input(8)
      else
         write(6,870)
      endif
      write(6,880) nx, n, d, lambda, ftol, stptol
      write(6,890) fnrm

c --------------------------------------------------------------------
c Call nitsol.
c --------------------------------------------------------------------

      call nitsol(n, x, fbratu, jacvbratu, ftol, stptol, 
     $     input, info, rwork, rpar, ipar, iterm, ddot, dnrm2)

c --------------------------------------------------------------------
c Write results.
c --------------------------------------------------------------------

      fnrm = dnrm2(n, rwork, 1)

      write(6,*) 
      write(6,900) iterm
      write(6,910) fnrm
      write(6,920) info(1)
      write(6,930) info(2)
      write(6,940) info(3)
      write(6,950) info(4)
      write(6,960) info(5)
      write(6,970) info(6)

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
     $/'ijacv  = 0 => finite-difference J*v',
     $/'         1 => analytic J*v',
     $/'-------------------------------------------------------------')
 810  format(' Type iplvl (informational printout level), ', 
     $     'and ipunit (printout unit):'
     $/ ' iplvl = 0 => no printout, '
     $/ '       = 1 => iteration number and F-norm, '
     $/ '       = 2 => ... + stats, step-norm, lin model norm,'
     $/ '       = 3 => ... + some Krylov method and backtrack info.'
     $/ '       = 4 => ... + more Krylov method and backtrack info.')
 820  format(' Solve generalized Bratu problem using Newton-', a8)
 821  format(t8, ' GMRES restart value:  ', i3 )
 830  format(' Use ', a8, ' for forcing term')
 831  format(t8,' alpha = ', 1pe8.2, ' gamma = ', 1pe8.2)
 832  format(t8,' (Constant value = ', 1pe8.2,')')
 840  format(' No preconditioner is used')
 850  format(' Preconditioner is a fast Poisson solver')
 860  format(' Use ',i2,'-th order finite differences for Jacobian')
 870  format(' Use analytic J*v evaluations')
 880  format(' Problem parameters:',/,
     $'     nx:', t10, i3, t15, '    n:', i9, t31, '      d:', 1pe9.2,
     $                                    t48, ' lambda:', 1pe8.2,/,
     $' ftol:', 1pe9.2,  t31, ' stptol:', 1pe9.2)
 890  format(' Initial f-norm:               ', 1pe9.3)
 900  format(' Termination flag iterm:       ', i9)
 910  format(' Final f-norm:                 ', t36, 1pe9.3)
 920  format(' No. function evaluations:     ', i9)
 930  format(' No. J*v evaluations:          ', i9) 
 940  format(' No. P(inverse)*v evaluations: ', i9)
 950  format(' No. linear iterations:        ', i9)
 960  format(' No. nonlinear iterations:     ', i9)
 970  format(' No. backtracks:               ', i9)
c --------------------------------------------------------------------
c End of main program. 
c --------------------------------------------------------------------
      end
