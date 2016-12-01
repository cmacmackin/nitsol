
      integer instep, newstep, krystat
      double precision avrate, fcurnrm
      common /nitinfo/ avrate, fcurnrm, instep, newstep, krystat

!
! If information on the current state of the nonlinear iteration is
! desired in a user-supplied subroutine (for example, deciding 
! whether to update a preconditioner), include this common block
! in the subroutine. The variables are as follows: 
!
!     instep - inexact Newton step number. 
!
!    newstep - set to 0 at the beginning of an inexact Newton step.
!              This may be checked in a user-supplied jacv to decide
!              whether to update the preconditioner.  If you test on
!              newstep .eq. 0 to determine whether to take some 
!              special action at the beginning of a nonlinear iteration, 
!              you must also set newstep to some nonzero value to
!              subsequently avoid taking that action unnecessarily. 
!
!    krystat - status of the Krylov iteration; same as itrmks (see 
!              the nitsol documentation). 
!
!    avrate  - average rate of convergence of the Krylov solver during
!              the previous inexact Newton step.  This may be checked
!              in a user-supplied jacv to decide when to update the
!              preconditioner.
!
!    fcurnrm - ||f(xcur)||. 
!
