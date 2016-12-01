
      double precision choice1_exp, choice2_exp, choice2_coef
      double precision eta_cutoff, etamax
      double precision thmin, thmax, etafixed

      common /nitparam/ choice1_exp, choice2_exp, choice2_coef,
     &                  eta_cutoff, etamax, thmin, thmax, etafixed

!
! nitparam contains some parameters that control the nonlinear
! iterations.  In some cases, the default values reflect prevailing  
! practice; in other cases, they are chosen to produce good 
! average-case behavior.  To change the default values, include this 
! common block in the main program and set the desired variables 
! according to the following:
!
!    choice1_exp -  parameter used in the update of the forcing term 
!                   eta when ieta = 0 (default).  This is the exponent
!                   for determining the etamin safeguard.  The default
!                   value is choice1_exp = (1+sqrt(5))/2.  A larger
!                   value will allow eta to decrease more rapidly,
!                   while a smaller value will result in a larger 
!                   value for the safeguard. 
!
!    choice2_exp  - parameter used in the update of the forcing term 
!                   eta when ieta = 2.  This is the exponent alpha 
!		    in the expression gamma*(||fcur||/||fprev||)**alpha; 
!		    it is also used to determine the etamin safeguard.  
!		    The default value is 2.0. Valid values are in the 
!		    range (1.0, 2.0].
!
!    choice2_coef - parameter used in the update of the forcing term eta 
!                   when ieta = 2.  This is the coefficient gamma used 
!		    in the expression gamma*(||fcur||/||fprev||)**alpha;
!                   it is also used to determine the etamin safeguard.
!                   The default value is 1.0. Valid values are in the 
!		    range (0.0, 1.0]. 
!
!    eta_cutoff   - parameter used to determine when to disable 
!                   safeguarding the update of the forcing term.  It
!                   only has meaning when ieta .ne. 3.  The default
!                   value is 0.1.  A value of 0.0 will enable 
!		    safeguarding always; a value of 1.0 will disable 
!		    safeguarding always. 
!
!    etamax       - parameter used to provide an upper bound on the 
!		    forcing terms when input(10) .ne. 3. This is 
!		    necessary to ensure convergence of the inexact Newton 
!		    iterates and is imposed whenever eta would otherwise 
!		    be too large. (An overly large eta can result from 
!		    the updating formulas when input(10) .ne. 3 or from 
!                   safeguarding when the previous forcing term has been 
!		    excessively increased during backtracking.) The 
!		    default value of etamax is 1.0 - 1.e-4.  When 
!		    backtracking occurs several times during a nonlinear 
!		    solve the forcing term can remain near etamax for several
!                   nonlinear steps and cause the nonlinear iterations
!                   to nearly stagnate.  In such cases a smaller value of 
!                   etamax may prevent this.  Valid values are in the 
!                   range (0.0, 1.0).
!
!    etafixed     - this is the user-supplied fixed eta when ieta = 3.
!                   The  default value is etafixed = 0.1.  Valid values
!                   are in the range (0.0,1.0).
!
!    thmin        - when backtracking occurs, this is the smallest
!                   reduction factor that will be applied to the current
!                   step in a single backtracking reduction.  The default
!                   value is 0.1.  Valid  values are in the range
!                   [0.0, thmax].
!
!    thmax        - when backtracking occurs, this is the largest
!                   reduction factor that will be applied to the current
!                   step in a single backtracking reduction.  The default
!                   value is 0.5.  Valid values are in the range
!                   [thmin, 1.0).
!
