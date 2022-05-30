!> @file
!! @brief ASWNS code
!! @details
!!    * this code is called by the users' program
!!    * this code is a modified version of XNS-v2 code (Bucciantini+Del Zanna [2011, A&A 528:A101], Pili+ [2014, MNRAS 439:3541])
!!    * modifications with respect to XNS-v2 are described in Camelio+ [2018, 2019, 2021]
!!    * code units: G = c = M☉ = kB = 1
!!
!! @author Giovanni Camelio, Niccolò Bucciantini, Antonio Pili, Luca Del Zanna
!! @date 2022
!! @copyright Giovanni Camelio, Niccolò Bucciantini, Antonio Pili, Luca Del Zanna.
!!            This code is released under the CC BY-NC-SA 4.0 license.


!> settings, type definitions, and functions for solving the rotating star
MODULE system
   IMPLICIT NONE

   ! **************** !
   ! *** SETTINGS *** !
   ! **************** !

   !> total number of radial grid points
   !! @details total = (nr/2) points in the inner grid + (nr-nr/2) in the outer grid
   INTEGER, PARAMETER :: nr = 6000
   !> number of angular mesh points (1 = 1D TOV solution)
   !! @note if @ref planar_symmetry == .false. then ieq corresponds to the equator only if nth is odd
   INTEGER, PARAMETER :: nth = 50
   !> number of Legendre polinomials for the expansion in theta (0 = 1D)
   !! @details must be 0 < mls < nth
   INTEGER, PARAMETER :: mls = 30
   !> determine the amount of logging
   LOGICAL, PARAMETER :: verbose = .true.
   !> end of the inner radial grid / begin of the outer radial grid
   DOUBLE PRECISION, PARAMETER :: rhalf = 30d0
   !> outer radius (end of the radial grid)
   DOUBLE PRECISION, PARAMETER :: rmax = 1000d0
   !> if planar_symmetry == .true., only even harmonics are used (for the metric)
   !! @details if planar_symmetry == .false. => nth should be odd (not checked)
   LOGICAL, PARAMETER :: planar_symmetry = .true.
   !> damping factor for the metric quantities in the convergence loop
   !! @details
   !!    * metric = qfactro*metric_new + (1-qfactor)*metric_old
   !!    * metric = (psi, lps*psi, btp).
   !!    * should be 0 < qfactor <= 1d0
   DOUBLE PRECISION, PARAMETER :: qfactor = .1d0
   !> criterion for the absolute convergence of the configuration
   !! @details MAX( |hden_new - hden_old| ) < abs_tol
   DOUBLE PRECISION, PARAMETER :: abs_tol = 1d-11
   !> number of iterations before start checking for convergence after hydro relaxation
   !! @details it is not a parameter because the user may want to modify it in its code
   INTEGER :: wait_iters= 50

   ! ******************************************* !
   ! *** PHYSICAL AND MATHEMATICAL CONSTANTS *** !
   ! ******************************************* !

   !> neutron mass [Msun]
   DOUBLE PRECISION, PARAMETER :: mn= 8.42335d-58
   !> nuclear saturation density [Msun^-2]
   DOUBLE PRECISION, PARAMETER :: rhon= 4.339d-4
   !> convertion from km to solar masses
   DOUBLE PRECISION, PARAMETER :: km_to_msol= 0.677218135d0
   !> greek pi
   DOUBLE PRECISION, PARAMETER :: pi=DACOS(-1d0)

   ! ********************** !
   ! *** PUBLIC GLOBALS *** !
   ! ********************** !

   !> radial grid (1 ghost zone per side)
   DOUBLE PRECISION :: r(0:nr+1)
   !> angular grid (1 ghost zone per side)
   DOUBLE PRECISION :: th(0:nth+1)

   ! *********************** !
   ! *** PRIVATE GLOBALS *** !
   ! *********************** !

   !> radial grid increments
   DOUBLE PRECISION, PRIVATE :: dr(1:nr+1)
   !> angular grid increments
   DOUBLE PRECISION, PRIVATE :: dth(1:nth+1)
   !> radial ODE coefficients for the main diagonal
   DOUBLE PRECISION, PRIVATE, DIMENSION(nr) :: dcp
   !> radial ODE coefficients for the lower diagonal
   DOUBLE PRECISION, PRIVATE, DIMENSION(nr-1) :: dlp
   !> radial ODE coefficients for the upper diagonal
   DOUBLE PRECISION, PRIVATE, DIMENSION(nr-1) :: dup
   !> radial differentiation coefficients
   DOUBLE PRECISION, PRIVATE, DIMENSION(nr) :: a1,a2,a3
   !> angular differentiation coefficients
   DOUBLE PRECISION, PRIVATE, DIMENSION(nth) :: b1,b2,b3
   !> Gauss quadrature weights
   DOUBLE PRECISION, PRIVATE, DIMENSION(nth) :: wgq
   !> spherical harmonics -- more or less
   !! @details
   !! * pn_init and pd_init contains the quadrature weights
   !! * pdx_init is divided by sin(th) to avoid problems close to the poles
   DOUBLE PRECISION, PRIVATE, DIMENSION(0:mls,nth) :: pn_init, pd_init, pnx_init, pdx_init

   ! *************** !
   ! *** TYPEDEF *** !
   ! *************** !

   !> stellar total quantities and profiles
   TYPE model_quantities
      !> gravitational mass
      DOUBLE PRECISION :: mg
      !> proper mass
      DOUBLE PRECISION :: mp
      !> rotational energy
      DOUBLE PRECISION :: tr
      !> angular momentum
      DOUBLE PRECISION :: j
      !> keplerian angular velocity
      DOUBLE PRECISION :: kep_omega
      !> equatorial angular velocity
      DOUBLE PRECISION :: eq_omega

      !> error flag
      !! @details
      !! * 0: relaxed
      !! * 1: initialized
      !! * 2: mass shedding
      !! * -1: not converged
      !! * -2: error
      !! * -3: not initialized
      INTEGER :: info

      !> the surface of the rotating star is between wsurf and wsurf+1
      INTEGER, DIMENSION(nth) :: wsurf

      !> index of the equator
      INTEGER :: ieq

      !> pressure profile
      DOUBLE PRECISION, DIMENSION(nth,nr) :: p
      !> enthalpy density profile
      DOUBLE PRECISION, DIMENSION(nth,nr) :: hden
      !> contravariant speed profile
      DOUBLE PRECISION, DIMENSION(nth,nr) :: vphi
      !> conformal factor profile
      DOUBLE PRECISION, DIMENSION(nth,nr) :: psi
      !> profile of third component of the shift beta^phi = minus ZAMO angular velocity -omega
      DOUBLE PRECISION, DIMENSION(nth,nr) :: btp
      !> lapse profile
      DOUBLE PRECISION, DIMENSION(nth,nr) :: lps
   END TYPE model_quantities

   ! ****************** !
   ! *** INTERFACES *** !
   ! ****************** !

   INTERFACE

      !> EOS for the TOV solver @ref tov
      !! @param[in] p pressure
      !! @returns enthalpy density
      PURE DOUBLE PRECISION FUNCTION fun_hden(p)
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(IN) :: p
      END FUNCTION fun_hden

      !> EOS for @ref compute_quantities
      !! @param[in] p pressure
      !! @param[in] hden enthalpy density
      !! @param[out] rho rest mass density
      !! @param[out] s entropy per baryon
      !! @param[out] t temperature
      PURE SUBROUTINE eos(p,hden,rho,s,t)
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(IN) :: hden,p
         DOUBLE PRECISION, INTENT(OUT) :: rho,s,t
      END SUBROUTINE eos
      
      !> solve the force balance equation
      !! @param[in] iter number of the current iteration in @ref hydroeq
      !! @param[in] q0 Euler potential in the center
      !! @param[in] alp2 lapse to the second power α²
      !! @param[in] btp shift component in the φ direction, equivalent to -ω (minus the ZAMO angular momentum)
      !! @param[in] rcyl2 cylindrical radius to the second power R²
      !! @param[inout] p pressure, input: first guess, output: result
      !! @param[inout] omg angular velocity, input: first guess, output: result
      !! @param[out] hden enthalpy density
      !! @param[out] info success flag (0 on success)
      PURE SUBROUTINE euler(iter,q0,alp2,btp,rcyl2,p,omg,hden,info)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: iter
         DOUBLE PRECISION, INTENT(IN) :: q0,alp2,btp,rcyl2
         DOUBLE PRECISION, INTENT(INOUT) :: p,omg
         DOUBLE PRECISION, INTENT(OUT) :: hden
         INTEGER, INTENT(OUT) :: info
      END SUBROUTINE euler

      !> subroutine definition for using lapack in a pure subroutine
      PURE SUBROUTINE dgtsv(n,nrhs,dl,d,du,b,ldb,info)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: n, nrhs, ldb
         DOUBLE PRECISION, INTENT(INOUT) :: dl(n-1),d(n),du(n-1),b(ldb,nrhs)
         INTEGER, INTENT(OUT) :: info
      END SUBROUTINE dgtsv

      !> dirty trick for printing on screen in a pure subroutine
      !! @param[in] text string that will be printed on screen
      PURE SUBROUTINE print_on_screen(text)
         IMPLICIT NONE
         CHARACTER(LEN=*), INTENT(IN) :: text
      END SUBROUTINE print_on_screen

   END INTERFACE

   PRIVATE dgtsv, print_on_screen
   
   SAVE

CONTAINS

   !> initialize ASWNS
   !! @details
   !! * consistency cheks
   !! * initialize grid
   !! * initialize Gaus-Legendre weights and knots
   !! * initialize Legendre polynomials
   !! * stellar metric and matter initialization with TOV solution
   !!
   !! @note the EOS should be already initialized
   !!
   !! @param[out] a stellar profile and quantities
   !! @param[in] p2hden barotropic EOS for TOV solver @ref tov
   !! @param[in] p0 central pressure of the TOV solution
   !! @param[out] info success flag (0 on success)
   SUBROUTINE ns_init(a,p2hden,p0,info)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: p0
      PROCEDURE(fun_hden) :: p2hden
      TYPE(model_quantities), INTENT(OUT) :: a
      INTEGER, INTENT(OUT) :: info
      LOGICAL :: correct_init_TOV
      DOUBLE PRECISION :: temp,dj
      DOUBLE PRECISION, DIMENSION(nr) :: hdentv,ptv,psitv,nutv
      INTEGER :: i,iz,isur,info_leg,j
      INTEGER, PARAMETER :: nrhalf = nr/2
      DOUBLE PRECISION, DIMENSION(0:mls), PARAMETER :: l=(/ (DBLE(i), i=0, mls, +1) /)
      DOUBLE PRECISION, DIMENSION(nth) :: xgq
      DOUBLE PRECISION, DIMENSION(nth*2 - 1) :: tmpx, tmpw

      ! ************************** !
      ! *** CONSISTENCY CHECKS *** !
      ! ************************** !

      info= -1
      a%info= -3

      IF (nr<=0) THEN
         PRINT*,"ERROR: must be 0 < nr"
         RETURN
      ENDIF
      IF (nth <= 0) THEN
         PRINT*,"ERROR: must be 0 < nth"
         RETURN
      ENDIF
      IF (rhalf>=rmax .OR. rhalf<=0d0) THEN
         PRINT*,"ERROR: must be 0 < rhalf < rmax"
         RETURN
      ENDIF
      IF (mls<=0 .OR. mls>=nth) THEN
         PRINT*,"ERROR: must be 0 < mls < nth"
         RETURN
      ENDIF
      IF (qfactor<=0 .OR. qfactor>1d0) THEN
         PRINT*,"ERROR: must be 0 < qfactor <= 1"
         RETURN
      ENDIF
      IF (abs_tol < 0d0) THEN
         PRINT*,"ERROR: must be abs_tol > 0"
         RETURN
      ENDIF

      ! ************ !
      ! *** GRID *** !
      ! ************ !

      ! inner radial grid
      DO i=1,nrhalf
         dr(i) = rhalf / DBLE(nrhalf)
      ENDDO

      ! outer radial grid
      temp = 2d0 * (rmax - rhalf - dr(nrhalf)*DBLE(nr-nrhalf)) / DBLE( (nr-nrhalf) * (nr-nrhalf+1) )
      DO i=nrhalf+1,nr+1
         dr(i) = dr(nrhalf) + temp * DBLE(i-nrhalf)
      ENDDO

      ! set the radial grid
      r(0) = - 0.5d0 * dr(1)
      DO i=1,nr+1
         r(i) = r(i-1) + dr(i)
      ENDDO

      IF(planar_symmetry) THEN
         a%ieq= nth
      ELSE
         a%ieq= nth/2 + 1 ! this works also when nth == 1
      ENDIF
      
      ! compute the Gauss-quadrature points XGQ and weights WGQ
      IF(planar_symmetry) THEN
         ! if planar_symemtry, xgq(nth)= 0 always
         CALL legendre_knots(2*nth-1,tmpx,tmpw,info_leg)
         xgq(:)= tmpx(1:nth)
         ! here I account for having only northen emisphere:
         wgq(1:nth-1)= 2*tmpw(1:nth-1)
         ! central point is always in x= 0 and is not doubled:
         wgq(nth)= tmpw(nth)
      ELSE
         CALL legendre_knots(nth,xgq,wgq,info_leg)
      ENDIF
      IF(info_leg /= 0) THEN
         info= -2
         PRINT*,"ERROR in legendre_knots!"
         RETURN
      ENDIF

      ! angular grid on the Legendre knots
      th(1:nth)=DACOS(xgq(1:nth))

      ! boundary constraints
      th(0)=-th(1)
      IF(planar_symmetry) THEN
         th(nth+1)= pi - th(nth-1)
      ELSE
         th(nth+1)= 2*pi - th(nth)
      ENDIF

      ! dtheta
      DO i=1,nth+1
         dth(i)=th(i)-th(i-1)
      ENDDO

      ! spherical harmonics on the quadrature points
      DO i=1,nth
         pn_init(0,i)= 1d0
         pn_init(1,i)= xgq(i)
         pd_init(0,i)= 0d0
         pd_init(1,i)= 1d0
         DO j=2, mls
            dj= DBLE(j)
            pn_init(j,i)= (2d0 - 1d0/dj)*xgq(i)*pn_init(j-1,i) - (1d0 - 1d0/dj)*pn_init(j-2,i)
            pd_init(j,i)= dj*pn_init(j-1,i) + xgq(i)*pd_init(j-1,i)
         ENDDO
         pnx_init(0:mls,i) = pn_init(0:mls,i)
         pdx_init(0:mls,i) = pd_init(0:mls,i)
         ! now I computed the Legendre polynomials and their derivatives

         ! spherical harmonics * Legendre weight
         pn_init(0:mls,i)=pn_init(0:mls,i)*wgq(i)*DSQRT( (2d0*l(0:mls) + 1d0) / (4d0*pi) )
         ! associated spherical harmonics m=1 * Legendre weight
         pd_init(0:mls,i)=-pd_init(0:mls,i)*wgq(i)*DSQRT(1d0-xgq(i)**2)*DSQRT( (2d0*l(0:mls) + 1d0) / (4d0*pi) )

         ! pnx = Legendre polynomial of degree l
         ! here I left out a factor sin(th), because it is simplified later by a division for sin(th)
         pdx_init(1:mls,i)=-DSQRT((2d0*l(1:mls)+1d0)/4d0/pi)*pdx_init(1:mls,i)
         ! spherical harmonic
         pnx_init(0:mls,i)=pnx_init(0:mls,i)*DSQRT( (2d0*l(0:mls) + 1d0) / (4d0*pi) )
      ENDDO

      ! *********************** !
      ! *** DIFFERENTIATION *** !
      ! *********************** !

      ! differentiation along r
      DO I=1,NR
         ! for derive
         A1(i)=-DR(I+1)/DR(I)/(DR(I)+DR(I+1))
         A2(i)= (DR(I+1)-DR(I))/DR(I)/DR(I+1)
         A3(i)= DR(I)/DR(I+1)/(DR(I)+DR(I+1))
         ! for DGTSV
         DCP(I)  = -2d0/DR(I)/DR(I+1) + 2d0*A2(i)/R(I)
         IF(i/=1)  DLP(I-1)  = 2d0/DR(I)/(DR(I)+DR(I+1)) + 2d0*A1(i)/R(I)
         IF(i/=nr) DUP(I)    = 2d0/DR(I+1)/(DR(I)+DR(I+1)) + 2d0*A3(i)/R(I)
      ENDDO

      ! differentiation along theta
      DO i=1,nth
         ! for derive
         B1(i)=-DTH(I+1)/DTH(I)/(DTH(I)+DTH(I+1))
         B2(i)= (DTH(I+1)-DTH(I))/DTH(I)/DTH(I+1)
         B3(i)= DTH(I)/DTH(I+1)/(DTH(I)+DTH(I+1))
      ENDDO

      ! ***************************** !
      ! *** INITIAL CONFIGURATION *** !
      ! ***************************** !
    
      ! initialize stellar quantities and atmosphere with TOV solution
      ! TOV initialization, in isotropic coordinates with the original code
      CALL tov(p0,p2hden,ptv,hdentv,nutv,psitv,isur,correct_init_tov)
      IF(correct_init_TOV.EQV..false.) THEN
         info= -3
         PRINT*,"ERROR in tovini!"
         RETURN
      ENDIF

      ! load TOV parameters
      DO iz=1,nr
         ! matter
         a%p(:,iz)   = ptv(iz)
         a%hden(:,iz)= hdentv(iz)
         ! metric
         a%psi(:,iz) = psitv(iz) ! conformal factor
         a%lps(:,iz) = DEXP(nutv(iz))
      ENDDO
      a%vphi(:,:)= 0d0
      a%btp(:,:) = 0d0
      ! the stellar surface is between isur and isur+1
      a%wsurf(:)=isur

      ! compute stellar quantities
      CALL compute_quantities(a)

      info= 0
      a%info= 1
   END SUBROUTINE ns_init


   !> rotating structure of a neutron star
   !!
   !! @details
   !!    * the rotating structure is determined throught iteration between the
   !!      metric solver and the force balance solver
   !!    * reference BdZ11 = Bucciantini & Del Zanna (2011), A&A 528:A101
   !!
   !! @param[inout] a stellar profiles and quantities
   !! @param[in] fun function for solving the force balance equation
   !! @param[out] info success flag (0 on success)
   PURE SUBROUTINE ns_solve(a,fun,info)
      !USE system, ONLY: nth,nr,pi,th,r,qfactor,model_quantities,verbose,abs_tol
      IMPLICIT NONE
      PROCEDURE(euler) :: fun
      TYPE(model_quantities), INTENT(INOUT) :: a
      INTEGER, INTENT(OUT) :: info
      DOUBLE PRECISION, DIMENSION(nth,nr) :: tmp,tmp2,dwphi_dr,dwphi_dth,r2,sin2,epp,psl
      INTEGER, PARAMETER :: maxloop = 5000
      INTEGER :: iloop, ix, iz, ihydro
      DOUBLE PRECISION :: max_old_minus_new
      CHARACTER(LEN=255) :: string

      DO ix=1,nth
         sin2(ix,:) = DSIN(th(ix))**2
      ENDDO
      DO iz=1,nr
         r2(:,iz) = r(iz)**2
      ENDDO

      psl= a%psi*a%lps
      ! cycle untill convergence
      DO iloop=1,maxloop

         ! ****************** !
         ! *** find W^phi *** !
         ! ****************** !

         tmp= a%vphi*a%vphi*r2*sin2*a%psi**4 ! tmp = v2 -- physical speed^2
         epp = a%hden/(1d0 - tmp) ! E + p = GLF^2 * (eps + p)
         ! solve phi-component of vector poisson
         CALL shiftphi( & ! BdZ11:Eq:27
            & 8d0*pi*epp*a%vphi*a%psi**10, & ! source
            & tmp) ! tmp = W^phi in coordinate basis

         ! **************** !
         ! *** find psi *** !
         ! **************** !

         CALL derive(tmp,dwphi_dr,dwphi_dth) ! tmp = W^phi in coordinate basis
         tmp = 2d0*sin2*(r2*dwphi_dr**2+dwphi_dth**2) ! tmp = curvc
         ! solve poisson equation
         tmp2= a%psi ! tmp2= old
         CALL scalar_laplacian( & ! BdZ11:Eq:28
            & -2d0*pi*a%psi**5*(epp-a%p) - tmp/(8d0*a%psi**7), & ! source -- tmp = curvc
            & a%psi)
         a%psi(:,:)=1d0+a%psi(:,:) ! conformal factor psi
         a%psi=qfactor*a%psi+(1d0-qfactor)*tmp2 ! relax the metric -- tmp2 = old

         ! **************************** !
         ! *** find psl = psi*alpha *** !
         ! **************************** !

         ! recompute some matter quantities after updating psi
         tmp2=a%vphi*a%vphi*r2*sin2*a%psi**4 ! tmp2 = v2 -- physical speed^2
         epp = a%hden/(1d0 - tmp2) ! E + p == GLF^2 * (eps + p)
         tmp = ( 2d0*pi*(epp*(1d0+2d0*tmp2) + 5d0*a%p)*a%psi**4 & ! tmp(LHS) = source; tmp2 = v2 -- physical speed^2
             &   + 7d0/8d0*tmp/a%psi**8 & ! tmp(RHS) = curvc
             & ) * psl
         ! Solve poisson equation
         tmp2= psl ! tmp2 = old
         CALL scalar_laplacian(tmp,psl) ! tmp = source
         psl= 1d0 + psl ! psl = alpha*psi
         psl=qfactor*psl + (1d0 - qfactor)*tmp2 ! relax the metric -- tmp2 = old
         a%lps= psl/a%psi

         ! *************************** !
         ! *** find btp = beta^phi *** !
         ! *************************** !

         CALL derive(a%lps/a%psi**6, tmp, tmp2) ! tmp= d_dr, tmp2= d_dth
         tmp = 16d0*pi*a%lps*a%psi**4*epp*a%vphi & ! tmp(LHS) = curvp
             & + 2d0*(dwphi_dr*tmp + dwphi_dth*tmp2/r2) ! tmp(RHS) = d_dr, tmp2= d_dth
         ! solve phi-component of vector poisson 
         tmp2= a%btp ! tmp2 = old
         CALL shiftphi(tmp,a%btp) ! tmp = curvp
         ! BEWARE: from here, btp = beta^phi
         a%btp=qfactor*a%btp+(1d0-qfactor)*tmp2 ! relax the metric -- tmp2 = old

         ! ************************************************************** !
         ! *** invert the first integral and find new matter and vphi *** !
         ! ************************************************************** !

         tmp= a%hden ! tmp = old
         ! Solve the hydrostatic equilibrium in the newly computed CFC metric
         CALL hydroeq(iloop,fun,a,ihydro)

         ! ************************* !
         ! *** check convergence *** !
         ! ************************* !

         IF(ihydro < 0) THEN
            IF(verbose) THEN
               WRITE(string,*)"ERROR in hydroeq. ihydro = ",ihydro
               CALL print_on_screen(string)
            ENDIF
            info= ihydro
            EXIT
         ENDIF

         max_old_minus_new= MAXVAL(DABS(tmp - a%hden)) ! tmp = old

         IF(verbose) THEN
            WRITE(string,'(A,I4,4(A,E10.4))')'         > ',iloop,&
               & ": max|hden - old hden| = ",max_old_minus_new
            CALL print_on_screen(string)
         ENDIF

         ! check for runaway solutions
         IF(a%hden(1,1) < abs_tol) THEN
            IF(verbose) THEN
               WRITE(string,*)"BEWARE! runaway! ihydro = ",ihydro
               CALL print_on_screen(string)
            ENDIF
            info= -2
            EXIT
         ENDIF

         IF(MAXVAL(a%wsurf) == 1) THEN
            IF(verbose) THEN
               WRITE(string,*)"BEWARE! no star!"
               CALL print_on_screen(string)
            ENDIF
            info= 3
            EXIT
         ENDIF

         ! check for convergence
         IF((iloop > wait_iters) .AND. (max_old_minus_new < abs_tol)) THEN
            info= 0
            EXIT
         ENDIF

         ! when there is mass shedding, it may happen that the convergence
         ! is very oscillatory, if any.
         ! It may happen that the code start ciclying between configurations
         ! with mass shedding -- no mass shedding, or between different
         ! mass-shedding configuration. In that case, I exit before
         ! reaching the convergence.
         IF((iloop > wait_iters) .AND. ihydro == 2) THEN
            info= 2
            EXIT
         ENDIF
      ENDDO

      ! convergence check
      IF(iloop == maxloop + 1) THEN
         IF(verbose) THEN
            WRITE(string,*)"ERROR: no convergence reached"
            CALL print_on_screen(string)
         ENDIF
         info= -1
      ENDIF

      a%info= info
      ! compute model quantities
      CALL compute_quantities(a)
   END SUBROUTINE ns_solve


   !> solve the force balance equation (Euler equation)
   !!
   !! @details inverts the first integral to determine the matter quantities,
   !! that is the pressure, the enthalpy density, and the contravariant velocity
   !!
   !! @param[in] iter number of the iteration in @ref ns_solve
   !! @param[in] fun user-provided solver of the Euler equation
   !! @param[inout] a stellar model
   !! @param[out] info success flag (0 on success)
   PURE SUBROUTINE hydroeq(iter,fun,a,info)
      !USE system, ONLY: nth,nr,r,th,model_quantities
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iter
      PROCEDURE(euler) :: fun
      TYPE(model_quantities), INTENT(INOUT) :: a
      INTEGER, INTENT(OUT) :: info
      DOUBLE PRECISION :: q0,alp,betaphi,alp2,rcyl2,omg,p,vphi,hden
      INTEGER :: ix,iz,ifun
      
      ! derive the central value of the metric via quadratic interpolation
      q0= -DLOG((a%lps(a%ieq,1)*R(2)**2 - a%lps(a%ieq,2)*R(1)**2)/(R(2)**2 - R(1)**2))

      ! Compute the hydro variables at all radii
      info= 0
      a%wsurf(:)= nr
      radial: DO iz=1,nr
         angular: DO ix=1,nth
            IF(iz > a%wsurf(ix)) CYCLE angular

            ! metric in r, theta
            alp= a%lps(ix,iz)
            alp2= alp*alp
            betaphi= a%btp(ix,iz)
            rcyl2= (r(iz)*DSIN(th(ix))*a%psi(ix,iz)*a%psi(ix,iz))**2

            ! first guess on p, Omega
            omg= alp*a%vphi(ix,iz) - betaphi
            p= a%p(ix,iz)

            ! solve the integer of motion equations
            CALL fun(iter,q0,alp2,betaphi,rcyl2,p,omg,hden,ifun)
            IF(ifun < 0) THEN
               info= ifun
               EXIT radial
            ENDIF

            ! first check: non superluminar motion
            vphi= (omg + betaphi)/alp
            IF(vphi*vphi*rcyl2 > 1d0) THEN
               info= -3
               EXIT radial
            ENDIF

            ! second check: p >= surface_c and meaningful hden
            IF (p <= 0d0 .OR. hden <= p .OR. ifun > 0) THEN
               ! this is not an error but means that you reached the surface
               a%wsurf(ix)= iz - 1
               ! the surface is between wsurf and wsurf + 1
               ! note that p, hden, vphi are already set to zero
               CYCLE angular
            ENDIF

            ! third check: mass shedding
            IF(ix == a%ieq .AND. omg > get_kep_omega(a,iz)) THEN
               info = 2
               EXIT radial
            ENDIF

            ! now we are sure we are inside the star
            a%p(ix,iz)= p
            a%hden(ix,iz)= hden
            a%vphi(ix,iz)= vphi 
         ENDDO angular ! cycle on radius coordinate
      ENDDO radial ! cycle on angle coordinate

      ! set to zero outside the surface/truncation radius
      DO ix=1, nth
         a%wsurf(ix)= MIN(a%wsurf(ix), iz-1)
         ! outside the surface the matter is zero
         a%p   (ix, a%wsurf(ix)+1:nr)= 0d0
         a%hden(ix, a%wsurf(ix)+1:nr)= 0d0
         a%vphi(ix, a%wsurf(ix)+1:nr)= 0d0
      ENDDO
   END SUBROUTINE hydroeq


   !> stellar quantities and profiles
   !! @param[inout] a stellar model
   PURE SUBROUTINE compute_quantities(a)
      !USE system, ONLY: r,th,dr,nth,pi,model_quantities,wgq
      IMPLICIT NONE
      TYPE(model_quantities), INTENT(INOUT) :: a
      DOUBLE PRECISION :: GLF,DET,V,OMGN,SINIX,temp
      INTEGER :: izeq,ix,iz

      ! ************** !
      ! *** GLOBAL *** !
      ! ************** !

      a%Mg=0d0 ! gravitational mass
      a%Mp=0d0 ! proper mass
      a%Tr=0d0  ! rotational energy
      a%J=0d0  ! angular momentum

      DO ix=1,nth
         sinix=DSIN(th(ix))

         DO iz=1,a%wsurf(ix)!+2
            ! preliminaries
            v=a%vphi(ix,iz)*r(iz)*a%psi(ix,iz)**2*sinix
            glf=1d0/DSQRT(1d0-v*v)
            omgn=a%lps(ix,iz)*a%vphi(ix,iz)-a%btp(ix,iz)
            det=a%psi(ix,iz)**6*r(iz)**2*(dr(iz)+dr(iz+1))/2d0*wgq(ix) ! dx = wgq = d(cos(th)) = -sin(th) dth
            temp=a%hden(ix,iz)*glf**2

            ! quantities
            ! compare with Gourgoulhon[2011, notes] Eq. (4.18)
            a%Mg=a%Mg+2d0*pi*(a%lps(ix,iz)*(2d0*a%p(ix,iz)+temp*(1d0+v*v)) &
                &-2d0*a%btp(ix,iz)*temp*a%psi(ix,iz)**2*v*sinix*r(iz))*det
            ! Friedman+Ipser+Parker[1986,ApJ] Eq. (22)
            a%Mp=a%Mp+2d0*pi*(a%hden(ix,iz)-a%p(ix,iz))*glf*det
            ! Friedman+Ipser+Parker[1986,ApJ] Eq. (20)
            a%Tr=a%Tr+pi*temp*v*omgn*r(iz)*sinix*det*a%psi(ix,iz)**2
            ! compare with Gourgoulhon[2010, notes] Eq (4.39)
            a%J=a%J+2d0*pi*temp*v*r(iz)*sinix*det*a%psi(ix,iz)**2
         ENDDO
      ENDDO
      ! For the binding energy, see Friedman+Ipser+Parker[1986,ApJ] Eq. (21)
      ! binding_energy = Mp + Tr - Mg

      ! ***************** !
      ! *** KEPLERIAN *** !
      ! ***************** !

      izeq=a%wsurf(a%ieq) ! index of the stellar border at the equator -- surface is between iz and iz+1
      a%eq_omega=a%lps(a%ieq,izeq)*a%vphi(a%ieq,izeq)-a%btp(a%ieq,izeq) ! omega @ equator
      a%kep_omega= get_kep_omega(a,izeq)
   END SUBROUTINE compute_quantities


   !> keplerian angular velocity
   !! @param[in] a stellar model
   !! @param[in] i angular grid index of the equator
   !! @results keplerian angular velocity
   PURE DOUBLE PRECISION FUNCTION get_kep_omega(a,i) RESULT(kep_omega)
      !USE system, ONLY: planar_symmetry, model_quantities,r
      IMPLICIT NONE
      TYPE(model_quantities), INTENT(IN) :: a
      INTEGER, INTENT(IN) :: i
      DOUBLE PRECISION, DIMENSION(-1:1) :: nufip,psifip
      DOUBLE PRECISION :: tmp,dnufip,dpsifip,domega

      ! reference: FIP86 = Friedman+Ipser+Parker [1986, ApJ]
      ! metric conversion: BDZ11:Eq:84 --> FIP86:Eq:3
      nuFIP=DLOG(a%lps(a%ieq,i-1:i+1)) ! nu in FIP86:Eq:3
      psiFIP=DLOG(a%psi(a%ieq,i-1:i+1)**2*r(i-1:i+1)) ! psi in FIP86:Eq:3
      ! metric derivatives
      dnuFIP  = nuFIP(+1)-nuFIP(-1)
      dpsiFIP = psiFIP(+1)-psiFIP(-1)
      domega  = - (a%btp(a%ieq,i+1) - a%btp(a%ieq,i-1)) ! omega=-beta^phi
      ! Keplerian velocity and angular velocity in GR: FIP86:Eq:23
      tmp=DEXP(2d0*(nuFIP(0)-psiFIP(0)))
      kep_omega = -a%btp(a%ieq,i) + 0.5d0*domega/dpsiFIP + &
                & DSIGN(DSQRT(dnuFIP/dpsiFIP*tmp + (0.5d0*domega/dpsiFIP)**2), a%eq_omega)! c=1
   END FUNCTION get_kep_omega


   !> total stellar rest mass, total stellar entropy, total rest mass of the disk
   !! @details see Gourgoulhon[2011, notes]
   !! @param[in] a stellar model
   !! @param[in] p_hden2rho_s EOS from rest mass density and enthalpy density
   !! @param[out] m0 total stellar rest mass
   !! @param[out] stot total stellar entropy
   !! @param[out] total rest mass of the disk
   PURE SUBROUTINE get_m0_stot(a,p_hden2rho_s,m0,stot,m0disk)
      !USE system, ONLY: model_quantities, nth, th, r, dr, wgq, pi, mn
      IMPLICIT NONE
      TYPE(model_quantities), INTENT(IN) :: a
      PROCEDURE(eos) :: p_hden2rho_s
      DOUBLE PRECISION, INTENT(OUT) :: m0, stot, m0disk
      CHARACTER(LEN=255) :: string
      DOUBLE PRECISION :: sinth, glf, det, rho, rcyl2, s, jisco, eisco, j_over_m, t
      INTEGER :: i,j

      m0= 0d0
      stot= 0d0
      m0disk= 0d0

      ! ISCO quantities

      CALL kerr_isco(a%mg, a%j, jisco, eisco)

      DO i= 1, nth
         sinth= DSIN(th(i))

         DO j= a%wsurf(i), 1, -1
            ! preliminaries
            rcyl2= (r(j) * a%psi(i,j)**2 * sinth)**2
            glf= 1d0 / DSQRT(1d0 - a%vphi(i,j) * a%vphi(i,j) * rcyl2)
            det= a%psi(i,j)**6 * r(j) * r(j) * 0.5d0 * (dr(j) + dr(j+1)) * wgq(i)

            ! get rho and s
            CALL p_hden2rho_s(a%p(i,j), a%hden(i,j), rho, s, t)
            ! check the solution (now this is not needed)
            IF(rho < 0d0 .OR. s < 0d0) THEN
               WRITE(string,*)"BEWARE! unphysical rho,s at ith, ir= ",i,j," rho= ",rho,"s= ",s
               CALL print_on_screen(string)
            ENDIF

            m0= m0 + 2d0 * pi * rho * glf * det
            stot= stot + 2d0 * pi * rho * s / mn * glf * det

            ! j/m, see Gourgoulhon[2011]:eq:3.85
            j_over_m= a%hden(i,j) / rho * a%vphi(i,j) * glf * rcyl2
            ! check disk formation condition
            IF(j_over_m > jisco) THEN
               ! Gourgoulhon[2011, notes] Eq. (4.5)
               m0disk= m0disk + 2d0 * pi * rho * glf * det
            ENDIF
         ENDDO
      ENDDO
   END SUBROUTINE get_m0_stot


   !> ISCO quantities for corotating orbits for a Kerr BH
   !! @details reference: BPT72 = Bardeen, Press & Teukolsky [1972, APJ]
   !! @note r is not the isotropic radius!
   !! @param[in] m gravitational mass
   !! @param[in] j angular momentum
   !! @param[out] jisco specific angular momentum at the ISCO
   !! @param[out] eisco specific energy at the ISCO
   PURE SUBROUTINE kerr_isco(m, j, jisco, eisco)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: m,& ! stellar gravitational mass
                                      j   ! stellar angular momentum
      DOUBLE PRECISION, INTENT(OUT) :: jisco,eisco ! specic angular momentum at the ISCO
      DOUBLE PRECISION :: l,& ! specific angular momentum at ISCO (l=L/m)
                          r,& ! ISCO areal radius
                          e   ! specific energy at ISCO (e=E/m)
      DOUBLE PRECISION :: a,z1,z2,denom

      a=j/m ! angular momentum per unit mass
      ! direct (corotating orbits)
      ! BPT72:eq:2.21
      z1=1d0+(1d0-(a/m)**2)**(1d0/3d0)*((1+a/m)**(1d0/3d0)+(1d0-a/m)**(1d0/3d0))
      z2=DSQRT(3d0*(a/m)**2+z1**2) ! here there WAS the error: z1->z1**2
      r=m*(3d0+z2-DSQRT((3d0-z1)*(3d0+z1+2d0*z2))) ! ISCO areal radius
      ! BPT72:eq:2.12
      denom=r**(0.75d0)*DSQRT(r**1.5d0-3d0*m*DSQRT(r)+2d0*a*DSQRT(m))
      e=(r**1.5d0-2d0*m*DSQRT(r)+a*DSQRT(M))/denom ! ISCO energy per unit mass
      ! BPT72:eq:2.13
      l=DSQRT(m)*(r*r-2d0*a*DSQRT(m*r)+a*a)/denom ! ISCO angular momentum per unit mass
      ! specific angular momentum at the ISCO
      jisco=l ! L/m
      eisco=e ! E/m
   END SUBROUTINE kerr_isco


   !> solve the laplacian of a scalar in spherical coordinates with a series expansion
   !! in Legendre polynomials in the angle and by matrix inversion in the radius
   !! @details
   !!    assumptions:
   !!    * axisymmetry
   !!    * function --> r^(-l(l+1))
   !!    * function even in r and theta
   !! @param[in] source source of the laplacian
   !! @param[out] fun solution of the laplacian
   PURE SUBROUTINE scalar_laplacian(source,fun)
      !USE system, ONLY: NTH,NR,MLS,R,DR,pn_init,pnx_init,dcp,dup,dlp,pi,planar_symmetry
      IMPLICIT NONE
      INTEGER :: I,INFO
      DOUBLE PRECISION,DIMENSION(nth,nr), INTENT(IN):: source
      DOUBLE PRECISION,DIMENSION(nth,nr), INTENT(OUT):: fun
      INTEGER :: ix,iz
      ! Array for the Legendre expansion
      DOUBLE PRECISION,DIMENSION(NR,0:MLS) :: TAB
      ! Array for the radial solution
      INTEGER :: IL
      DOUBLE PRECISION :: A1,A2,A4,A5
      DOUBLE PRECISION,DIMENSION(NR) :: DC1
      DOUBLE PRECISION,DIMENSION(NR-1) :: DL1,DU1
      INTEGER :: step

      IF(planar_symmetry) THEN
         step= 2
      ELSE
         step= 1
      ENDIF

      DO IZ=1,NR     
         ! Evaluate the coefficient of the Legendre poly. expansion
         DO I=0,MLS,step
            ! Compute the matrix of legendre coefficients
            TAB(IZ,i)=2d0*pi*DOT_PRODUCT(pn_init(I,1:nth),source(1:nth,iz))
         END DO
      END DO
     
      DO IL=0,MLS,step
         dc1(:) = dcp(:) - DBLE(il*(il+1)) / r(1:nr)**2

         ! BC at inner radius ( Origin )
         A1= 2./DR(1)/(DR(1)+DR(1+1))
         A4=-DR(1+1)/DR(1)/(DR(1)+DR(1+1))
         IF(MOD(il,2) == 0) THEN
            DC1(1)=dc1(1)+(A1+2.*A4/R(1))
         ELSE
            DC1(1)=dc1(1)-(A1+2.*A4/R(1))
         ENDIF
         
         ! BC at outer radius (multipoles must dacay as R^(L+1))
         A2= 2./DR(NR+1)/(DR(NR)+DR(NR+1))
         A5= DR(NR)/DR(NR+1)/(DR(NR)+DR(NR+1))
         dc1(nr)=dc1(nr)+(a2+2.*a5/r(nr))*(r(nr)/(r(nr)+dr(nr+1)))**(il+1)
        
         ! solve the system
         DL1=DLP
         DU1=DUP
         CALL dgtsv(NR,1,DL1,DC1,DU1,TAB(1:NR,IL),NR,INFO)
      ENDDO

      ! Compute the matrix of Legendre Polyn in the grid points
      DO IZ=1,NR
         DO IX=1,NTH
            fun(IX,IZ)=DOT_PRODUCT(pnx_init(0:mls:step,ix),TAB(IZ,0:MLS:step))
         ENDDO
      ENDDO
   END SUBROUTINE scalar_laplacian


   !> solve the phi component of the vector Laplacian equation with a series expansion
   !! in Legendre polynomials in the angle and by matrix inversion in the radius
   !! @details
   !!    assumptions:
   !!    * axisymmetry
   !!    * circularity
   !!    * input/output in coordinate basis
   !! @note
   !!    the subroutine actually works in orthonormal basis, but it makes the
   !!    conversions coordinate basis --> orthonormal basis automatically
   !! @param[in] source source of the vector Laplacian equation
   !! @param[out] btp phi component of the shift, that is, minus ZAMO angular velocity
   PURE SUBROUTINE SHIFTPHI(source,btp)
      !USE system, ONLY: MLS,NTH,NR,PI,R,DR,pd_init,pdx_init,th,dcp,dup,dlp,planar_symmetry
      IMPLICIT NONE
      INTEGER :: I,INFO
      DOUBLE PRECISION,DIMENSION(NTH,NR), INTENT(IN):: source
      DOUBLE PRECISION,DIMENSION(NTH,NR), INTENT(OUT):: btp
      INTEGER :: ix,iz
      ! Array for the Legendre expansion
      DOUBLE PRECISION,DIMENSION(NR,0:MLS) :: TAB
      ! Array for the radial solution
      INTEGER :: IL
      DOUBLE PRECISION :: A1,A2,A4,A5
      DOUBLE PRECISION,DIMENSION(NR) :: DC1
      DOUBLE PRECISION,DIMENSION(NR-1) :: DL1,DU1
      INTEGER :: step

      IF(planar_symmetry) THEN
         step= 2
      ELSE
         step= 1
      ENDIF

      ! Initialize the source in the domain
      DO IZ=1,NR
         DO IX=1,NTH
            ! from coordinate basis --> orthornormal basis (gio:2018/03/26)
            ! BEWARE! to save memory I am using btp as temp storage of the source!
            btp(IX,IZ) = source(IX,IZ)*R(IZ)*SIN(TH(IX)) ! here btp = source
         ENDDO
      ENDDO   

      DO IZ=1,NR
         ! Evaluate the coefficient of the Spheric Harm. expansion (L=0 => 0)
         DO I=1,MLS,step
            TAB(IZ,i)=2d0*pi/DBLE(I*(I+1))*DOT_PRODUCT(pd_init(I,1:nth),btp(1:nth,iz)) ! here btp = source
         ENDDO
      ENDDO
      ! BEWARE: in the next use btp will be re-assigned to the real output profile!

      ! Parity imposes that there is not L=0 term (we work on coord quantity)
      DO IL=1,MLS,step

         dc1(:) = dcp(:) - DBLE(il*(il+1)) / r(1:nr)**2

         ! BC at inner radius (origin)
         A1= 2./DR(1)/(DR(1)+DR(1+1))
         A4=-DR(1+1)/DR(1)/(DR(1)+DR(1+1))
         IF(MOD(il+1,2) == 0) THEN
            DC1(1)=dc1(1)-(A1+2.*A4/R(1))
         ELSE
            DC1(1)=dc1(1)+(A1+2.*A4/R(1))
         ENDIF
         ! note that the source is a vector field therefore the parity is different from the scalar case
         
         ! BC at outer radius: multipoles must decay as R^(L+1)
         A2= 2./DR(NR+1)/(DR(NR)+DR(NR+1))
         A5= DR(NR)/DR(NR+1)/(DR(NR)+DR(NR+1))
         dc1(nr)=dc1(nr)+(A2+2.*A5/R(NR))*(R(NR)/(R(NR)+DR(NR+1)))**(IL+1)
        
         DL1=DLP
         DU1=DUP
         CALL DGTSV(NR,1,DL1,DC1,DU1,TAB(1:NR,IL),NR,INFO)
      END DO

      ! Compute the matrix of Legendre Polyn in the grid points
      DO IZ=1,NR
         DO IX=1,NTH
            btp(IX,IZ)=DOT_PRODUCT(pdx_init(1:mls:step,ix),tab(IZ,1:MLS:step))
            ! orthonormal basis --> coordinate basis
            ! one avoids to divide by sin(theta) noting that
            ! the associated spherical harmonics already contain sin(theta)
            btp(ix,iz)=btp(ix,iz) / r(iz)
         ENDDO
      ENDDO
   END SUBROUTINE  SHIFTPHI


   !> differentiate a quantity along r and th
   !! @details
   !!    assumptions:
   !!    * even parity at the poles,
   !!    * axisymmetry,
   !!    * smoothness at the outer radial boundary
   !! @param[in] fun quantity to derive
   !! @param[out] der_r differentiation along the radial coordinate
   !! @param[out] der_th differentiation along the angular coordinate
   PURE SUBROUTINE derive(fun,der_r,der_th)
      !USE system, ONLY: a1,a2,a3,b1,b2,b3,nr,nth,dr,planar_symmetry
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(nth,nr), INTENT(IN) :: fun
      DOUBLE PRECISION, DIMENSION(nth,nr), INTENT(OUT) :: der_r, der_th
      INTEGER :: i

      ! ************************ !
      ! *** radial direction *** !
      ! ************************ !

      ! use axisymmetry at the center (see a1)
      IF(planar_symmetry) THEN
         der_r(:,1) = A1(1)*fun(:,1) + A3(1)*fun(:,2) + A2(1)*fun(:,1)
      ELSE
         der_r(1:nth,1) = A1(1)*fun(nth:1:-1,1) + A3(1)*fun(:,2) + A2(1)*fun(:,1)
      ENDIF

      ! internal points
      DO I=2,NR-1
         der_r(:,i) = A1(i)*fun(:,I-1) + A3(i)*fun(:,I+1) + A2(i)*fun(:,I)
      ENDDO

      ! use smoothness at the outer boundary (see a3)
      der_r(:,nr) = A1(nr)*fun(:,nr-1) + A3(nr)*&
                  & (  fun(:,nr) + DR(NR+1)/DR(NR)*(fun(:,NR)-fun(:,NR-1))  ) &
                  & + A2(nr)*fun(:,nr)

      ! ************************* !
      ! *** angular direction *** !
      ! ************************* !

      ! even at the north pole (see b1)
      der_th(1,:) = B1(1)*fun(1,:) + B3(1)*fun(2,:) + B2(1)*fun(1,:)

      ! internal points
      DO I=2,NTH-1
         der_th(i,:) = B1(i)*fun(I-1,:) + B3(i)*fun(I+1,:) + B2(i)*fun(I,:)
      ENDDO

      IF(planar_symmetry) THEN
         ! even around the equator
         der_th(nth,:) = B1(nth)*fun(nth-1,:) + B3(nth)*fun(nth-1,:) + B2(nth)*fun(nth,:)
      ELSE
         ! even at the south pole (see b3)
         der_th(nth,:) = B1(nth)*fun(nth-1,:) + B3(nth)*fun(nth,:) + B2(nth)*fun(nth,:)
      ENDIF
   END SUBROUTINE derive


   !> TOV solver in isotropic coordinates
   !! @copyright Giovanni Camelio 2019/06/26 Stockholm university
   !! @details
   !!    for theory, look at:
   !!    * Gourgoulhon 2011, arxiv:1003.5015v2, Sec. 3.3
   !!    * website (as today, 2019/06/26):
   !!      ion.uwinnipeg.ca/~vincent/4500.6-001/Cosmology/IsotropicCoordinates.htm
   !!    * my notes
   !! @note
   !!    inside the star you cannot directly obtain the Schwartzchild
   !!    radius from the isotropic radius as in references above but you have to integrate it.
   !!    The integration is particularly sensitive to the initial Schw. radius.
   !!    For this reason you need to find the initial Schw. radius with a Newton-Rapson.
   !!    This is similar to the TOV initializer of the original XNS-v2.
   !! @param[in] p0 central pressure
   !! @param[in] p2hden barotropic EOS
   !! @param[out] p pressure profile
   !! @param[out] hden enthalpy density profile
   !! @param[out] nu gravitational potential profile
   !! @param[out] psi profile of the conformal factor
   !! @param[out] surf position of the surface (the surface is between surf and surf+1)
   !! @param[out] correct_init success flag (true on success)
   PURE SUBROUTINE tov(p0,p2hden,p,hden,nu,psi,surf,correct_init)
      !USE system, ONLY: nr, r, pi, verbose
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: p0 ! central pressure
      DOUBLE PRECISION, DIMENSION(nr), INTENT(OUT) :: p, hden, nu, psi
      PROCEDURE(fun_hden) :: p2hden
      INTEGER, INTENT(OUT) :: surf
      LOGICAL, INTENT(OUT) :: correct_init
      CHARACTER(LEN=255) :: string
      INTEGER, PARAMETER :: maxit= 20 ! original: 20
      DOUBLE PRECISION, PARAMETER :: tol= 1.5d-15 ! original: 1.5d-15
      DOUBLE PRECISION :: y(4),dy(4),ynew(4),m,rsch,x,f,fold,xold,xnew
      INTEGER :: i,j

      ! Newton-Raphson on the initial Schwarzchild radius
      x= 1d0
      DO j=1, maxit
         ! *very* approximate initial conditions for schwartzchild radius
         y(1)= r(1)*x ! schwartzchild radius
         y(2)= 0d0    ! gravitational mass
         y(3)= 0d0    ! log(lapse)
         y(4)= p0     ! pressure

         ! stellar interior
         DO i=1, nr-1
            rsch= y(1)                ! Schwartzchild radius
            psi(i)= DSQRT(y(1)/r(i))  ! conformal factor
            m= y(2)                   ! gravitational mass
            nu(i)= y(3)               ! log(lapse)
            p(i)= y(4)                ! pressure
            hden(i)= p2hden(p(i))     ! enthalpy density

            ! derivatives with respect to the isotropic radius
            dy(1)= DSQRT(rsch*rsch - 2d0*m*rsch)/r(i) ! drsch/dr
            dy(2)= 4d0*pi*rsch*rsch*(hden(i) - p(i))*dy(1) ! dm/dr
            dy(3)= (m + 4d0*pi*rsch*rsch*rsch*p(i))/(rsch*rsch - 2d0*m*rsch)*dy(1) ! dnu/dr
            dy(4)= -hden(i)*dy(3) ! dp/dr

            ! a simple Euler integrator is enough
            ynew(:)= y(:) + (r(i+1) - r(i))*dy(:)
            ! if the pressure is not positive, you reached the stellar surface
            IF(ynew(4) <= 0d0) EXIT
            y= ynew
         ENDDO
         
         ! update quantities for Newton-Rapson
         f= m - 2d0*(psi(i) - 1d0)*r(i)
         IF(j == 1) THEN
            xnew= 2d0*x
         ELSE
            xnew= x - f*(x - xold)/(f - fold)
         ENDIF
         IF(verbose) THEN
            WRITE(string,*)"TOV: cycle=",j,"surf=",i,"M=",m,"x=",x,"dx=",xnew-x
            CALL print_on_screen(string)
         ENDIF

         ! check result
         IF(DABS(xnew - x) < tol) EXIT

         xold= x
         fold= f
         x= xnew
      ENDDO

      ! check integration
      IF(i >= nr) THEN ! original: i + 10 > nr/2 - 1, later i > nr/2
         IF(verbose) THEN
            WRITE(string,*)"ERROR in TOV: neutron star too large, surf index=",i
            CALL print_on_screen(string)
         ENDIF
         correct_init= .false.
      ELSE
         IF(i > nr/2 .AND. verbose) THEN
            WRITE(string,*)"BEWARE in TOV: neutron star too large, surf index=",i
            CALL print_on_screen(string)
         ENDIF
         IF(j == maxit + 1) THEN
            WRITE(string,*)"BEWARE in TOV: no convergence"
            CALL print_on_screen(string)
            ! I do not set correct_init to false because this is anyway a tentative initial condition
         ENDIF
         correct_init= .true.
      ENDIF

      ! surface is between i and i + 1
      surf= i

      ! adjust the lapse in the stellar interior
      nu(:surf)= nu(:surf) - nu(surf) + DLOG((1d0 - 0.5d0*m/r(surf))/(1d0 + 0.5d0*m/r(surf))) 
      
      ! stellar exterior
      p(surf+1:nr)= 0d0
      hden(surf+1:nr)= 0d0
      nu(surf+1:nr)= DLOG((1d0 - 0.5d0*m/r(surf+1:nr))/(1d0 + 0.5d0*m/r(surf+1:nr)))
      psi(surf+1:nr)= 1d0 + 0.5d0*m/r(surf+1:nr)
   END SUBROUTINE tov

   !> compute the Legendre knots and weights
   !! @details
   !!    references:
   !!    * math.stackexchange.com/questions/12160/roots-of-legendre-polynomial
   !!    * Francesco Tricomi (1950), https://doi.org/10.1007/BF02428258
   !! @param[in] n number of the knots
   !! @param[in] z position of the knots
   !! @param[in] w weight of the knots
   !! @param[out] info success flag (0 on success)
   PURE SUBROUTINE legendre_knots(n,z,w,info)
      !USE system, ONLY: pi
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      DOUBLE PRECISION, INTENT(OUT) :: z(n), w(n)
      INTEGER, INTENT(OUT) :: info
      INTEGER, PARAMETER :: maxit= 100
      DOUBLE PRECISION, PARAMETER :: tol= 1d-15 ! tol=0 ==> machine precision
      DOUBLE PRECISION :: p0, pn, pd, dk, dn, x, tmp
      INTEGER :: i,j,k

      IF(n < 1) THEN
         info= -1
         RETURN
      ENDIF

      info= 0
      dn= DBLE(n)
      DO i=1,(n + 1)/2
         ! first guess on the i-th root
         ! this approximation due to Francesco Tricomi (1950) is O(n^-4),
         ! cf. Eqs. (2) & (13) of https://doi.org/10.1007/BF02428258
         x= (1d0 - 0.125d0/(dn*dn)*(1d0 - 1d0/dn))*DCOS(pi*DBLE(4*i - 1)/DBLE(4*n + 2))
         ! now find the i-th zero with Newton-Rapson
         DO j=1, maxit
            ! compute Pn(x) and Pn'(x)
            p0= 1d0
            pn= x
            pd= 1d0
            DO k=2,n
               dk= DBLE(k)
               tmp= ((2d0*dk - 1d0)*x*pn - (dk - 1d0)*p0)/dk
               pd= dk*pn + x*pd ! <== valid also for |x| == 1
               ! pd= dn*(x*pn - p0)/(x*x - 1d0) ! <== valid if |x| /= 1
               p0= pn
               pn= tmp
            ENDDO
            IF(pd == 0d0 .OR. x == 0d0) EXIT
            ! Newton-Raphson
            tmp= x - pn/pd
            IF(DABS(x - tmp) <= DABS(x)*tol) EXIT
            x= tmp
         ENDDO
         IF(j == maxit+1) info= -2
         z(i)= x
         w(i)= 2d0/((1 - x*x)*pd*pd)
         ! symmetry
         z(n+1-i)= -x
         w(n+1-i)= w(i)
      ENDDO

      ! check result
      tmp= 1d0
      DO i=1, n
         ! all roots found ==> strictly descendent order
         IF(z(i) >= tmp) info= -3
         tmp= z(i)
         ! pd == 0 somewhere in the procedure ==> NaN
         IF(w(i) /= w(i)) info= -4
      ENDDO
   END SUBROUTINE legendre_knots

END MODULE system


!> for printing on screen from a pure function
!! @details
!!    this is a dirty trick in order to have side effects from pure functions.
!!    for this reason it is outside the module
!! @param[in] text string to be printed on screen
SUBROUTINE print_on_screen(text)
   IMPLICIT NONE
   CHARACTER(LEN=*), INTENT(IN) :: text

   WRITE(6,*) TRIM(text)
END SUBROUTINE print_on_screen
