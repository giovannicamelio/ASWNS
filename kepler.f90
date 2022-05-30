!> @file
!! @brief example of usage of the ASWNS code
!! @details
!!    * EOS: ɛ = ρ + k1·ρ^Γ
!!    * cold rigidly rotating star
!!    * model: Q = Q0 + Fcal + H + b·H·Fcal
!!    * if the user want to change the EOS, s/he should modify only this file
!!    * code units: G = c = M☉ = kB = 1
!! @author Giovanni Camelio
!! @date 2022
!! @copyright Giovanni Camelio.
!!            This code is released under the CC BY-NC-SA 4.0 license.


!> example of user module for dealing with the EOS and the parametric model
!! @details it contains the parameters of the model and the EOS,
!!          and the functions to solve the EOS and the force balance equation (Euler)
MODULE kepler
   USE system, ONLY: rhon, km_to_msol
   IMPLICIT NONE

   ! *********************** !
   ! *** OUTPUT SETTINGS *** !
   ! *********************** !

   !> wether to print on screen a brief summary of the results or not
   LOGICAL, PARAMETER :: verbose= .true.
   !> starting density of the Keplerian the curve
   DOUBLE PRECISION, PARAMETER :: rhostart= 4d0*rhon
   !> final density of the Keplerian the curve
   DOUBLE PRECISION, PARAMETER :: rhoend= 7d0*rhon
   !> file where to print the Keplerian curve
   CHARACTER(*), PARAMETER :: filename= "kepler.dat"

   ! ******************** !
   ! *** EOS SETTINGS *** !
   ! ******************** !

   !> EOS exponent of the cold component
   DOUBLE PRECISION, PARAMETER :: gam= 3d0
   !> EOS proportionality constant of the cold component
   DOUBLE PRECISION, PARAMETER :: k1= 5d4

   ! ********************** !
   ! *** MODEL SETTINGS *** !
   ! ********************** !

   !> angular velocity in the center (will be set in the main)
   DOUBLE PRECISION :: omg0= 0d0
   !> central rest mass density (will be set in the main)
   DOUBLE PRECISION :: rho0= rhostart

   ! ********************** !
   ! *** OTHER SETTINGS *** !
   ! ********************** !

   !> relaxation iterations in the force balance equation (Euler) solver
   !! @details used in subroutine @ref solve_euler, must be 0 <= relax_iters < @ref wait_iters
   INTEGER, PARAMETER :: relax_iters= 0

   SAVE

CONTAINS

   !> with this subroutine the user can initialize the EOS and the barotropic model
   !! @details
   !!    * it should be called before the EOS is needed
   !!    * at the moment it just check that the EOS and the model are physical
   !! @param[out] info success flag (success = 0)
   SUBROUTINE user_initialization(info)
      USE system, ONLY: wait_iters
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: info

      info= -3

      ! check settings
      IF(k1 <= 0d0) THEN
         PRINT*,"ERROR: must be k1 > 0"
         RETURN
      ENDIF
      IF(gam <= 1d0) THEN
         PRINT*,"ERROR: must be gam > 1"
         RETURN
      ENDIF
      IF(relax_iters < 0 .OR. relax_iters >= wait_iters) THEN
         PRINT*,"ERROR: must be 0 <=  relax_iters < wait_iters"
         RETURN
      ENDIF

      ! initialize the EOS

      ! initialize the barotropic profile
      
      ! other user initializations
      
      info= 0

   END SUBROUTINE user_initialization


   !> print the stellar results and profiles.
   !!
   !! @details
   !!    * if @ref verbose is true, print on screen a short summary of the results.
   !!    * if @ref logfile is empty, print on logfile an extended summary of the results.
   !!    * if @ref binfile is empty, print on binfile the stellar profiles in binary format.
   !!
   !! @param[in] a stellar profiles
   !! @param[in] fun function to solve the EOS given the pressure and the enthalpy density
   !! @param[out] m0 total stellar baryonic mass
   !! @param[out] stot total stellar entropy
   !! @param[out] m0disk total baryonic mass of the disk
   SUBROUTINE print_quantities(a,fun,m0,stot,m0disk)
      USE system, system_verbose => verbose
      IMPLICIT NONE
      TYPE(model_quantities), INTENT(IN) :: a
      PROCEDURE(eos) :: fun
      DOUBLE PRECISION, INTENT(OUT) :: m0,stot,m0disk
      INTEGER :: ireq,irpol

      ireq= a%wsurf(a%ieq) ! index of the stellar border at the equator -- surface is between ireq and ireq+1
      irpol= a%wsurf(1)  ! index of the stellar border at the pole -- surface is between irpol and irpol+1
      CALL get_m0_stot(a,fun,m0,stot,m0disk)

      ! print on screen
      IF(verbose) THEN
         WRITE(6,*)""
         WRITE(6,*)'gravitivational mass = ',a%Mg
         WRITE(6,*)"rest mass            = ",m0
         WRITE(6,*)"total entropy        = ",stot
         WRITE(6,*)'angular momentum     = ',a%J
         WRITE(6,*)'equatorial radius    = ',r(ireq)," = ",r(ireq)/km_to_msol,' km'
         WRITE(6,*)'radii ratio          = ',r(irpol)/r(ireq)
         WRITE(6,*)'equatorial Omega     = ',a%eq_omega
         WRITE(6,*)"keplerian Omega      = ",a%kep_omega
         WRITE(6,*)"disk mass            = ",m0disk
         WRITE(6,*)"central rho          = ",rho0
         WRITE(6,*)"axial Omega          = ",omg0
         WRITE(6,*)"error flag           = ",a%info
         WRITE(6,*)""
      ENDIF
   END SUBROUTINE print_quantities


   !> enthalpy density from pressure for the cold EOS
   !! @details this function is passed to the TOV solver in order to initialize the star
   !! @param[in] p pressure
   !! @result enthalpy density
   PURE DOUBLE PRECISION FUNCTION p2hden(p) RESULT(hden)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: p
      DOUBLE PRECISION :: tmp

      tmp= p/((gam - 1d0)*k1)
      hden= tmp**(1d0/gam) + k1*tmp + p
   END FUNCTION p2hden


   !> EOS from pressure and enthalpy density
   !! @details this subroutine is passed to print_quantities and get_m0_stot
   !!    to compute additional quantities non necessary for solving the stellar structure
   !! @param[in] p pressure
   !! @param[in] hden enthalpy density
   !! @param[out] rho rest mass density
   !! @param[out] s entropy
   !! @param[out] t temperature
   PURE SUBROUTINE p_hden2eos(p,hden,rho,s,t)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: hden,p
      DOUBLE PRECISION, INTENT(OUT) :: rho,s,t

      rho= (p / k1 / (gam - 1d0))**(1d0/gam)
      s= 0d0
      t= 0d0
   END SUBROUTINE p_hden2eos


   !> EOS from rest mass density
   !! @param[in] rho rest mass density
   !! @result p pressure
   PURE DOUBLE PRECISION FUNCTION rho2p(rho) RESULT(p)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: rho
      
      p= k1*(gam - 1d0)*rho**gam
   END FUNCTION rho2p


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
   PURE SUBROUTINE solve_euler(iter,q0,alp2,btp,rcyl2,p,omg,hden,info)
      ! uses maxit and tol from user module
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iter
      DOUBLE PRECISION, INTENT(IN) :: q0,alp2,btp,rcyl2
      DOUBLE PRECISION, INTENT(INOUT) :: p,omg
      DOUBLE PRECISION, INTENT(OUT) :: hden
      INTEGER, INTENT(OUT) :: info
      DOUBLE PRECISION :: heat,tmp1,tmp2,h,h0,rho

      info= -1

      ! the target central omg is slowly increased at each iteration to improve the convergence
      omg= omg0*MIN(1d0, DBLE(iter)/DBLE(relax_iters + 1))
      ! compute Q(Omega,r,theta)
      tmp1= omg + btp
      tmp2= alp2 - rcyl2*tmp1*tmp1
      ! check velocity < c
      IF(tmp2 <= 0d0) RETURN

      heat= -0.5d0*DLOG(tmp2) - q0
      h0= 1d0 + gam*k1*rho0**(gam - 1d0)
      h= h0*DEXP(heat)
      IF(h <= 1d0) THEN
         p= 0d0
         hden= 0d0
      ELSE
         rho= ((h - 1d0)/(gam*k1))**(1d0/(gam - 1d0))
         p= (gam - 1d0)*k1*rho**gam
         hden= rho*h
      ENDIF

      info= 0
   END SUBROUTINE solve_euler

END MODULE kepler


!> find the Keplerian configuration corresponding to a given central density
!!
!! @details it assumes that the first configuration is already computed
!!
!! @param[in] increment step increment in the angular velocity
!! @param[inout] a stellar model
SUBROUTINE find_kep(increment,a)
   USE system, ONLY: model_quantities, ns_solve
   USE kepler, ONLY: print_quantities, p_hden2eos, user_initialization, solve_euler, verbose, omg0
   IMPLICIT NONE
   DOUBLE PRECISION, INTENT(IN) :: increment
   TYPE(model_quantities), INTENT(INOUT) :: a
   INTEGER, PARAMETER :: maxiterations= 100
   TYPE(model_quantities) :: b
   DOUBLE PRECISION :: bomg0,m0,stot,m0disk
   INTEGER :: i, info

   ! save the first configuration
   b= a
   bomg0= omg0

   DO i=1, maxiterations
      ! choose the new parameters
      omg0= omg0 + increment
      ! I have to load the barotropic thermal profile
      CALL user_initialization(info)

      ! compute the configuration and the stellar quantities
      CALL ns_solve(a,solve_euler,info)

      ! logging
      IF(verbose) THEN
         PRINT*, ""
         PRINT*, "## find_kep: iteration= ", i
         CALL print_quantities(a,p_hden2eos,m0,stot,m0disk)
      ENDIF

      ! if the new configuration is above the keplerian limit OR it is wrong,
      ! return the older (almost keplerian) configuration
      IF(info /= 0) THEN
         a= b
         omg0= bomg0
         EXIT
      ENDIF
      ! now the configuration is under the keplerian limit, AND it is correct

      ! save the new conf
      b= a
      bomg0= omg0
   ENDDO
END SUBROUTINE find_kep


!> keplerian curve for a cold, rigidly rotating neutron star
PROGRAM keplerian_curve
   USE system, ONLY: model_quantities, r, ns_init, ns_solve
   USE kepler, ONLY: verbose, p2hden, solve_euler, p_hden2eos, &
      & user_initialization, print_quantities, filename, rhostart, rhoend, rho2p, omg0, rho0
   IMPLICIT NONE
   ! maximal number of steps in the keplerian curve
   INTEGER, PARAMETER :: imax= 1000
   ! differentiating step for the angular velocity
   DOUBLE PRECISION, PARAMETER :: deltaomg= 1d-4
   ! differentiating step for the central rest mass density
   ! and distance between two densities in the curve
   DOUBLE PRECISION, PARAMETER :: deltarho= 1d-5
   DOUBLE PRECISION :: m0,stot,m0disk,dmdrho,bj,bm,rcirc,aomg0
   TYPE(model_quantities) :: a,b
   INTEGER :: i,ireq,info

   ! init stellar parameters
   rho0= rhostart
   CALL user_initialization(info)
   IF(info /= 0) RETURN

   ! initialize the non rotating star
   CALL ns_init(a,p2hden,rho2p(rho0),info)
   IF(info /= 0) RETURN

   ! relax the non rotating configuration
   omg0= 0d0
   CALL ns_solve(a,solve_euler,info)
   IF(info /= 0) RETURN

   OPEN(FILE=filename,UNIT=8)

   ! find the keplerian line
   DO i=1,imax
      IF(verbose) PRINT*,"# keplerian_curve: iteration= ", i

      ! find the Keplerian configuration
      IF(i == 1) CALL find_kep(1d-2,a)
      CALL find_kep(1d-3,a)
      CALL find_kep(1d-4,a)
      ! in Camelio+2018, the minimum increment was 1d-4
      CALL find_kep(1d-5,a)
      aomg0= omg0

      ! differentiate w/ to Omega0
      b= a
      omg0= aomg0 - deltaomg
      CALL ns_solve(b,solve_euler,info)
      IF(info /= 0) EXIT
      bj= b%j
      bm= b%mg

      ! differentiate w/ to rho0
      b= a
      omg0= aomg0
      rho0= rho0 + deltarho
      CALL ns_solve(b,solve_euler,info)
      IF(info /= 0) EXIT

      ! stability Camelio+(2018):Eq:14
      dmdrho= (b%mg - a%mg - (b%j - a%j)*(a%mg - bm)/(a%j - bj))/deltarho

      ! print result
      CALL print_quantities(a,p_hden2eos,m0,stot,m0disk)
      ireq= a%wsurf(a%ieq)
      rcirc= r(ireq)*a%psi(a%ieq,ireq)**2 ! circumferential radius
      WRITE(8,*)rho0,omg0,a%info,a%mg,m0,a%mp,stot,a%j,a%tr,a%eq_omega,a%kep_omega, &
               & r(ireq),rcirc,r(a%wsurf(1)),m0disk,dmdrho

      ! increase central density
      IF(rho0 > rhoend) EXIT
      a= b
   ENDDO

   CLOSE(8)
END PROGRAM keplerian_curve
