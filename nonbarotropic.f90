!> @file
!! @brief example of usage of the ASWNS code
!! @details
!!    * EOS: ɛ = ρ + k1·ρ^Γ + k2·s²·ρ^Γth
!!    * thermal barotropic law: s = k3·ρ^[(Γ - Γth)/2]
!!    * rotational barotropic law: F = -R0²(Ω - Ω0)/2
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
MODULE user
   USE system, ONLY: rhon, km_to_msol
   IMPLICIT NONE

   ! *********************** !
   ! *** OUTPUT SETTINGS *** !
   ! *********************** !

   !> wether to print on screen a brief summary of the results or not
   LOGICAL, PARAMETER :: verbose= .true.
   !> file where to print the log of the parameters and the extended summary of the results
   !! @details if empty, does not print the log
   CHARACTER(*), PARAMETER :: logfile= "star.log"
   !> file where to print the stellar profile in binary format
   !! @details if empty, does not print the profiles
   CHARACTER(*), PARAMETER :: binfile= "star.out"

   ! ******************** !
   ! *** EOS SETTINGS *** !
   ! ******************** !

   !> EOS exponent of the cold component
   DOUBLE PRECISION, PARAMETER :: gam= 3d0
   !> EOS proportionality constant of the cold component
   DOUBLE PRECISION, PARAMETER :: k1= 5d4
   !> EOS exponent of the thermal component
   DOUBLE PRECISION, PARAMETER :: gamth= 1.75d0
   !> EOS proportionality constant of the thermal component
   DOUBLE PRECISION, PARAMETER :: k2= 1.5d0

   ! ********************** !
   ! *** MODEL SETTINGS *** !
   ! ********************** !

   !> angular velocity in the center
   DOUBLE PRECISION, PARAMETER :: omg0= 0.035d0
   !> baroclinic parameter
   DOUBLE PRECISION, PARAMETER :: bvalue= -2d0
   !> central rest mass density
   DOUBLE PRECISION, PARAMETER :: rho0= 4d0*rhon
   !> entropy proportionality constant of the barotropic law
   DOUBLE PRECISION, PARAMETER :: k3= 2d0/rho0**(0.5d0*(gam - gamth))
   !> inverse scale radius of the differential rotation sigma = 1/R0
   !! @details if sigma <= 0, the star is in rigid rotation
   DOUBLE PRECISION, PARAMETER :: sigma= 1d0/(15d0*km_to_msol)

   ! ********************** !
   ! *** OTHER SETTINGS *** !
   ! ********************** !

   !> max iterations for the Newton-Raphson
   !! @details used in subroutines @ref newton_raphson and @ref solve_euler, must be maxit > 0
   INTEGER, PARAMETER :: maxit= 1000
   !> relative tolerance for the Newton-Raphson
   !! @details used in subroutines @ref newton_raphson and @ref solve_euler, must be tol > 0d0
   DOUBLE PRECISION, PARAMETER :: tol= 1d-15
   !> relaxation iterations in the force balance equation (Euler) solver
   !! @details used in subroutine @ref solve_euler, must be 0 <= relax_iters < @ref wait_iters
   INTEGER, PARAMETER :: relax_iters= 0

   ! *************** !
   ! *** GLOBALS *** !
   ! *************** !

   !> critical density for the EOS inversion from (p, hden) to (rho, s)
   !> @details will be set in @ref user_initialization
   DOUBLE PRECISION :: rhocrit= 0d0
   !> maximal value of the function needed for the EOS inversion from (p, hden) to (rho, s)
   !> @details will be set in @ref user_initialization
   DOUBLE PRECISION :: funmax= HUGE(1d0)

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
      IF(k2 <= 0d0) THEN
         PRINT*,"ERROR: must be k2 > 0"
         RETURN
      ENDIF
      IF(k3 < 0d0) THEN
         PRINT*, "ERROR: must be k3 >= 0"
         RETURN
      ENDIF
      IF(gam <= 1d0) THEN
         PRINT*,"ERROR: must be gam > 1"
         RETURN
      ENDIF
      IF(gamth <= 1d0) THEN
         PRINT*,"ERROR: must be gamth > 1"
         RETURN
      ENDIF
      IF(relax_iters < 0 .OR. relax_iters >= wait_iters) THEN
         PRINT*,"ERROR: must be 0 <=  relax_iters < wait_iters"
         RETURN
      ENDIF
      IF(maxit <= 0) THEN
         PRINT*,"ERROR: must be maxit > 0"
         RETURN
      ENDIF
      IF(tol <= 0d0) THEN
         PRINT*,"ERROR: must be tol > 0"
         RETURN
      ENDIF

      ! initialize the EOS

      IF(gam > gamth) THEN
         ! I need to use rhocrit as tmp to avoid (stupid) compilation errors (with gfortran)
         rhocrit= gam - gamth
         rhocrit= ((gamth - 1d0)/(rhocrit*k1*gam))**(1d0/(gam - 1d0))
         funmax= (gamth - 1d0)*rhocrit + (gamth - gam)*k1*rhocrit**gam
         IF(rho0 > rhocrit) PRINT*,"WARNING: rho0 > rhocrit, you cannot invert from (p, hden) to (rho, s)"
      ENDIF

      ! initialize the barotropic profile
      
      ! other user initializations
      
      info= 0

   END SUBROUTINE user_initialization


   !> print the stellar results and profiles.
   !! @details
   !!    * if @ref verbose is true, print on screen a short summary of the results.
   !!    * if @ref logfile is empty, print on logfile an extended summary of the results.
   !!    * if @ref binfile is empty, print on binfile the stellar profiles in binary format.
   !! @param[in] a stellar profiles
   !! @param[in] fun function to solve the EOS given the pressure and the enthalpy density
   SUBROUTINE print_quantities(a,fun)
      USE system
      IMPLICIT NONE
      TYPE(model_quantities), INTENT(IN) :: a
      PROCEDURE(eos) :: fun
      CHARACTER(len=8) :: date
      CHARACTER(len=10) :: time
      INTEGER :: ix,iz,reclen,mr,ireq,irpol
      DOUBLE PRECISION :: rho(nth,nr),s(nth,nr),t(nth,nr),rcirc,m0,stot,m0disk

      ireq= a%wsurf(a%ieq) ! index of the stellar border at the equator -- surface is between ireq and ireq+1
      irpol= a%wsurf(1)  ! index of the stellar border at the pole -- surface is between irpol and irpol+1
      rcirc= r(ireq)*a%psi(a%ieq,ireq)**2 ! circumferential radius
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
         WRITE(6,*)"central Omega        = ",omg0
         WRITE(6,*)"error flag           = ",a%info
         WRITE(6,*)""
      ENDIF
      
      ! print the log
      IF(logfile/="") THEN
         ! Write Log file with Run Parameters
         OPEN(8,FILE=logfile)
         WRITE(8,*)"ASWNS log"
         WRITE(8,*)"version: ASWNS-v1"
         CALL DATE_AND_TIME(date,time)
         WRITE(8,*)"date: ",date
         WRITE(8,*)"time: ",time
         WRITE(8,*)"c = G = Msol = kB = 1"
         WRITE(8,*)""
         WRITE(8,*)"ASWNS SETTINGS:"
         WRITE(8,*)' nr   = ', NR   
         WRITE(8,*)' rhalf= ', rhalf
         WRITE(8,*)' rmax = ', RMAX
         WRITE(8,*)' nth  = ', NTH   
         WRITE(8,*)' mls  = ', MLS
         WRITE(8,*)" planar_symmetry = ",planar_symmetry
         WRITE(8,*)' qfactor  = ',QFACTOR
         WRITE(8,*)' abs. tol.= ',abs_tol
         WRITE(8,*)' wait iter=', wait_iters
         WRITE(8,*)""
         WRITE(8,*)"EOS SETTINGS:"
         WRITE(8,*)" k1       = ",k1
         WRITE(8,*)" K        = ",k1*(gam - 1d0)
         WRITE(8,*)" Gamma    = ",gam
         WRITE(8,*)" kth      = ",k2
         WRITE(8,*)" Gamma_th = ",gamth
         WRITE(8,*)""
         WRITE(8,*)"MODEL SETTINGS:"
         WRITE(8,*)" central rho             = ",rho0
         WRITE(8,*)" central Omega           = ",omg0
         WRITE(8,*)" inverse radial scale    = ",sigma
         WRITE(8,*)" thermal law scale       = ",k3
         WRITE(8,*)" baroclinic constant     = ",bvalue
         WRITE(8,*)''
         WRITE(8,*)'OUTPUT PARAMETERS'
         WRITE(8,*)' gravitational mass     = ',a%mg
         WRITE(8,*)" rest mass              = ",m0
         WRITE(8,*)" total entropy          = ",stot
         WRITE(8,*)' proper mass            = ',a%mp
         WRITE(8,*)' rotational energy      = ',a%tr
         WRITE(8,*)' angular momentum       = ',a%j
         WRITE(8,*)' equatorial radius      = ',r(ireq)
         WRITE(8,*)' polar radius           = ',r(irpol)
         WRITE(8,*)' radii ratio            = ',r(irpol)/r(ireq)
         WRITE(8,*)' circumferential radius = ',rcirc
         WRITE(8,*)' equatorial Omega       = ',a%eq_omega
         WRITE(8,*)" keplerian Omega        = ",a%kep_omega
         WRITE(8,*)" disk mass              = ",m0disk
         WRITE(8,*)" error flag             = ",a%info
         CLOSE(8)
      ENDIF
      
      ! print the profiles
      IF(binfile/="") THEN
         ! compute rest mass density rho, entropy per baryon s, temperature T
         rho(:,:)= 0d0
         s(:,:)= 0d0
         t(:,:)= 0d0
         DO ix=1, nth
            DO iz= a%wsurf(ix),1,-1
               CALL fun(a%p(ix,iz),a%hden(ix,iz),rho(ix,iz),s(ix,iz),t(ix,iz))
            ENDDO
         ENDDO
         ! print the unformatted file
         mr= MAXVAL(a%wsurf)
         INQUIRE(IOLENGTH=reclen)nth,nr,a%ieq,mr,r(1:nr),th(1:nth),a%wsurf,&
                                &rho(:,:mr),a%p(:,:mr),a%hden(:,:mr),a%vphi(:,:mr),&
                                &-a%btp,a%psi,a%lps,s(:,:mr),t(:,:mr)
         OPEN(8,FILE=binfile,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=reclen)
         ! 1:nth   2:nr   3:itheq   4:mr  :: int
         ! 5:r(nr)     6:th(nth)          :: double(:)
         ! 7:surf                         :: int(nth)
         ! 8:rho  9:p  10:hden  11:v^phi  :: double(nth,mr)
         ! 12:omega    13:psi    14:alpha :: double(nth,nr)
         ! 15:s        16:t               :: double(nth,mr)
         WRITE(8,REC=1)nth,nr,a%ieq,mr,r(1:nr),th(1:nth),a%wsurf,&
                      &rho(:,:mr),a%p(:,:mr),a%hden(:,:mr),a%vphi(:,:mr),&
                      &-a%btp,a%psi,a%lps,s(:,:mr),t(:,:mr)
         CLOSE(8)
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
   !!    (see Appendix B of Camelio+2019 for details)
   !! @param[in] p pressure
   !! @param[in] hden enthalpy density
   !! @param[out] rho rest mass density
   !! @param[out] s entropy
   !! @param[out] t temperature
   PURE SUBROUTINE p_hden2eos(p,hden,rho,s,t)
      USE system, ONLY: mn, rhon
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: hden,p
      DOUBLE PRECISION, INTENT(OUT) :: rho,s,t
      DOUBLE PRECISION :: f, tmp1, tmp2
      INTEGER :: info

      ! cold EOS
      IF(k3 <= 0d0) THEN
         rho= (p / k1 / (gam - 1d0))**(1d0/gam)
         s= 0d0
         t= 0d0
         RETURN
      ENDIF

      f= (gamth - 1d0)*hden - gamth*p

      IF(f <= 0d0) THEN
         ! BEWARE: I am discarding the high density solution even if it is physically valid in case of f <=0 and gam > gamth
         rho= 0d0
      ELSE IF(gam == gamth) THEN
         ! I need to invert a linear eq to find rho
         rho= f/(gamth - 1d0)
      ELSE IF(gam < gamth) THEN
         ! NOTE: I am not separating gam < gamth in subcases like I am doing for gam > gamth to avoid complicating the code
         ! solve A·ρ + B·ρ^Γ + C = 0
         rho= rhon ! first guess
         CALL newton_raphson( &  ! use Newton-Raphson
            & rho, &             ! output
            & (gamth - 1d0), &   ! coef. of ρ
            & (gamth - gam)*k1, &! coef. of ρ^Γ
            & gam, &             ! Γ
            & -f, &              ! constant term (C)
            & 0d0, &             ! left search boundary (min ρ)
            & HUGE(1d0), &       ! right search boundary (max ρ)
            & info)              ! status flag: 0 on success
         IF(info /= 0) rho= 0d0
      ! now it is gam > gamth and f > 0
      ELSE IF(f > funmax) THEN
         rho= 0d0
      ELSE IF(gam == 2d0) THEN
         ! I take the low density branch of the quadratic equation
         tmp1= gamth - 1d0
         tmp2= k1*(gamth - gam)
         ! note: the determinant cannot be < 0 because f <= fmax
         rho= 0.5d0*(DSQRT(MAX(0d0, tmp1*tmp1 + 4d0*tmp2*f)) - tmp1)/tmp2
      ! NOTE: I am commenting the next case because the inversion with the cubic for gam=1.5
      !    is less robust than Newton-Raphson, for a reason I have no time to identify.
      !    However, it works for the majority of the points and somebody may wish to debug it,
      !    so I leave the comment for this somebody or for the future me
      !ELSE IF(gam == 1.5d0) THEN
      !   ! I have to solve a cubic to find sqrt(rho)
      !   CALL solve_cubic((gamth - gam)*k1, & ! coeff of x^3
      !                   & gamth - 1d0,      & ! coeff of x^2
      !                   & 0d0,              & ! coeff of x
      !                   & -f,               & ! constant term
      !                   & tmp1, tmp2, rho)         ! roots
      !   rho= rho*rho
      ELSE IF(gam == 3d0) THEN
         ! I have to solve a cubic to find rho
         CALL solve_cubic((gamth - gam)*k1, & ! coeff of x^3
                         & 0d0,              & ! coeff of x^2
                         & gamth - 1d0,      & ! coeff of x
                         & -f,               & ! constant term
                         & tmp1, tmp2, rho)         ! roots
         ! BEWARE: taking the 3rd root is valid only when rho < rhocrit!
         !         when rho > rhocrit the correct root is the 1st
         !         the 2nd root is negative
      ELSE
         ! solve A·ρ + B·ρ^Γ + C = 0
         rho= 0.5d0*rhocrit ! first guess
         CALL newton_raphson( &  ! use Newton-Raphson
            & rho, &             ! output
            & (gamth - 1d0), &   ! coef. of ρ
            & (gamth - gam)*k1, &! coef. of ρ^Γ
            & gam, &             ! Γ
            & -f, &              ! constant term (C)
            & 0d0, &             ! left search boundary (min ρ)
            & rhocrit, &         ! right search boundary (max ρ)
            & info)              ! status flag: 0 on success
         IF(info /= 0) rho= 0d0
      ENDIF

      ! sanitize rho and compute entropy s and temperature t
      IF(rho <= 0d0) THEN
         rho= 0d0
         s= 0d0
         t= 0d0
      ELSE
         s= DSQRT(MAX(0d0, (p - (gam - 1d0)*k1*rho**gam)/((gamth - 1d0)*k2*rho**gamth)))
         t= 2d0*mn*k2*s*rho**(gamth - 1d0)
      ENDIF

   CONTAINS

      !> find the 3 real roots of: a x^3 + b x^2 + c x + d = 0
      !! @details it uses the Cardano's formula
      !! @param[in] a coefficient of x^3
      !! @param[in] b coefficient of x^2
      !! @param[in] c coefficient of x
      !! @param[in] d constant term
      !! @param[out] x1 first root
      !! @param[out] x2 second root
      !! @param[out] x3 third root
      PURE SUBROUTINE solve_cubic(a,b,c,d,x1,x2,x3)
      USE system, ONLY: pi
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: a,b,c,d
      DOUBLE PRECISION, INTENT(OUT) :: x1,x2,x3
      DOUBLE PRECISION :: p,q,r,tmp,tmp1,tmp2,angle
      DOUBLE PRECISION, PARAMETER :: oo3= 1d0/3d0
      DOUBLE PRECISION, PARAMETER :: pi23= pi*2d0/3d0
      
      p= -b/(3d0*a)
      q= p*p*p + (b*c - 3d0*a*d)/(6d0*a*a)
      r= c/(3d0*a)
      tmp= r - p*p
      tmp= q*q + tmp*tmp*tmp
      IF(tmp >= 0d0) THEN
         tmp= DSQRT(tmp)
         tmp1= q + tmp
         tmp1= DSIGN(DABS(tmp1)**oo3, tmp1)
         tmp2= q - tmp
         tmp2= DSIGN(DABS(tmp2)**oo3, tmp2)
         x1= tmp1 + tmp2 + p
         x2= x1
         x3= x1
      ELSE
         tmp= -tmp
         angle= DATAN2(DSQRT(tmp),q)/3d0
         tmp= 2d0*DSQRT(q*q + tmp)**oo3
         x1= tmp*DCOS(angle) + p
         x2= tmp*DCOS(angle + pi23) + p
         x3= tmp*DCOS(angle + pi23*2d0) + p
      ENDIF
      END SUBROUTINE solve_cubic

      !> solve a·ρ + b·ρ^Γ + c = 0 with Newton-Raphson
      !! @param[inout] x input: first guess, output: root
      !! @param[in] a coefficient of ρ
      !! @param[in] b coefficient of ρ^Γ
      !! @param[in] c constant term
      !! @param[in] x0 left bracket of the root
      !! @param[in] x1 right bracket of the root
      !! @param[out] info success flag (0 on success)
      PURE SUBROUTINE newton_raphson(x,a,b,g,c,x0,x1,info)
      ! uses maxit and tol from user module
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(INOUT) :: x
      DOUBLE PRECISION, INTENT(IN) :: a,b,g,c,x0,x1
      INTEGER, INTENT(OUT) :: info
      DOUBLE PRECISION :: f,df,dx
      INTEGER :: i

      IF(x < x0 .OR. x > x1) x= 0.5d0*(x0 + x1)
      DO i=1,maxit
         f= a*x + b*x**g + c
         df= a + g*b*x**(g - 1d0)
         dx= -f/df
         ! error
         IF(dx /= dx) THEN
            info= -1
            RETURN
         ENDIF
         ! point found
         IF(DABS(dx) < tol*DABS(x)) THEN
            info= 0
            RETURN
         ENDIF
         x= MIN(x1, MAX(x0, x + dx))
      ENDDO
      ! not converged
      info= 1
      END SUBROUTINE newton_raphson

   END SUBROUTINE p_hden2eos


   !> EOS from density and entropy
   !! @param[out] p pressure
   !! @param[out] hden enthalpy density
   !! @param[out] t temperature
   !! @param[in] rho rest mass density
   !! @param[in] s entropy per baryon
   PURE SUBROUTINE rho_s2eos(p,hden,rho,s,t)
      USE system, ONLY: mn
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: rho,s
      DOUBLE PRECISION, INTENT(OUT) :: p,hden,t
      
      t= 2d0*mn*k2*s*rho**(gamth - 1d0)
      p= k1*(gam - 1d0)*rho**gam + k2*s*s*(gamth - 1d0)*rho**gamth
      hden= rho + k1*gam*rho**gam + k2*s*s*gamth*rho**gamth
   END SUBROUTINE rho_s2eos


   !> effective barotropic law
   !! @param[in] rho rest mass density
   !! @result entropy per baryon
   PURE DOUBLE PRECISION FUNCTION rho2s(rho) RESULT(s)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: rho

      s= k3*rho**(0.5d0*(gam - gamth)) 
   END FUNCTION rho2s


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
      ! BEWARE: h and rho are the physical ones only when bvalue = 0!
      DOUBLE PRECISION :: heat,tmp1,tmp2,h,h0,rho,r02,tomg0
      DOUBLE PRECISION :: f,df,q,fcal,dfcal,ddfcal,eq,deq,domg
      INTEGER :: i

      info= -1
      
      ! ********************************** !
      ! *** determine angular velocity *** !
      ! ********************************** !

      ! the target central omg is slowly increased at each iteration to improve the convergence
      tomg0= omg0*MIN(1d0, DBLE(iter)/DBLE(relax_iters + 1))

      ! if it is not rigidly rotating
      IF(sigma > 0d0) THEN
         r02= 1d0/(sigma*sigma)
         ! solve with Newton-Raphson
         DO i=1, maxit
            ! compute F(Omega,r,theta)
            tmp1= omg + btp
            tmp2= alp2 - rcyl2*tmp1*tmp1
            ! check velocity < c
            IF(tmp2 <= 0d0) RETURN
            f= rcyl2*tmp1/tmp2
            df= rcyl2*(alp2 + rcyl2*tmp1*tmp1)/(tmp2*tmp2)
            ! compute Q(Omega,r,theta)
            q= -0.5d0*DLOG(tmp2)
            ! compute F(Omega)
            tmp1= omg - tomg0
            fcal= -0.5d0*tmp1*tmp1*r02
            dfcal= -tmp1*r02
            ddfcal= -r02
            ! compute equation and its derivate
            eq= f*(1d0 + bvalue*fcal) - dfcal*(1d0 + bvalue*(q - q0))
            deq= df*(1d0 + bvalue*fcal) - ddfcal*(1d0 + bvalue*(q - q0))
            ! compute new solution
            domg= -eq/deq
            ! check for NaN
            IF(domg /= domg) RETURN
            ! check convergence
            IF(DABS(domg)<tol*DABS(omg)) EXIT
            omg= omg + domg
         ENDDO
         IF(i == maxit + 1) RETURN
      ! if it is rigidly rotating
      ELSE
         omg= tomg0
         ! fcal is 0 for omg = tomg0, EFFECTIVELY accounting for rigid rotation
         fcal= 0d0
         ! compute Q(Omega,r,theta)
         tmp1= omg + btp
         tmp2= alp2 - rcyl2*tmp1*tmp1
         ! check velocity < c
         IF(tmp2 <= 0d0) RETURN
         q= -0.5d0*DLOG(tmp2)
      ENDIF

      ! *********************************************** !
      ! *** determine pressure and enthalpy density *** !
      ! *********************************************** !

      heat= (q - q0 - fcal)/(1d0 + bvalue*fcal)
      tmp1= (gam - 1d0)*k1 + (gamth - 1d0)*k2*k3*k3
      tmp2= gam*k1 + gamth*k2*k3*k3
      h0= 1d0 + tmp2*rho0**(gam - 1d0)
      ! BEWARE: h is the physical specific enthalpy only when bvalue = 0!
      h= h0*DEXP(heat*(gam - 1d0)/gam*tmp2/tmp1)
      IF(h <= 1d0) THEN
         p= 0d0
         hden= 0d0
      ELSE
         ! BEWARE: rho is the physical density only when bvalue = 0!
         rho= ((h - 1d0)/tmp2)**(1d0/(gam - 1d0))
         ! while p is always physical
         p= tmp1*rho**gam
         ! hden is obtained from derivative of the potential
         hden= rho*h/(1d0 + bvalue*fcal)
      ENDIF

      info= 0
   END SUBROUTINE solve_euler

END MODULE user


!> solve the stellar structure of a baroclinic neutron star
!! (the convective model of Camelio et al. (2019), PRD)
PROGRAM baroclinic_neutron_star
   USE system, ONLY: model_quantities,ns_init,ns_solve,fun_hden,eos
   USE user, ONLY: rho0,p2hden,p_hden2eos,solve_euler,rho2s,rho_s2eos,print_quantities,user_initialization
   IMPLICIT NONE
   TYPE(model_quantities) :: a
   DOUBLE PRECISION :: t,hden,p0
   INTEGER :: info

   ! initialize parameters, code settings, and thermal law
   CALL user_initialization(info)

   ! determine the stellar structure
   IF(info == 0) THEN
      ! get the central pressure
      CALL rho_s2eos(p0, hden, rho0, rho2s(rho0), t)
      ! initialize the non rotating star
      CALL ns_init(a, p2hden, p0, info)
      ! solve the rotating star and print result
      IF(info == 0) THEN
         ! solve the rotating star
         CALL ns_solve(a, solve_euler, info)
         ! print the results
         CALL print_quantities(a,p_hden2eos)
      ELSE
         PRINT*,"ERROR! bad initialization. Nothing to do"
      ENDIF
   ELSE
      a%info= info
      PRINT*,"ERROR in initialization!"
   ENDIF
END PROGRAM baroclinic_neutron_star
