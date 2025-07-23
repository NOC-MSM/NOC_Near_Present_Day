MODULE tide_mod
   !!======================================================================
   !!                       ***  MODULE  tide_mod  ***
   !! Compute nodal modulations corrections and pulsations
   !!======================================================================
   !! History :  1.0  !  2007  (O. Le Galloudec)  Original code
   !!                 !  2019  (S. Mueller)
   !!----------------------------------------------------------------------
   !!
   !! ** Reference :
   !!       S58) Schureman, P. (1958): Manual of Harmonic Analysis and
   !!            Prediction of Tides (Revised (1940) Edition (Reprinted 1958
   !!            with corrections). Reprinted June 2001). U.S. Department of
   !!            Commerce, Coast and Geodetic Survey Special Publication
   !!            No. 98. Washington DC, United States Government Printing
   !!            Office. 317 pp. DOI: 10.25607/OBP-155.
   !!----------------------------------------------------------------------

   USE oce, ONLY : ssh        ! sea-surface height
   USE par_oce                ! ocean parameters
   USE phycst, ONLY : rpi, rad, rday
   USE daymod, ONLY : ndt05   ! half-length of time step
   USE in_out_manager         ! I/O units
   USE iom                    ! xIOs server

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tide_init
   PUBLIC   tide_update           ! called by stp
   PUBLIC   tide_init_harmonics   ! called internally and by module diaharm
   PUBLIC   upd_tide              ! called in dynspg_... modules
   PUBLIC tide_init_diss
   INTEGER, PUBLIC, PARAMETER ::   jpmax_harmo = 64   !: maximum number of harmonic components

   TYPE ::    tide
      CHARACTER(LEN=4) ::   cname_tide = ''
      REAL(wp)         ::   equitide
      INTEGER          ::   nt, ns, nh, np, np1, shift
      INTEGER          ::   nksi, nnu0, nnu1, nnu2, R
      INTEGER          ::   nformula
   END TYPE tide

   TYPE(tide), DIMENSION(:), POINTER ::   tide_components   !: Array of selected tidal component parameters

   TYPE, PUBLIC ::   tide_harmonic       !:   Oscillation parameters of harmonic tidal components
      CHARACTER(LEN=4) ::   cname_tide   !    Name of component
      REAL(wp)         ::   equitide     !    Amplitude of equilibrium tide
      REAL(wp)         ::   f            !    Node factor
      REAL(wp)         ::   omega        !    Angular velocity
      REAL(wp)         ::   v0           !    Initial phase at prime meridian
      REAL(wp)         ::   u            !    Phase correction
   END type tide_harmonic

   TYPE(tide_harmonic), PUBLIC, DIMENSION(:), POINTER ::   tide_harmonics   !: Oscillation parameters of selected tidal components

   LOGICAL , PUBLIC ::   ln_tide             !:
   LOGICAL , PUBLIC ::   ln_tide_pot         !:
   INTEGER          ::   nn_tide_var         !  Variant of tidal parameter set and tide-potential computation
   LOGICAL          ::   ln_tide_dia         !  Enable tidal diagnostic output
   LOGICAL          ::   ln_read_load        !:
   LOGICAL , PUBLIC ::   ln_scal_load        !:
   LOGICAL , PUBLIC ::   ln_tide_ramp        !:
   INTEGER , PUBLIC ::   nb_harmo            !: Number of active tidal components
   REAL(wp), PUBLIC ::   rn_tide_ramp_dt     !:
   REAL(wp), PUBLIC ::   rn_scal_load        !:
   CHARACTER(lc), PUBLIC ::   cn_tide_load   !: 
   REAL(wp)         ::   rn_tide_gamma       ! Tidal tilt factor

   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   pot_astro            !: tidal potential
   REAL(wp),         ALLOCATABLE, DIMENSION(:,:)   ::   pot_astro_comp       !  tidal-potential component
   REAL(wp),         ALLOCATABLE, DIMENSION(:,:,:) ::   amp_pot, phi_pot
   REAL(wp),         ALLOCATABLE, DIMENSION(:,:,:) ::   amp_load, phi_load

   REAL(wp) ::   rn_tide_ramp_t   !   Elapsed time in seconds

   REAL(wp) ::   sh_T, sh_s, sh_h, sh_p, sh_p1             ! astronomic angles
   REAL(wp) ::   sh_xi, sh_nu, sh_nuprim, sh_nusec, sh_R   !
   REAL(wp) ::   sh_I, sh_x1ra, sh_N                       !

   ! Longitudes on 1 Jan 1900, 00h and angular velocities (units of deg and
   ! deg/h, respectively. The values of these module variables have been copied
   ! from subroutine astronomic_angle of the version of this module used in
   ! release version 4.0 of NEMO.
   REAL(wp) ::   rlon00_N  =  259.1560564_wp               ! Longitude of ascending lunar node
   REAL(wp) ::   romega_N  = -.0022064139_wp
   REAL(wp) ::   rlon00_T  =  180.0_wp                     ! Mean solar angle (GMT)
   REAL(wp) ::   romega_T  =  15.0_wp
   REAL(wp) ::   rlon00_h  =  280.1895014_wp               ! Mean solar Longitude
   REAL(wp) ::   romega_h  =  .0410686387_wp
   REAL(wp) ::   rlon00_s  =  277.0256206_wp               ! Mean lunar Longitude
   REAL(wp) ::   romega_s  =  .549016532_wp
   REAL(wp) ::   rlon00_p1 =  281.2208569_wp               ! Longitude of solar perigee
   REAL(wp) ::   romega_p1 =  .000001961_wp
   REAL(wp) ::   rlon00_p  =  334.3837214_wp               ! Longitude of lunar perigee
   REAL(wp) ::   romega_p  =  .004641834_wp
   ! Values of cos(i)*cos(epsilon), rcice, and sin(incl)*sin(epsilon), rsise,
   ! where i is the inclination of the orbit of the Moon w.r.t. the ecliptic and
   ! epsilon the obliquity of the ecliptic on 1 January 1900, 00h. The values of
   ! these module variables have been copied from subroutine astronomic_angle
   ! (computation of the cosine of inclination of orbit of Moon to the celestial
   ! equator) of the version of this module used in release version 4.0 of NEMO.
   REAL(wp) ::   rcice    = 0.913694997_wp
   REAL(wp) ::   rsise    = 0.035692561_wp
   ! Coefficients used to compute sh_xi and sh_nu in subroutine astronomic_angle
   ! according to two equations given in the explanation of Table 6 of S58
   REAL(wp) ::   rxinu1, rxinu2
   ! davbyr : Switch for internal wave drag on barotropic currents.
   LOGICAL , PUBLIC :: ln_int_wave_drag, ln_calc_tdiss !:
   CHARACTER(lc), PUBLIC ::   cn_int_wave_drag   !:
   CHARACTER(lc), PUBLIC ::   cn_h_rough        !:
   REAL(wp), PUBLIC      ::   rn_kappa_tdiss,tdiss_mindepth   !: Parameters for IW dissipation
   ! davbyr: Tide drag
   REAL(wp), PUBLIC,  ALLOCATABLE, SAVE, DIMENSION(:,:) ::   tdiss, h2rough,N2mean
   ! END davbyr
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: tide_mod.F90 13286 2020-07-09 15:48:29Z smasson $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
!! * Substitutions
#  include "do_loop_substitute.h90"

CONTAINS

   SUBROUTINE tide_init
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE tide_init  ***
      !!----------------------------------------------------------------------      
      INTEGER  :: ji, jk
      CHARACTER(LEN=4), DIMENSION(jpmax_harmo) :: sn_tide_cnames ! Names of selected tidal components
      INTEGER  ::   ios                 ! Local integer output status for namelist read
      !
      NAMELIST/nam_tide/ln_tide, nn_tide_var, ln_tide_dia, ln_tide_pot, rn_tide_gamma, &
         &              ln_scal_load, ln_read_load, cn_tide_load,         &
         &              ln_tide_ramp, rn_scal_load, rn_tide_ramp_dt,      &
         &              sn_tide_cnames,                                    &
         &              ln_int_wave_drag,ln_calc_tdiss, cn_int_wave_drag, &
         &              cn_h_rough,rn_kappa_tdiss,tdiss_mindepth
      !!----------------------------------------------------------------------
      !
      ! Initialise all array elements of sn_tide_cnames, as some of them
      ! typically do not appear in namelist_ref or namelist_cfg
      sn_tide_cnames(:) = ''
      ! Read Namelist nam_tide
      READ  ( numnam_ref, nam_tide, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nam_tide in reference namelist' )
      !
      READ  ( numnam_cfg, nam_tide, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'nam_tide in configuration namelist' )
      IF(lwm) WRITE ( numond, nam_tide )
      !
      IF( ln_tide ) THEN
         IF (lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'tide_init : Initialization of the tidal components'
            WRITE(numout,*) '~~~~~~~~~ '
            WRITE(numout,*) '   Namelist nam_tide'
            WRITE(numout,*) '      Use tidal components                       ln_tide         = ', ln_tide
            WRITE(numout,*) '         Variant (1: default; 0: legacy option)  nn_tide_var     = ', nn_tide_var
            WRITE(numout,*) '         Tidal diagnostic output                 ln_tide_dia     = ', ln_tide_dia
            WRITE(numout,*) '         Apply astronomical potential            ln_tide_pot     = ', ln_tide_pot
            WRITE(numout,*) '         Tidal tilt factor                       rn_tide_gamma   = ', rn_tide_gamma
            WRITE(numout,*) '         Use scalar approx. for load potential   ln_scal_load    = ', ln_scal_load
            WRITE(numout,*) '         Read load potential from file           ln_read_load    = ', ln_read_load
            WRITE(numout,*) '         Apply ramp on tides at startup          ln_tide_ramp    = ', ln_tide_ramp
            WRITE(numout,*) '         Fraction of SSH used in scal. approx.   rn_scal_load    = ', rn_scal_load
            WRITE(numout,*) '         Duration (days) of ramp                 rn_tide_ramp_dt = ', rn_tide_ramp_dt
            WRITE(numout,*) '         Internal Wave Drag Parameterization     ln_int_wave_drag = ', ln_int_wave_drag
            IF(ln_int_wave_drag) WRITE(numout,*) 'cn_int_wave_drag = ', cn_int_wave_drag
         ENDIF
            IF( ln_int_wave_drag )THEN !davbyr: Allocation and read tdiss
               ALLOCATE( tdiss(jpi, jpj) )
               tdiss(:,:) = 0.0_wp
               IF ( ln_calc_tdiss ) THEN
                ALLOCATE( h2rough(jpi, jpj) )
                ALLOCATE( N2mean(jpi, jpj) )
                h2rough(:,:)     = 0.0_wp
                N2mean(:, :) = 0.0_wp

               ENDIF
               CALL tide_init_diss
            ENDIF
      ELSE
         rn_scal_load = 0._wp 

         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tide_init : tidal components not used (ln_tide = F)'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~ '
         RETURN
      ENDIF
      !
      IF( ln_read_load.AND.(.NOT.ln_tide_pot) ) &
          &   CALL ctl_stop('ln_read_load requires ln_tide_pot')
      IF( ln_scal_load.AND.(.NOT.ln_tide_pot) ) &
          &   CALL ctl_stop('ln_scal_load requires ln_tide_pot')
      IF( ln_scal_load.AND.ln_read_load ) &
          &   CALL ctl_stop('Choose between ln_scal_load and ln_read_load')
      IF( ln_tide_ramp.AND.((nitend-nit000+1)*rn_Dt/rday < rn_tide_ramp_dt) )   &
         &   CALL ctl_stop('rn_tide_ramp_dt must be lower than run duration')
      IF( ln_tide_ramp.AND.(rn_tide_ramp_dt<0.) ) &
         &   CALL ctl_stop('rn_tide_ramp_dt must be positive')
      !
      ! Compute coefficients which are used in subroutine astronomic_angle to
      ! compute sh_xi and sh_nu according to two equations given in the
      ! explanation of Table 6 of S58
      rxinu1 = COS( 0.5_wp * ( ABS( ACOS( rcice + rsise ) ) ) ) / COS( 0.5_wp * ( ACOS( rcice - rsise ) ) )
      rxinu2 = SIN( 0.5_wp * ( ABS( ACOS( rcice + rsise ) ) ) ) / SIN( 0.5_wp * ( ACOS( rcice - rsise ) ) )
      !
      ! Initialise array used to store tidal oscillation parameters (frequency,
      ! amplitude, phase); also retrieve and store array of information about
      ! selected tidal components
      CALL tide_init_harmonics(sn_tide_cnames, tide_harmonics, tide_components)
      !
      ! Number of active tidal components
      nb_harmo = size(tide_components)
      !       
      ! Ensure that tidal components have been set in namelist_cfg
      IF( nb_harmo == 0 )   CALL ctl_stop( 'tide_init : No tidal components set in nam_tide' )
      !
      IF (.NOT.ln_scal_load ) rn_scal_load = 0._wp
      !
      ALLOCATE( amp_pot(jpi,jpj,nb_harmo),                      &
           &      phi_pot(jpi,jpj,nb_harmo), pot_astro(jpi,jpj)   )
      IF( ln_tide_dia ) ALLOCATE( pot_astro_comp(jpi,jpj) )
      IF( ln_read_load ) THEN
         ALLOCATE( amp_load(jpi,jpj,nb_harmo), phi_load(jpi,jpj,nb_harmo) )
         CALL tide_init_load
         amp_pot(:,:,:) = amp_load(:,:,:)
         phi_pot(:,:,:) = phi_load(:,:,:)
      ELSE		
         amp_pot(:,:,:) = 0._wp
         phi_pot(:,:,:) = 0._wp	 
      ENDIF
      !
   END SUBROUTINE tide_init


   SUBROUTINE tide_init_components(pcnames, ptide_comp)
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE tide_init_components  ***
      !!
      !! Returns pointer to array of variables of type 'tide' that contain
      !! information about the selected tidal components
      !! ----------------------------------------------------------------------
      CHARACTER(LEN=4),              DIMENSION(jpmax_harmo), INTENT(in)  ::   pcnames             ! Names of selected components
      TYPE(tide),       POINTER,     DIMENSION(:),           INTENT(out) ::   ptide_comp          ! Selected components
      INTEGER,          ALLOCATABLE, DIMENSION(:)                        ::   icomppos            ! Indices of selected components
      INTEGER                                                            ::   icomp, jk, jj, ji   ! Miscellaneous integers
      LOGICAL                                                            ::   llmatch             ! Local variables used for
      INTEGER                                                            ::   ic1, ic2            !    string comparison
      TYPE(tide),       POINTER,     DIMENSION(:)                        ::   tide_components     ! All available components
      
      ! Populate local array with information about all available tidal
      ! components
      !
      ! Note, here 'tide_components' locally overrides the global module
      ! variable of the same name to enable the use of the global name in the
      ! include file that contains the initialisation of elements of array
      ! 'tide_components'
      ALLOCATE(tide_components(jpmax_harmo), icomppos(jpmax_harmo))
      ! Initialise array of indices of the selected componenents
      icomppos(:) = 0
      ! Include tidal component parameters for all available components
      IF (nn_tide_var < 1) THEN
#define TIDE_VAR_0
#include "tide.h90"
#undef TIDE_VAR_0
      ELSE
#include "tide.h90"
      END IF
      ! Identify the selected components that are availble
      icomp = 0
      DO jk = 1, jpmax_harmo
         IF (TRIM(pcnames(jk)) /= '') THEN
            DO jj = 1, jpmax_harmo
               ! Find matches between selected and available constituents
               ! (ignore capitalisation unless legacy variant has been selected)
               IF (nn_tide_var < 1) THEN
                  llmatch = (TRIM(pcnames(jk)) == TRIM(tide_components(jj)%cname_tide))
               ELSE
                  llmatch = .TRUE.
                  ji = MAX(LEN_TRIM(pcnames(jk)), LEN_TRIM(tide_components(jj)%cname_tide))
                  DO WHILE (llmatch.AND.(ji > 0))
                     ic1 = IACHAR(pcnames(jk)(ji:ji))
                     IF ((ic1 >= 97).AND.(ic1 <= 122)) ic1 = ic1 - 32
                     ic2 = IACHAR(tide_components(jj)%cname_tide(ji:ji))
                     IF ((ic2 >= 97).AND.(ic2 <= 122)) ic2 = ic2 - 32
                     llmatch = (ic1 == ic2)
                     ji = ji - 1
                  END DO
               END IF
               IF (llmatch) THEN
                  ! Count and record the match
                  icomp = icomp + 1
                  icomppos(icomp) = jj
                  ! Set the capitalisation of the tidal constituent identifier
                  ! as specified in the namelist
                  tide_components(jj)%cname_tide = pcnames(jk)
                  IF (lwp) WRITE(numout, '(10X,"Tidal component #",I2.2,36X,"= ",A4)') icomp, tide_components(jj)%cname_tide
                  EXIT
               END IF
            END DO
            IF ((lwp).AND.(jj > jpmax_harmo)) WRITE(numout, '(10X,"Tidal component ",A4," is not available!")') pcnames(jk)
         END IF
      END DO
      
      ! Allocate and populate reduced list of components
      ALLOCATE(ptide_comp(icomp))
      DO jk = 1, icomp
         ptide_comp(jk) = tide_components(icomppos(jk))
      END DO
      
      ! Release local array of available components and list of selected
      ! components
      DEALLOCATE(tide_components, icomppos)
      
   END SUBROUTINE tide_init_components


   SUBROUTINE tide_init_harmonics(pcnames, ptide_harmo, ptide_comp)
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE tide_init_harmonics  ***
      !!
      !! Returns pointer to array of variables of type 'tide_harmonics' that
      !! contain oscillation parameters of the selected harmonic tidal
      !! components
      !! ----------------------------------------------------------------------
      CHARACTER(LEN=4),             DIMENSION(jpmax_harmo), INTENT(in) ::   pcnames     ! Names of selected components
      TYPE(tide_harmonic), POINTER, DIMENSION(:)                       ::   ptide_harmo ! Oscillation parameters of tidal components
      TYPE(tide),          POINTER, DIMENSION(:), OPTIONAL             ::   ptide_comp  ! Selected components
      TYPE(tide),          POINTER, DIMENSION(:)                       ::   ztcomp      ! Selected components

      ! Retrieve information about selected tidal components
      ! If requested, prepare tidal component array for returning
      IF (PRESENT(ptide_comp)) THEN
         CALL tide_init_components(pcnames, ptide_comp)
         ztcomp => ptide_comp
      ELSE
         CALL tide_init_components(pcnames, ztcomp)
      END IF

      ! Allocate and populate array of oscillation parameters
      ALLOCATE(ptide_harmo(size(ztcomp)))
      ptide_harmo(:)%cname_tide = ztcomp(:)%cname_tide
      ptide_harmo(:)%equitide = ztcomp(:)%equitide
      CALL tide_harmo(ztcomp, ptide_harmo)

   END SUBROUTINE tide_init_harmonics


   SUBROUTINE tide_init_potential
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE tide_init_potential  ***
      !!
      !! ** Reference :
      !!      CT71) Cartwright, D. E. and Tayler, R. J. (1971): New computations of
      !!            the Tide-generating Potential. Geophys. J. R. astr. Soc. 23,
      !!            pp. 45-74. DOI: 10.1111/j.1365-246X.1971.tb01803.x
      !!
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zcons, ztmp1, ztmp2, zlat, zlon, ztmp, zamp, zcs   ! local scalar
      !!----------------------------------------------------------------------

      IF( ln_read_load ) THEN
	 amp_pot(:,:,:) = amp_load(:,:,:)
	 phi_pot(:,:,:) = phi_load(:,:,:)
      ELSE		
         amp_pot(:,:,:) = 0._wp
         phi_pot(:,:,:) = 0._wp	 
      ENDIF
      DO jk = 1, nb_harmo
         zcons = rn_tide_gamma * tide_components(jk)%equitide * tide_harmonics(jk)%f
         DO ji = 1, jpi
            DO jj = 1, jpj
               ztmp1 =  tide_harmonics(jk)%f * amp_pot(ji,jj,jk) * COS( phi_pot(ji,jj,jk) &
                  &                                                   + tide_harmonics(jk)%v0 + tide_harmonics(jk)%u )
               ztmp2 = -tide_harmonics(jk)%f * amp_pot(ji,jj,jk) * SIN( phi_pot(ji,jj,jk) &
                  &                                                   + tide_harmonics(jk)%v0 + tide_harmonics(jk)%u )
               zlat = gphit(ji,jj)*rad !! latitude en radian
               zlon = glamt(ji,jj)*rad !! longitude en radian
               ztmp = tide_harmonics(jk)%v0 + tide_harmonics(jk)%u + tide_components(jk)%nt * zlon
               ! le potentiel est composé des effets des astres:
               SELECT CASE( tide_components(jk)%nt )
               CASE( 0 )                                                  !   long-periodic tidal constituents (included unless
                  zcs = zcons * ( 0.5_wp - 1.5_wp * SIN( zlat )**2 )      !   compatibility with original formulation is requested)
                  IF ( nn_tide_var < 1 ) zcs = 0.0_wp
               CASE( 1 )                                                  !   diurnal tidal constituents
                  zcs = zcons * SIN( 2.0_wp*zlat )
               CASE( 2 )                                                  !   semi-diurnal tidal constituents
                  zcs = zcons * COS( zlat )**2
               CASE( 3 )                                                  !   Terdiurnal tidal constituents; the colatitude-dependent
                  zcs = zcons * COS( zlat )**3                            !   factor is sin(theta)^3 (Table 2 of CT71)
               CASE DEFAULT                                               !   constituents of higher frequency are not included
                  zcs = 0.0_wp
               END SELECT
               ztmp1 = ztmp1 + zcs * COS( ztmp )
               ztmp2 = ztmp2 - zcs * SIN( ztmp )
               zamp = SQRT( ztmp1*ztmp1 + ztmp2*ztmp2 )
               amp_pot(ji,jj,jk) = zamp
               phi_pot(ji,jj,jk) = ATAN2( -ztmp2 / MAX( 1.e-10_wp , zamp ) ,   &
                  &                        ztmp1 / MAX( 1.e-10_wp,  zamp )   )
            END DO
         END DO
      END DO
      !
   END SUBROUTINE tide_init_potential


   SUBROUTINE tide_init_load
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE tide_init_load  ***
      !!----------------------------------------------------------------------
      INTEGER :: inum                 ! Logical unit of input file
      INTEGER :: ji, jj, itide        ! dummy loop indices
      REAL(wp), DIMENSION(jpi,jpj) ::   ztr, zti   !: workspace to read in tidal harmonics data 
      !!----------------------------------------------------------------------
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'tide_init_load : Initialization of load potential from file'
         WRITE(numout,*) '~~~~~~~~~~~~~~ '
      ENDIF
      !
      CALL iom_open ( cn_tide_load , inum )
      !
      DO itide = 1, nb_harmo
         CALL iom_get  ( inum, jpdom_global,TRIM(tide_components(itide)%cname_tide)//'_z1', ztr(:,:) )
         CALL iom_get  ( inum, jpdom_global,TRIM(tide_components(itide)%cname_tide)//'_z2', zti(:,:) )
         !
         DO ji=1,jpi
            DO jj=1,jpj
               amp_load(ji,jj,itide) =  SQRT( ztr(ji,jj)**2. + zti(ji,jj)**2. )
               phi_load(ji,jj,itide) = ATAN2(-zti(ji,jj), ztr(ji,jj) )
            END DO
         END DO
         !
      END DO
      CALL iom_close( inum )
      !
   END SUBROUTINE tide_init_load


   SUBROUTINE tide_update( kt )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE tide_update  ***
      !!----------------------------------------------------------------------      
      INTEGER, INTENT( in ) ::   kt     ! ocean time-step
      INTEGER               ::   jk     ! dummy loop index
      !!----------------------------------------------------------------------
      
      IF( nsec_day == NINT(0.5_wp * rn_Dt) .OR. kt == nit000 ) THEN      ! start a new day
         !
         CALL tide_harmo(tide_components, tide_harmonics, ndt05) ! Update oscillation parameters of tidal components for start of current day
         !
         !
         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'tide_update : Update of the components and (re)Init. the potential at kt=', kt
            WRITE(numout,*) '~~~~~~~~~~~ '
            DO jk = 1, nb_harmo
               WRITE(numout,*) tide_harmonics(jk)%cname_tide, tide_harmonics(jk)%u, &
                  &            tide_harmonics(jk)%f,tide_harmonics(jk)%v0, tide_harmonics(jk)%omega
            END DO
         ENDIF
         !
         IF( ln_tide_pot )   CALL tide_init_potential
         !
         rn_tide_ramp_t = (kt - nit000)*rn_Dt !   Elapsed time in seconds
      ENDIF
      !
   END SUBROUTINE tide_update


   SUBROUTINE tide_harmo( ptide_comp, ptide_harmo, psec_day )
      !
      TYPE(tide),          DIMENSION(:), POINTER ::   ptide_comp   ! Array of selected tidal component parameters
      TYPE(tide_harmonic), DIMENSION(:), POINTER ::   ptide_harmo  ! Oscillation parameters of selected tidal components
      INTEGER, OPTIONAL ::   psec_day                              ! Number of seconds since the start of the current day
      !
      IF (PRESENT(psec_day)) THEN 
         CALL astronomic_angle(psec_day)
      ELSE
         CALL astronomic_angle(nsec_day)
      END IF
      CALL tide_pulse( ptide_comp, ptide_harmo )
      CALL tide_vuf(   ptide_comp, ptide_harmo )
      !
   END SUBROUTINE tide_harmo


   SUBROUTINE astronomic_angle(psec_day)
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE astronomic_angle  ***
      !!                      
      !! ** Purpose : Compute astronomic angles
      !!----------------------------------------------------------------------
      INTEGER  ::   psec_day !   Number of seconds from midnight
      REAL(wp) ::   zp, zq, zt2, zs2, ztgI2, zP1, ztgn2, zat1, zat2
      REAL(wp) ::   zqy , zsy, zday, zdj, zhfrac, zt
      !!----------------------------------------------------------------------
      !
      ! Computation of the time from 1 Jan 1900, 00h in years
      zqy = AINT( (nyear - 1901.0_wp) / 4.0_wp )
      zsy = nyear - 1900.0_wp
      !
      zdj  = dayjul( nyear, nmonth, nday )
      zday = zdj + zqy - 1.0_wp
      !
      zhfrac = psec_day / 3600.0_wp
      !
      zt = zsy * 365.0_wp * 24.0_wp + zday * 24.0_wp + zhfrac
      !
      ! Longitude of ascending lunar node
      sh_N  = ( rlon00_N  + romega_N * zt     ) * rad
      sh_N = MOD( sh_N,  2*rpi )
      ! Mean solar angle (Greenwhich time)
      sh_T  = ( rlon00_T  + romega_T * zhfrac ) * rad
      ! Mean solar Longitude
      sh_h  = ( rlon00_h  + romega_h * zt     ) * rad
      sh_h = MOD( sh_h,  2*rpi )
      ! Mean lunar Longitude
      sh_s  = ( rlon00_s  + romega_s * zt     ) * rad
      sh_s = MOD( sh_s,  2*rpi )
      ! Longitude of solar perigee
      sh_p1 = ( rlon00_p1 + romega_p1 * zt    ) * rad
      sh_p1= MOD( sh_p1, 2*rpi )
      ! Longitude of lunar perigee
      sh_p  = ( rlon00_p  + romega_p  * zt    ) * rad
      sh_p = MOD( sh_p,  2*rpi )
      !
      ! Inclination of the orbit of the moon w.r.t. the celestial equator, see
      ! explanation of Table 6 of S58
      sh_I = ACOS( rcice - rsise * COS( sh_N ) )
      !
      ! Computation of sh_xi and sh_nu, see explanation of Table 6 of S58
      ztgn2 = TAN( sh_N / 2.0_wp )
      zat1 = ATAN( rxinu1 * ztgn2 )
      zat2 = ATAN( rxinu2 * ztgn2 )
      sh_xi = sh_N - zat1 - zat2
      IF( sh_N > rpi ) sh_xi = sh_xi - 2.0_wp * rpi
      sh_nu = zat1 - zat2
      !
      ! Computation of sh_x1ra, sh_R, sh_nuprim, and sh_nusec used for tidal
      ! constituents L2, K1, and K2
      !
      ! Computation of sh_x1ra and sh_R (Equations 204, 213, and 214 of S58)
      ztgI2 = tan( sh_I / 2.0_wp )
      zP1   = sh_p - sh_xi
      zt2   = ztgI2 * ztgI2
      sh_x1ra = SQRT( 1.0 - 12.0 * zt2 * COS( 2.0_wp * zP1 ) + 36.0_wp * zt2 * zt2 )
      zp = SIN( 2.0_wp * zP1 )
      zq = 1.0_wp / ( 6.0_wp * zt2 ) - COS( 2.0_wp * zP1 )
      sh_R = ATAN( zp / zq )
      !
      ! Computation of sh_nuprim (Equation 224 of S58)
      zp = SIN( 2.0_wp * sh_I ) * SIN( sh_nu )
      zq = SIN( 2.0_wp * sh_I ) * COS( sh_nu ) + 0.3347_wp
      sh_nuprim = ATAN( zp / zq )
      !
      ! Computation of sh_nusec  (Equation 232 of S58)
      zs2 = SIN( sh_I ) * SIN( sh_I )
      zp  = zs2 * SIN( 2.0_wp * sh_nu )
      zq  = zs2 * COS( 2.0_wp * sh_nu ) + 0.0727_wp
      sh_nusec = 0.5_wp * ATAN( zp / zq )
      !
   END SUBROUTINE astronomic_angle


   SUBROUTINE tide_pulse( ptide_comp, ptide_harmo )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE tide_pulse  ***
      !!                      
      !! ** Purpose : Compute tidal frequencies
      !!----------------------------------------------------------------------
      TYPE(tide),          DIMENSION(:), POINTER ::   ptide_comp   ! Array of selected tidal component parameters
      TYPE(tide_harmonic), DIMENSION(:), POINTER ::   ptide_harmo  ! Oscillation parameters of selected tidal components
      !
      INTEGER  ::   jh
      REAL(wp) ::   zscale
      !!----------------------------------------------------------------------
      !
      zscale =  rad / 3600.0_wp
      !
      DO jh = 1, size(ptide_harmo)
         ptide_harmo(jh)%omega = (  romega_T * ptide_comp( jh )%nT   &
            &                     + romega_s * ptide_comp( jh )%ns   &
            &                     + romega_h * ptide_comp( jh )%nh   &
            &                     + romega_p * ptide_comp( jh )%np   &
            &                     + romega_p1* ptide_comp( jh )%np1  ) * zscale
      END DO
      !
   END SUBROUTINE tide_pulse


   SUBROUTINE tide_vuf( ptide_comp, ptide_harmo )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE tide_vuf  ***
      !!                      
      !! ** Purpose : Compute nodal modulation corrections
      !!
      !! ** Outputs : vt: Phase of tidal potential relative to Greenwich (radians)
      !!              ut: Phase correction u due to nodal motion (radians)
      !!              ft: Nodal correction factor
      !!----------------------------------------------------------------------
      TYPE(tide),          DIMENSION(:), POINTER ::   ptide_comp   ! Array of selected tidal component parameters
      TYPE(tide_harmonic), DIMENSION(:), POINTER ::   ptide_harmo  ! Oscillation parameters of selected tidal components
      !
      INTEGER ::   jh   ! dummy loop index
      !!----------------------------------------------------------------------
      !
      DO jh = 1, size(ptide_harmo)
         !  Phase of the tidal potential relative to the Greenwhich 
         !  meridian (e.g. the position of the fictuous celestial body). Units are radian:
         ptide_harmo(jh)%v0 = sh_T * ptide_comp( jh )%nT    &
            &                   + sh_s * ptide_comp( jh )%ns    &
            &                   + sh_h * ptide_comp( jh )%nh    &
            &                   + sh_p * ptide_comp( jh )%np    &
            &                   + sh_p1* ptide_comp( jh )%np1   &
            &                   +        ptide_comp( jh )%shift * rad
         !
         !  Phase correction u due to nodal motion. Units are radian:
         ptide_harmo(jh)%u = sh_xi     * ptide_comp( jh )%nksi   &
            &                  + sh_nu     * ptide_comp( jh )%nnu0   &
            &                  + sh_nuprim * ptide_comp( jh )%nnu1   &
            &                  + sh_nusec  * ptide_comp( jh )%nnu2   &
            &                  + sh_R      * ptide_comp( jh )%R

         !  Nodal correction factor:
         ptide_harmo(jh)%f = nodal_factort( ptide_comp( jh )%nformula )
      END DO
      !
   END SUBROUTINE tide_vuf


   RECURSIVE FUNCTION nodal_factort( kformula ) RESULT( zf )
      !!----------------------------------------------------------------------
      !!                  ***  FUNCTION nodal_factort  ***
      !!
      !! ** Purpose : Compute amplitude correction factors due to nodal motion
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kformula
      !
      REAL(wp)            ::   zf
      REAL(wp)            ::   zs, zf1, zf2
      CHARACTER(LEN=3)    ::   clformula
      !!----------------------------------------------------------------------
      !
      SELECT CASE( kformula )
      !
      CASE( 0 )                  ! Formula 0, solar waves
         zf = 1.0
         !
      CASE( 1 )                  ! Formula 1, compound waves (78 x 78)
         zf=nodal_factort( 78 )
         zf = zf * zf
         !
      CASE ( 4 )                 ! Formula 4,  compound waves (78 x 235) 
         zf1 = nodal_factort( 78 )
         zf  = nodal_factort(235)
         zf  = zf1 * zf
         !
      CASE( 18 )                 ! Formula 18,  compound waves (78 x 78 x 78 )
         zf1 = nodal_factort( 78 )
         zf  = zf1 * zf1 * zf1
         !
      CASE( 20 )                 ! Formula 20, compound waves ( 78 x 78 x 78 x 78 )
         zf1 = nodal_factort( 78 )
         zf  = zf1 * zf1 * zf1 * zf1
         !
      CASE( 73 )                 ! Formula 73 of S58
         zs = SIN( sh_I )
         zf = ( 2.0_wp / 3.0_wp - zs * zs ) / 0.5021_wp
         !
      CASE( 74 )                 ! Formula 74 of S58
         zs = SIN(sh_I)
         zf = zs * zs / 0.1578_wp
         !
      CASE( 75 )                 ! Formula 75 of S58
         zs = COS( sh_I / 2.0_wp )
         zf = SIN( sh_I ) * zs * zs / 0.3800_wp
         !
      CASE( 76 )                 ! Formula 76 of S58
         zf = SIN( 2.0_wp * sh_I ) / 0.7214_wp
         !
      CASE( 78 )                 ! Formula 78 of S58
         zs = COS( sh_I/2 )
         zf = zs * zs * zs * zs / 0.9154_wp
         !
      CASE( 149 )                ! Formula 149 of S58
         zs = COS( sh_I/2 )
         zf = zs * zs * zs * zs * zs * zs / 0.8758_wp
         !
      CASE( 215 )                ! Formula 215 of S58 with typo correction (0.9154 instead of 0.9145)
         zs = COS( sh_I/2 )
         zf = zs * zs * zs * zs / 0.9154_wp * sh_x1ra
         !
      CASE( 227 )                ! Formula 227 of S58
         zs = SIN( 2.0_wp * sh_I )
         zf = SQRT( 0.8965_wp * zs * zs + 0.6001_wp * zs * COS( sh_nu ) + 0.1006_wp )
         !
      CASE ( 235 )               ! Formula 235 of S58
         zs = SIN( sh_I )
         zf = SQRT( 19.0444_wp * zs * zs * zs * zs + 2.7702_wp * zs * zs * cos( 2.0_wp * sh_nu ) + 0.0981_wp )
         !
      CASE DEFAULT
         WRITE( clformula, '(I3)' ) kformula
         CALL ctl_stop('nodal_factort: formula ' // clformula // ' is not available')
      END SELECT
      !
   END FUNCTION nodal_factort


   FUNCTION dayjul( kyr, kmonth, kday )
      !!----------------------------------------------------------------------
      !!                   ***  FUNCTION dayjul  ***
      !!
      !! Purpose :  compute the Julian day
      !!----------------------------------------------------------------------
      INTEGER,INTENT(in)    ::   kyr, kmonth, kday
      !
      INTEGER,DIMENSION(12) ::   idayt = (/ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 /)
      INTEGER,DIMENSION(12) ::   idays
      INTEGER               ::   inc, ji, zyq
      REAL(wp)              ::   dayjul
      !!----------------------------------------------------------------------
      !
      idays(1) = 0
      idays(2) = 31
      inc = 0.0_wp
      zyq = MOD( kyr - 1900 , 4 )
      IF( zyq == 0 ) inc = 1
      DO ji = 3, 12
         idays(ji) = idayt(ji) + inc
      END DO
      dayjul = REAL( idays(kmonth) + kday, KIND=wp )
      !
   END FUNCTION dayjul


   SUBROUTINE upd_tide(pdelta, Kmm)
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE upd_tide  ***
      !!
      !! ** Purpose :   provide at each time step the astronomical potential
      !!
      !! ** Method  :   computed from pulsation and amplitude of all tide components
      !!
      !! ** Action  :   pot_astro   actronomical potential
      !!----------------------------------------------------------------------      
      REAL(wp), INTENT(in)          ::   pdelta      ! Temporal offset in seconds
      INTEGER, INTENT(IN)           ::   Kmm         ! Time level index
      INTEGER                       ::   jk          ! Dummy loop index
      REAL(wp)                      ::   zt, zramp   ! Local scalars
      REAL(wp), DIMENSION(nb_harmo) ::   zwt         ! Temporary array
      !!----------------------------------------------------------------------      
      !
      zwt(:) = tide_harmonics(:)%omega * pdelta
      !
      IF( ln_tide_ramp ) THEN         ! linear increase if asked
         zt = rn_tide_ramp_t + pdelta
         zramp = MIN(  MAX( zt / (rn_tide_ramp_dt*rday) , 0._wp ) , 1._wp  )
      ENDIF
      !
      pot_astro(:,:) = 0._wp          ! update tidal potential (sum of all harmonics)
      DO jk = 1, nb_harmo
         IF ( .NOT. ln_tide_dia ) THEN
            pot_astro(:,:) = pot_astro(:,:) + amp_pot(:,:,jk) * COS( zwt(jk) + phi_pot(:,:,jk) )
         ELSE
            pot_astro_comp(:,:) = amp_pot(:,:,jk) * COS( zwt(jk) + phi_pot(:,:,jk) )
            pot_astro(:,:) = pot_astro(:,:) + pot_astro_comp(:,:)
            IF ( iom_use( "tide_pot_" // TRIM( tide_harmonics(jk)%cname_tide ) ) ) THEN   ! Output tidal potential (incl. load potential)
               IF ( ln_tide_ramp ) pot_astro_comp(:,:) = zramp * pot_astro_comp(:,:)
               CALL iom_put( "tide_pot_" // TRIM( tide_harmonics(jk)%cname_tide ), pot_astro_comp(:,:) )
            END IF
         END IF
      END DO
      !
      IF ( ln_tide_ramp ) pot_astro(:,:) = zramp * pot_astro(:,:)
      !
      IF( ln_tide_dia ) THEN          ! Output total tidal potential (incl. load potential)
         IF ( iom_use( "tide_pot" ) ) CALL iom_put( "tide_pot", pot_astro(:,:) + rn_scal_load * ssh(:,:,Kmm) )
      END IF
      !
   END SUBROUTINE upd_tide

   SUBROUTINE tide_init_diss !davbyr: subroutine
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE tide_init_diss  ***
      !!----------------------------------------------------------------------
      INTEGER :: inum                 ! Logical unit of input file
      INTEGER :: ji, jj,jidbg,jjdbg    ! dummy loop indices

      !!----------------------------------------------------------------------
      IF (ln_calc_tdiss) THEN
        IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'tide_init_diss : Read roughness (h) from file:',cn_h_rough
         WRITE(numout,*) '~~~~~~~~~~~~~~ '
        ENDIF

        CALL iom_open(cn_h_rough, inum)
        CALL iom_get (inum, jpdom_auto, 'h', h2rough(:,:))
        CALL iom_close(inum)
        CALL iom_close( inum )
        ! mask hrough in shllow water

       DO_2D( 0, 0, 0, 0 )
          h2rough(ji,jj) = h2rough(ji,jj) * h2rough(ji,jj) ! read in h, need h^2
          IF ( gdepw_0(ji,jj,mbkt(ji,jj)+1) < tdiss_mindepth ) THEN
           h2rough(ji,jj) = 0.0
          ENDIF
       END_2D
        CALL lbc_lnk( 'tide_init_diss', h2rough, 'T', 1._wp )

        !!!!!!!
      ELSE
        IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'tide_init_diss : Read tidal dissipation from file:',cn_int_wave_drag
         WRITE(numout,*) '~~~~~~~~~~~~~~ '
        ENDIF

        CALL iom_open(cn_int_wave_drag, inum)
        CALL iom_get (inum, jpdom_auto, 'tdiss', tdiss(:,:))
        CALL iom_close(inum)
        CALL iom_close( inum )

      ENDIF

   END SUBROUTINE tide_init_diss

   !!======================================================================
END MODULE tide_mod
