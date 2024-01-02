MODULE icethd_zdf_BL99
   !!======================================================================
   !!                       ***  MODULE icethd_zdf_BL99 ***
   !!   sea-ice: vertical heat diffusion in sea ice (computation of temperatures)
   !!======================================================================
   !! History :       !  2003-02  (M. Vancoppenolle) original 1D code
   !!                 !  2005-06  (M. Vancoppenolle) 3d version
   !!            4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!  ice_thd_zdf_BL99 : vertical diffusion computation
   !!----------------------------------------------------------------------
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants (ocean directory)
   USE ice            ! sea-ice: variables
   USE ice1D          ! sea-ice: thermodynamics variables
   USE icevar         ! sea-ice: operations
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_thd_zdf_BL99   ! called by icethd_zdf

   !!----------------------------------------------------------------------
   !! NEMO/ICE 4.0 , NEMO Consortium (2018)
   !! $Id: icethd_zdf_bl99.F90 14072 2020-12-04 07:48:38Z laurent $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_thd_zdf_BL99( k_cnd )
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE ice_thd_zdf_BL99  ***
      !!
      !! ** Purpose : computes the time evolution of snow and sea-ice temperature
      !!              profiles, using the original Bitz and Lipscomb (1999) algorithm
      !!
      !! ** Method  : solves the heat equation diffusion with a Neumann boundary
      !!              condition at the surface and a Dirichlet one at the bottom.
      !!              Solar radiation is partially absorbed into the ice.
      !!              The specific heat and thermal conductivities depend on ice
      !!              salinity and temperature to take into account brine pocket
      !!              melting. The numerical scheme is an iterative Crank-Nicolson
      !!              on a non-uniform multilayer grid in the ice and snow system.
      !!
      !!           The successive steps of this routine are
      !!           1.  initialization of ice-snow layers thicknesses
      !!           2.  Internal absorbed and transmitted radiation
      !!           Then iterative procedure begins
      !!           3.  Thermal conductivity
      !!           4.  Kappa factors
      !!           5.  specific heat in the ice
      !!           6.  eta factors
      !!           7.  surface flux computation
      !!           8.  tridiagonal system terms
      !!           9.  solving the tridiagonal system with Gauss elimination
      !!           Iterative procedure ends according to a criterion on evolution
      !!           of temperature
      !!           10. Fluxes at the interfaces
      !!
      !! ** Inputs / Ouputs : (global commons)
      !!           surface temperature              : t_su_1d
      !!           ice/snow temperatures            : t_i_1d, t_s_1d
      !!           ice salinities                   : sz_i_1d
      !!           number of layers in the ice/snow : nlay_i, nlay_s
      !!           total ice/snow thickness         : h_i_1d, h_s_1d
      !!-------------------------------------------------------------------
      INTEGER, INTENT(in) ::   k_cnd     ! conduction flux (off, on, emulated)
      !
      INTEGER ::   ji, jk         ! spatial loop index
      INTEGER ::   jm             ! current reference number of equation
      INTEGER ::   jm_mint, jm_maxt
      INTEGER ::   iconv          ! number of iterations in iterative procedure
      INTEGER ::   iconv_max = 50 ! max number of iterations in iterative procedure
      !
      INTEGER, DIMENSION(jpij) ::   jm_min    ! reference number of top equation
      INTEGER, DIMENSION(jpij) ::   jm_max    ! reference number of bottom equation

      LOGICAL, DIMENSION(jpij) ::   l_T_converged   ! true when T converges (per grid point)
      !
      REAL(wp) ::   zg1s      =  2._wp        ! for the tridiagonal system
      REAL(wp) ::   zg1       =  2._wp        !
      REAL(wp) ::   zgamma    =  18009._wp    ! for specific heat
      REAL(wp) ::   zbeta     =  0.117_wp     ! for thermal conductivity (could be 0.13)
      REAL(wp) ::   zkimin    =  0.10_wp      ! minimum ice thermal conductivity
      REAL(wp) ::   ztsu_err  =  1.e-5_wp     ! range around which t_su is considered at 0C
      REAL(wp) ::   zdti_bnd  =  1.e-4_wp     ! maximal authorized error on temperature
      REAL(wp) ::   zhs_ssl   =  0.03_wp      ! surface scattering layer in the snow
      REAL(wp) ::   zhi_ssl   =  0.10_wp      ! surface scattering layer in the ice
      REAL(wp) ::   zh_min    =  1.e-3_wp     ! minimum ice/snow thickness for conduction
      REAL(wp) ::   ztmelts                   ! ice melting temperature
      REAL(wp) ::   zdti_max                  ! current maximal error on temperature
      REAL(wp) ::   zcpi                      ! Ice specific heat
      REAL(wp) ::   zhfx_err, zdq             ! diag errors on heat
      !
      REAL(wp), DIMENSION(jpij) ::   zraext_s     ! extinction coefficient of radiation in the snow
      REAL(wp), DIMENSION(jpij) ::   ztsub        ! surface temperature at previous iteration
      REAL(wp), DIMENSION(jpij) ::   zh_i, z1_h_i ! ice layer thickness
      REAL(wp), DIMENSION(jpij) ::   zh_s, z1_h_s ! snow layer thickness
      REAL(wp), DIMENSION(jpij) ::   zqns_ice_b   ! solar radiation absorbed at the surface
      REAL(wp), DIMENSION(jpij) ::   zfnet        ! surface flux function
      REAL(wp), DIMENSION(jpij) ::   zdqns_ice_b  ! derivative of the surface flux function
      !
      REAL(wp), DIMENSION(jpij       )   ::   ztsuold     ! Old surface temperature in the ice
      REAL(wp), DIMENSION(jpij,nlay_i)   ::   ztiold      ! Old temperature in the ice
      REAL(wp), DIMENSION(jpij,nlay_s)   ::   ztsold      ! Old temperature in the snow
      REAL(wp), DIMENSION(jpij,nlay_i)   ::   ztib        ! Temporary temperature in the ice to check the convergence
      REAL(wp), DIMENSION(jpij,nlay_s)   ::   ztsb        ! Temporary temperature in the snow to check the convergence
      REAL(wp), DIMENSION(jpij,0:nlay_i) ::   ztcond_i    ! Ice thermal conductivity
      REAL(wp), DIMENSION(jpij,0:nlay_i) ::   ztcond_i_cp ! copy
      REAL(wp), DIMENSION(jpij,0:nlay_i) ::   zradtr_i    ! Radiation transmitted through the ice
      REAL(wp), DIMENSION(jpij,0:nlay_i) ::   zradab_i    ! Radiation absorbed in the ice
      REAL(wp), DIMENSION(jpij,0:nlay_i) ::   zkappa_i    ! Kappa factor in the ice
      REAL(wp), DIMENSION(jpij,0:nlay_i) ::   zeta_i      ! Eta factor in the ice
      REAL(wp), DIMENSION(jpij,0:nlay_s) ::   zradtr_s    ! Radiation transmited through the snow
      REAL(wp), DIMENSION(jpij,0:nlay_s) ::   zradab_s    ! Radiation absorbed in the snow
      REAL(wp), DIMENSION(jpij,0:nlay_s) ::   zkappa_s    ! Kappa factor in the snow
      REAL(wp), DIMENSION(jpij,0:nlay_s) ::   zeta_s      ! Eta factor in the snow
      REAL(wp), DIMENSION(jpij)          ::   zkappa_comb ! Combined snow and ice surface conductivity
      REAL(wp), DIMENSION(jpij)          ::   zq_ini      ! diag errors on heat
      REAL(wp), DIMENSION(jpij)          ::   zghe        ! G(he), th. conduct enhancement factor, mono-cat
      REAL(wp), DIMENSION(jpij)          ::   za_s_fra    ! ice fraction covered by snow
      REAL(wp), DIMENSION(jpij)          ::   isnow       ! snow presence (1) or not (0)
      REAL(wp), DIMENSION(jpij)          ::   isnow_comb  ! snow presence for met-office
      REAL(wp), DIMENSION(jpij,nlay_i+nlay_s+1)   ::   zindterm    ! 'Ind'ependent term
      REAL(wp), DIMENSION(jpij,nlay_i+nlay_s+1)   ::   zindtbis    ! Temporary 'ind'ependent term
      REAL(wp), DIMENSION(jpij,nlay_i+nlay_s+1)   ::   zdiagbis    ! Temporary 'dia'gonal term
      REAL(wp), DIMENSION(jpij,nlay_i+nlay_s+1,3) ::   ztrid       ! Tridiagonal system terms
      !
      ! Mono-category
      REAL(wp) ::   zepsilon   ! determines thres. above which computation of G(h) is done
      REAL(wp) ::   zhe        ! dummy factor
      REAL(wp) ::   zcnd_i     ! mean sea ice thermal conductivity
      !!------------------------------------------------------------------

      ! --- diag error on heat diffusion - PART 1 --- !
      DO ji = 1, npti
         zq_ini(ji) = ( SUM( e_i_1d(ji,1:nlay_i) ) * h_i_1d(ji) * r1_nlay_i +  &
            &           SUM( e_s_1d(ji,1:nlay_s) ) * h_s_1d(ji) * r1_nlay_s )
      END DO

      ! calculate ice fraction covered by snow for radiation
      CALL ice_var_snwfra( h_s_1d(1:npti), za_s_fra(1:npti) )

      !------------------
      ! 1) Initialization
      !------------------
      !
      ! extinction radiation in the snow
      IF    ( nn_qtrice == 0 ) THEN   ! constant
         zraext_s(1:npti) = rn_kappa_s
      ELSEIF( nn_qtrice == 1 ) THEN   ! depends on melting/freezing conditions
         WHERE( t_su_1d(1:npti) < rt0 )   ;   zraext_s(1:npti) = rn_kappa_sdry   ! no surface melting
         ELSEWHERE                        ;   zraext_s(1:npti) = rn_kappa_smlt   !    surface melting
         END WHERE
      ENDIF
      !
      ! thicknesses
      DO ji = 1, npti
         ! ice thickness
         IF( h_i_1d(ji) > 0._wp ) THEN
            zh_i  (ji) = MAX( zh_min , h_i_1d(ji) ) * r1_nlay_i ! set a minimum thickness for conduction
            z1_h_i(ji) = 1._wp / zh_i(ji)                       !       it must be very small
         ELSE
            zh_i  (ji) = 0._wp
            z1_h_i(ji) = 0._wp
         ENDIF
         ! snow thickness
         IF( h_s_1d(ji) > 0._wp ) THEN
            zh_s  (ji) = MAX( zh_min , h_s_1d(ji) ) * r1_nlay_s ! set a minimum thickness for conduction
            z1_h_s(ji) = 1._wp / zh_s(ji)                       !       it must be very small
            isnow (ji) = 1._wp
         ELSE
            zh_s  (ji) = 0._wp
            z1_h_s(ji) = 0._wp
            isnow (ji) = 0._wp
         ENDIF
         ! for Met-Office
         IF( h_s_1d(ji) < zh_min ) THEN
            isnow_comb(ji) = h_s_1d(ji) / zh_min
         ELSE
            isnow_comb(ji) = 1._wp
         ENDIF
      END DO
      ! clem: we should apply correction on snow thickness to take into account snow fraction
      !       it must be a distribution, so it is a bit complicated
      !
      ! Store initial temperatures and non solar heat fluxes
      IF( k_cnd == np_cnd_OFF .OR. k_cnd == np_cnd_EMU ) THEN
         ztsub      (1:npti) = t_su_1d(1:npti)                          ! surface temperature at iteration n-1
         ztsuold    (1:npti) = t_su_1d(1:npti)                          ! surface temperature initial value
         t_su_1d    (1:npti) = MIN( t_su_1d(1:npti), rt0 - ztsu_err )   ! required to leave the choice between melting or not
         zdqns_ice_b(1:npti) = dqns_ice_1d(1:npti)                      ! derivative of incoming nonsolar flux
         zqns_ice_b (1:npti) = qns_ice_1d(1:npti)                       ! store previous qns_ice_1d value
         !
      ENDIF
      !
      ztsold (1:npti,:) = t_s_1d(1:npti,:)   ! Old snow temperature
      ztiold (1:npti,:) = t_i_1d(1:npti,:)   ! Old ice temperature

      !-------------
      ! 2) Radiation
      !-------------
      ! --- Transmission/absorption of solar radiation in the ice --- !
      zradtr_s(1:npti,0) = qtr_ice_top_1d(1:npti)
      DO jk = 1, nlay_s
         DO ji = 1, npti
            !                             ! radiation transmitted below the layer-th snow layer
            zradtr_s(ji,jk) = zradtr_s(ji,0) * EXP( - zraext_s(ji) * MAX( 0._wp, zh_s(ji) * REAL(jk) - zhs_ssl ) )
            !                             ! radiation absorbed by the layer-th snow layer
            zradab_s(ji,jk) = zradtr_s(ji,jk-1) - zradtr_s(ji,jk)
         END DO
      END DO
      !
      zradtr_i(1:npti,0) = zradtr_s(1:npti,nlay_s) * za_s_fra(1:npti) + qtr_ice_top_1d(1:npti) * ( 1._wp - za_s_fra(1:npti) )
      DO jk = 1, nlay_i
         DO ji = 1, npti
            !                             ! radiation transmitted below the layer-th ice layer
            zradtr_i(ji,jk) =           za_s_fra(ji)   * zradtr_s(ji,nlay_s)                       &   ! part covered by snow
               &                                       * EXP( - rn_kappa_i * MAX( 0._wp, zh_i(ji) * REAL(jk) - zh_min  ) ) &
               &            + ( 1._wp - za_s_fra(ji) ) * qtr_ice_top_1d(ji)                        &   ! part snow free
               &                                       * EXP( - rn_kappa_i * MAX( 0._wp, zh_i(ji) * REAL(jk) - zhi_ssl ) )
            !                             ! radiation absorbed by the layer-th ice layer
            zradab_i(ji,jk) = zradtr_i(ji,jk-1) - zradtr_i(ji,jk)
         END DO
      END DO
      !
      qtr_ice_bot_1d(1:npti) = zradtr_i(1:npti,nlay_i)   ! record radiation transmitted below the ice
      !
      iconv = 0          ! number of iterations
      !
      l_T_converged(:) = .FALSE.
      ! Convergence calculated until all sub-domain grid points have converged
      ! Calculations keep going for all grid points until sub-domain convergence (vectorisation optimisation)
      ! but values are not taken into account (results independant of MPI partitioning)
      !
      !                                                                            !============================!
      DO WHILE ( ( .NOT. ALL (l_T_converged(1:npti)) ) .AND. iconv < iconv_max )   ! Iterative procedure begins !
         !                                                                         !============================!
         iconv = iconv + 1
         !
         ztib(1:npti,:) = t_i_1d(1:npti,:)
         ztsb(1:npti,:) = t_s_1d(1:npti,:)
         !
         !--------------------------------
         ! 3) Sea ice thermal conductivity
         !--------------------------------
         IF( ln_cndi_U64 ) THEN         !-- Untersteiner (1964) formula: k = k0 + beta.S/T
            !
            DO ji = 1, npti
               ztcond_i_cp(ji,0)      = rcnd_i + zbeta * sz_i_1d(ji,1)      / MIN( -epsi10, t_i_1d(ji,1) - rt0 )
               ztcond_i_cp(ji,nlay_i) = rcnd_i + zbeta * sz_i_1d(ji,nlay_i) / MIN( -epsi10, t_bo_1d(ji)  - rt0 )
            END DO
            DO jk = 1, nlay_i-1
               DO ji = 1, npti
                  ztcond_i_cp(ji,jk) = rcnd_i + zbeta * 0.5_wp * ( sz_i_1d(ji,jk) + sz_i_1d(ji,jk+1) ) /  &
                     &                    MIN( -epsi10, 0.5_wp * (  t_i_1d(ji,jk) +  t_i_1d(ji,jk+1) ) - rt0 )
               END DO
            END DO
            !
         ELSEIF( ln_cndi_P07 ) THEN     !-- Pringle et al formula: k = k0 + beta1.S/T - beta2.T
            !
            DO ji = 1, npti
               ztcond_i_cp(ji,0)      = rcnd_i + 0.09_wp  *  sz_i_1d(ji,1)      / MIN( -epsi10, t_i_1d(ji,1) - rt0 )  &
                  &                            - 0.011_wp * ( t_i_1d(ji,1) - rt0 )
               ztcond_i_cp(ji,nlay_i) = rcnd_i + 0.09_wp  *  sz_i_1d(ji,nlay_i) / MIN( -epsi10, t_bo_1d(ji)  - rt0 )  &
                  &                            - 0.011_wp * ( t_bo_1d(ji) - rt0 )
            END DO
            DO jk = 1, nlay_i-1
               DO ji = 1, npti
                  ztcond_i_cp(ji,jk) = rcnd_i + 0.09_wp  *   0.5_wp * ( sz_i_1d(ji,jk) + sz_i_1d(ji,jk+1) ) /       &
                     &                         MIN( -epsi10, 0.5_wp * (  t_i_1d(ji,jk) +  t_i_1d(ji,jk+1) ) - rt0 ) &
                     &                        - 0.011_wp * ( 0.5_wp * (  t_i_1d(ji,jk) +  t_i_1d(ji,jk+1) ) - rt0 )
               END DO
            END DO
            !
         ENDIF

         ! Variable used after iterations
         ! Value must be frozen after convergence for MPP independance reason
         DO ji = 1, npti
            IF ( .NOT. l_T_converged(ji) ) &
               ztcond_i(ji,:) = MAX( zkimin, ztcond_i_cp(ji,:) )
         END DO
         !
         !--- G(he) : enhancement of thermal conductivity in mono-category case
         ! Computation of effective thermal conductivity G(h)
         ! Used in mono-category case only to simulate an ITD implicitly
         ! Fichefet and Morales Maqueda, JGR 1997
         zghe(1:npti) = 1._wp
         !
         IF( ln_virtual_itd ) THEN
            !
            zepsilon = 0.1_wp
            DO ji = 1, npti
               zcnd_i = SUM( ztcond_i(ji,:) ) / REAL( nlay_i+1, wp )                                ! Mean sea ice thermal conductivity
               zhe = ( rn_cnd_s * h_i_1d(ji) + zcnd_i * h_s_1d(ji) ) / ( rn_cnd_s + zcnd_i )        ! Effective thickness he (zhe)
               IF( zhe >=  zepsilon * 0.5_wp * EXP(1._wp) )  &
                  &   zghe(ji) = MIN( 2._wp, 0.5_wp * ( 1._wp + LOG( 2._wp * zhe / zepsilon ) ) )   ! G(he)
            END DO
            !
         ENDIF
         !
         !-----------------
         ! 4) kappa factors
         !-----------------
         !--- Snow
         ! Variable used after iterations
         ! Value must be frozen after convergence for MPP independance reason
         DO jk = 0, nlay_s-1
            DO ji = 1, npti
               IF ( .NOT. l_T_converged(ji) ) &
                  zkappa_s(ji,jk) = zghe(ji) * rn_cnd_s * z1_h_s(ji)
            END DO
         END DO
         DO ji = 1, npti   ! Snow-ice interface
            IF ( .NOT. l_T_converged(ji) ) &
               zkappa_s(ji,nlay_s) = isnow(ji) * zghe(ji) * rn_cnd_s * ztcond_i(ji,0) &
                  &                            / ( 0.5_wp * ( ztcond_i(ji,0) * zh_s(ji) + rn_cnd_s * zh_i(ji) ) )
         END DO

         !--- Ice
         ! Variable used after iterations
         ! Value must be frozen after convergence for MPP independance reason
         DO jk = 0, nlay_i
            DO ji = 1, npti
               IF ( .NOT. l_T_converged(ji) ) &
                  zkappa_i(ji,jk) = zghe(ji) * ztcond_i(ji,jk) * z1_h_i(ji)
            END DO
         END DO
         DO ji = 1, npti   ! Snow-ice interface
            IF ( .NOT. l_T_converged(ji) ) THEN
               ! Calculate combined surface snow and ice conductivity to pass through the coupler (met-office)
               zkappa_comb(ji) = isnow_comb(ji) * zkappa_s(ji,0) + ( 1._wp - isnow_comb(ji) ) * zkappa_i(ji,0)
               ! If there is snow then use the same snow-ice interface conductivity for the top layer of ice
               IF( h_s_1d(ji) > 0._wp )   zkappa_i(ji,0) = zkappa_s(ji,nlay_s)
           ENDIF
         END DO
         !
         !--------------------------------------
         ! 5) Sea ice specific heat, eta factors
         !--------------------------------------
         DO jk = 1, nlay_i
            DO ji = 1, npti
               zcpi = rcpi + zgamma * sz_i_1d(ji,jk) / MAX( ( t_i_1d(ji,jk) - rt0 ) * ( ztiold(ji,jk) - rt0 ), epsi10 )
               zeta_i(ji,jk) = rDt_ice * r1_rhoi * z1_h_i(ji) / zcpi
            END DO
         END DO

         DO jk = 1, nlay_s
            DO ji = 1, npti
               zeta_s(ji,jk) = rDt_ice * r1_rhos * r1_rcpi * z1_h_s(ji)
            END DO
         END DO
         !
         !----------------------------------------!
         !                                        !
         !   Conduction flux is off or emulated   !
         !                                        !
         !----------------------------------------!
         !
         IF( k_cnd == np_cnd_OFF .OR. k_cnd == np_cnd_EMU ) THEN
            !
            ! ==> The original BL99 temperature computation is used
            !       (with qsr_ice, qns_ice and dqns_ice as inputs)
            !
            !----------------------------
            ! 6) surface flux computation
            !----------------------------
            ! update of the non solar flux according to the update in T_su
            DO ji = 1, npti
               ! Variable used after iterations
               ! Value must be frozen after convergence for MPP independance reason
               IF ( .NOT. l_T_converged(ji) ) &
                  qns_ice_1d(ji) = qns_ice_1d(ji) + dqns_ice_1d(ji) * ( t_su_1d(ji) - ztsub(ji) )
            END DO

            DO ji = 1, npti
               zfnet(ji) = qsr_ice_1d(ji) - qtr_ice_top_1d(ji) + qns_ice_1d(ji) ! net heat flux = net - transmitted solar + non solar
            END DO
            !
            !----------------------------
            ! 7) tridiagonal system terms
            !----------------------------
            ! layer denotes the number of the layer in the snow or in the ice
            ! jm denotes the reference number of the equation in the tridiagonal
            ! system, terms of tridiagonal system are indexed as following :
            ! 1 is subdiagonal term, 2 is diagonal and 3 is superdiagonal one

            ! ice interior terms (top equation has the same form as the others)
            ztrid   (1:npti,:,:) = 0._wp
            zindterm(1:npti,:)   = 0._wp
            zindtbis(1:npti,:)   = 0._wp
            zdiagbis(1:npti,:)   = 0._wp

            DO jm = nlay_s + 2, nlay_s + nlay_i
               DO ji = 1, npti
                  jk = jm - nlay_s - 1
                  ztrid   (ji,jm,1) =       - zeta_i(ji,jk) *   zkappa_i(ji,jk-1)
                  ztrid   (ji,jm,2) = 1._wp + zeta_i(ji,jk) * ( zkappa_i(ji,jk-1) + zkappa_i(ji,jk) )
                  ztrid   (ji,jm,3) =       - zeta_i(ji,jk) *                       zkappa_i(ji,jk)
                  zindterm(ji,jm)   = ztiold(ji,jk) + zeta_i(ji,jk) * zradab_i(ji,jk)
               END DO
            END DO

            jm =  nlay_s + nlay_i + 1
            DO ji = 1, npti
               ! ice bottom term
               ztrid   (ji,jm,1) =       - zeta_i(ji,nlay_i) *   zkappa_i(ji,nlay_i-1)
               ztrid   (ji,jm,2) = 1._wp + zeta_i(ji,nlay_i) * ( zkappa_i(ji,nlay_i-1) + zkappa_i(ji,nlay_i) * zg1 )
               ztrid   (ji,jm,3) = 0._wp
               zindterm(ji,jm)   = ztiold(ji,nlay_i) + zeta_i(ji,nlay_i) *  &
                  &              ( zradab_i(ji,nlay_i) + zkappa_i(ji,nlay_i) * zg1 * t_bo_1d(ji) )
            END DO

            DO ji = 1, npti
               !                               !---------------------!
               IF( h_s_1d(ji) > 0._wp ) THEN   !  snow-covered cells !
                  !                            !---------------------!
                  ! snow interior terms (bottom equation has the same form as the others)
                  DO jm = 3, nlay_s + 1
                     jk = jm - 1
                     ztrid   (ji,jm,1) =       - zeta_s(ji,jk) *   zkappa_s(ji,jk-1)
                     ztrid   (ji,jm,2) = 1._wp + zeta_s(ji,jk) * ( zkappa_s(ji,jk-1) + zkappa_s(ji,jk) )
                     ztrid   (ji,jm,3) =       - zeta_s(ji,jk) *                       zkappa_s(ji,jk)
                     zindterm(ji,jm)   = ztsold(ji,jk) + zeta_s(ji,jk) * zradab_s(ji,jk)
                  END DO

                  ! case of only one layer in the ice (ice equation is altered)
                  IF( nlay_i == 1 ) THEN
                     ztrid   (ji,nlay_s+2,3) = 0._wp
                     zindterm(ji,nlay_s+2)   = zindterm(ji,nlay_s+2) + zeta_i(ji,1) * zkappa_i(ji,1) * t_bo_1d(ji)
                  ENDIF

                  IF( t_su_1d(ji) < rt0 ) THEN   !--  case 1 : no surface melting

                     jm_min(ji) = 1
                     jm_max(ji) = nlay_i + nlay_s + 1

                     ! surface equation
                     ztrid   (ji,1,1) = 0._wp
                     ztrid   (ji,1,2) = zdqns_ice_b(ji) - zg1s * zkappa_s(ji,0)
                     ztrid   (ji,1,3) =                   zg1s * zkappa_s(ji,0)
                     zindterm(ji,1)   = zdqns_ice_b(ji) * t_su_1d(ji) - zfnet(ji)

                     ! first layer of snow equation
                     ztrid   (ji,2,1) =       - zeta_s(ji,1) *                    zkappa_s(ji,0) * zg1s
                     ztrid   (ji,2,2) = 1._wp + zeta_s(ji,1) * ( zkappa_s(ji,1) + zkappa_s(ji,0) * zg1s )
                     ztrid   (ji,2,3) =       - zeta_s(ji,1) *   zkappa_s(ji,1)
                     zindterm(ji,2)   = ztsold(ji,1) + zeta_s(ji,1) * zradab_s(ji,1)

                  ELSE                            !--  case 2 : surface is melting
                     !
                     jm_min(ji) = 2
                     jm_max(ji) = nlay_i + nlay_s + 1

                     ! first layer of snow equation
                     ztrid   (ji,2,1) = 0._wp
                     ztrid   (ji,2,2) = 1._wp + zeta_s(ji,1) * ( zkappa_s(ji,1) + zkappa_s(ji,0) * zg1s )
                     ztrid   (ji,2,3) =       - zeta_s(ji,1) *   zkappa_s(ji,1)
                     zindterm(ji,2)   = ztsold(ji,1) + zeta_s(ji,1) * ( zradab_s(ji,1) + zkappa_s(ji,0) * zg1s * t_su_1d(ji) )
                  ENDIF
                  !                            !---------------------!
               ELSE                            ! cells without snow  !
                  !                            !---------------------!
                  !
                  IF( t_su_1d(ji) < rt0 ) THEN   !--  case 1 : no surface melting
                     !
                     jm_min(ji) = nlay_s + 1
                     jm_max(ji) = nlay_i + nlay_s + 1

                     ! surface equation
                     ztrid   (ji,jm_min(ji),1) = 0._wp
                     ztrid   (ji,jm_min(ji),2) = zdqns_ice_b(ji) - zkappa_i(ji,0) * zg1
                     ztrid   (ji,jm_min(ji),3) =                   zkappa_i(ji,0) * zg1
                     zindterm(ji,jm_min(ji))   = zdqns_ice_b(ji) * t_su_1d(ji) - zfnet(ji)

                     ! first layer of ice equation
                     ztrid   (ji,jm_min(ji)+1,1) =       - zeta_i(ji,1) *                    zkappa_i(ji,0) * zg1
                     ztrid   (ji,jm_min(ji)+1,2) = 1._wp + zeta_i(ji,1) * ( zkappa_i(ji,1) + zkappa_i(ji,0) * zg1 )
                     ztrid   (ji,jm_min(ji)+1,3) =       - zeta_i(ji,1) *   zkappa_i(ji,1)
                     zindterm(ji,jm_min(ji)+1)   = ztiold(ji,1) + zeta_i(ji,1) * zradab_i(ji,1)

                     ! case of only one layer in the ice (surface & ice equations are altered)
                     IF( nlay_i == 1 ) THEN
                        ztrid   (ji,jm_min(ji),1)   = 0._wp
                        ztrid   (ji,jm_min(ji),2)   = zdqns_ice_b(ji)      -   zkappa_i(ji,0) * 2._wp
                        ztrid   (ji,jm_min(ji),3)   =                          zkappa_i(ji,0) * 2._wp
                        ztrid   (ji,jm_min(ji)+1,1) =       - zeta_i(ji,1) *   zkappa_i(ji,0) * 2._wp
                        ztrid   (ji,jm_min(ji)+1,2) = 1._wp + zeta_i(ji,1) * ( zkappa_i(ji,0) * 2._wp + zkappa_i(ji,1) )
                        ztrid   (ji,jm_min(ji)+1,3) = 0._wp
                        zindterm(ji,jm_min(ji)+1)   = ztiold(ji,1) + zeta_i(ji,1) * (zradab_i(ji,1) + zkappa_i(ji,1) * t_bo_1d(ji))
                     ENDIF

                  ELSE                            !--  case 2 : surface is melting

                     jm_min(ji) = nlay_s + 2
                     jm_max(ji) = nlay_i + nlay_s + 1

                     ! first layer of ice equation
                     ztrid   (ji,jm_min(ji),1) = 0._wp
                     ztrid   (ji,jm_min(ji),2) = 1._wp + zeta_i(ji,1) * ( zkappa_i(ji,1) + zkappa_i(ji,0) * zg1 )
                     ztrid   (ji,jm_min(ji),3) =       - zeta_i(ji,1) *   zkappa_i(ji,1)
                     zindterm(ji,jm_min(ji))   = ztiold(ji,1) + zeta_i(ji,1) * (zradab_i(ji,1) + zkappa_i(ji,0) * zg1 * t_su_1d(ji))

                     ! case of only one layer in the ice (surface & ice equations are altered)
                     IF( nlay_i == 1 ) THEN
                        ztrid   (ji,jm_min(ji),1) = 0._wp
                        ztrid   (ji,jm_min(ji),2) = 1._wp + zeta_i(ji,1) * ( zkappa_i(ji,0) * 2._wp + zkappa_i(ji,1) )
                        ztrid   (ji,jm_min(ji),3) = 0._wp
                        zindterm(ji,jm_min(ji))   = ztiold(ji,1) + zeta_i(ji,1) * ( zradab_i(ji,1) + zkappa_i(ji,1) * t_bo_1d(ji) ) &
                           &                      + t_su_1d(ji) * zeta_i(ji,1) * zkappa_i(ji,0) * 2._wp
                     ENDIF

                  ENDIF
               ENDIF
               !
               zindtbis(ji,jm_min(ji)) = zindterm(ji,jm_min(ji))
               zdiagbis(ji,jm_min(ji)) = ztrid   (ji,jm_min(ji),2)
               !
            END DO
            !
            !------------------------------
            ! 8) tridiagonal system solving
            !------------------------------
            ! Solve the tridiagonal system with Gauss elimination method.
            ! Thomas algorithm, from Computational fluid Dynamics, J.D. ANDERSON, McGraw-Hill 1984
!!$            jm_maxt = 0
!!$            jm_mint = nlay_i+5
!!$            DO ji = 1, npti
!!$               jm_mint = MIN(jm_min(ji),jm_mint)
!!$               jm_maxt = MAX(jm_max(ji),jm_maxt)
!!$            END DO
!!$            !!clem SNWLAY => check why LIM1D does not get this loop. Is nlay_i+5 correct?
!!$
!!$            DO jk = jm_mint+1, jm_maxt
!!$               DO ji = 1, npti
!!$                  jm = MIN(MAX(jm_min(ji)+1,jk),jm_max(ji))
!!$                  zdiagbis(ji,jm) = ztrid   (ji,jm,2) - ztrid(ji,jm,1) * ztrid   (ji,jm-1,3) / zdiagbis(ji,jm-1)
!!$                  zindtbis(ji,jm) = zindterm(ji,jm  ) - ztrid(ji,jm,1) * zindtbis(ji,jm-1  ) / zdiagbis(ji,jm-1)
!!$               END DO
!!$            END DO
            ! clem: maybe one should find a way to reverse this loop for mpi performance
            DO ji = 1, npti
               jm_mint = jm_min(ji)
               jm_maxt = jm_max(ji)
               DO jm = jm_mint+1, jm_maxt
                  zdiagbis(ji,jm) = ztrid   (ji,jm,2) - ztrid(ji,jm,1) * ztrid   (ji,jm-1,3) / zdiagbis(ji,jm-1)
                  zindtbis(ji,jm) = zindterm(ji,jm  ) - ztrid(ji,jm,1) * zindtbis(ji,jm-1  ) / zdiagbis(ji,jm-1)
               END DO
            END DO

            ! ice temperatures
            DO ji = 1, npti
               ! Variable used after iterations
               ! Value must be frozen after convergence for MPP independance reason
               IF ( .NOT. l_T_converged(ji) ) &
                  t_i_1d(ji,nlay_i) = zindtbis(ji,jm_max(ji)) / zdiagbis(ji,jm_max(ji))
            END DO

            DO jm = nlay_i + nlay_s, nlay_s + 2, -1
               DO ji = 1, npti
                  jk = jm - nlay_s - 1
                  IF ( .NOT. l_T_converged(ji) ) &
                     t_i_1d(ji,jk) = ( zindtbis(ji,jm) - ztrid(ji,jm,3) * t_i_1d(ji,jk+1) ) / zdiagbis(ji,jm)
               END DO
            END DO

            ! snow temperatures
            DO ji = 1, npti
               ! Variables used after iterations
               ! Value must be frozen after convergence for MPP independance reason
               IF ( .NOT. l_T_converged(ji) .AND. h_s_1d(ji) > 0._wp ) &
                  &   t_s_1d(ji,nlay_s) = ( zindtbis(ji,nlay_s+1) - ztrid(ji,nlay_s+1,3) * t_i_1d(ji,1) ) / zdiagbis(ji,nlay_s+1)
            END DO
            !!clem SNWLAY
            DO jm = nlay_s, 2, -1
               DO ji = 1, npti
                  jk = jm - 1
                  IF ( .NOT. l_T_converged(ji) .AND. h_s_1d(ji) > 0._wp ) &
                     &   t_s_1d(ji,jk) = ( zindtbis(ji,jm) - ztrid(ji,jm,3) * t_s_1d(ji,jk+1) ) / zdiagbis(ji,jm)
               END DO
            END DO

            ! surface temperature
            DO ji = 1, npti
               IF( .NOT. l_T_converged(ji) ) THEN
                  ztsub(ji) = t_su_1d(ji)
                  IF( t_su_1d(ji) < rt0 ) THEN
                     t_su_1d(ji) = ( zindtbis(ji,jm_min(ji)) - ztrid(ji,jm_min(ji),3) *  &
                        &          ( isnow(ji) * t_s_1d(ji,1) + ( 1._wp - isnow(ji) ) * t_i_1d(ji,1) ) ) / zdiagbis(ji,jm_min(ji))
                  ENDIF
               ENDIF
            END DO
            !
            !--------------------------------------------------------------
            ! 9) Has the scheme converged?, end of the iterative procedure
            !--------------------------------------------------------------
            ! check that nowhere it has started to melt
            ! zdti_max is a measure of error, it has to be under zdti_bnd

            DO ji = 1, npti

               zdti_max = 0._wp

               IF ( .NOT. l_T_converged(ji) ) THEN

                  t_su_1d(ji) = MAX( MIN( t_su_1d(ji) , rt0 ) , rt0 - 100._wp )
                  zdti_max    = MAX( zdti_max, ABS( t_su_1d(ji) - ztsub(ji) ) )

                  IF( h_s_1d(ji) > 0._wp ) THEN
                     DO jk = 1, nlay_s
                        t_s_1d(ji,jk) = MAX( MIN( t_s_1d(ji,jk), rt0 ), rt0 - 100._wp )
                        zdti_max      = MAX ( zdti_max , ABS( t_s_1d(ji,jk) - ztsb(ji,jk) ) )
                     END DO
                  ENDIF

                  DO jk = 1, nlay_i
                     ztmelts       = -rTmlt * sz_i_1d(ji,jk) + rt0
                     t_i_1d(ji,jk) =  MAX( MIN( t_i_1d(ji,jk), ztmelts ), rt0 - 100._wp )
                     zdti_max      =  MAX( zdti_max, ABS( t_i_1d(ji,jk) - ztib(ji,jk) ) )
                  END DO

                  ! convergence test
                  IF( ln_zdf_chkcvg ) THEN
                     tice_cvgerr_1d(ji) = zdti_max
                     tice_cvgstp_1d(ji) = REAL(iconv)
                  ENDIF

                  IF( zdti_max < zdti_bnd )   l_T_converged(ji) = .TRUE.

               ENDIF

            END DO

         !----------------------------------------!
         !                                        !
         !      Conduction flux is on             !
         !                                        !
         !----------------------------------------!
         !
         ELSEIF( k_cnd == np_cnd_ON ) THEN
            !
            ! ==> we use a modified BL99 solver with conduction flux (qcn_ice) as forcing term
            !
            !----------------------------
            ! 7) tridiagonal system terms
            !----------------------------
            ! layer denotes the number of the layer in the snow or in the ice
            ! jm denotes the reference number of the equation in the tridiagonal
            ! system, terms of tridiagonal system are indexed as following :
            ! 1 is subdiagonal term, 2 is diagonal and 3 is superdiagonal one

            ! ice interior terms (top equation has the same form as the others)
            ztrid   (1:npti,:,:) = 0._wp
            zindterm(1:npti,:)   = 0._wp
            zindtbis(1:npti,:)   = 0._wp
            zdiagbis(1:npti,:)   = 0._wp

            DO jm = nlay_s + 2, nlay_s + nlay_i
               DO ji = 1, npti
                  jk = jm - nlay_s - 1
                  ztrid   (ji,jm,1) =       - zeta_i(ji,jk) *   zkappa_i(ji,jk-1)
                  ztrid   (ji,jm,2) = 1._wp + zeta_i(ji,jk) * ( zkappa_i(ji,jk-1) + zkappa_i(ji,jk) )
                  ztrid   (ji,jm,3) =       - zeta_i(ji,jk) *                       zkappa_i(ji,jk)
                  zindterm(ji,jm)   = ztiold(ji,jk) + zeta_i(ji,jk) * zradab_i(ji,jk)
               END DO
            ENDDO

            jm =  nlay_s + nlay_i + 1
            DO ji = 1, npti
               ! ice bottom term
               ztrid   (ji,jm,1) =       - zeta_i(ji,nlay_i) *   zkappa_i(ji,nlay_i-1)
               ztrid   (ji,jm,2) = 1._wp + zeta_i(ji,nlay_i) * ( zkappa_i(ji,nlay_i-1) + zkappa_i(ji,nlay_i) * zg1 )
               ztrid   (ji,jm,3) = 0._wp
               zindterm(ji,jm)   = ztiold(ji,nlay_i) + zeta_i(ji,nlay_i) *  &
                  &              ( zradab_i(ji,nlay_i) + zkappa_i(ji,nlay_i) * zg1 * t_bo_1d(ji) )
            ENDDO

            DO ji = 1, npti
               !                               !---------------------!
               IF( h_s_1d(ji) > 0._wp ) THEN   !  snow-covered cells !
                  !                            !---------------------!
                  ! snow interior terms (bottom equation has the same form as the others)
                  DO jm = 3, nlay_s + 1
                     jk = jm - 1
                     ztrid   (ji,jm,1) =       - zeta_s(ji,jk) *   zkappa_s(ji,jk-1)
                     ztrid   (ji,jm,2) = 1._wp + zeta_s(ji,jk) * ( zkappa_s(ji,jk-1) + zkappa_s(ji,jk) )
                     ztrid   (ji,jm,3) =       - zeta_s(ji,jk) *                       zkappa_s(ji,jk)
                     zindterm(ji,jm)   = ztsold(ji,jk) + zeta_s(ji,jk) * zradab_s(ji,jk)
                  END DO

                  ! case of only one layer in the ice (ice equation is altered)
                  IF ( nlay_i == 1 ) THEN
                     ztrid   (ji,nlay_s+2,3) = 0._wp
                     zindterm(ji,nlay_s+2)   = zindterm(ji,nlay_s+2) + zeta_i(ji,1) * zkappa_i(ji,1) * t_bo_1d(ji)
                  ENDIF

                  jm_min(ji) = 2
                  jm_max(ji) = nlay_i + nlay_s + 1

                  ! first layer of snow equation
                  ztrid   (ji,2,1) = 0._wp
                  ztrid   (ji,2,2) = 1._wp + zeta_s(ji,1) * zkappa_s(ji,1)
                  ztrid   (ji,2,3) =       - zeta_s(ji,1) * zkappa_s(ji,1)
                  zindterm(ji,2)   = ztsold(ji,1) + zeta_s(ji,1) * ( zradab_s(ji,1) + qcn_ice_1d(ji) )

                  !                            !---------------------!
               ELSE                            ! cells without snow  !
                  !                            !---------------------!
                  jm_min(ji) = nlay_s + 2
                  jm_max(ji) = nlay_i + nlay_s + 1

                  ! first layer of ice equation
                  ztrid   (ji,jm_min(ji),1) = 0._wp
                  ztrid   (ji,jm_min(ji),2) = 1._wp + zeta_i(ji,1) * zkappa_i(ji,1)
                  ztrid   (ji,jm_min(ji),3) =       - zeta_i(ji,1) * zkappa_i(ji,1)
                  zindterm(ji,jm_min(ji))   = ztiold(ji,1) + zeta_i(ji,1) * ( zradab_i(ji,1) + qcn_ice_1d(ji) )

                  ! case of only one layer in the ice (surface & ice equations are altered)
                  IF( nlay_i == 1 ) THEN
                     ztrid   (ji,jm_min(ji),1) = 0._wp
                     ztrid   (ji,jm_min(ji),2) = 1._wp + zeta_i(ji,1) * zkappa_i(ji,1)
                     ztrid   (ji,jm_min(ji),3) = 0._wp
                     zindterm(ji,jm_min(ji))   = ztiold(ji,1) + zeta_i(ji,1) *  &
                        &                                     ( zradab_i(ji,1) + zkappa_i(ji,1) * t_bo_1d(ji) + qcn_ice_1d(ji) )
                  ENDIF

               ENDIF
               !
               zindtbis(ji,jm_min(ji)) = zindterm(ji,jm_min(ji))
               zdiagbis(ji,jm_min(ji)) = ztrid   (ji,jm_min(ji),2)
               !
            END DO
            !
            !------------------------------
            ! 8) tridiagonal system solving
            !------------------------------
            ! Solve the tridiagonal system with Gauss elimination method.
            ! Thomas algorithm, from Computational fluid Dynamics, J.D. ANDERSON, McGraw-Hill 1984
!!$            jm_maxt = 0
!!$            jm_mint = nlay_i+5
!!$            DO ji = 1, npti
!!$               jm_mint = MIN(jm_min(ji),jm_mint)
!!$               jm_maxt = MAX(jm_max(ji),jm_maxt)
!!$            END DO
!!$
!!$            DO jk = jm_mint+1, jm_maxt
!!$               DO ji = 1, npti
!!$                  jm = MIN(MAX(jm_min(ji)+1,jk),jm_max(ji))
!!$                  zdiagbis(ji,jm) = ztrid   (ji,jm,2) - ztrid(ji,jm,1) * ztrid   (ji,jm-1,3) / zdiagbis(ji,jm-1)
!!$                  zindtbis(ji,jm) = zindterm(ji,jm)   - ztrid(ji,jm,1) * zindtbis(ji,jm-1)   / zdiagbis(ji,jm-1)
!!$               END DO
!!$            END DO
            ! clem: maybe one should find a way to reverse this loop for mpi performance
            DO ji = 1, npti
               jm_mint = jm_min(ji)
               jm_maxt = jm_max(ji)
               DO jm = jm_mint+1, jm_maxt
                  zdiagbis(ji,jm) = ztrid   (ji,jm,2) - ztrid(ji,jm,1) * ztrid   (ji,jm-1,3) / zdiagbis(ji,jm-1)
                  zindtbis(ji,jm) = zindterm(ji,jm  ) - ztrid(ji,jm,1) * zindtbis(ji,jm-1  ) / zdiagbis(ji,jm-1)
               END DO
            END DO

            ! ice temperatures
            DO ji = 1, npti
               ! Variable used after iterations
               ! Value must be frozen after convergence for MPP independance reason
               IF ( .NOT. l_T_converged(ji) ) &
                  t_i_1d(ji,nlay_i) = zindtbis(ji,jm_max(ji)) / zdiagbis(ji,jm_max(ji))
            END DO

            DO jm = nlay_i + nlay_s, nlay_s + 2, -1
               DO ji = 1, npti
                  IF ( .NOT. l_T_converged(ji) ) THEN
                     jk = jm - nlay_s - 1
                     t_i_1d(ji,jk) = ( zindtbis(ji,jm) - ztrid(ji,jm,3) * t_i_1d(ji,jk+1) ) / zdiagbis(ji,jm)
                  ENDIF
               END DO
            END DO

            ! snow temperatures
            DO ji = 1, npti
               ! Variables used after iterations
               ! Value must be frozen after convergence for MPP independance reason
               IF ( .NOT. l_T_converged(ji) .AND. h_s_1d(ji) > 0._wp ) &
                  &   t_s_1d(ji,nlay_s) = ( zindtbis(ji,nlay_s+1) - ztrid(ji,nlay_s+1,3) * t_i_1d(ji,1) ) / zdiagbis(ji,nlay_s+1)
            END DO
            !!clem SNWLAY
            DO jm = nlay_s, 2, -1
               DO ji = 1, npti
                  jk = jm - 1
                  IF ( .NOT. l_T_converged(ji) .AND. h_s_1d(ji) > 0._wp ) &
                     &   t_s_1d(ji,jk) = ( zindtbis(ji,jm) - ztrid(ji,jm,3) * t_s_1d(ji,jk+1) ) / zdiagbis(ji,jm)
               END DO
            END DO
            !
            !--------------------------------------------------------------
            ! 9) Has the scheme converged?, end of the iterative procedure
            !--------------------------------------------------------------
            ! check that nowhere it has started to melt
            ! zdti_max is a measure of error, it has to be under zdti_bnd

            DO ji = 1, npti

               zdti_max = 0._wp

               IF ( .NOT. l_T_converged(ji) ) THEN

                  IF( h_s_1d(ji) > 0._wp ) THEN
                     DO jk = 1, nlay_s
                        t_s_1d(ji,jk) = MAX( MIN( t_s_1d(ji,jk), rt0 ), rt0 - 100._wp )
                        zdti_max      = MAX ( zdti_max , ABS( t_s_1d(ji,jk) - ztsb(ji,jk) ) )
                     END DO
                  ENDIF

                  DO jk = 1, nlay_i
                     ztmelts       = -rTmlt * sz_i_1d(ji,jk) + rt0
                     t_i_1d(ji,jk) =  MAX( MIN( t_i_1d(ji,jk), ztmelts ), rt0 - 100._wp )
                     zdti_max      =  MAX ( zdti_max, ABS( t_i_1d(ji,jk) - ztib(ji,jk) ) )
                  END DO

                  ! convergence test
                  IF( ln_zdf_chkcvg ) THEN
                     tice_cvgerr_1d(ji) = zdti_max
                     tice_cvgstp_1d(ji) = REAL(iconv)
                  ENDIF

                  IF( zdti_max < zdti_bnd )   l_T_converged(ji) = .TRUE.

               ENDIF

            END DO

         ENDIF ! k_cnd

      END DO  ! End of the do while iterative procedure
      !
      !-----------------------------
      ! 10) Fluxes at the interfaces
      !-----------------------------
      !
      ! --- calculate conduction fluxes (positive downward)
      !     bottom ice conduction flux
      DO ji = 1, npti
         qcn_ice_bot_1d(ji) = - zkappa_i(ji,nlay_i) * zg1 * ( t_bo_1d(ji ) - t_i_1d (ji,nlay_i) )
      END DO
      !     surface ice conduction flux
      IF( k_cnd == np_cnd_OFF .OR. k_cnd == np_cnd_EMU ) THEN
         !
         DO ji = 1, npti
            qcn_ice_top_1d(ji) = -           isnow(ji)   * zkappa_s(ji,0) * zg1s * ( t_s_1d(ji,1) - t_su_1d(ji) ) &
               &                 - ( 1._wp - isnow(ji) ) * zkappa_i(ji,0) * zg1  * ( t_i_1d(ji,1) - t_su_1d(ji) )
         END DO
         !
      ELSEIF( k_cnd == np_cnd_ON ) THEN
         !
         DO ji = 1, npti
            qcn_ice_top_1d(ji) = qcn_ice_1d(ji)
         END DO
         !
      ENDIF
      !     surface ice temperature
      IF( k_cnd == np_cnd_ON .AND. ln_cndemulate ) THEN
         !
         DO ji = 1, npti
            t_su_1d(ji) = ( qcn_ice_top_1d(ji) +          isnow(ji)   * zkappa_s(ji,0) * zg1s * t_s_1d(ji,1) + &
               &                                ( 1._wp - isnow(ji) ) * zkappa_i(ji,0) * zg1  * t_i_1d(ji,1) ) &
               &          / MAX( epsi10, isnow(ji) * zkappa_s(ji,0) * zg1s + ( 1._wp - isnow(ji) ) * zkappa_i(ji,0) * zg1 )
            t_su_1d(ji) = MAX( MIN( t_su_1d(ji), rt0 ), rt0 - 100._wp )  ! cap t_su
         END DO
         !
      ENDIF
      !
      ! --- Diagnose the heat loss due to changing non-solar / conduction flux --- !
      !
      IF( k_cnd == np_cnd_OFF .OR. k_cnd == np_cnd_EMU ) THEN
         !
         DO ji = 1, npti
            hfx_err_dif_1d(ji) = hfx_err_dif_1d(ji) - ( qns_ice_1d(ji) - zqns_ice_b(ji) ) * a_i_1d(ji)
         END DO
         !
      ENDIF
      !
      ! --- Diagnose the heat loss due to non-fully converged temperature solution (should not be above 10-4 W-m2) --- !
      !
      IF( k_cnd == np_cnd_OFF .OR. k_cnd == np_cnd_ON ) THEN

         CALL ice_var_enthalpy

         ! zhfx_err = correction on the diagnosed heat flux due to non-convergence of the algorithm used to solve heat equation
         DO ji = 1, npti
            zdq = - zq_ini(ji) + ( SUM( e_i_1d(ji,1:nlay_i) ) * h_i_1d(ji) * r1_nlay_i +  &
               &                   SUM( e_s_1d(ji,1:nlay_s) ) * h_s_1d(ji) * r1_nlay_s )

            IF( k_cnd == np_cnd_OFF ) THEN

               IF( t_su_1d(ji) < rt0 ) THEN  ! case T_su < 0degC
                  zhfx_err = ( qns_ice_1d(ji)     + qsr_ice_1d(ji)     - zradtr_i(ji,nlay_i) - qcn_ice_bot_1d(ji)  &
                     &       + zdq * r1_Dt_ice ) * a_i_1d(ji)
               ELSE                          ! case T_su = 0degC
                  zhfx_err = ( qcn_ice_top_1d(ji) + qtr_ice_top_1d(ji) - zradtr_i(ji,nlay_i) - qcn_ice_bot_1d(ji)  &
                     &       + zdq * r1_Dt_ice ) * a_i_1d(ji)
               ENDIF

            ELSEIF( k_cnd == np_cnd_ON ) THEN

               zhfx_err    = ( qcn_ice_top_1d(ji) + qtr_ice_top_1d(ji) - zradtr_i(ji,nlay_i) - qcn_ice_bot_1d(ji)  &
                  &          + zdq * r1_Dt_ice ) * a_i_1d(ji)

            ENDIF
            !
            ! total heat sink to be sent to the ocean
            hfx_err_dif_1d(ji) = hfx_err_dif_1d(ji) + zhfx_err
            !
            ! hfx_dif = Heat flux diagnostic of sensible heat used to warm/cool ice in W.m-2
            hfx_dif_1d(ji) = hfx_dif_1d(ji) - zdq * r1_Dt_ice * a_i_1d(ji)
            !
         END DO
         !
      ENDIF
      !
      !--------------------------------------------------------------------
      ! 11) reset inner snow and ice temperatures, update conduction fluxes
      !--------------------------------------------------------------------
      ! effective conductivity and 1st layer temperature (needed by Met Office)
      ! this is a conductivity at mid-layer, hence the factor 2
      DO ji = 1, npti
         IF( h_i_1d(ji) >= zhi_ssl ) THEN
            cnd_ice_1d(ji) = 2._wp * zkappa_comb(ji)
            !!cnd_ice_1d(ji) = 2._wp * zkappa_i(ji,0)
         ELSE
            cnd_ice_1d(ji) = 2._wp * ztcond_i(ji,0) / zhi_ssl ! cnd_ice is capped by: cond_i/zhi_ssl
         ENDIF
         t1_ice_1d(ji) = isnow_comb(ji) * t_s_1d(ji,1) + ( 1._wp - isnow_comb(ji) ) * t_i_1d(ji,1)
      END DO
      !
      IF( k_cnd == np_cnd_EMU ) THEN
         ! Restore temperatures to their initial values
         t_s_1d    (1:npti,:) = ztsold        (1:npti,:)
         t_i_1d    (1:npti,:) = ztiold        (1:npti,:)
         qcn_ice_1d(1:npti)   = qcn_ice_top_1d(1:npti)
      ENDIF
      !
      ! --- SIMIP diagnostics
      !
      DO ji = 1, npti
         !--- Snow-ice interfacial temperature (diagnostic SIMIP)
         IF( h_s_1d(ji) >= zhs_ssl ) THEN
            t_si_1d(ji) = (   rn_cnd_s       * h_i_1d(ji) * r1_nlay_i * t_s_1d(ji,nlay_s)   &
               &            + ztcond_i(ji,1) * h_s_1d(ji) * r1_nlay_s * t_i_1d(ji,1)      ) &
               &          / ( rn_cnd_s       * h_i_1d(ji) * r1_nlay_i &
               &            + ztcond_i(ji,1) * h_s_1d(ji) * r1_nlay_s )
         ELSE
            t_si_1d(ji) = t_su_1d(ji)
         ENDIF
      END DO
      !
   END SUBROUTINE ice_thd_zdf_BL99

#else
   !!----------------------------------------------------------------------
   !!   Default option       Dummy Module             No SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE icethd_zdf_BL99
