MODULE icestp
   !!======================================================================
   !!                       ***  MODULE  icestp  ***
   !! sea ice : Master routine for all the sea ice model
   !!=====================================================================
   !!
   !! The sea ice model SI3 (Sea Ice modelling Integrated Initiative),
   !!                        aka Sea Ice cube for its nickname
   !!
   !!    is originally based on LIM3, developed in Louvain-la-Neuve by:
   !!       * Martin Vancoppenolle (UCL-ASTR, Belgium)
   !!       * Sylvain Bouillon (UCL-ASTR, Belgium)
   !!       * Miguel Angel Morales Maqueda (NOC-L, UK)
   !!      thanks to valuable earlier work by
   !!       * Thierry Fichefet
   !!       * Hugues Goosse
   !!      thanks also to the following persons who contributed
   !!       * Gurvan Madec, Claude Talandier, Christian Ethe (LOCEAN, France)
   !!       * Xavier Fettweis (UCL-ASTR), Ralph Timmermann (AWI, Germany)
   !!       * Bill Lipscomb (LANL), Cecilia Bitz (UWa) and Elisabeth Hunke (LANL), USA.
   !!
   !! SI3 has been made possible by a handful of persons who met as working group
   !!      (from France, Belgium, UK and Italy)
   !!    * Clement Rousset, Martin Vancoppenolle & Gurvan Madec (LOCEAN, France)
   !!    * Matthieu Chevalier & David Salas (Meteo France, France)
   !!    * Gilles Garric (Mercator Ocean, France)
   !!    * Thierry Fichefet & Francois Massonnet (UCL, Belgium)
   !!    * Ed Blockley & Jeff Ridley (Met Office, UK)
   !!    * Danny Feltham & David Schroeder (CPOM, UK)
   !!    * Yevgeny Aksenov (NOC, UK)
   !!    * Paul Holland (BAS, UK)
   !!    * Dorotea Iovino (CMCC, Italy)
   !!======================================================================
   !! History :  4.0  !  2018     (C. Rousset)      Original code SI3
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_stp       : sea-ice model time-stepping and update ocean SBC over ice-covered area
   !!   ice_init      : initialize sea-ice
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE c1d            ! 1D vertical configuration
   USE ice            ! sea-ice: variables
   USE ice1D          ! sea-ice: thermodynamical 1D variables
   !
   USE phycst         ! Define parameters for the routines
   USE eosbn2         ! equation of state
   USE sbc_oce        ! Surface boundary condition: ocean fields
   USE sbc_ice        ! Surface boundary condition: ice   fields
   !
   USE icesbc         ! sea-ice: Surface boundary conditions
   USE icedyn         ! sea-ice: dynamics
   USE icethd         ! sea-ice: thermodynamics
   USE iceupdate      ! sea-ice: sea surface boundary condition update
   USE icedia         ! sea-ice: budget diagnostics
   USE icewri         ! sea-ice: outputs
   USE icerst         ! sea-ice: restarts
   USE icevar         ! sea-ice: operations
   USE icectl         ! sea-ice: control
   USE iceistate      ! sea-ice: initial state
   USE iceitd         ! sea-ice: remapping thickness distribution
   USE icealb         ! sea-ice: albedo
   !
   USE bdy_oce , ONLY : ln_bdy   ! flag for bdy
   USE bdyice         ! unstructured open boundary data for sea-ice
# if defined key_agrif
   USE agrif_ice
   USE agrif_ice_interp
# endif
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)
   USE timing         ! Timing
   USE prtctl         ! Print control

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_stp    ! called by sbcmod.F90
   PUBLIC   ice_init   ! called by sbcmod.F90

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/ICE 4.0 , NEMO Consortium (2018)
   !! $Id: icestp.F90 15023 2021-06-18 14:35:25Z gsamson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_stp( kt, Kbb, Kmm, ksbc )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE ice_stp  ***
      !!
      !! ** Purpose :   sea-ice model time-stepping and update ocean surface
      !!              boundary condition over ice-covered area
      !!
      !! ** Method  :   ice model time stepping
      !!              - call the ice dynamics routine
      !!              - call the ice advection/diffusion routine
      !!              - call the ice thermodynamics routine
      !!              - call the routine that computes mass and
      !!                heat fluxes at the ice/ocean interface
      !!              - save the outputs
      !!              - save the outputs for restart when necessary
      !!
      !! ** Action  : - time evolution of the LIM sea-ice model
      !!              - update all sbc variables below sea-ice:
      !!                utau, vtau, taum, wndm, qns , qsr, emp , sfx
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt       ! ocean time step
      INTEGER, INTENT(in) ::   Kbb, Kmm ! ocean time level indices
      INTEGER, INTENT(in) ::   ksbc     ! flux formulation (user defined, bulk, or Pure Coupled)
      !
      INTEGER ::   jl   ! dummy loop index
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('icestp')
      !
      !                                      !-----------------------!
      IF( MOD( kt-1, nn_fsbc ) == 0 ) THEN   ! --- Ice time step --- !
         !                                   !-----------------------!
         !
         kt_ice = kt                              ! -- Ice model time step
         !
         u_oce(:,:) = ssu_m(:,:)                  ! -- mean surface ocean current
         v_oce(:,:) = ssv_m(:,:)
         !
         CALL eos_fzp( sss_m(:,:) , t_bo(:,:) )   ! -- freezing temperature [Kelvin] (set to rt0 over land)
         t_bo(:,:) = ( t_bo(:,:) + rt0 ) * tmask(:,:,1) + rt0 * ( 1._wp - tmask(:,:,1) )
         !
         !                          !==  AGRIF Parent to Child  ==!
#if defined key_agrif
         !                              ! nbstep_ice ranges from 1 to the nb of child ocean steps inside one parent ice step
         IF( .NOT. Agrif_Root() )       nbstep_ice = MOD( nbstep_ice, Agrif_irhot() * Agrif_Parent(nn_fsbc) / nn_fsbc ) + 1
         !                              ! these calls must remain here for restartability purposes
                                        CALL agrif_interp_ice( 'T' )
                                        CALL agrif_interp_ice( 'U' )
                                        CALL agrif_interp_ice( 'V' )
#endif
                                        CALL store_fields             ! Store now ice values
         !
         !------------------------------------------------!
         ! --- Dynamical coupling with the atmosphere --- !
         !------------------------------------------------!
         ! It provides the following fields used in sea ice model:
         !    utau_ice, vtau_ice = surface ice stress [N/m2]
         !------------------------------------------------!
                                        CALL ice_sbc_tau( kt, ksbc, utau_ice, vtau_ice )
         !-------------------------------------!
         ! --- ice dynamics and advection  --- !
         !-------------------------------------!
                                        CALL diag_set0                ! set diag of mass, heat and salt fluxes to 0
                                        CALL ice_rst_opn( kt )        ! Open Ice restart file (if necessary)
         !
         IF( ln_icedyn .AND. .NOT.ln_c1d )   &
            &                           CALL ice_dyn( kt, Kmm )       ! -- Ice dynamics
         !
                                        CALL diag_trends( 1 )         ! record dyn trends
         !
         !                          !==  lateral boundary conditions  ==!
         IF( ln_icethd .AND. ln_bdy )   CALL bdy_ice( kt )            ! -- bdy ice thermo
         !
         !                          !==  previous lead fraction and ice volume for flux calculations
                                        CALL ice_var_glo2eqv          ! h_i and h_s for ice albedo calculation
                                        CALL ice_var_agg(1)           ! at_i for coupling
                                        CALL store_fields             ! Store now ice values
         !
         !------------------------------------------------------!
         ! --- Thermodynamical coupling with the atmosphere --- !
         !------------------------------------------------------!
         ! It provides the following fields used in sea ice model:
         !    emp_oce , emp_ice    = E-P over ocean and sea ice                    [Kg/m2/s]
         !    sprecip              = solid precipitation                           [Kg/m2/s]
         !    evap_ice             = sublimation                                   [Kg/m2/s]
         !    qsr_tot , qns_tot    = solar & non solar heat flux (total)           [W/m2]
         !    qsr_ice , qns_ice    = solar & non solar heat flux over ice          [W/m2]
         !    dqns_ice             = non solar  heat sensistivity                  [W/m2]
         !    qemp_oce, qemp_ice,  = sensible heat (associated with evap & precip) [W/m2]
         !    qprec_ice, qevap_ice
         !------------------------------------------------------!
                                        CALL ice_sbc_flx( kt, ksbc )
         !----------------------------!
         ! --- ice thermodynamics --- !
         !----------------------------!
         IF( ln_icethd )                CALL ice_thd( kt )            ! -- Ice thermodynamics
         !
                                        CALL diag_trends( 2 )         ! record thermo trends
                                        CALL ice_var_glo2eqv          ! necessary calls (at least for coupling)
                                        CALL ice_var_agg( 2 )         ! necessary calls (at least for coupling)
         !
                                        CALL ice_update_flx( kt )     ! -- Update ocean surface mass, heat and salt fluxes
         !
         IF( ln_icediahsb )             CALL ice_dia( kt )            ! -- Diagnostics outputs
         !
         IF( ln_icediachk )             CALL ice_drift_wri( kt )      ! -- Diagnostics outputs for conservation
         !
                                        CALL ice_wri( kt )            ! -- Ice outputs
         !
         IF( lrst_ice )                 CALL ice_rst_write( kt )      ! -- Ice restart file
         !
         IF( ln_icectl )                CALL ice_ctl( kt )            ! -- Control checks
         !
      ENDIF   ! End sea-ice time step only

      !-------------------------!
      ! --- Ocean time step --- !
      !-------------------------!
      CALL ice_update_tau( kt, uu(:,:,1,Kbb), vv(:,:,1,Kbb) )         ! -- update surface ocean stresses
!!gm   remark, the ocean-ice stress is not saved in ice diag call above .....  find a solution!!!
      !
      IF( ln_timing )   CALL timing_stop('icestp')
      !
   END SUBROUTINE ice_stp


   SUBROUTINE ice_init( Kbb, Kmm, Kaa )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ice_init  ***
      !!
      !! ** purpose :   Initialize sea-ice parameters
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) :: Kbb, Kmm, Kaa
      !
      INTEGER ::   ierr
      !!----------------------------------------------------------------------
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'Sea Ice Model: SI3 (Sea Ice modelling Integrated Initiative)'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~'
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'ice_init: Arrays allocation & Initialization of all routines & init state'
      IF(lwp) WRITE(numout,*) '~~~~~~~~'
      !
      !                                ! Load the reference and configuration namelist files and open namelist output file
      CALL load_nml( numnam_ice_ref, 'namelist_ice_ref',    numout, lwm )
      CALL load_nml( numnam_ice_cfg, 'namelist_ice_cfg',    numout, lwm )
      IF(lwm) CALL ctl_opn( numoni , 'output.namelist.ice', 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp, 1 )
      !
      CALL par_init                ! set some ice run parameters
      !
#if defined key_agrif
      CALL Agrif_Declare_Var_ice  !  "      "   "   "      "  Sea ice
#endif
      !
      !                                ! Allocate the ice arrays (sbc_ice already allocated in sbc_init)
      ierr =        ice_alloc        ()      ! ice variables
      ierr = ierr + sbc_ice_alloc    ()      ! surface boundary conditions
      ierr = ierr + ice1D_alloc      ()      ! thermodynamics
      !
      CALL mpp_sum( 'icestp', ierr )
      IF( ierr /= 0 )   CALL ctl_stop('STOP', 'ice_init : unable to allocate ice arrays')
      !
      !                                ! set max concentration in both hemispheres
      WHERE( gphit(:,:) > 0._wp )   ;   rn_amax_2d(:,:) = rn_amax_n  ! NH
      ELSEWHERE                     ;   rn_amax_2d(:,:) = rn_amax_s  ! SH
      END WHERE
      !
      CALL diag_set0                   ! set diag of mass, heat and salt fluxes to 0: needed for Agrif child grids
      !
      CALL ice_itd_init                ! ice thickness distribution initialization
      !
      CALL ice_thd_init                ! set ice thermodynics parameters (clem: important to call it first for melt ponds)
      !
      CALL ice_sbc_init                ! set ice-ocean and ice-atm. coupling parameters
      !
      CALL ice_istate_init             ! Initial sea-ice state
      IF ( ln_rstart .OR. nn_iceini_file == 2 ) THEN
         CALL ice_rst_read( Kbb, Kmm, Kaa )         ! start from a restart file
      ELSE
         CALL ice_istate( nit000, Kbb, Kmm, Kaa )   ! start from rest or read a file
      ENDIF
      CALL ice_var_glo2eqv
      CALL ice_var_agg(2)
      !
      CALL ice_dyn_init                ! set ice dynamics parameters
      !
      CALL ice_update_init             ! ice surface boundary condition
      !
      CALL ice_alb_init                ! ice surface albedo
      !
      CALL ice_dia_init                ! initialization for diags
      !
      CALL ice_drift_init              ! initialization for diags of conservation
      !
      fr_i  (:,:)   = at_i(:,:)        ! initialisation of sea-ice fraction
      tn_ice(:,:,:) = t_su(:,:,:)      ! initialisation of surface temp for coupled simu
      !
      IF( ln_rstart )  THEN
         CALL iom_close( numrir )  ! close input ice restart file
         IF(lrxios) CALL iom_context_finalize(      cr_icerst_cxt         )
      ENDIF
      !
   END SUBROUTINE ice_init


   SUBROUTINE par_init
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE par_init ***
      !!
      !! ** Purpose :   Definition generic parameters for ice model
      !!
      !! ** Method  :   Read namelist and check the parameter
      !!                values called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist nampar
      !!-------------------------------------------------------------------
      INTEGER  ::   ios                 ! Local integer
      !!
      NAMELIST/nampar/ jpl, nlay_i, nlay_s, ln_virtual_itd, ln_icedyn, ln_icethd, rn_amax_n, rn_amax_s,  &
         &             cn_icerst_in, cn_icerst_indir, cn_icerst_out, cn_icerst_outdir
      !!-------------------------------------------------------------------
      !
      READ  ( numnam_ice_ref, nampar, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nampar in reference namelist' )
      READ  ( numnam_ice_cfg, nampar, IOSTAT = ios, ERR = 902 )
902   IF( ios > 0 )   CALL ctl_nam ( ios , 'nampar in configuration namelist' )
      IF(lwm) WRITE( numoni, nampar )
      !
      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) '   par_init: ice parameters shared among all the routines'
         WRITE(numout,*) '   ~~~~~~~~'
         WRITE(numout,*) '      Namelist nampar: '
         WRITE(numout,*) '         number of ice  categories                           jpl       = ', jpl
         WRITE(numout,*) '         number of ice  layers                               nlay_i    = ', nlay_i
         WRITE(numout,*) '         number of snow layers                               nlay_s    = ', nlay_s
         WRITE(numout,*) '         virtual ITD param for jpl=1 (T) or not (F)     ln_virtual_itd = ', ln_virtual_itd
         WRITE(numout,*) '         Ice dynamics       (T) or not (F)                   ln_icedyn = ', ln_icedyn
         WRITE(numout,*) '         Ice thermodynamics (T) or not (F)                   ln_icethd = ', ln_icethd
         WRITE(numout,*) '         maximum ice concentration for NH                              = ', rn_amax_n
         WRITE(numout,*) '         maximum ice concentration for SH                              = ', rn_amax_s
      ENDIF
      !                                        !--- change max ice concentration for roundoff errors
      rn_amax_n = MIN( rn_amax_n, 1._wp - epsi10 )
      rn_amax_s = MIN( rn_amax_s, 1._wp - epsi10 )
      !                                        !--- check consistency
      IF ( jpl > 1 .AND. ln_virtual_itd ) THEN
         ln_virtual_itd = .FALSE.
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   ln_virtual_itd forced to false as jpl>1, no need with multiple categories to emulate them'
      ENDIF
      !
      IF( ln_cpl .AND. nn_cats_cpl /= 1 .AND. nn_cats_cpl /= jpl ) THEN
         CALL ctl_stop( 'STOP', 'par_init: in coupled mode, nn_cats_cpl should be either 1 or jpl' )
      ENDIF
      !
      rDt_ice   = REAL(nn_fsbc) * rn_Dt          !--- sea-ice timestep and its inverse
      r1_Dt_ice = 1._wp / rDt_ice
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '      ice timestep rDt_ice = nn_fsbc*rn_Dt = ', rDt_ice
      !
      r1_nlay_i = 1._wp / REAL( nlay_i, wp )   !--- inverse of nlay_i and nlay_s
      r1_nlay_s = 1._wp / REAL( nlay_s, wp )
      !
   END SUBROUTINE par_init


   SUBROUTINE store_fields
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE store_fields  ***
      !!
      !! ** purpose :  store ice variables at "before" time step
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jl      ! dummy loop index
      !!----------------------------------------------------------------------
      !
      a_i_b (:,:,:)   = a_i (:,:,:)     ! ice area
      v_i_b (:,:,:)   = v_i (:,:,:)     ! ice volume
      v_s_b (:,:,:)   = v_s (:,:,:)     ! snow volume
      v_ip_b(:,:,:)   = v_ip(:,:,:)     ! pond volume
      v_il_b(:,:,:)   = v_il(:,:,:)     ! pond lid volume
      sv_i_b(:,:,:)   = sv_i(:,:,:)     ! salt content
      e_s_b (:,:,:,:) = e_s (:,:,:,:)   ! snow thermal energy
      e_i_b (:,:,:,:) = e_i (:,:,:,:)   ! ice thermal energy
      WHERE( a_i_b(:,:,:) >= epsi20 )
         h_i_b(:,:,:) = v_i_b(:,:,:) / a_i_b(:,:,:)   ! ice thickness
         h_s_b(:,:,:) = v_s_b(:,:,:) / a_i_b(:,:,:)   ! snw thickness
      ELSEWHERE
         h_i_b(:,:,:) = 0._wp
         h_s_b(:,:,:) = 0._wp
      END WHERE
      !
      ! ice velocities & total concentration
      at_i_b(:,:)  = SUM( a_i_b(:,:,:), dim=3 )
      u_ice_b(:,:) = u_ice(:,:)
      v_ice_b(:,:) = v_ice(:,:)
      !
   END SUBROUTINE store_fields


   SUBROUTINE diag_set0
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE diag_set0  ***
      !!
      !! ** purpose :  set ice-ocean and ice-atm. fluxes to zeros at the beggining
      !!               of the time step
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jl      ! dummy loop index
      !!----------------------------------------------------------------------

      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )   ! needed for (at least) diag_adv_mass -> to be removed 
         sfx    (ji,jj) = 0._wp   ;
         sfx_bri(ji,jj) = 0._wp   ;   sfx_lam(ji,jj) = 0._wp
         sfx_sni(ji,jj) = 0._wp   ;   sfx_opw(ji,jj) = 0._wp
         sfx_bog(ji,jj) = 0._wp   ;   sfx_dyn(ji,jj) = 0._wp
         sfx_bom(ji,jj) = 0._wp   ;   sfx_sum(ji,jj) = 0._wp
         sfx_res(ji,jj) = 0._wp   ;   sfx_sub(ji,jj) = 0._wp
         !
         wfx_snw(ji,jj) = 0._wp   ;   wfx_ice(ji,jj) = 0._wp
         wfx_sni(ji,jj) = 0._wp   ;   wfx_opw(ji,jj) = 0._wp
         wfx_bog(ji,jj) = 0._wp   ;   wfx_dyn(ji,jj) = 0._wp
         wfx_bom(ji,jj) = 0._wp   ;   wfx_sum(ji,jj) = 0._wp
         wfx_res(ji,jj) = 0._wp   ;   wfx_sub(ji,jj) = 0._wp
         wfx_spr(ji,jj) = 0._wp   ;   wfx_lam(ji,jj) = 0._wp
         wfx_snw_dyn(ji,jj) = 0._wp ; wfx_snw_sum(ji,jj) = 0._wp
         wfx_snw_sub(ji,jj) = 0._wp ; wfx_ice_sub(ji,jj) = 0._wp
         wfx_snw_sni(ji,jj) = 0._wp
         wfx_pnd(ji,jj) = 0._wp

         hfx_thd(ji,jj) = 0._wp   ;
         hfx_snw(ji,jj) = 0._wp   ;   hfx_opw(ji,jj) = 0._wp
         hfx_bog(ji,jj) = 0._wp   ;   hfx_dyn(ji,jj) = 0._wp
         hfx_bom(ji,jj) = 0._wp   ;   hfx_sum(ji,jj) = 0._wp
         hfx_res(ji,jj) = 0._wp   ;   hfx_sub(ji,jj) = 0._wp
         hfx_spr(ji,jj) = 0._wp   ;   hfx_dif(ji,jj) = 0._wp
         hfx_err_dif(ji,jj) = 0._wp
         wfx_err_sub(ji,jj) = 0._wp
         !
         diag_heat(ji,jj) = 0._wp ;   diag_sice(ji,jj) = 0._wp
         diag_vice(ji,jj) = 0._wp ;   diag_vsnw(ji,jj) = 0._wp
         diag_aice(ji,jj) = 0._wp ;   diag_vpnd(ji,jj) = 0._wp

         tau_icebfr (ji,jj) = 0._wp   ! landfast ice param only (clem: important to keep the init here)
         qsb_ice_bot(ji,jj) = 0._wp   ! (needed if ln_icethd=F)

         fhld(ji,jj) = 0._wp   ! needed if ln_icethd=F

         ! for control checks (ln_icediachk)
         diag_trp_vi(ji,jj) = 0._wp   ;   diag_trp_vs(ji,jj) = 0._wp
         diag_trp_ei(ji,jj) = 0._wp   ;   diag_trp_es(ji,jj) = 0._wp
         diag_trp_sv(ji,jj) = 0._wp
         !
         diag_adv_mass(ji,jj) = 0._wp
         diag_adv_salt(ji,jj) = 0._wp
         diag_adv_heat(ji,jj) = 0._wp
      END_2D

      DO jl = 1, jpl
         DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
            ! SIMIP diagnostics
            t_si       (ji,jj,jl) = rt0     ! temp at the ice-snow interface
            qcn_ice_bot(ji,jj,jl) = 0._wp
            qcn_ice_top(ji,jj,jl) = 0._wp   ! conductive fluxes
            cnd_ice    (ji,jj,jl) = 0._wp   ! effective conductivity at the top of ice/snow (ln_cndflx=T)
            qcn_ice    (ji,jj,jl) = 0._wp   ! conductive flux (ln_cndflx=T & ln_cndemule=T)
            qtr_ice_bot(ji,jj,jl) = 0._wp   ! part of solar radiation transmitted through the ice needed at least for outputs
            qml_ice    (ji,jj,jl) = 0._wp   ! surface melt heat flux
            ! Melt pond surface melt diagnostics (mv - more efficient: grouped into one water volume flux)
            dh_i_sum_2d(ji,jj,jl) = 0._wp
            dh_s_mlt_2d(ji,jj,jl) = 0._wp
         END_2D
      ENDDO

   END SUBROUTINE diag_set0


   SUBROUTINE diag_trends( kn )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE diag_trends  ***
      !!
      !! ** purpose : diagnostics of the trends. Used for conservation purposes
      !!              and outputs
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kn    ! 1 = after dyn ; 2 = after thermo
      !!----------------------------------------------------------------------
      !
      ! --- trends of heat, salt, mass (used for conservation controls)
      IF( ln_icediachk .OR. iom_use('hfxdhc') ) THEN
         !
         diag_heat(:,:) = diag_heat(:,:) &
            &             - SUM(SUM( e_i (:,:,1:nlay_i,:) - e_i_b (:,:,1:nlay_i,:), dim=4 ), dim=3 ) * r1_Dt_ice &
            &             - SUM(SUM( e_s (:,:,1:nlay_s,:) - e_s_b (:,:,1:nlay_s,:), dim=4 ), dim=3 ) * r1_Dt_ice
         diag_sice(:,:) = diag_sice(:,:) &
            &             + SUM(     sv_i(:,:,:)          - sv_i_b(:,:,:)                  , dim=3 ) * r1_Dt_ice * rhoi
         diag_vice(:,:) = diag_vice(:,:) &
            &             + SUM(     v_i (:,:,:)          - v_i_b (:,:,:)                  , dim=3 ) * r1_Dt_ice * rhoi
         diag_vsnw(:,:) = diag_vsnw(:,:) &
            &             + SUM(     v_s (:,:,:)          - v_s_b (:,:,:)                  , dim=3 ) * r1_Dt_ice * rhos
         diag_vpnd(:,:) = diag_vpnd(:,:) &
            &             + SUM(     v_ip + v_il          - v_ip_b - v_il_b                , dim=3 ) * r1_Dt_ice * rhow
         !
         IF( kn == 2 )    CALL iom_put ( 'hfxdhc' , diag_heat )   ! output of heat trend
         !
      ENDIF
      !
      ! --- trends of concentration (used for simip outputs)
      IF( iom_use('afxdyn') .OR. iom_use('afxthd') .OR. iom_use('afxtot') ) THEN
         !
         diag_aice(:,:) = diag_aice(:,:) + SUM( a_i(:,:,:) - a_i_b(:,:,:), dim=3 ) * r1_Dt_ice
         !
         IF( kn == 1 )   CALL iom_put( 'afxdyn' , diag_aice )                                           ! dyn trend
         IF( kn == 2 )   CALL iom_put( 'afxthd' , SUM( a_i(:,:,:) - a_i_b(:,:,:), dim=3 ) * r1_Dt_ice ) ! thermo trend
         IF( kn == 2 )   CALL iom_put( 'afxtot' , diag_aice )                                           ! total trend
         !
      ENDIF
      !
   END SUBROUTINE diag_trends

#else
   !!----------------------------------------------------------------------
   !!   Default option           Dummy module         NO SI3 sea-ice model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE ice_stp ( kt, ksbc )     ! Dummy routine
      INTEGER, INTENT(in) ::   kt, ksbc
      WRITE(*,*) 'ice_stp: You should not have seen this print! error?', kt
   END SUBROUTINE ice_stp
   SUBROUTINE ice_init                 ! Dummy routine
   END SUBROUTINE ice_init
#endif

   !!======================================================================
END MODULE icestp
