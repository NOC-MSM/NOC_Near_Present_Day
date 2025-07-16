MODULE dynspg_ts

   !! Includes ROMS wd scheme with diagnostic outputs ; puu(:,:,:,Kmm) and puu(:,:,:,Krhs) updates are commented out ! 

   !!======================================================================
   !!                   ***  MODULE  dynspg_ts  ***
   !! Ocean dynamics:  surface pressure gradient trend, split-explicit scheme
   !!======================================================================
   !! History :   1.0  ! 2004-12  (L. Bessieres, G. Madec)  Original code
   !!              -   ! 2005-11  (V. Garnier, G. Madec)  optimization
   !!              -   ! 2006-08  (S. Masson)  distributed restart using iom
   !!             2.0  ! 2007-07  (D. Storkey) calls to BDY routines
   !!              -   ! 2008-01  (R. Benshila)  change averaging method
   !!             3.2  ! 2009-07  (R. Benshila, G. Madec) Complete revisit associated to vvl reactivation
   !!             3.3  ! 2010-09  (D. Storkey, E. O'Dea) update for BDY for Shelf configurations
   !!             3.3  ! 2011-03  (R. Benshila, R. Hordoir, P. Oddo) update calculation of ub_b
   !!             3.5  ! 2013-07  (J. Chanut) Switch to Forward-backward time stepping
   !!             3.6  ! 2013-11  (A. Coward) Update for z-tilde compatibility
   !!             3.7  ! 2015-11  (J. Chanut) free surface simplification
   !!              -   ! 2016-12  (G. Madec, E. Clementi) update for Stoke-Drift divergence
   !!             4.0  ! 2017-05  (G. Madec)  drag coef. defined at t-point (zdfdrg.F90)
   !!---------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_spg_ts     : compute surface pressure gradient trend using a time-splitting scheme 
   !!   dyn_spg_ts_init: initialisation of the time-splitting scheme
   !!   ts_wgt         : set time-splitting weights for temporal averaging (or not)
   !!   ts_rst         : read/write time-splitting fields in restart file
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! surface boundary condition: ocean
   USE isf_oce         ! ice shelf variable (fwfisf)
   USE zdf_oce         ! vertical physics: variables
   USE zdfdrg          ! vertical physics: top/bottom drag coef.
   USE sbcapr          ! surface boundary condition: atmospheric pressure
   USE dynadv    , ONLY: ln_dynadv_vec
   USE dynvor          ! vortivity scheme indicators
   USE phycst          ! physical constants
   USE dynvor          ! vorticity term
   USE wet_dry         ! wetting/drying flux limter
   USE bdy_oce         ! open boundary
   USE bdyvol          ! open boundary volume conservation
   USE bdytides        ! open boundary condition data
   USE bdydyn2d        ! open boundary conditions on barotropic variables
   USE tide_mod        !
   USE sbcwave         ! surface wave
#if defined key_agrif
   USE agrif_oce_interp ! agrif
   USE agrif_oce
#endif
#if defined key_asminc   
   USE asminc          ! Assimilation increment
#endif
   !
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distributed memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE prtctl          ! Print control
   USE iom             ! IOM library
   USE restart         ! only for lrst_oce
   USE trd_oce        ! trends: ocean variables
   USE trddyn         ! trend manager: dynamics

   USE iom   ! to remove
!jth IWD
   USE tide_mod
!
   IMPLICIT NONE
   PRIVATE

   PUBLIC dyn_spg_ts        ! called by dyn_spg 
   PUBLIC dyn_spg_ts_init   !    -    - dyn_spg_init

   !! Time filtered arrays at baroclinic time step:
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   un_adv , vn_adv   !: Advection vel. at "now" barocl. step
   !
   INTEGER, SAVE :: icycle      ! Number of barotropic sub-steps for each internal step nn_e <= 2.5 nn_e
   REAL(dp),SAVE :: rDt_e       ! Barotropic time step
   !
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:)    :: wgtbtp2! 1st & 2nd weights used in time filtering of barotropic fields
   REAL(dp), ALLOCATABLE, SAVE, DIMENSION(:)    :: wgtbtp1! 1st & 2nd weights used in time filtering of barotropic fields
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   zwz                ! ff_f/h at F points
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ftnw, ftne         ! triad of coriolis parameter
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ftsw, ftse         ! (only used with een vorticity scheme)
   REAL(wp), ALLOCATABLE,       DIMENSION(:,:) ::   zua_tmp, zva_tmp     ! RDP: temporary save of 2d momentum for WAD and BDY
!jth_IWD
   REAL(wp),PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:)   :: un_detide,vn_detide ! arrays for de-tided barotropic currents
   REAL(wp),PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:)   :: un_25,vn_25 !arrays for 25hr running means
!
   REAL(wp) ::   r1_12 = 1._wp / 12._wp   ! local ratios
   REAL(wp) ::   r1_8  = 0.125_wp         !
   REAL(wp) ::   r1_4  = 0.25_wp          !
   REAL(wp) ::   r1_2  = 0.5_wp           !

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "single_precision_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: dynspg_ts.F90 15489 2021-11-10 09:18:39Z jchanut $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION dyn_spg_ts_alloc()
      !!----------------------------------------------------------------------
      !!                  ***  routine dyn_spg_ts_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER :: ierr(3)
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !
      ALLOCATE( wgtbtp1(3*nn_e), wgtbtp2(3*nn_e), zwz(jpi,jpj), STAT=ierr(1) )
      IF( ln_dynvor_een .OR. ln_dynvor_eeT )   &
         &     ALLOCATE( ftnw(jpi,jpj) , ftne(jpi,jpj) , ftsw(jpi,jpj) , ftse(jpi,jpj), STAT=ierr(2)   )
         !
      ALLOCATE( un_adv(jpi,jpj), vn_adv(jpi,jpj)                    , STAT=ierr(3) )
      !
      dyn_spg_ts_alloc = MAXVAL( ierr(:) )
      !
      CALL mpp_sum( 'dynspg_ts', dyn_spg_ts_alloc )
      IF( dyn_spg_ts_alloc /= 0 )   CALL ctl_stop( 'STOP', 'dyn_spg_ts_alloc: failed to allocate arrays' )
      !
   END FUNCTION dyn_spg_ts_alloc


   SUBROUTINE dyn_spg_ts( kt, Kbb, Kmm, Krhs, puu, pvv, pssh, puu_b, pvv_b, Kaa )
      !!----------------------------------------------------------------------
      !!
      !! ** Purpose : - Compute the now trend due to the explicit time stepping
      !!              of the quasi-linear barotropic system, and add it to the
      !!              general momentum trend. 
      !!
      !! ** Method  : - split-explicit schem (time splitting) :
      !!      Barotropic variables are advanced from internal time steps
      !!      "n"   to "n+1" if ln_bt_fw=T
      !!      or from 
      !!      "n-1" to "n+1" if ln_bt_fw=F
      !!      thanks to a generalized forward-backward time stepping (see ref. below).
      !!
      !! ** Action :
      !!      -Update the filtered free surface at step "n+1"      : pssh(:,:,Kaa)
      !!      -Update filtered barotropic velocities at step "n+1" : puu_b(:,:,:,Kaa), vv_b(:,:,:,Kaa)
      !!      -Compute barotropic advective fluxes at step "n"     : un_adv, vn_adv
      !!      These are used to advect tracers and are compliant with discrete
      !!      continuity equation taken at the baroclinic time steps. This 
      !!      ensures tracers conservation.
      !!      - (puu(:,:,:,Krhs), pvv(:,:,:,Krhs)) momentum trend updated with barotropic component.
      !!
      !! References : Shchepetkin and McWilliams, Ocean Modelling, 2005. 
      !!---------------------------------------------------------------------
      INTEGER                             , INTENT( in )  ::  kt                  ! ocean time-step index
      INTEGER                             , INTENT( in )  ::  Kbb, Kmm, Krhs, Kaa ! ocean time level indices
      REAL(dp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) ::  puu, pvv            ! ocean velocities and RHS of momentum equation
      REAL(wp), DIMENSION(jpi,jpj,jpt)    , INTENT(inout)  :: puu_b, pvv_b! SSH and barotropic velocities at main time levels
      REAL(dp), DIMENSION(jpi,jpj,jpt)    , INTENT(inout)  :: pssh! SSH and barotropic velocities at main time levels
      !
      INTEGER  ::   ji, jj, jk, jn        ! dummy loop indices
      LOGICAL  ::   ll_fw_start           ! =T : forward integration 
      LOGICAL  ::   ll_init               ! =T : special startup of 2d equations
      INTEGER  ::   noffset               ! local integers  : time offset for bdy update
      REAL(wp) ::   r1_Dt_b, z1_hu, z1_hv          ! local scalars
      REAL(wp)  :: za0, za2, za3!   -      -
      REAL(dp)  :: za1!   -      -
      REAL(wp) ::   zztmp, zldg               !   -      -
      REAL(wp)  :: zhu_bck, zhv_bck!   -      -
      REAL(dp)  :: zhdiv!   -      -
      REAL(wp) ::   zun_save, zvn_save              !   -      -
      REAL(wp), DIMENSION(jpi,jpj) :: zu_trd, zu_frc, zu_spg
      REAL(wp), DIMENSION(jpi,jpj) :: zv_trd, zv_frc, zv_spg
      REAL(wp), DIMENSION(jpi,jpj) :: zsshu_a, zhup2_e, zhtp2_e
      REAL(wp), DIMENSION(jpi,jpj) :: zsshv_a, zhvp2_e, zsshp2_e
      REAL(wp), DIMENSION(jpi,jpj) :: zCdU_u, zCdU_v   ! top/bottom stress at u- & v-points
      REAL(wp), DIMENSION(jpi,jpj)  :: zhV! fluxes
      REAL(dp), DIMENSION(jpi,jpj)  :: zhU! fluxes
!!st#if defined key_qco 
!!st      REAL(wp), DIMENSION(jpi, jpj, jpk) :: ze3u, ze3v
!!st#endif
      !
      REAL(wp) ::   zwdramp                     ! local scalar - only used if ln_wd_dl = .True. 

      INTEGER  :: iwdg, jwdg, kwdg   ! short-hand values for the indices of the output point

      REAL(wp) ::   zepsilon, zgamma            !   -      -
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: ztwdmask, zuwdmask, zvwdmask ! ROMS wetting and drying masks at t,u,v points
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zuwdav2, zvwdav2    ! averages over the sub-steps of zuwdmask and zvwdmask
      REAL(wp) ::   zt0substep !   Time of day at the beginning of the time substep
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zspgtrdu, zspgtrdv, zpvotrdu, zpvotrdv  ! SPG and PVO trends (if l_trddyn)
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: ztautrdu, ztautrdv, zbfrtrdu, zbfrtrdv  ! TAU and BFR trends (if l_trddyn)
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: ztfrtrdu, ztfrtrdv, ztottrdu, ztottrdv  ! TFR and TOT trends (if l_trddyn)
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zatmtrdu, zatmtrdv, zbdytrdu, zbdytrdv  ! ATM and BDY trends (if l_trddyn) 
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: ztidepotgradu, ztidepotgradv ! temp variables for tide potential gradients
      !!----------------------------------------------------------------------
      !
      !                                         !* Allocate temporary arrays
      IF( ln_wd_dl ) ALLOCATE( ztwdmask(jpi,jpj), zuwdmask(jpi,jpj), zvwdmask(jpi,jpj), zuwdav2(jpi,jpj), zvwdav2(jpi,jpj))
      !
      IF( l_trddyn ) THEN
          ALLOCATE( zspgtrdu(jpi,jpj), zspgtrdv(jpi,jpj), zpvotrdu(jpi,jpj), zpvotrdv(jpi,jpj), &
         &          ztautrdu(jpi,jpj), ztautrdv(jpi,jpj), zbfrtrdu(jpi,jpj), zbfrtrdv(jpi,jpj), &
         &          ztottrdu(jpi,jpj), ztottrdv(jpi,jpj) )
          IF( ln_isfcav .OR. ln_drgice_imp ) THEN
             ALLOCATE ( ztfrtrdu(jpi,jpj), ztfrtrdv(jpi,jpj) )
             ztfrtrdu(:,:) = 0._wp
             ztfrtrdv(:,:) = 0._wp
          ENDIF
          IF( ln_apr_dyn ) THEN
             ALLOCATE( zatmtrdu(jpi,jpj), zatmtrdv(jpi,jpj) )
          ENDIF
          IF( ln_bdy ) THEN
             ALLOCATE( zbdytrdu(jpi,jpj), zbdytrdv(jpi,jpj) )
             zbdytrdu(:,:) = 0._wp
             zbdytrdv(:,:) = 0._wp
          ENDIF
          zspgtrdu(:,:) = 0._wp
          zspgtrdv(:,:) = 0._wp
          zpvotrdu(:,:) = 0._wp
          zpvotrdv(:,:) = 0._wp
      ENDIF
      !
      zu_trd(:,:) = 0._wp
      zv_trd(:,:) = 0._wp
      zu_spg(:,:) = 0._wp
      zv_spg(:,:) = 0._wp
      !
      zwdramp = r_rn_wdmin1               ! simplest ramp 
!     zwdramp = 1._wp / (rn_wdmin2 - rn_wdmin1) ! more general ramp
      !                                         ! inverse of baroclinic time step 
      r1_Dt_b = 1._wp / rDt 
      !
      ll_init     = ln_bt_av                    ! if no time averaging, then no specific restart 
      ll_fw_start = .FALSE.
      !                                         ! time offset in steps for bdy data update
      IF( .NOT.ln_bt_fw ) THEN   ;   noffset = - nn_e
      ELSE                       ;   noffset =   0 
      ENDIF
      !
      IF( kt == nit000 ) THEN                   !* initialisation
         !
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_spg_ts : surface pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~   free surface with time splitting'
         IF(lwp) WRITE(numout,*)
         !
         IF( l_1st_euler )   ll_init=.TRUE.
         !
         IF( ln_bt_fw .OR. l_1st_euler ) THEN
            ll_fw_start =.TRUE.
            noffset     = 0
         ELSE
            ll_fw_start =.FALSE.
         ENDIF
         !                    ! Set averaging weights and cycle length:
         CALL ts_wgt( ln_bt_av, ll_fw_start, icycle, wgtbtp1, wgtbtp2 )
         !
      ELSEIF( kt == nit000 + 1 ) THEN           !* initialisation 2nd time-step
         !
         IF( .NOT.ln_bt_fw ) THEN
            ! If we did an Euler timestep on the first timestep we need to reset ll_fw_start
            ! and the averaging weights. We don't have an easy way of telling whether we did
            ! an Euler timestep on the first timestep (because l_1st_euler is reset to .false.
            ! at the end of the first timestep) so just do this in all cases. 
            ll_fw_start = .FALSE.
            CALL ts_wgt( ln_bt_av, ll_fw_start, icycle, wgtbtp1, wgtbtp2 )
         ENDIF
         !
      ENDIF
      !
      ! -----------------------------------------------------------------------------
      !  Phase 1 : Coupling between general trend and barotropic estimates (1st step)
      ! -----------------------------------------------------------------------------
      !      
      !
      !                                   !=  zu_frc =  1/H e3*d/dt(Ua)  =!  (Vertical mean of Ua, the 3D trends)
      !                                   !  ---------------------------  !
#if defined key_qco 
      zu_frc(:,:) = SUM( e3u_0(:,:,:  ) * puu(:,:,:,Krhs) * umask(:,:,:), DIM=3 ) * r1_hu_0(:,:)
      zv_frc(:,:) = SUM( e3v_0(:,:,:  ) * pvv(:,:,:,Krhs) * vmask(:,:,:), DIM=3 ) * r1_hv_0(:,:)
#else
      zu_frc(:,:) = SUM( e3u(:,:,:,Kmm) * puu(:,:,:,Krhs) * umask(:,:,:), DIM=3 ) * r1_hu(:,:,Kmm)
      zv_frc(:,:) = SUM( e3v(:,:,:,Kmm) * pvv(:,:,:,Krhs) * vmask(:,:,:), DIM=3 ) * r1_hv(:,:,Kmm)
#endif 
      !
      !
      !                                   !=  U(Krhs) => baroclinic trend  =!   (remove its vertical mean)
      DO jk = 1, jpkm1                    !  -----------------------------  !
         puu(:,:,jk,Krhs) = ( puu(:,:,jk,Krhs) - zu_frc(:,:) ) * umask(:,:,jk)
         pvv(:,:,jk,Krhs) = ( pvv(:,:,jk,Krhs) - zv_frc(:,:) ) * vmask(:,:,jk)
      END DO
      
!!gm  Question here when removing the Vertically integrated trends, we remove the vertically integrated NL trends on momentum....
!!gm  Is it correct to do so ?   I think so...
      
      !                                   !=  remove 2D Coriolis trend  =!
      !                                   !  --------------------------  !
      !
      IF( kt == nit000 .OR. .NOT. ln_linssh )   CALL dyn_cor_2D_init( Kmm )   ! Set zwz, the barotropic Coriolis force coefficient
      !                      ! recompute zwz = f/depth  at every time step for (.NOT.ln_linssh) as the water colomn height changes
      !
      zhU(:,:) = puu_b(:,:,Kmm) * hu(:,:,Kmm) * e2u(:,:)        ! now fluxes 
      zhV(:,:) = pvv_b(:,:,Kmm) * hv(:,:,Kmm) * e1v(:,:)        ! NB: FULL domain : put a value in last row and column
      !
      CALL dyn_cor_2d( ht(:,:), hu(:,:,Kmm), hv(:,:,Kmm), puu_b(:,:,Kmm), pvv_b(:,:,Kmm), zhU, zhV,  &   ! <<== in
         &                                                                          zu_trd, zv_trd   )   ! ==>> out
      IF( l_trddyn ) THEN
         ! send correction to baroclinic planetary vorticity trend to trd_dyn
         CALL trd_dyn( zu_trd, zv_trd, jpdyn_pvo_corr, kt, Kmm)
      ENDIF
      !
      !
      DO_2D( 0, 0, 0, 0 )                          ! Remove coriolis term (and possibly spg) from barotropic trend
         zu_frc(ji,jj) = zu_frc(ji,jj) - zu_trd(ji,jj) * ssumask(ji,jj)
         zv_frc(ji,jj) = zv_frc(ji,jj) - zv_trd(ji,jj) * ssvmask(ji,jj)
      END_2D
      !
      !                                   !=  Add bottom stress contribution from baroclinic velocities  =!
      !                                   !  -----------------------------------------------------------  !
      IF( l_trddyn ) THEN
         !
         ! Output constant forcing terms (excluding top and bottom stresses) as diagnostics.
         CALL trd_dyn( zu_frc, zv_frc, jpdyn_frc2d, kt, Kmm)
         !
         CALL dyn_drg_init(kt, Kbb, Kmm, puu, pvv, puu_b, pvv_b, &     ! also provide the barotropic drag coefficients
              &             zu_frc, zv_frc, zCdU_u, zCdU_v,   &
              &             ztfrtrdu, ztfrtrdv, zbfrtrdu, zbfrtrdv )
         !
      ELSE				     !* remove baroclinic drag AND provide the barotropic drag coefficients
         CALL dyn_drg_init(kt, Kbb, Kmm, puu, pvv, puu_b, pvv_b, zu_frc, zv_frc, zCdU_u, zCdU_v )
      ENDIF
      !
      !                                   !=  Add atmospheric pressure forcing  =!
      !                                   !  ----------------------------------  !
      IF( ln_apr_dyn ) THEN
         IF( l_trddyn ) THEN
            ! initialise fields for atmospheric pressure trends
            zatmtrdu(:,:) = zu_frc(:,:)
            zatmtrdv(:,:) = zv_frc(:,:)
         ENDIF
         IF( ln_bt_fw ) THEN                          ! FORWARD integration: use kt+1/2 pressure (NOW+1/2)
            DO_2D( 0, 0, 0, 0 )
               zu_frc(ji,jj) = zu_frc(ji,jj) + grav * (  ssh_ib (ji+1,jj  ) - ssh_ib (ji,jj) ) * r1_e1u(ji,jj)
               zv_frc(ji,jj) = zv_frc(ji,jj) + grav * (  ssh_ib (ji  ,jj+1) - ssh_ib (ji,jj) ) * r1_e2v(ji,jj)
            END_2D
         ELSE                                         ! CENTRED integration: use kt-1/2 + kt+1/2 pressure (NOW)
            zztmp = grav * r1_2
            DO_2D( 0, 0, 0, 0 )
               zu_frc(ji,jj) = zu_frc(ji,jj) + zztmp * (  ssh_ib (ji+1,jj  ) - ssh_ib (ji,jj)  &
                    &                                   + ssh_ibb(ji+1,jj  ) - ssh_ibb(ji,jj)  ) * r1_e1u(ji,jj)
               zv_frc(ji,jj) = zv_frc(ji,jj) + zztmp * (  ssh_ib (ji  ,jj+1) - ssh_ib (ji,jj)  &
                    &                                   + ssh_ibb(ji  ,jj+1) - ssh_ibb(ji,jj)  ) * r1_e2v(ji,jj)
            END_2D
         ENDIF
         IF( l_trddyn ) THEN
            ! atmospheric pressure trend diagnostic
            zatmtrdu(:,:) = zu_frc(:,:) - zatmtrdu(:,:)
            zatmtrdv(:,:) = zv_frc(:,:) - zatmtrdv(:,:)
            CALL trd_dyn( zatmtrdu, zatmtrdv, jpdyn_atm, kt)
         ENDIF
      ENDIF
      !
      !                                   !=  Add wind forcing  =!
      !                                   !  ------------------  !
      IF( l_trddyn ) THEN
         IF( nn_ice == 2 ) THEN  !  Calculate and output the (partial) ice-ocean stress if using SI3.
            ztautrdu(:,:) = 0._wp ; ztautrdv(:,:) = 0._wp
            IF( ln_bt_fw ) THEN                        
               DO_2D( 0, 0, 0, 0 )
                  ztautrdu(ji,jj) =  r1_rho0 * uiceoc(ji,jj) * r1_hu(ji,jj,Kmm)
                  ztautrdv(ji,jj) =  r1_rho0 * viceoc(ji,jj) * r1_hv(ji,jj,Kmm)
               END_2D
            ELSE
               zztmp = r1_rho0 * r1_2
               DO_2D( 0, 0, 0, 0 )
                  ztautrdu(ji,jj) =  zztmp * ( uiceoc_b(ji,jj) + uiceoc(ji,jj) ) * r1_hu(ji,jj,Kmm)
                  ztautrdv(ji,jj) =  zztmp * ( viceoc_b(ji,jj) + viceoc(ji,jj) ) * r1_hv(ji,jj,Kmm)
               END_2D
            ENDIF
            IF( ln_isfcav.OR.ln_drgice_imp ) THEN
               ! Save this part of the ice-ocean drag as the first installment of top friction
               CALL trd_dyn( ztautrdu, ztautrdv, jpdyn_iceoc2d, kt, Kmm)
            ELSE
               ! In this case this is the whole top friction trend
               CALL trd_dyn( ztautrdu, ztautrdv, jpdyn_tfr, kt, Kmm)
            ENDIF
         ENDIF
         ! initialise fields for wind stress trends
         ztautrdu(:,:) = zu_frc(:,:)
         ztautrdv(:,:) = zv_frc(:,:)
      ENDIF
      !
      IF( ln_bt_fw ) THEN
         DO_2D( 0, 0, 0, 0 )
            zu_frc(ji,jj) =  zu_frc(ji,jj) + r1_rho0 * utau(ji,jj) * r1_hu(ji,jj,Kmm)
            zv_frc(ji,jj) =  zv_frc(ji,jj) + r1_rho0 * vtau(ji,jj) * r1_hv(ji,jj,Kmm)
         END_2D
      ELSE
         zztmp = r1_rho0 * r1_2
         DO_2D( 0, 0, 0, 0 )
            zu_frc(ji,jj) =  zu_frc(ji,jj) + zztmp * ( utau_b(ji,jj) + utau(ji,jj) ) * r1_hu(ji,jj,Kmm)
            zv_frc(ji,jj) =  zv_frc(ji,jj) + zztmp * ( vtau_b(ji,jj) + vtau(ji,jj) ) * r1_hv(ji,jj,Kmm)
         END_2D
      ENDIF  
      !
      IF( l_trddyn ) THEN
         ! wind stress trend diagnostic
         ztautrdu(:,:) = zu_frc(:,:) - ztautrdu(:,:)
         ztautrdv(:,:) = zv_frc(:,:) - ztautrdv(:,:) 
         CALL trd_dyn( ztautrdu, ztautrdv, jpdyn_tau2d, kt, Kmm)
      ENDIF
      !
      !              !----------------!
      !              !==  sssh_frc  ==!   Right-Hand-Side of the barotropic ssh equation   (over the FULL domain)
      !              !----------------!
      !                                   !=  Net water flux forcing applied to a water column  =!
      !                                   ! ---------------------------------------------------  !
      IF (ln_bt_fw) THEN                          ! FORWARD integration: use kt+1/2 fluxes (NOW+1/2)
         ssh_frc(:,:) = r1_rho0 * ( emp(:,:) - rnf(:,:) - fwfisf_cav(:,:) - fwfisf_par(:,:) )
      ELSE                                        ! CENTRED integration: use kt-1/2 + kt+1/2 fluxes (NOW)
         zztmp = r1_rho0 * r1_2
         ssh_frc(:,:) = zztmp * (   emp(:,:)        + emp_b(:,:)          &
            &                     - rnf(:,:)        - rnf_b(:,:)          &
            &                     - fwfisf_cav(:,:) - fwfisf_cav_b(:,:)   &
            &                     - fwfisf_par(:,:) - fwfisf_par_b(:,:)   )
      ENDIF
      !                                   !=  Add Stokes drift divergence  =!   (if exist)
      IF( ln_sdw ) THEN                   !  -----------------------------  !
         ssh_frc(:,:) = ssh_frc(:,:) + div_sd(:,:)
      ENDIF
      !
      !                                         ! ice sheet coupling
      IF ( ln_isf .AND. ln_isfcpl ) THEN
         !
         ! ice sheet coupling
         IF( ln_rstart .AND. kt == nit000 ) THEN
            ssh_frc(:,:) = ssh_frc(:,:) + risfcpl_ssh(:,:)
         END IF
         !
         ! conservation option
         IF( ln_isfcpl_cons ) THEN
            ssh_frc(:,:) = ssh_frc(:,:) + risfcpl_cons_ssh(:,:)
         END IF
         !
      END IF
      !
#if defined key_asminc
      !                                   !=  Add the IAU weighted SSH increment  =!
      !                                   !  ------------------------------------  !
      IF( lk_asminc .AND. ln_sshinc .AND. ln_asmiau ) THEN
         ssh_frc(:,:) = ssh_frc(:,:) - ssh_iau(:,:)
      ENDIF
#endif
      !                                   != Fill boundary data arrays for AGRIF
      !                                   ! ------------------------------------
#if defined key_agrif
         IF( .NOT.Agrif_Root() ) CALL agrif_dta_ts( kt )
#endif
      !
      ! -----------------------------------------------------------------------
      !  Phase 2 : Integration of the barotropic equations 
      ! -----------------------------------------------------------------------
      !
      !                                             ! ==================== !
      !                                             !    Initialisations   !
      !                                             ! ==================== !  
      ! Initialize barotropic variables:      
      IF( ll_init )THEN
         sshbb_e(:,:) = 0._wp
         ubb_e  (:,:) = 0._wp
         vbb_e  (:,:) = 0._wp
         sshb_e (:,:) = 0._wp
         ub_e   (:,:) = 0._wp
         vb_e   (:,:) = 0._wp
      ENDIF
      !
      IF( ln_linssh ) THEN    ! mid-step ocean depth is fixed (hup2_e=hu_n=hu_0)
         zhup2_e(:,:) = hu_0(:,:)
         zhvp2_e(:,:) = hv_0(:,:)
         zhtp2_e(:,:) = ht_0(:,:)
      ENDIF
      !
      IF( ln_bt_fw ) THEN                 ! FORWARD integration: start from NOW fields                    
         sshn_e(:,:) =    pssh (:,:,Kmm)            
         un_e  (:,:) =    puu_b(:,:,Kmm)            
         vn_e  (:,:) =    pvv_b(:,:,Kmm)
         !
         hu_e  (:,:) =    hu(:,:,Kmm)       
         hv_e  (:,:) =    hv(:,:,Kmm) 
         hur_e (:,:) = r1_hu(:,:,Kmm)    
         hvr_e (:,:) = r1_hv(:,:,Kmm)
      ELSE                                ! CENTRED integration: start from BEFORE fields
         sshn_e(:,:) =    pssh (:,:,Kbb)
         un_e  (:,:) =    puu_b(:,:,Kbb)         
         vn_e  (:,:) =    pvv_b(:,:,Kbb)
         !
         hu_e  (:,:) =    hu(:,:,Kbb)       
         hv_e  (:,:) =    hv(:,:,Kbb) 
         hur_e (:,:) = r1_hu(:,:,Kbb)    
         hvr_e (:,:) = r1_hv(:,:,Kbb)
      ENDIF
      !
      ! Initialize sums:
      puu_b (:,:,Kaa) = 0._wp       ! After barotropic velocities (or transport if flux form)          
      pvv_b (:,:,Kaa) = 0._wp
      pssh  (:,:,Kaa) = 0._wp       ! Sum for after averaged sea level
      un_adv(:,:)     = 0._wp       ! Sum for now transport issued from ts loop
      vn_adv(:,:)     = 0._wp
      !
      IF( ln_wd_dl ) THEN
         zuwdmask(:,:) = 0._wp  ! set to zero for definiteness (not sure this is necessary) 
         zvwdmask(:,:) = 0._wp  ! 
         zuwdav2 (:,:) = 0._wp 
         zvwdav2 (:,:) = 0._wp   
      END IF 

      !                                             ! ==================== !
      DO jn = 1, icycle                             !  sub-time-step loop  !
         !                                          ! ==================== !
         !
         l_full_nf_update = jn == icycle   ! false: disable full North fold update (performances) for jn = 1 to icycle-1
         !
         !                    !==  Update the forcing ==! (BDY and tides)
         !
         IF( ln_bdy      .AND. ln_tide )   CALL bdy_dta_tides( kt, kit=jn, pt_offset= REAL(noffset+1,wp) )
         ! Update tide potential at the beginning of current time substep
         IF( ln_tide_pot .AND. ln_tide ) THEN
            zt0substep = REAL(nsec_day, wp) - 0.5_wp*rn_Dt + (jn + noffset - 1) * rn_Dt / REAL(nn_e, wp)
            CALL upd_tide(zt0substep, Kmm)
         END IF
         !
         !                    !==  extrapolation at mid-step  ==!   (jn+1/2)
         !
         !                       !* Set extrapolation coefficients for predictor step:
         IF ((jn<3).AND.ll_init) THEN      ! Forward           
           za1 = 1._wp                                          
           za2 = 0._wp                        
           za3 = 0._wp                        
         ELSE                              ! AB3-AM4 Coefficients: bet=0.281105 
           za1 =  1.781105_wp              ! za1 =   3/2 +   bet
           za2 = -1.06221_wp               ! za2 = -(1/2 + 2*bet)
           za3 =  0.281105_wp              ! za3 = bet
         ENDIF
         !
         !                       !* Extrapolate barotropic velocities at mid-step (jn+1/2)
         !--        m+1/2               m                m-1           m-2       --!
         !--       u      = (3/2+beta) u   -(1/2+2beta) u      + beta u          --!
         !-------------------------------------------------------------------------!
         ua_e(:,:) = za1 * un_e(:,:) + za2 * ub_e(:,:) + za3 * ubb_e(:,:)
         va_e(:,:) = za1 * vn_e(:,:) + za2 * vb_e(:,:) + za3 * vbb_e(:,:)

         IF( .NOT.ln_linssh ) THEN                        !* Update ocean depth (variable volume case only)
            !                                             !  ------------------
            ! Extrapolate Sea Level at step jit+0.5:
            !--         m+1/2                 m                  m-1             m-2       --!
            !--      ssh      = (3/2+beta) ssh   -(1/2+2beta) ssh      + beta ssh          --!
            !--------------------------------------------------------------------------------!
            zsshp2_e(:,:) = za1 * sshn_e(:,:)  + za2 * sshb_e(:,:) + za3 * sshbb_e(:,:)
            
            ! set wetting & drying mask at tracer points for this barotropic mid-step
            IF( ln_wd_dl )   CALL wad_tmsk( zsshp2_e, ztwdmask )
            !
            !                          ! ocean t-depth at mid-step
            zhtp2_e(:,:) = ht_0(:,:) + zsshp2_e(:,:)
            !
            !                          ! ocean u- and v-depth at mid-step   (separate DO-loops remove the need of a lbc_lnk)
#if defined key_qcoTest_FluxForm
            !                                ! 'key_qcoTest_FluxForm' : simple ssh average
            DO_2D( 1, 0, 1, 1 )   ! not jpi-column
               zhup2_e(ji,jj) = hu_0(ji,jj) + r1_2 * (  zsshp2_e(ji,jj) + zsshp2_e(ji+1,jj  )  ) * ssumask(ji,jj)
            END_2D
            DO_2D( 1, 1, 1, 0 )
               zhvp2_e(ji,jj) = hv_0(ji,jj) + r1_2 * (  zsshp2_e(ji,jj) + zsshp2_e(ji  ,jj+1)  ) * ssvmask(ji,jj)
            END_2D
#else
            !                                ! no 'key_qcoTest_FluxForm' : surface weighted ssh average
            DO_2D( 1, 0, 1, 1 )   ! not jpi-column
               zhup2_e(ji,jj) = hu_0(ji,jj) + r1_2 * r1_e1e2u(ji,jj)                        &
                    &                              * (  e1e2t(ji  ,jj) * zsshp2_e(ji  ,jj)  &
                    &                                 + e1e2t(ji+1,jj) * zsshp2_e(ji+1,jj)  ) * ssumask(ji,jj)
            END_2D
            DO_2D( 1, 1, 1, 0 )   ! not jpj-row
               zhvp2_e(ji,jj) = hv_0(ji,jj) + r1_2 * r1_e1e2v(ji,jj)                        &
                    &                              * (  e1e2t(ji,jj  ) * zsshp2_e(ji,jj  )  &
                    &                                 + e1e2t(ji,jj+1) * zsshp2_e(ji,jj+1)  ) * ssvmask(ji,jj)
            END_2D
#endif               
            !
         ENDIF
         !
         !                    !==  after SSH  ==!   (jn+1)
         !
         !                             ! update (ua_e,va_e) to enforce volume conservation at open boundaries
         !                             ! values of zhup2_e and zhvp2_e on the halo are not needed in bdy_vol2d
         IF( ln_bdy .AND. ln_vol ) THEN
            IF (l_trddyn) THEN
               ALLOCATE( zua_tmp(jpi,jpj), zva_tmp(jpi,jpj) )
               zua_tmp(:,:) = ua_e(:,:)
               zva_tmp(:,:) = va_e(:,:)
            ENDIF
            CALL bdy_vol2d( kt, jn, ua_e, va_e, zhup2_e, zhvp2_e )
            IF (l_trddyn) THEN
               zbdytrdu(:,:) = zbdytrdu(:,:) + za2 * ( ua_e(:,:) - zua_tmp(:,:) ) / rDt_e
               zbdytrdv(:,:) = zbdytrdv(:,:) + za2 * ( va_e(:,:) - zva_tmp(:,:) ) / rDt_e
               DEALLOCATE( zua_tmp, zva_tmp )
            ENDIF
         ENDIF
         !      
         !                             ! resulting flux at mid-step (not over the full domain)
         DO_2D( 1, 0, 1, 1 )   ! not jpi-column
            zhU(ji,jj) = e2u(ji,jj) * ua_e(ji,jj) * zhup2_e(ji,jj)
         END_2D
         DO_2D( 1, 1, 1, 0 )   ! not jpj-row
            zhV(ji,jj) = e1v(ji,jj) * va_e(ji,jj) * zhvp2_e(ji,jj)
         END_2D
         !
#if defined key_agrif
         ! Set fluxes during predictor step to ensure volume conservation
         IF( ln_bt_fw )   CALL agrif_dyn_ts_flux( jn, zhU, zhV )
#endif
         IF( ln_wd_dl ) THEN           ! un_e and vn_e are set to zero at faces where 
            !                          ! the direction of the flow is from dry cells
            CALL wad_Umsk( ztwdmask, zhU, zhV, un_e, vn_e, zuwdmask, zvwdmask )   ! not jpi colomn for U, not jpj row for V
            !
         ENDIF    
         !
         !
         !     Compute Sea Level at step jit+1
         !--           m+1        m                               m+1/2          --!
         !--        ssh    =  ssh   - delta_t' * [ frc + div( flux      ) ]      --!
         !-------------------------------------------------------------------------!
         DO_2D( 0, 0, 0, 0 )
            zhdiv = (   zhU(ji,jj) - zhU(ji-1,jj) + zhV(ji,jj) - zhV(ji,jj-1)   ) * r1_e1e2t(ji,jj)
            ssha_e(ji,jj) = (  sshn_e(ji,jj) - rDt_e * ( ssh_frc(ji,jj) + zhdiv )  ) * ssmask(ji,jj)
         END_2D
         !
         CALL lbc_lnk( 'dynspg_ts', ssha_e, 'T', 1._dp,  zhU, 'U', -1._dp)
         CALL lbc_lnk( 'dynspg_ts',  zhV, 'V', -1._wp )
         !
         ! Duplicate sea level across open boundaries (this is only cosmetic if linssh=T)
         IF( ln_bdy )   CALL bdy_ssh( ssha_e )
#if defined key_agrif
         CALL agrif_ssh_ts( jn )
#endif
         !
         !                             ! Sum over sub-time-steps to compute advective velocities
         za2 = wgtbtp2(jn)             ! zhU, zhV hold fluxes extrapolated at jn+0.5
         un_adv(:,:) = un_adv(:,:) + za2 * zhU(:,:) * r1_e2u(:,:)
         vn_adv(:,:) = vn_adv(:,:) + za2 * zhV(:,:) * r1_e1v(:,:)
         ! sum over sub-time-steps to decide which baroclinic velocities to set to zero (zuwdav2 is only used when ln_wd_dl_bc=True) 
         IF ( ln_wd_dl_bc ) THEN
            DO_2D( 1, 0, 1, 1 )   ! not jpi-column
               zuwdav2(ji,jj) = zuwdav2(ji,jj) + za2 * zuwdmask(ji,jj)
            END_2D
            DO_2D( 1, 1, 1, 0 )   ! not jpj-row
               zvwdav2(ji,jj) = zvwdav2(ji,jj) + za2 * zvwdmask(ji,jj)
            END_2D
         END IF
         !
         !  
         ! Sea Surface Height at u-,v-points (vvl case only)
         IF( .NOT.ln_linssh ) THEN
#if defined key_qcoTest_FluxForm
            !                                ! 'key_qcoTest_FluxForm' : simple ssh average
            DO_2D( 1, 0, 1, 1 )
               zsshu_a(ji,jj) = r1_2 * (  ssha_e(ji,jj) + ssha_e(ji+1,jj  )  ) * ssumask(ji,jj)
            END_2D
            DO_2D( 1, 1, 1, 0 )
               zsshv_a(ji,jj) = r1_2 * (  ssha_e(ji,jj) + ssha_e(ji  ,jj+1)  ) * ssvmask(ji,jj)
            END_2D
#else
            DO_2D( 0, 0, 0, 0 )
               zsshu_a(ji,jj) = r1_2 * r1_e1e2u(ji,jj) * ( e1e2t(ji  ,jj  ) * ssha_e(ji  ,jj  )   &
                  &                                      + e1e2t(ji+1,jj  ) * ssha_e(ji+1,jj  ) ) * ssumask(ji,jj)
               zsshv_a(ji,jj) = r1_2 * r1_e1e2v(ji,jj) * ( e1e2t(ji  ,jj  ) * ssha_e(ji  ,jj  )   &
                  &                                      + e1e2t(ji  ,jj+1) * ssha_e(ji  ,jj+1) ) * ssvmask(ji,jj)
            END_2D
#endif
         ENDIF
         !         
         ! Half-step back interpolation of SSH for surface pressure computation at step jit+1/2
         !--            m+1/2           m+1              m               m-1              m-2     --!
         !--        ssh'    =  za0 * ssh     +  za1 * ssh   +  za2 * ssh      +  za3 * ssh        --!
         !------------------------------------------------------------------------------------------!
         CALL ts_bck_interp( jn, ll_init, za0, za1, za2, za3 )   ! coeficients of the interpolation
         zsshp2_e(:,:) = za0 *  ssha_e(:,:) + za1 *  sshn_e (:,:)   &
            &          + za2 *  sshb_e(:,:) + za3 *  sshbb_e(:,:)
         !
         !                             ! Surface pressure gradient
         zldg = ( 1._wp - rn_scal_load ) * grav    ! local factor
         DO_2D( 0, 0, 0, 0 )
            zu_spg(ji,jj) = - zldg * ( zsshp2_e(ji+1,jj) - zsshp2_e(ji,jj) ) * r1_e1u(ji,jj)
            zv_spg(ji,jj) = - zldg * ( zsshp2_e(ji,jj+1) - zsshp2_e(ji,jj) ) * r1_e2v(ji,jj)
         END_2D
         !
         IF( l_trddyn ) THEN
            za2 = wgtbtp2(jn)
            zspgtrdu(:,:) = zspgtrdu(:,:) + za2 * zu_spg(:,:) * ssumask(:,:)
            zspgtrdv(:,:) = zspgtrdv(:,:) + za2 * zv_spg(:,:) * ssvmask(:,:)
         ENDIF
         !
         ! Add Coriolis trend:
         ! zwz array below or triads normally depend on sea level with ln_linssh=F and should be updated
         ! at each time step. We however keep them constant here for optimization.
         ! Recall that zhU and zhV hold fluxes at jn+0.5 (extrapolated not backward interpolated)
         CALL dyn_cor_2d( zhtp2_e, zhup2_e, zhvp2_e, ua_e, va_e, zhU, zhV,    zu_trd, zv_trd   )
         !
         IF( l_trddyn ) THEN
            za2 = wgtbtp2(jn)
            zpvotrdu(:,:) = zpvotrdu(:,:) + za2 * zu_trd(:,:) * ssumask(:,:)
            zpvotrdv(:,:) = zpvotrdv(:,:) + za2 * zv_trd(:,:) * ssvmask(:,:)
         ENDIF
         ! Add tidal astronomical forcing if defined
         IF ( ln_tide .AND. ln_tide_pot ) THEN
            IF(l_trddyn) THEN
               ALLOCATE( ztidepotgradu(jpi,jpj), ztidepotgradv(jpi,jpj) )            
               ztidepotgradu(:,:) = zu_trd(:,:)
               ztidepotgradv(:,:) = zv_trd(:,:) 
            ENDIF
            !
            DO_2D( 0, 0, 0, 0 )
               zu_trd(ji,jj) = zu_trd(ji,jj) + grav * ( pot_astro(ji+1,jj) - pot_astro(ji,jj) ) * r1_e1u(ji,jj)
               zv_trd(ji,jj) = zv_trd(ji,jj) + grav * ( pot_astro(ji,jj+1) - pot_astro(ji,jj) ) * r1_e2v(ji,jj)
            END_2D
            !
            IF( l_trddyn ) THEN
               za2 = wgtbtp2(jn)
               ztidepotgradu(:,:) = zu_trd(:,:) - ztidepotgradu(:,:)
               ztidepotgradv(:,:) = zv_trd(:,:) - ztidepotgradv(:,:)
               zspgtrdu(:,:) = zspgtrdu + za2 * ztidepotgradu(:,:) * ssumask(:,:)
               zspgtrdv(:,:) = zspgtrdv + za2 * ztidepotgradv(:,:) * ssvmask(:,:)
               DEALLOCATE( ztidepotgradu, ztidepotgradv )
            ENDIF 
         ENDIF 
         !
         ! Add bottom stresses:
!jth do implicitly instead
         IF ( .NOT. ll_wd ) THEN ! Revert to explicit for bit comparison tests in non wad runs
            DO_2D( 0, 0, 0, 0 )
               zu_trd(ji,jj) = zu_trd(ji,jj) + zCdU_u(ji,jj) * un_e(ji,jj) * hur_e(ji,jj)
               zv_trd(ji,jj) = zv_trd(ji,jj) + zCdU_v(ji,jj) * vn_e(ji,jj) * hvr_e(ji,jj)
            END_2D
            IF( l_trddyn ) THEN
               za2 = wgtbtp2(jn)
               IF( ln_isfcav.OR.ln_drgice_imp ) THEN          ! top+bottom friction (ocean cavities)
                  DO_2D( 0, 0, 0, 0 )
                     ztfrtrdu(ji,jj) = ztfrtrdu(ji,jj) + za2 * 0.5_wp*( rCdU_top(ji+1,jj)+rCdU_top(ji,jj)) * un_e(ji,jj) * hur_e(ji,jj)
                     ztfrtrdv(ji,jj) = ztfrtrdv(ji,jj) + za2 * 0.5_wp*( rCdU_top(ji,jj+1)+rCdU_top(ji,jj)) * vn_e(ji,jj) * hvr_e(ji,jj)
                  END_2D
               ENDIF
               DO_2D( 0, 0, 0, 0 )
                  zbfrtrdu(ji,jj) = zbfrtrdu(ji,jj) + za2 * 0.5_wp*( rCdU_bot(ji+1,jj)+rCdU_bot(ji,jj)) * un_e(ji,jj) * hur_e(ji,jj)
                  zbfrtrdv(ji,jj) = zbfrtrdv(ji,jj) + za2 * 0.5_wp*( rCdU_bot(ji,jj+1)+rCdU_bot(ji,jj)) * vn_e(ji,jj) * hvr_e(ji,jj)
               END_2D
            ENDIF
         ENDIF
!jth_IWD Add internal wave drag
        IF ( ln_int_wave_drag .AND. .NOT. ll_wd ) THEN
            DO_2D( 0, 1, 0, 1 )
                  zu_trd(ji,jj) = zu_trd(ji,jj) - r1_2*( tdiss(ji+1,jj)+ tdiss(ji,jj))  * ( un_e(ji,jj) - un_detide(ji,jj) ) * hur_e(ji,jj)
                  zv_trd(ji,jj) = zv_trd(ji,jj) - r1_2*( tdiss(ji,jj+1)+ tdiss(ji,jj))  * ( vn_e(ji,jj) - vn_detide(ji,jj) ) * hvr_e(ji,jj)
            END_2D
        ENDIF
!!
         !
         ! Set next velocities:
         !     Compute barotropic speeds at step jit+1    (h : total height of the water colomn)
         !--                              VECTOR FORM
         !--   m+1                 m               /                                                       m+1/2           \    --!
         !--  u     =             u   + delta_t' * \         (1-r)*g * grad_x( ssh') -         f * k vect u      +     frc /    --!
         !--                                                                                                                    --!
         !--                             FLUX FORM                                                                              --!
         !--  m+1   __1__  /  m    m               /  m+1/2                             m+1/2              m+1/2    n      \ \  --!
         !-- u    =   m+1 |  h  * u   + delta_t' * \ h     * (1-r)*g * grad_x( ssh') - h     * f * k vect u      + h * frc /  | --!
         !--         h     \                                                                                                 /  --!
         !------------------------------------------------------------------------------------------------------------------------!
         IF( ln_dynadv_vec .OR. ln_linssh ) THEN      !* Vector form
            DO_2D( 0, 0, 0, 0 )
               ua_e(ji,jj) = (                                 un_e(ji,jj)   & 
                         &     + rDt_e * (                   zu_spg(ji,jj)   &
                         &                                 + zu_trd(ji,jj)   &
                         &                                 + zu_frc(ji,jj) ) & 
                         &   ) * ssumask(ji,jj)

               va_e(ji,jj) = (                                 vn_e(ji,jj)   &
                         &     + rDt_e * (                   zv_spg(ji,jj)   &
                         &                                 + zv_trd(ji,jj)   &
                         &                                 + zv_frc(ji,jj) ) &
                         &   ) * ssvmask(ji,jj)
            END_2D
            !
         ELSE                           !* Flux form
            DO_2D( 0, 0, 0, 0 )
               !                    ! hu_e, hv_e hold depth at jn,  zhup2_e, zhvp2_e hold extrapolated depth at jn+1/2
               !                    ! backward interpolated depth used in spg terms at jn+1/2
#if defined key_qcoTest_FluxForm
            !                                ! 'key_qcoTest_FluxForm' : simple ssh average
               zhu_bck = hu_0(ji,jj) + r1_2 * (  zsshp2_e(ji,jj) + zsshp2_e(ji+1,jj  )  ) * ssumask(ji,jj)
               zhv_bck = hv_0(ji,jj) + r1_2 * (  zsshp2_e(ji,jj) + zsshp2_e(ji  ,jj+1)  ) * ssvmask(ji,jj)
#else
               zhu_bck = hu_0(ji,jj) + r1_2*r1_e1e2u(ji,jj) * (  e1e2t(ji  ,jj) * zsshp2_e(ji  ,jj)    &
                    &                                          + e1e2t(ji+1,jj) * zsshp2_e(ji+1,jj)  ) * ssumask(ji,jj)
               zhv_bck = hv_0(ji,jj) + r1_2*r1_e1e2v(ji,jj) * (  e1e2t(ji,jj  ) * zsshp2_e(ji,jj  )    &
                    &                                          + e1e2t(ji,jj+1) * zsshp2_e(ji,jj+1)  ) * ssvmask(ji,jj)
#endif
               !                    ! inverse depth at jn+1
               z1_hu = ssumask(ji,jj) / ( hu_0(ji,jj) + zsshu_a(ji,jj) + 1._wp - ssumask(ji,jj) )
               z1_hv = ssvmask(ji,jj) / ( hv_0(ji,jj) + zsshv_a(ji,jj) + 1._wp - ssvmask(ji,jj) )
               !
               ua_e(ji,jj) = (               hu_e  (ji,jj) *   un_e (ji,jj)      & 
                    &            + rDt_e * (  zhu_bck        * zu_spg (ji,jj)  &   !
                    &                       + zhup2_e(ji,jj) * zu_trd (ji,jj)  &   !
                    &                       +  hu(ji,jj,Kmm) * zu_frc (ji,jj)  )   ) * z1_hu
               !
               va_e(ji,jj) = (               hv_e  (ji,jj) *   vn_e (ji,jj)      &
                    &            + rDt_e * (  zhv_bck        * zv_spg (ji,jj)  &   !
                    &                       + zhvp2_e(ji,jj) * zv_trd (ji,jj)  &   !
                    &                       +  hv(ji,jj,Kmm) * zv_frc (ji,jj)  )   ) * z1_hv
            END_2D
         ENDIF
         IF ( ll_wd ) THEN ! RDP implicit drag
            IF (l_trddyn) THEN
               ALLOCATE( zua_tmp(jpi,jpj), zva_tmp(jpi,jpj) )
               zua_tmp(:,:) = ua_e(:,:)
               zva_tmp(:,:) = va_e(:,:)
            ENDIF
            DO_2D( 0, 0, 0, 0 )
               ua_e(ji,jj) =  ua_e(ji,jj) / ( 1._wp - rDt_e * zCdU_u(ji,jj) * hur_e(ji,jj) )
               va_e(ji,jj) =  va_e(ji,jj) / ( 1._wp - rDt_e * zCdU_v(ji,jj) * hvr_e(ji,jj) )
            END_2D
            IF (l_trddyn) THEN ! get bottom drag diagnositc (TODO: add ice drag)
               za2 = wgtbtp2(jn)
               zbfrtrdu(:,:) = zbfrtrdu(:,:) + za2 * (ua_e(:,:) - zua_tmp(:,:)) / rDt_e
               zbfrtrdv(:,:) = zbfrtrdu(:,:) + za2 * (va_e(:,:) - zva_tmp(:,:)) / rDt_e
               DEALLOCATE( zua_tmp, zva_tmp )
            ENDIF
         ENDIF
       
         IF( .NOT.ln_linssh ) THEN !* Update ocean depth (variable volume case only)
            DO_2D( 0, 0, 0, 0 )
               hu_e (ji,jj) =    hu_0(ji,jj) + zsshu_a(ji,jj)
               hur_e(ji,jj) = ssumask(ji,jj) / (  hu_e(ji,jj) + 1._wp - ssumask(ji,jj)  )
               hv_e (ji,jj) =    hv_0(ji,jj) + zsshv_a(ji,jj)
               hvr_e(ji,jj) = ssvmask(ji,jj) / (  hv_e(ji,jj) + 1._wp - ssvmask(ji,jj)  )
            END_2D
         ENDIF
         !
         IF( .NOT.ln_linssh ) THEN   !* Update ocean depth (variable volume case only)
            CALL lbc_lnk( 'dynspg_ts', ua_e , 'U', -1._wp, va_e , 'V', -1._wp  &
                 &                   , hu_e , 'U',  1._wp, hv_e , 'V',  1._wp  &
                 &                   , hur_e, 'U',  1._wp, hvr_e, 'V',  1._wp  )
         ELSE
            CALL lbc_lnk( 'dynspg_ts', ua_e , 'U', -1._wp, va_e , 'V', -1._wp  )
         ENDIF
         !                                                 ! open boundaries
         IF( ln_bdy ) THEN
            IF (l_trddyn) THEN
               ALLOCATE( zua_tmp(jpi,jpj), zva_tmp(jpi,jpj) )
               zua_tmp(:,:) = ua_e(:,:)
               zva_tmp(:,:) = va_e(:,:)
            ENDIF
            CALL bdy_dyn2d( jn, ua_e, va_e, un_e, vn_e, hur_e, hvr_e, ssha_e )
            IF (l_trddyn) THEN
               za2 = wgtbtp2(jn)
               zbdytrdu(:,:) = zbdytrdu(:,:) + za2 * ( ua_e(:,:) - zua_tmp(:,:) ) / rDt_e
               zbdytrdv(:,:) = zbdytrdv(:,:) + za2 * ( va_e(:,:) - zva_tmp(:,:) ) / rDt_e
               DEALLOCATE( zua_tmp, zva_tmp )
            ENDIF
         ENDIF
#if defined key_agrif                                                           
         CALL agrif_dyn_ts( jn )  ! Agrif
#endif
         !                                             !* Swap
         !                                             !  ----
         ubb_e  (:,:) = ub_e  (:,:)
         ub_e   (:,:) = un_e  (:,:)
         un_e   (:,:) = ua_e  (:,:)
         !
         vbb_e  (:,:) = vb_e  (:,:)
         vb_e   (:,:) = vn_e  (:,:)
         vn_e   (:,:) = va_e  (:,:)
         !
         sshbb_e(:,:) = sshb_e(:,:)
         sshb_e (:,:) = sshn_e(:,:)
         sshn_e (:,:) = ssha_e(:,:)

         !                                             !* Sum over whole bt loop
         !                                             !  ----------------------
         za1 = wgtbtp1(jn)                                    
         IF( ln_dynadv_vec .OR. ln_linssh ) THEN    ! Sum velocities
            puu_b  (:,:,Kaa) = puu_b  (:,:,Kaa) + za1 * ua_e  (:,:) 
            pvv_b  (:,:,Kaa) = pvv_b  (:,:,Kaa) + za1 * va_e  (:,:) 
         ELSE                                       ! Sum transports
            IF ( .NOT.ln_wd_dl ) THEN  
               puu_b  (:,:,Kaa) = puu_b  (:,:,Kaa) + za1 * ua_e  (:,:) * hu_e (:,:)
               pvv_b  (:,:,Kaa) = pvv_b  (:,:,Kaa) + za1 * va_e  (:,:) * hv_e (:,:)
            ELSE 
               puu_b  (:,:,Kaa) = puu_b  (:,:,Kaa) + za1 * ua_e  (:,:) * hu_e (:,:) * zuwdmask(:,:)
               pvv_b  (:,:,Kaa) = pvv_b  (:,:,Kaa) + za1 * va_e  (:,:) * hv_e (:,:) * zvwdmask(:,:)
            END IF 
         ENDIF
         !                                          ! Sum sea level
         pssh(:,:,Kaa) = pssh(:,:,Kaa) + za1 * ssha_e(:,:)

         !                                                 ! ==================== !
      END DO                                               !        end loop      !
      !                                                    ! ==================== !
      ! -----------------------------------------------------------------------------
      ! Phase 3. update the general trend with the barotropic trend
      ! -----------------------------------------------------------------------------
      !
      ! Set advection velocity correction:
      IF (ln_bt_fw) THEN
         IF( .NOT.( kt == nit000 .AND. l_1st_euler ) ) THEN
            DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
               zun_save = un_adv(ji,jj)
               zvn_save = vn_adv(ji,jj)
               !                          ! apply the previously computed correction 
               un_adv(ji,jj) = r1_2 * ( ub2_b(ji,jj) + zun_save - rn_atfp * un_bf(ji,jj) )
               vn_adv(ji,jj) = r1_2 * ( vb2_b(ji,jj) + zvn_save - rn_atfp * vn_bf(ji,jj) )
               !                          ! Update corrective fluxes for next time step
               un_bf(ji,jj)  = rn_atfp * un_bf(ji,jj) + ( zun_save - ub2_b(ji,jj) )
               vn_bf(ji,jj)  = rn_atfp * vn_bf(ji,jj) + ( zvn_save - vb2_b(ji,jj) )
               !                          ! Save integrated transport for next computation
               ub2_b(ji,jj) = zun_save
               vb2_b(ji,jj) = zvn_save
            END_2D
         ELSE
            un_bf(:,:) = 0._wp            ! corrective fluxes for next time step set to zero
            vn_bf(:,:) = 0._wp
            ub2_b(:,:) = un_adv(:,:)      ! Save integrated transport for next computation
            vb2_b(:,:) = vn_adv(:,:)
         END IF
      ENDIF


      !
      ! Update barotropic trend:
      IF( ln_dynadv_vec .OR. ln_linssh ) THEN
         DO jk=1,jpkm1
            puu(:,:,jk,Krhs) = puu(:,:,jk,Krhs) + ( puu_b(:,:,Kaa) - puu_b(:,:,Kbb) ) * r1_Dt_b
            pvv(:,:,jk,Krhs) = pvv(:,:,jk,Krhs) + ( pvv_b(:,:,Kaa) - pvv_b(:,:,Kbb) ) * r1_Dt_b
         END DO
         IF( l_trddyn ) THEN
            ztottrdu(:,:) = ( puu_b(:,:,Kaa) - puu_b(:,:,Kbb) ) * r1_Dt_b
            ztottrdv(:,:) = ( pvv_b(:,:,Kaa) - pvv_b(:,:,Kbb) ) * r1_Dt_b
            CALL trd_dyn( ztottrdu, ztottrdv, jpdyn_tot, kt)
         ENDIF
      ELSE
         ! At this stage, pssh(:,:,:,Krhs) has been corrected: compute new depths at velocity points
#if defined key_qcoTest_FluxForm
         !                                ! 'key_qcoTest_FluxForm' : simple ssh average
         DO_2D( 1, 0, 1, 0 )
            zsshu_a(ji,jj) = r1_2 * ( pssh(ji,jj,Kaa) + pssh(ji+1,jj  ,Kaa) ) * ssumask(ji,jj)
            zsshv_a(ji,jj) = r1_2 * ( pssh(ji,jj,Kaa) + pssh(ji  ,jj+1,Kaa) ) * ssvmask(ji,jj)
         END_2D
#else
         DO_2D( 1, 0, 1, 0 )
            zsshu_a(ji,jj) = r1_2 * r1_e1e2u(ji,jj) * ( e1e2t(ji  ,jj) * pssh(ji  ,jj,Kaa)   &
               &                                      + e1e2t(ji+1,jj) * pssh(ji+1,jj,Kaa) ) * ssumask(ji,jj)
            zsshv_a(ji,jj) = r1_2 * r1_e1e2v(ji,jj) * ( e1e2t(ji,jj  ) * pssh(ji,jj  ,Kaa)   &
               &                                      + e1e2t(ji,jj+1) * pssh(ji,jj+1,Kaa) ) * ssvmask(ji,jj)
         END_2D
#endif   
         CALL lbc_lnk( 'dynspg_ts', zsshu_a, 'U', 1._wp, zsshv_a, 'V', 1._wp ) ! Boundary conditions
         !
         DO jk=1,jpkm1
            puu(:,:,jk,Krhs) = puu(:,:,jk,Krhs) + r1_hu(:,:,Kmm)   &
               &             * ( puu_b(:,:,Kaa) - puu_b(:,:,Kbb) * hu(:,:,Kbb) ) * r1_Dt_b
            pvv(:,:,jk,Krhs) = pvv(:,:,jk,Krhs) + r1_hv(:,:,Kmm)   &
               &             * ( pvv_b(:,:,Kaa) - pvv_b(:,:,Kbb) * hv(:,:,Kbb) ) * r1_Dt_b
         END DO
         ! Save barotropic velocities not transport:
         puu_b(:,:,Kaa) =  puu_b(:,:,Kaa) / ( hu_0(:,:) + zsshu_a(:,:) + 1._wp - ssumask(:,:) )
         pvv_b(:,:,Kaa) =  pvv_b(:,:,Kaa) / ( hv_0(:,:) + zsshv_a(:,:) + 1._wp - ssvmask(:,:) )
         IF( l_trddyn ) THEN
            ztottrdu(:,:) = r1_hu(:,:,Kmm) * ( puu_b(:,:,Kaa) - puu_b(:,:,Kbb) * hu(:,:,Kbb) ) * r1_Dt_b
            ztottrdv(:,:) = r1_hv(:,:,Kmm) * ( pvv_b(:,:,Kaa) - pvv_b(:,:,Kbb) * hv(:,:,Kbb) ) * r1_Dt_b
            CALL trd_dyn( ztottrdu, ztottrdv, jpdyn_tot, kt)
         ENDIF
      ENDIF


      ! Correct velocities so that the barotropic velocity equals (un_adv, vn_adv) (in all cases)  
      DO jk = 1, jpkm1
         puu(:,:,jk,Kmm) = ( puu(:,:,jk,Kmm) + un_adv(:,:)*r1_hu(:,:,Kmm) - puu_b(:,:,Kmm) ) * umask(:,:,jk)
         pvv(:,:,jk,Kmm) = ( pvv(:,:,jk,Kmm) + vn_adv(:,:)*r1_hv(:,:,Kmm) - pvv_b(:,:,Kmm) ) * vmask(:,:,jk)
      END DO

      IF ( ln_wd_dl .and. ln_wd_dl_bc) THEN 
         DO jk = 1, jpkm1
            puu(:,:,jk,Kmm) = ( un_adv(:,:)*r1_hu(:,:,Kmm) &
                       & + zuwdav2(:,:)*(puu(:,:,jk,Kmm) - un_adv(:,:)*r1_hu(:,:,Kmm)) ) * umask(:,:,jk) 
            pvv(:,:,jk,Kmm) = ( vn_adv(:,:)*r1_hv(:,:,Kmm) & 
                       & + zvwdav2(:,:)*(pvv(:,:,jk,Kmm) - vn_adv(:,:)*r1_hv(:,:,Kmm)) ) * vmask(:,:,jk)  
         END DO
      END IF 

      
      CALL iom_put(  "ubar", un_adv(:,:)*r1_hu(:,:,Kmm) )    ! barotropic i-current
      CALL iom_put(  "vbar", vn_adv(:,:)*r1_hv(:,:,Kmm) )    ! barotropic i-current
      !
#if defined key_agrif
      ! Save time integrated fluxes during child grid integration
      ! (used to update coarse grid transports at next time step)
      !
      IF( .NOT.Agrif_Root() .AND. ln_bt_fw .AND. ln_agrif_2way ) THEN
         IF( Agrif_NbStepint() == 0 ) THEN
            ub2_i_b(:,:) = 0._wp
            vb2_i_b(:,:) = 0._wp
         END IF
         !
         za1 = 1._wp / REAL(Agrif_rhot(), wp)
         ub2_i_b(:,:) = ub2_i_b(:,:) + za1 * ub2_b(:,:)
         vb2_i_b(:,:) = vb2_i_b(:,:) + za1 * vb2_b(:,:)
      ENDIF
#endif      
      !                                   !* write time-spliting arrays in the restart
      IF( lrst_oce .AND.ln_bt_fw )   CALL ts_rst( kt, 'WRITE' )
      !
      IF( ln_wd_dl )   DEALLOCATE( ztwdmask, zuwdmask, zvwdmask, zuwdav2, zvwdav2 )
      !
      IF( l_trddyn ) THEN
         CALL trd_dyn( zspgtrdu, zspgtrdv, jpdyn_spg, kt, Kmm)
         CALL trd_dyn( zpvotrdu, zpvotrdv, jpdyn_pvo, kt, Kmm)
         CALL trd_dyn( zbfrtrdu, zbfrtrdv, jpdyn_bfr, kt, Kmm)
         DEALLOCATE( zspgtrdu, zspgtrdv, zpvotrdu, zpvotrdv, ztautrdu, ztautrdv, zbfrtrdu, zbfrtrdv )
         IF ( ln_bdy ) THEN
            CALL trd_dyn( zbdytrdu, zbdytrdv, jpdyn_bdy2d, kt )
            DEALLOCATE( zbdytrdu, zbdytrdv )
         ENDIF
         IF( ln_isfcav.OR.ln_drgice_imp ) THEN          ! top+bottom friction (ocean cavities)
            CALL trd_dyn( ztfrtrdu, ztfrtrdv, jpdyn_tfr, kt, Kmm)
            DEALLOCATE( ztfrtrdu, ztfrtrdv )
         ENDIF
      ENDIF
      !
      CALL iom_put( "baro_u" , puu_b(:,:,Kmm)*ssumask(:,:) )  ! Barotropic  U Velocity
      CALL iom_put( "baro_v" , pvv_b(:,:,Kmm)*ssvmask(:,:) )  ! Barotropic  V Velocity
   END SUBROUTINE dyn_spg_ts


   SUBROUTINE ts_wgt( ll_av, ll_fw, jpit, zwgt1, zwgt2)
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE ts_wgt  ***
      !!
      !! ** Purpose : Set time-splitting weights for temporal averaging (or not)
      !!----------------------------------------------------------------------
      LOGICAL, INTENT(in) ::   ll_av      ! temporal averaging=.true.
      LOGICAL, INTENT(in) ::   ll_fw      ! forward time splitting =.true.
      INTEGER, INTENT(inout) :: jpit      ! cycle length    
      REAL(wp), DIMENSION(3*nn_e), INTENT(inout)  :: zwgt2 ! Primary weights
      REAL(dp), DIMENSION(3*nn_e), INTENT(inout)  :: zwgt1 ! Secondary weights
      
      INTEGER ::  jic, jn, ji                      ! temporary integers
      REAL(wp)  :: za2
      REAL(dp)  :: za1
      !!----------------------------------------------------------------------

      zwgt1(:) = 0._wp
      zwgt2(:) = 0._wp

      ! Set time index when averaged value is requested
      IF (ll_fw) THEN 
         jic = nn_e
      ELSE
         jic = 2 * nn_e
      ENDIF

      ! Set primary weights:
      IF (ll_av) THEN
           ! Define simple boxcar window for primary weights 
           ! (width = nn_e, centered around jic)     
         SELECT CASE ( nn_bt_flt )
              CASE( 0 )  ! No averaging
                 zwgt1(jic) = 1._wp
                 jpit = jic

              CASE( 1 )  ! Boxcar, width = nn_e
                 DO jn = 1, 3*nn_e
                    za1 = ABS(float(jn-jic))/float(nn_e) 
                    IF (za1 < 0.5_wp) THEN
                      zwgt1(jn) = 1._wp
                      jpit = jn
                    ENDIF
                 ENDDO

              CASE( 2 )  ! Boxcar, width = 2 * nn_e
                 DO jn = 1, 3*nn_e
                    za1 = ABS(float(jn-jic))/float(nn_e) 
                    IF (za1 < 1._wp) THEN
                      zwgt1(jn) = 1._wp
                      jpit = jn
                    ENDIF
                 ENDDO
              CASE DEFAULT   ;   CALL ctl_stop( 'unrecognised value for nn_bt_flt' )
         END SELECT

      ELSE ! No time averaging
         zwgt1(jic) = 1._wp
         jpit = jic
      ENDIF
    
      ! Set secondary weights
      DO jn = 1, jpit
        DO ji = jn, jpit
             zwgt2(jn) = zwgt2(jn) + zwgt1(ji)
        END DO
      END DO

      ! Normalize weigths:
      za1 = 1._wp / SUM(zwgt1(1:jpit))
      za2 = 1._wp / SUM(zwgt2(1:jpit))
      DO jn = 1, jpit
        zwgt1(jn) = zwgt1(jn) * za1
        zwgt2(jn) = zwgt2(jn) * za2
      END DO
      !
   END SUBROUTINE ts_wgt


   SUBROUTINE ts_rst( kt, cdrw )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE ts_rst  ***
      !!
      !! ** Purpose : Read or write time-splitting arrays in restart file
      !!----------------------------------------------------------------------
      INTEGER         , INTENT(in) ::   kt     ! ocean time-step
      CHARACTER(len=*), INTENT(in) ::   cdrw   ! "READ"/"WRITE" flag
      !!----------------------------------------------------------------------
      !
      IF( TRIM(cdrw) == 'READ' ) THEN        ! Read/initialise 
         !                                   ! ---------------
         IF( ln_rstart .AND. ln_bt_fw .AND. .NOT.l_1st_euler ) THEN    !* Read the restart file
            CALL iom_get( numror, jpdom_auto, 'ub2_b'  , ub2_b  (:,:), cd_type = 'U', psgn = -1._wp )   
            CALL iom_get( numror, jpdom_auto, 'vb2_b'  , vb2_b  (:,:), cd_type = 'V', psgn = -1._wp ) 
            CALL iom_get( numror, jpdom_auto, 'un_bf'  , un_bf  (:,:), cd_type = 'U', psgn = -1._wp )   
            CALL iom_get( numror, jpdom_auto, 'vn_bf'  , vn_bf  (:,:), cd_type = 'V', psgn = -1._wp ) 
            IF( .NOT.ln_bt_av ) THEN
               CALL iom_get( numror, jpdom_auto, 'sshbb_e'  , sshbb_e(:,:), cd_type = 'T', psgn =  1._wp )   
               CALL iom_get( numror, jpdom_auto, 'ubb_e'    ,   ubb_e(:,:), cd_type = 'U', psgn = -1._wp )   
               CALL iom_get( numror, jpdom_auto, 'vbb_e'    ,   vbb_e(:,:), cd_type = 'V', psgn = -1._wp )
               CALL iom_get( numror, jpdom_auto, 'sshb_e'   ,  sshb_e(:,:), cd_type = 'T', psgn =  1._wp ) 
               CALL iom_get( numror, jpdom_auto, 'ub_e'     ,    ub_e(:,:), cd_type = 'U', psgn = -1._wp )   
               CALL iom_get( numror, jpdom_auto, 'vb_e'     ,    vb_e(:,:), cd_type = 'V', psgn = -1._wp )
            ENDIF
!jth IWD
            IF (ln_int_wave_drag) THEN
                CALL iom_get( numror, jpdom_auto, 'un_detide'  , un_detide  (:,:), cd_type = 'U', psgn = -1._wp )
                CALL iom_get( numror, jpdom_auto, 'vn_detide'  , vn_detide  (:,:), cd_type = 'V', psgn = -1._wp )
                CALL iom_get( numror, jpdom_auto, 'un_25'      , un_25  (:,:,:), cd_type = 'U', psgn = -1._wp )
                CALL iom_get( numror, jpdom_auto, 'vn_25'      , vn_25  (:,:,:), cd_type = 'V', psgn = -1._wp )
            END IF
#if defined key_agrif
            ! Read time integrated fluxes
            IF ( .NOT.Agrif_Root() ) THEN
               CALL iom_get( numror, jpdom_auto, 'ub2_i_b'  , ub2_i_b(:,:), cd_type = 'U', psgn = -1._wp )   
               CALL iom_get( numror, jpdom_auto, 'vb2_i_b'  , vb2_i_b(:,:), cd_type = 'V', psgn = -1._wp )
            ELSE
               ub2_i_b(:,:) = 0._wp   ;   vb2_i_b(:,:) = 0._wp   ! used in the 1st update of agrif
            ENDIF
#endif
         ELSE                                   !* Start from rest
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) '   ==>>>   start from rest: set barotropic values to 0'
            ub2_b  (:,:) = 0._wp   ;   vb2_b  (:,:) = 0._wp   ! used in the 1st interpol of agrif
            un_adv (:,:) = 0._wp   ;   vn_adv (:,:) = 0._wp   ! used in the 1st interpol of agrif
            un_bf  (:,:) = 0._wp   ;   vn_bf  (:,:) = 0._wp   ! used in the 1st update   of agrif
#if defined key_agrif
            ub2_i_b(:,:) = 0._wp   ;   vb2_i_b(:,:) = 0._wp   ! used in the 1st update of agrif
#endif
         ENDIF
         !
      ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN   ! Create restart file
         !                                   ! -------------------
         IF(lwp) WRITE(numout,*) '---- ts_rst ----'
         CALL iom_rstput( kt, nitrst, numrow, 'ub2_b'   , ub2_b  (:,:) )
         CALL iom_rstput( kt, nitrst, numrow, 'vb2_b'   , vb2_b  (:,:) )
         CALL iom_rstput( kt, nitrst, numrow, 'un_bf'   , un_bf  (:,:) )
         CALL iom_rstput( kt, nitrst, numrow, 'vn_bf'   , vn_bf  (:,:) )
         !
         IF (.NOT.ln_bt_av) THEN
            CALL iom_rstput( kt, nitrst, numrow, 'sshbb_e'  , sshbb_e(:,:) ) 
            CALL iom_rstput( kt, nitrst, numrow, 'ubb_e'    ,   ubb_e(:,:) )
            CALL iom_rstput( kt, nitrst, numrow, 'vbb_e'    ,   vbb_e(:,:) )
            CALL iom_rstput( kt, nitrst, numrow, 'sshb_e'   ,  sshb_e(:,:) )
            CALL iom_rstput( kt, nitrst, numrow, 'ub_e'     ,    ub_e(:,:) )
            CALL iom_rstput( kt, nitrst, numrow, 'vb_e'     ,    vb_e(:,:) )
         ENDIF
!jth IWD
         IF (ln_int_wave_drag) THEN
            CALL iom_rstput( kt, nitrst, numrow, 'un_detide'    ,   un_detide(:,:) )
            CALL iom_rstput( kt, nitrst, numrow, 'vn_detide'    ,   vn_detide(:,:) )
            CALL iom_rstput( kt, nitrst, numrow, 'un_25'    ,       un_25(:,:,:) )
            CALL iom_rstput( kt, nitrst, numrow, 'vn_25'    ,       vn_25(:,:,:) )
         END IF
#if defined key_agrif
         ! Save time integrated fluxes
         IF ( .NOT.Agrif_Root() ) THEN
            CALL iom_rstput( kt, nitrst, numrow, 'ub2_i_b'  , ub2_i_b(:,:) )
            CALL iom_rstput( kt, nitrst, numrow, 'vb2_i_b'  , vb2_i_b(:,:) )
         ENDIF
#endif
      ENDIF
      !
   END SUBROUTINE ts_rst


   SUBROUTINE dyn_spg_ts_init
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE dyn_spg_ts_init  ***
      !!
      !! ** Purpose : Set time splitting options
      !!----------------------------------------------------------------------
      INTEGER  ::   ji ,jj              ! dummy loop indices
      REAL(wp) ::   zxr2, zyr2, zcmax   ! local scalar
      REAL(wp), DIMENSION(jpi,jpj) ::   zcu
!jth_IWD
      INTEGER  ::   ierr
      !!----------------------------------------------------------------------
      !
      ! Max courant number for ext. grav. waves
      !
      DO_2D( 0, 0, 0, 0 )
         zxr2 = r1_e1t(ji,jj) * r1_e1t(ji,jj)
         zyr2 = r1_e2t(ji,jj) * r1_e2t(ji,jj)
         zcu(ji,jj) = SQRT( grav * MAX(ht_0(ji,jj),0._wp) * (zxr2 + zyr2) )
      END_2D
      !
#if defined key_agrif
     ! Discard points that are not stepped by 2d mode:
     zcu(Nis0:Nie0,Njs0:Nje0) = zcu(Nis0:Nie0,Njs0:Nje0) &
                              & * (1._wp - tmask_upd(Nis0:Nie0,Njs0:Nje0))
#endif
      !
      zcmax = MAXVAL( zcu(Nis0:Nie0,Njs0:Nje0) )
      CALL mpp_max( 'dynspg_ts', zcmax )

      ! Estimate number of iterations to satisfy a max courant number= rn_bt_cmax
      IF( ln_bt_auto )   nn_e = CEILING( rn_Dt / rn_bt_cmax * zcmax)
      
      rDt_e = rn_Dt / REAL( nn_e , dp )
      zcmax = zcmax * rDt_e
      ! Print results
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'dyn_spg_ts_init : split-explicit free surface'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~'
      IF( ln_bt_auto ) THEN
         IF(lwp) WRITE(numout,*) '     ln_ts_auto =.true. Automatically set nn_e '
         IF(lwp) WRITE(numout,*) '     Max. courant number allowed: ', rn_bt_cmax
      ELSE
         IF(lwp) WRITE(numout,*) '     ln_ts_auto=.false.: Use nn_e in namelist   nn_e = ', nn_e
      ENDIF

      IF(ln_bt_av) THEN
         IF(lwp) WRITE(numout,*) '     ln_bt_av =.true.  ==> Time averaging over nn_e time steps is on '
      ELSE
         IF(lwp) WRITE(numout,*) '     ln_bt_av =.false. => No time averaging of barotropic variables '
      ENDIF
      !
      !
      IF(ln_bt_fw) THEN
         IF(lwp) WRITE(numout,*) '     ln_bt_fw=.true.  => Forward integration of barotropic variables '
      ELSE
         IF(lwp) WRITE(numout,*) '     ln_bt_fw =.false.=> Centred integration of barotropic variables '
      ENDIF
      !
#if defined key_agrif
      ! Restrict the use of Agrif to the forward case only
      IF( .NOT.ln_bt_fw .AND. .NOT.Agrif_Root() )   CALL ctl_stop( 'AGRIF not implemented if ln_bt_fw=.FALSE.' )
#endif
      !
      IF(lwp) WRITE(numout,*)    '     Time filter choice, nn_bt_flt: ', nn_bt_flt
      SELECT CASE ( nn_bt_flt )
         CASE( 0 )      ;   IF(lwp) WRITE(numout,*) '           Dirac'
         CASE( 1 )      ;   IF(lwp) WRITE(numout,*) '           Boxcar: width = nn_e'
         CASE( 2 )      ;   IF(lwp) WRITE(numout,*) '           Boxcar: width = 2*nn_e' 
         CASE DEFAULT   ;   CALL ctl_stop( 'unrecognised value for nn_bt_flt: should 0,1, or 2' )
      END SELECT
      !
      IF(lwp) WRITE(numout,*) ' '
      IF(lwp) WRITE(numout,*) '     nn_e = ', nn_e
      IF(lwp) WRITE(numout,*) '     Barotropic time step [s] is :', rDt_e
      IF(lwp) WRITE(numout,*) '     Maximum Courant number is   :', zcmax
      !
      IF(lwp) WRITE(numout,*)    '     Time diffusion parameter rn_bt_alpha: ', rn_bt_alpha
      IF ((ln_bt_av.AND.nn_bt_flt/=0).AND.(rn_bt_alpha>0._wp)) THEN
         CALL ctl_stop( 'dynspg_ts ERROR: if rn_bt_alpha > 0, remove temporal averaging' )
      ENDIF
      !
      IF( .NOT.ln_bt_av .AND. .NOT.ln_bt_fw ) THEN
         CALL ctl_stop( 'dynspg_ts ERROR: No time averaging => only forward integration is possible' )
      ENDIF
      IF( zcmax>0.9_wp ) THEN
         CALL ctl_stop( 'dynspg_ts ERROR: Maximum Courant number is greater than 0.9: Inc. nn_e !' )          
      ENDIF
      !
      !                             ! Allocate time-splitting arrays
      IF( dyn_spg_ts_alloc() /= 0    )   CALL ctl_stop('STOP', 'dyn_spg_init: failed to allocate dynspg_ts  arrays' )

!jth_IWD
!     allocate arrays and initialise
      IF (ln_int_wave_drag) THEN
         ALLOCATE( un_detide(jpi,jpj), vn_detide(jpi,jpj),un_25(jpi,jpj,25), vn_25(jpi,jpj,25), STAT=ierr )
         IF( ierr /= 0 )   CALL ctl_stop( 'dynspg_ts_init IWD: failed to allocate arrays' )
         un_detide(:,:) = 0._wp
         vn_detide(:,:) = 0._wp
         un_25(:,:,:) = 0._wp
         vn_25(:,:,:) = 0._wp
      ENDIF
      !
      !                             ! read restart when needed
      CALL ts_rst( nit000, 'READ' )
      !
   END SUBROUTINE dyn_spg_ts_init

   
   SUBROUTINE dyn_cor_2D_init( Kmm )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE dyn_cor_2D_init  ***
      !!
      !! ** Purpose : Set time splitting options
      !! Set arrays to remove/compute coriolis trend.
      !! Do it once during initialization if volume is fixed, else at each long time step.
      !! Note that these arrays are also used during barotropic loop. These are however frozen
      !! although they should be updated in the variable volume case. Not a big approximation.
      !! To remove this approximation, copy lines below inside barotropic loop
      !! and update depths at T- points (ht) at each barotropic time step
      !!
      !! Compute zwz = f / ( height of the water colomn )
      !!----------------------------------------------------------------------
      INTEGER,  INTENT(in)         ::  Kmm  ! Time index
      INTEGER  ::   ji ,jj, jk              ! dummy loop indices
      REAL(wp) ::   z1_ht
      !!----------------------------------------------------------------------
      !
      SELECT CASE( nvor_scheme )
      CASE( np_EEN, np_ENE, np_ENS , np_MIX )   !=  schemes using the same e3f definition
         SELECT CASE( nn_e3f_typ )                  !* ff_f/e3 at F-point
         CASE ( 0 )                                   ! original formulation  (masked averaging of e3t divided by 4)
            DO_2D( 0, 0, 0, 0 )
               zwz(ji,jj) = ( ht(ji,jj+1) + ht(ji+1,jj+1)   &
                    &       + ht(ji,jj  ) + ht(ji+1,jj  ) ) * 0.25_wp  
               IF( zwz(ji,jj) /= 0._wp )   zwz(ji,jj) = ff_f(ji,jj) / zwz(ji,jj)
            END_2D
         CASE ( 1 )                                   ! new formulation  (masked averaging of e3t divided by the sum of mask)
            DO_2D( 0, 0, 0, 0 )
               zwz(ji,jj) =     (    ht(ji,jj+1) +     ht(ji+1,jj+1)      &
                    &            +   ht(ji,jj  ) +     ht(ji+1,jj  )  )   &
                    &    / ( MAX(ssmask(ji,jj+1) + ssmask(ji+1,jj+1)      &
                    &          + ssmask(ji,jj  ) + ssmask(ji+1,jj  ) , 1._wp  )   )
               IF( zwz(ji,jj) /= 0._wp )   zwz(ji,jj) = ff_f(ji,jj) / zwz(ji,jj)
            END_2D
         END SELECT
         CALL lbc_lnk( 'dynspg_ts', zwz, 'F', 1._wp )
      END SELECT
      !
      SELECT CASE( nvor_scheme )
      CASE( np_EEN )
         !
         ftne(1,:) = 0._wp   ;   ftnw(1,:) = 0._wp   ;   ftse(1,:) = 0._wp   ;   ftsw(1,:) = 0._wp
         DO_2D( 0, 1, 0, 1 )
            ftne(ji,jj) = zwz(ji-1,jj  ) + zwz(ji  ,jj  ) + zwz(ji  ,jj-1)
            ftnw(ji,jj) = zwz(ji-1,jj-1) + zwz(ji-1,jj  ) + zwz(ji  ,jj  )
            ftse(ji,jj) = zwz(ji  ,jj  ) + zwz(ji  ,jj-1) + zwz(ji-1,jj-1)
            ftsw(ji,jj) = zwz(ji  ,jj-1) + zwz(ji-1,jj-1) + zwz(ji-1,jj  )
         END_2D
         !
      CASE( np_EET )                            != EEN scheme using e3t energy conserving scheme
         ftne(1,:) = 0._wp   ;   ftnw(1,:) = 0._wp   ;   ftse(1,:) = 0._wp   ;   ftsw(1,:) = 0._wp
         DO_2D( 0, 1, 0, 1 )
            z1_ht = ssmask(ji,jj) / ( ht(ji,jj) + 1._wp - ssmask(ji,jj) )
            ftne(ji,jj) = ( ff_f(ji-1,jj  ) + ff_f(ji  ,jj  ) + ff_f(ji  ,jj-1) ) * z1_ht
            ftnw(ji,jj) = ( ff_f(ji-1,jj-1) + ff_f(ji-1,jj  ) + ff_f(ji  ,jj  ) ) * z1_ht
            ftse(ji,jj) = ( ff_f(ji  ,jj  ) + ff_f(ji  ,jj-1) + ff_f(ji-1,jj-1) ) * z1_ht
            ftsw(ji,jj) = ( ff_f(ji  ,jj-1) + ff_f(ji-1,jj-1) + ff_f(ji-1,jj  ) ) * z1_ht
         END_2D
         !
      END SELECT
      
   END SUBROUTINE dyn_cor_2d_init


   SUBROUTINE dyn_cor_2d( pht, phu, phv, punb, pvnb, zhU, zhV,    zu_trd, zv_trd   )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE dyn_cor_2d  ***
      !!
      !! ** Purpose : Compute u and v coriolis trends
      !!----------------------------------------------------------------------
      INTEGER  ::   ji ,jj                             ! dummy loop indices
      REAL(wp) ::   zx1, zx2, zy1, zy2, z1_hu, z1_hv   !   -      -
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   )  :: pht, phu, phv, punb, pvnb, zhV
      REAL(dp), DIMENSION(jpi,jpj), INTENT(in   )  :: zhU
      REAL(wp), DIMENSION(jpi,jpj), INTENT(  out) :: zu_trd, zv_trd
      !!----------------------------------------------------------------------
      SELECT CASE( nvor_scheme )
      CASE( np_ENT )                ! enstrophy conserving scheme (f-point)
         DO_2D( 0, 0, 0, 0 )
            z1_hu = ssumask(ji,jj) / ( phu(ji,jj) + 1._wp - ssumask(ji,jj) )
            z1_hv = ssvmask(ji,jj) / ( phv(ji,jj) + 1._wp - ssvmask(ji,jj) )
            zu_trd(ji,jj) = + r1_4 * r1_e1e2u(ji,jj) * z1_hu                    &
               &               * (  e1e2t(ji+1,jj)*pht(ji+1,jj)*ff_t(ji+1,jj) * ( pvnb(ji+1,jj) + pvnb(ji+1,jj-1) )   &
               &                  + e1e2t(ji  ,jj)*pht(ji  ,jj)*ff_t(ji  ,jj) * ( pvnb(ji  ,jj) + pvnb(ji  ,jj-1) )   )
               !
            zv_trd(ji,jj) = - r1_4 * r1_e1e2v(ji,jj) * z1_hv                    &
               &               * (  e1e2t(ji,jj+1)*pht(ji,jj+1)*ff_t(ji,jj+1) * ( punb(ji,jj+1) + punb(ji-1,jj+1) )   & 
               &                  + e1e2t(ji,jj  )*pht(ji,jj  )*ff_t(ji,jj  ) * ( punb(ji,jj  ) + punb(ji-1,jj  ) )   ) 
         END_2D
         !         
      CASE( np_ENE , np_MIX )        ! energy conserving scheme (t-point) ENE or MIX
         DO_2D( 0, 0, 0, 0 )
            zy1 = ( zhV(ji,jj-1) + zhV(ji+1,jj-1) ) * r1_e1u(ji,jj)
            zy2 = ( zhV(ji,jj  ) + zhV(ji+1,jj  ) ) * r1_e1u(ji,jj)
            zx1 = ( zhU(ji-1,jj) + zhU(ji-1,jj+1) ) * r1_e2v(ji,jj)
            zx2 = ( zhU(ji  ,jj) + zhU(ji  ,jj+1) ) * r1_e2v(ji,jj)
            ! energy conserving formulation for planetary vorticity term
            zu_trd(ji,jj) =   r1_4 * ( zwz(ji  ,jj-1) * zy1 + zwz(ji,jj) * zy2 )
            zv_trd(ji,jj) = - r1_4 * ( zwz(ji-1,jj  ) * zx1 + zwz(ji,jj) * zx2 )
         END_2D
         !
      CASE( np_ENS )                ! enstrophy conserving scheme (f-point)
         DO_2D( 0, 0, 0, 0 )
            zy1 =   r1_8 * ( zhV(ji  ,jj-1) + zhV(ji+1,jj-1) &
              &            + zhV(ji  ,jj  ) + zhV(ji+1,jj  ) ) * r1_e1u(ji,jj)
            zx1 = - r1_8 * ( zhU(ji-1,jj  ) + zhU(ji-1,jj+1) &
              &            + zhU(ji  ,jj  ) + zhU(ji  ,jj+1) ) * r1_e2v(ji,jj)
            zu_trd(ji,jj)  = zy1 * ( zwz(ji  ,jj-1) + zwz(ji,jj) )
            zv_trd(ji,jj)  = zx1 * ( zwz(ji-1,jj  ) + zwz(ji,jj) )
         END_2D
         !
      CASE( np_EET , np_EEN )      ! energy & enstrophy scheme (using e3t or e3f)         
         DO_2D( 0, 0, 0, 0 )
            zu_trd(ji,jj) = + r1_12 * r1_e1u(ji,jj) * (  ftne(ji,jj  ) * zhV(ji  ,jj  ) &
             &                                         + ftnw(ji+1,jj) * zhV(ji+1,jj  ) &
             &                                         + ftse(ji,jj  ) * zhV(ji  ,jj-1) &
             &                                         + ftsw(ji+1,jj) * zhV(ji+1,jj-1) )
            zv_trd(ji,jj) = - r1_12 * r1_e2v(ji,jj) * (  ftsw(ji,jj+1) * zhU(ji-1,jj+1) &
             &                                         + ftse(ji,jj+1) * zhU(ji  ,jj+1) &
             &                                         + ftnw(ji,jj  ) * zhU(ji-1,jj  ) &
             &                                         + ftne(ji,jj  ) * zhU(ji  ,jj  ) )
         END_2D
         !
      END SELECT
      !
   END SUBROUTINE dyn_cor_2D


   SUBROUTINE wad_tmsk( pssh, ptmsk )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE wad_tmsk  ***
      !!                    
      !! ** Purpose :   set wetting & drying mask at tracer points 
      !!              for the current barotropic sub-step 
      !!
      !! ** Method  :   ??? 
      !!
      !! ** Action  :  ptmsk : wetting & drying t-mask
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pssh    !
      REAL(wp), DIMENSION(jpi,jpj), INTENT(  out) ::   ptmsk   !
      !
      INTEGER  ::   ji, jj   ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      IF( ln_wd_dl_rmp ) THEN     
         DO_2D( 1, 1, 1, 1 )
            IF    ( pssh(ji,jj) + ht_0(ji,jj) >  2._wp * rn_wdmin1 ) THEN 
               !           IF    ( pssh(ji,jj) + ht_0(ji,jj) >          rn_wdmin2 ) THEN 
               ptmsk(ji,jj) = 1._wp
            ELSEIF( pssh(ji,jj) + ht_0(ji,jj) >          rn_wdmin1 ) THEN
               ptmsk(ji,jj) = TANH( 50._wp*( ( pssh(ji,jj) + ht_0(ji,jj) -  rn_wdmin1 )*r_rn_wdmin1) )
            ELSE 
               ptmsk(ji,jj) = 0._wp
            ENDIF
         END_2D
      ELSE  
         DO_2D( 1, 1, 1, 1 )
            IF ( pssh(ji,jj) + ht_0(ji,jj) >  rn_wdmin1 ) THEN   ;   ptmsk(ji,jj) = 1._wp
            ELSE                                                 ;   ptmsk(ji,jj) = 0._wp
            ENDIF
         END_2D
      ENDIF
      !
   END SUBROUTINE wad_tmsk


   SUBROUTINE wad_Umsk( pTmsk, phU, phV, pu, pv, pUmsk, pVmsk )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE wad_Umsk  ***
      !!                    
      !! ** Purpose :   set wetting & drying mask at tracer points 
      !!              for the current barotropic sub-step 
      !!
      !! ** Method  :   ??? 
      !!
      !! ** Action  :  ptmsk : wetting & drying t-mask
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pTmsk              ! W & D t-mask
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout)  :: phV, pu, pv! ocean velocities and transports
      REAL(dp), DIMENSION(jpi,jpj), INTENT(inout)  :: phU! ocean velocities and transports
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   pUmsk, pVmsk       ! W & D u- and v-mask
      !
      INTEGER  ::   ji, jj   ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      DO_2D( 1, 0, 1, 1 )   ! not jpi-column
         IF ( phU(ji,jj) > 0._wp ) THEN   ;   pUmsk(ji,jj) = pTmsk(ji  ,jj) 
         ELSE                             ;   pUmsk(ji,jj) = pTmsk(ji+1,jj)  
         ENDIF
         phU(ji,jj) = pUmsk(ji,jj)*phU(ji,jj)
         pu (ji,jj) = pUmsk(ji,jj)*pu (ji,jj)
      END_2D
      !
      DO_2D( 1, 1, 1, 0 )   ! not jpj-row
         IF ( phV(ji,jj) > 0._wp ) THEN   ;   pVmsk(ji,jj) = pTmsk(ji,jj  )
         ELSE                             ;   pVmsk(ji,jj) = pTmsk(ji,jj+1)  
         ENDIF
         phV(ji,jj) = pVmsk(ji,jj)*phV(ji,jj) 
         pv (ji,jj) = pVmsk(ji,jj)*pv (ji,jj)
      END_2D
      !
   END SUBROUTINE wad_Umsk


   SUBROUTINE dyn_drg_init(kt, Kbb, Kmm, puu, pvv, puu_b, pvv_b, pu_RHSi, pv_RHSi, pCdU_u, pCdU_v, &
              &             ptfrtrdu, ptfrtrdv, pbfrtrdu, pbfrtrdv )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_drg_init  ***
      !!                    
      !! ** Purpose : - add the baroclinic top/bottom drag contribution to 
      !!              the baroclinic part of the barotropic RHS
      !!              - compute the barotropic drag coefficients
      !!
      !! ** Method  :   computation done over the INNER domain only 
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT(in   ) ::  Kbb, Kmm           ! ocean time level indices
      REAL(dp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(in   ) ::  puu, pvv           ! ocean velocities and RHS of momentum equation
      REAL(wp), DIMENSION(jpi,jpj,jpt)    , INTENT(in   ) ::  puu_b, pvv_b       ! barotropic velocities at main time levels
      REAL(wp), DIMENSION(jpi,jpj)        , INTENT(inout) ::  pu_RHSi, pv_RHSi   ! baroclinic part of the barotropic RHS
      REAL(wp), DIMENSION(jpi,jpj)        , INTENT(  out) ::  pCdU_u , pCdU_v    ! barotropic drag coefficients
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout), OPTIONAL :: ptfrtrdu, ptfrtrdv ! top friction trends
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout), OPTIONAL :: pbfrtrdu, pbfrtrdv ! bottom friction trends
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      INTEGER  ::   ikbu, ikbv, iktu, iktv
      REAL(wp) ::   zztmp
      REAL(wp), DIMENSION(jpi,jpj) ::   zu_i, zv_i
!!jth_IWD
      INTEGER, INTENT(in)  ::   kt   ! ocean time-step index

      REAL(wp), SAVE ::   r1_25 = 0.04_wp   ! =1/25
      INTEGER                          ::   i_steps                           ! no of timesteps per hour
      INTEGER                          ::   i_is_25hrs                           ! flag if 1st 25hrs have passed (1=yes)
!!
      !!----------------------------------------------------------------------
      !
      !                    !==  Set the barotropic drag coef.  ==!
      !
      IF( ln_isfcav.OR.ln_drgice_imp ) THEN          ! top+bottom friction (ocean cavities)
         
         DO_2D( 0, 0, 0, 0 )
            pCdU_u(ji,jj) = r1_2*( rCdU_bot(ji+1,jj)+rCdU_bot(ji,jj) + rCdU_top(ji+1,jj)+rCdU_top(ji,jj) )
            pCdU_v(ji,jj) = r1_2*( rCdU_bot(ji,jj+1)+rCdU_bot(ji,jj) + rCdU_top(ji,jj+1)+rCdU_top(ji,jj) )
         END_2D
      ELSE                          ! bottom friction only
         DO_2D( 0, 0, 0, 0 )
            pCdU_u(ji,jj) = r1_2*( rCdU_bot(ji+1,jj)+rCdU_bot(ji,jj) )
            pCdU_v(ji,jj) = r1_2*( rCdU_bot(ji,jj+1)+rCdU_bot(ji,jj) )
         END_2D
      ENDIF
      !
      !                    !==  BOTTOM stress contribution from baroclinic velocities  ==!
      !
      IF( l_trddyn ) THEN
           pbfrtrdu(:,:) = pu_RHSi(:,:)
           pbfrtrdv(:,:) = pv_RHSi(:,:)
      ENDIF
      !
      IF( ln_bt_fw ) THEN                 ! FORWARD integration: use NOW bottom baroclinic velocities
         
         DO_2D( 0, 0, 0, 0 )
            ikbu = mbku(ji,jj)       
            ikbv = mbkv(ji,jj)    
            zu_i(ji,jj) = puu(ji,jj,ikbu,Kmm) - puu_b(ji,jj,Kmm)
            zv_i(ji,jj) = pvv(ji,jj,ikbv,Kmm) - pvv_b(ji,jj,Kmm)
         END_2D
      ELSE                                ! CENTRED integration: use BEFORE bottom baroclinic velocities
         
         DO_2D( 0, 0, 0, 0 )
            ikbu = mbku(ji,jj)       
            ikbv = mbkv(ji,jj)    
            zu_i(ji,jj) = puu(ji,jj,ikbu,Kbb) - puu_b(ji,jj,Kbb)
            zv_i(ji,jj) = pvv(ji,jj,ikbv,Kbb) - pvv_b(ji,jj,Kbb)
         END_2D
      ENDIF
      !
      ! use "unclipped" drag (even if explicit friction is used in 3D calculation)
         
      DO_2D( 0, 0, 0, 0 )
         pu_RHSi(ji,jj) = pu_RHSi(ji,jj) + r1_hu(ji,jj,Kmm) * r1_2*( rCdU_bot(ji+1,jj)+rCdU_bot(ji,jj) ) * zu_i(ji,jj)
         pv_RHSi(ji,jj) = pv_RHSi(ji,jj) + r1_hv(ji,jj,Kmm) * r1_2*( rCdU_bot(ji,jj+1)+rCdU_bot(ji,jj) ) * zv_i(ji,jj)
      END_2D
      !
      IF( l_trddyn ) THEN
         pbfrtrdu(:,:) = pu_RHSi(:,:) - pbfrtrdu(:,:) 
         pbfrtrdv(:,:) = pv_RHSi(:,:) - pbfrtrdv(:,:) 
      ENDIF
      !                    !==  TOP stress contribution from baroclinic velocities  ==!   (no W/D case)
      !
      IF( ln_isfcav.OR.ln_drgice_imp ) THEN
         !
         IF( l_trddyn ) THEN
            ptfrtrdu(:,:) = pu_RHSi(:,:)
            ptfrtrdv(:,:) = pv_RHSi(:,:)
         ENDIF
         !
         IF( ln_bt_fw ) THEN                ! FORWARD integration: use NOW top baroclinic velocity
            
            DO_2D( 0, 0, 0, 0 )
               iktu = miku(ji,jj)
               iktv = mikv(ji,jj)
               zu_i(ji,jj) = puu(ji,jj,iktu,Kmm) - puu_b(ji,jj,Kmm)
               zv_i(ji,jj) = pvv(ji,jj,iktv,Kmm) - pvv_b(ji,jj,Kmm)
            END_2D
         ELSE                                ! CENTRED integration: use BEFORE top baroclinic velocity
            
            DO_2D( 0, 0, 0, 0 )
               iktu = miku(ji,jj)
               iktv = mikv(ji,jj)
               zu_i(ji,jj) = puu(ji,jj,iktu,Kbb) - puu_b(ji,jj,Kbb)
               zv_i(ji,jj) = pvv(ji,jj,iktv,Kbb) - pvv_b(ji,jj,Kbb)
            END_2D
         ENDIF
         !
         !                    ! use "unclipped" top drag (even if explicit friction is used in 3D calculation)
         
         DO_2D( 0, 0, 0, 0 )
            pu_RHSi(ji,jj) = pu_RHSi(ji,jj) + r1_hu(ji,jj,Kmm) * r1_2*( rCdU_top(ji+1,jj)+rCdU_top(ji,jj) ) * zu_i(ji,jj)
            pv_RHSi(ji,jj) = pv_RHSi(ji,jj) + r1_hv(ji,jj,Kmm) * r1_2*( rCdU_top(ji,jj+1)+rCdU_top(ji,jj) ) * zv_i(ji,jj)
         END_2D
         !
         IF( l_trddyn ) THEN
            ptfrtrdu(:,:) = pu_RHSi(:,:) - ptfrtrdu(:,:) 
            ptfrtrdv(:,:) = pv_RHSi(:,:) - ptfrtrdv(:,:)
         ENDIF
         !
      ENDIF
      !
!jth IWD!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Calculate detided barotropic velocity for use with internal wave drag
      IF ( ln_int_wave_drag ) THEN
       IF( MOD( 3600,NINT(rdt) ) == 0 ) THEN
          i_steps = 3600/NINT(rdt)
       ELSE
          CALL ctl_stop('STOP', "dynspg_ts int_wave_drag : timestep must give MOD(3600,rdt) = 0 otherwise no hourly values are possible")
       ENDIF

!cw_IWD start
! The following loop may no longer be necessary, as un_detide, vn_detide are
! now initialised to zero and are only updated if at least 25 hours after the
! initialisation time or if restarted and stored variables read in.
!cw_IWD end

! is it first hour in a initialised (i.e. non-restarted) run?  If so, set un_detide, vn_detide to zero.
       IF (kt < i_steps) THEN
          un_detide(:,:) = 0.0 ! un_b(:,:)
          vn_detide(:,:) = 0.0 ! vn_b(:,:)
       ENDIF

!cw_IWD start
! Doing restarts with un_25, un_detide, ... now, so amend following line to:
!       IF( MOD( kt, i_steps ) == 0  .AND. kt /= nn_it000 ) THEN
       IF( MOD( kt, i_steps ) == 0  ) THEN ! XXXXXXXXX HOURLY LOOP START XXXXXX
!cw_IWD end

! if at least 25 hrs after the initialisation time (from kt=1), then set
! i_is_25hrs =1 else 0.
                IF (kt >= 25 * i_steps) THEN
                  i_is_25hrs = 1 !start accumulating hourly values into window
                ELSE
                  i_is_25hrs = 0 !do not accumulate yet
                ENDIF

!
         DO_2D( 0, 0, 0, 0 )
        ! Shift arrays forward by 1 hr
                        un_25(ji,jj,1:24) = un_25(ji,jj,2:25)
                        vn_25(ji,jj,1:24) = vn_25(ji,jj,2:25)
        ! fill last value
                        un_25(ji,jj,25) = puu_b(ji,jj,Kmm)
                        vn_25(ji,jj,25) = pvv_b(ji,jj,Kmm)
        ! set un_detide if more than 25hrs since start of run, otherwise set to un_b
                        un_detide(ji,jj) = r1_25 * sum(un_25(ji,jj,:)) * i_is_25hrs + ( 1 - i_is_25hrs ) * 0.0 !un_b(ji,jj)
                        vn_detide(ji,jj) = r1_25 * sum(vn_25(ji,jj,:)) * i_is_25hrs + ( 1 - i_is_25hrs ) * 0.0 !vn_b(ji,jj)
        END_2D

                 IF ( ln_calc_tdiss ) THEN
!!!! Simmons et al. Ocean Modelling (2004)
!jth              N2mean(:,:) =  SUM( e3w_n(:,:,:) * rn2(:,:,:) * wmask(:,:,:) , DIM=3 ) / SUM( e3w_n(:,:,:) * wmask(:,:,:) , DIM=3 )  * tmask(:,:,1)
!jth avoid divide by zero
                 N2mean(:,:) = 0.0
                 DO jk = 1, jpkm1
                    N2mean(:,:)= N2mean(:,:) + e3w(:,:,jk,Kmm) * rn2(:,:,jk) * wmask(:,:,jk)
                 END DO
                 N2mean(:,:)= N2mean(:,:) / (ht(:,:) + 1._wp - tmask(:,:,1) ) *  tmask(:,:,1)
                 !N2mean(:,:) =  SUM( e3w(:,:,:,Kmm) * rn2(:,:,:) * wmask(:,:,:) , DIM=3 ) / (ht(:,:) + 1._wp - tmask(:,:,1) ) * tmask(:,:,1)

         DO_2D( 1, 0, 1, 0 )  !need to check  original is DO jj = 2, jpj
                      tdiss(ji,jj) = 0.5 * rn_kappa_tdiss * h2rough(ji,jj)*sqrt(max(N2mean(ji,jj),0.0))
         END_2D
                 END IF
!!!!!!!!!!!!!!!!!!!!
                 !jidbg= 477
                 !jjdbg= 170

                 !ji=jidbg-nimpp+1
                 !jj=jjdbg-njmpp+1
                 !if (ji>0 .and. ji < jpi .and. jj > 0 .and. jj < jpj)  then
                 !write(68,'(1I5,4e12.4)') kt,tdiss(ji,jj),h2rough(ji,jj),sqrt(max(N2mean(ji,jj),0.0)), wmask(ji,jj,imk)
                 !write(68,*) kt,(ht_n(ji,jj) + 1._wp - tmask(ji,jj,1) ),SUM( e3w_n(ji,jj,:) * wmask(ji,jj,:) , DIM=1 )

                 !    flush(68)
                 !endif


                ENDIF
      ENDIF              ! XXXXXXXXXXX END OF HOURLY LOOP XXXXXX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   END SUBROUTINE dyn_drg_init

   SUBROUTINE ts_bck_interp( jn, ll_init,       &   ! <== in
      &                      za0, za1, za2, za3 )   ! ==> out
      !!----------------------------------------------------------------------
      INTEGER ,INTENT(in   ) ::   jn                   ! index of sub time step
      LOGICAL ,INTENT(in   ) ::   ll_init              !
      REAL(wp),INTENT(  out)  :: za0, za2, za3! Half-step back interpolation coefficient
      REAL(dp),INTENT(  out)  :: za1! Half-step back interpolation coefficient
      !
      REAL(wp) ::   zepsilon, zgamma                   !   -      -
      !!----------------------------------------------------------------------
      !                             ! set Half-step back interpolation coefficient
      IF    ( jn==1 .AND. ll_init ) THEN   !* Forward-backward
         za0 = 1._wp                        
         za1 = 0._wp                           
         za2 = 0._wp
         za3 = 0._wp
      ELSEIF( jn==2 .AND. ll_init ) THEN   !* AB2-AM3 Coefficients; bet=0 ; gam=-1/6 ; eps=1/12
         za0 = 1.0833333333333_wp                 ! za0 = 1-gam-eps
         za1 =-0.1666666666666_wp                 ! za1 = gam
         za2 = 0.0833333333333_wp                 ! za2 = eps
         za3 = 0._wp              
      ELSE                                 !* AB3-AM4 Coefficients; bet=0.281105 ; eps=0.013 ; gam=0.0880 
         IF( rn_bt_alpha == 0._wp ) THEN      ! Time diffusion  
            za0 = 0.614_wp                        ! za0 = 1/2 +   gam + 2*eps
            za1 = 0.285_wp                        ! za1 = 1/2 - 2*gam - 3*eps
            za2 = 0.088_wp                        ! za2 = gam
            za3 = 0.013_wp                        ! za3 = eps
         ELSE                                 ! no time diffusion
            zepsilon = 0.00976186_wp - 0.13451357_wp * rn_bt_alpha
            zgamma   = 0.08344500_wp - 0.51358400_wp * rn_bt_alpha
            za0 = 0.5_wp + zgamma + 2._wp * rn_bt_alpha + 2._wp * zepsilon
            za1 = 1._wp - za0 - zgamma - zepsilon
            za2 = zgamma
            za3 = zepsilon
         ENDIF 
      ENDIF
   END SUBROUTINE ts_bck_interp


   !!======================================================================
END MODULE dynspg_ts
