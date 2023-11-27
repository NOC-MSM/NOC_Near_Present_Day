MODULE iceupdate
   !!======================================================================
   !!                       ***  MODULE iceupdate   ***
   !!  Sea-ice :   computation of the flux at the sea ice/ocean interface
   !!======================================================================
   !! History :  4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_update_alloc : allocate the iceupdate arrays
   !!   ice_update_init  : initialisation
   !!   ice_update_flx   : updates mass, heat and salt fluxes at the ocean surface
   !!   ice_update_tau   : update i- and j-stresses, and its modulus at the ocean surface
   !!----------------------------------------------------------------------
   USE phycst         ! physical constants
   USE dom_oce        ! ocean domain
   USE ice            ! sea-ice: variables
   USE sbc_ice        ! Surface boundary condition: ice   fields
   USE sbc_oce        ! Surface boundary condition: ocean fields
   USE sbccpl         ! Surface boundary condition: coupled interface
   USE icealb         ! sea-ice: albedo parameters
   USE traqsr         ! add penetration of solar flux in the calculation of heat budget
   USE icectl         ! sea-ice: control prints
   USE zdfdrg  , ONLY : ln_drgice_imp
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)
   USE lbclnk         ! lateral boundary conditions (or mpp links)
   USE timing         ! Timing
   USE oce
   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_update_init   ! called by ice_init
   PUBLIC   ice_update_flx    ! called by ice_stp
   PUBLIC   ice_update_tau    ! called by ice_stp

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   utau_oce, vtau_oce   ! air-ocean surface i- & j-stress     [N/m2]
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   tmod_io              ! modulus of the ice-ocean velocity   [m/s]

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/ICE 4.0 , NEMO Consortium (2018)
   !! $Id: iceupdate.F90 15385 2021-10-15 13:52:48Z clem $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION ice_update_alloc()
      !!-------------------------------------------------------------------
      !!             ***  ROUTINE ice_update_alloc ***
      !!-------------------------------------------------------------------
      ALLOCATE( utau_oce(jpi,jpj), vtau_oce(jpi,jpj), tmod_io(jpi,jpj), STAT=ice_update_alloc )
      !
      CALL mpp_sum( 'iceupdate', ice_update_alloc )
      IF( ice_update_alloc /= 0 )   CALL ctl_stop( 'STOP', 'ice_update_alloc: failed to allocate arrays' )
      !
   END FUNCTION ice_update_alloc


   SUBROUTINE ice_update_flx( kt )
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE ice_update_flx ***
      !!
      !! ** Purpose :   Update the surface ocean boundary condition for heat
      !!                salt and mass over areas where sea-ice is non-zero
      !!
      !! ** Action  : - computes the heat and freshwater/salt fluxes
      !!                at the ice-ocean interface.
      !!              - Update the ocean sbc
      !!
      !! ** Outputs : - qsr     : sea heat flux:     solar
      !!              - qns     : sea heat flux: non solar
      !!              - emp     : freshwater budget: volume flux
      !!              - sfx     : salt flux
      !!              - fr_i    : ice fraction
      !!              - tn_ice  : sea-ice surface temperature
      !!              - alb_ice : sea-ice albedo (recomputed only for coupled mode)
      !!
      !! References : Goosse, H. et al. 1996, Bul. Soc. Roy. Sc. Liege, 65, 87-90.
      !!              Tartinville et al. 2001 Ocean Modelling, 3, 95-108.
      !!              These refs are now obsolete since everything has been revised
      !!              The ref should be Rousset et al., 2015
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! number of iteration
      !
      INTEGER  ::   ji, jj, jl, jk   ! dummy loop indices
      REAL(wp) ::   zqsr             ! New solar flux received by the ocean
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   z2d                  ! 2D workspace
      !!---------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('iceupdate')

      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*)'ice_update_flx: update fluxes (mass, salt and heat) at the ice-ocean interface'
         WRITE(numout,*)'~~~~~~~~~~~~~~'
      ENDIF

      ! Net heat flux on top of the ice-ocean (W.m-2)
      !----------------------------------------------
      IF( ln_cndflx ) THEN   ! ice-atm interface = conduction (and melting) fluxes
         qt_atm_oi(:,:) = ( 1._wp - at_i_b(:,:) ) * ( qns_oce(:,:) + qsr_oce(:,:) ) + qemp_oce(:,:) + &
            &             SUM( a_i_b * ( qcn_ice + qml_ice + qtr_ice_top ), dim=3 ) + qemp_ice(:,:)
      ELSE                   ! ice-atm interface = solar and non-solar fluxes
         qt_atm_oi(:,:) = qns_tot(:,:) + qsr_tot(:,:) 
      ENDIF

      ! --- case we bypass ice thermodynamics --- !
      IF( .NOT. ln_icethd ) THEN   ! we suppose ice is impermeable => ocean is isolated from atmosphere
         qt_atm_oi  (:,:)   = ( 1._wp - at_i_b(:,:) ) * ( qns_oce(:,:) + qsr_oce(:,:) ) + qemp_oce(:,:)
         qt_oce_ai  (:,:)   = ( 1._wp - at_i_b(:,:) ) *   qns_oce(:,:)                  + qemp_oce(:,:)
         emp_ice    (:,:)   = 0._wp
         qemp_ice   (:,:)   = 0._wp
         qevap_ice  (:,:,:) = 0._wp
      ENDIF

      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )

         ! Solar heat flux reaching the ocean (max) = zqsr (W.m-2)
         !---------------------------------------------------
         IF( ln_cndflx ) THEN   ! ice-atm interface = conduction (and melting) fluxes
            zqsr = ( 1._wp - at_i_b(ji,jj) ) * qsr_oce(ji,jj) + SUM( a_i_b (ji,jj,:) * qtr_ice_bot(ji,jj,:) )
         ELSE                   ! ice-atm interface = solar and non-solar fluxes
            zqsr = qsr_tot(ji,jj) - SUM( a_i_b(ji,jj,:) * ( qsr_ice(ji,jj,:) - qtr_ice_bot(ji,jj,:) ) )
         ENDIF

         ! Total heat flux reaching the ocean = qt_oce_ai (W.m-2)
         !---------------------------------------------------
         IF( ln_icethd ) THEN
            qt_oce_ai(ji,jj) = qt_atm_oi(ji,jj) - hfx_sum(ji,jj) - hfx_bom(ji,jj) - hfx_bog(ji,jj) &
               &                                - hfx_dif(ji,jj) - hfx_opw(ji,jj) - hfx_snw(ji,jj) &
               &                                + hfx_thd(ji,jj) + hfx_dyn(ji,jj) + hfx_res(ji,jj) &
               &                                + hfx_sub(ji,jj) - SUM( qevap_ice(ji,jj,:) * a_i_b(ji,jj,:) ) + hfx_spr(ji,jj)
         ENDIF
         
         ! New qsr and qns used to compute the oceanic heat flux at the next time step
         !----------------------------------------------------------------------------
         ! if warming and some ice remains, then we suppose that the whole solar flux has been consumed to melt the ice
         ! else ( cooling or no ice left ), then we suppose that     no    solar flux has been consumed
         !
         IF( fhld(ji,jj) > 0._wp .AND. at_i(ji,jj) > 0._wp ) THEN   !-- warming and some ice remains
            !                                        solar flux transmitted thru the 1st level of the ocean (i.e. not used by sea-ice)
            qsr(ji,jj) = ( 1._wp - at_i_b(ji,jj) ) * qsr_oce(ji,jj) * ( 1._wp - frq_m(ji,jj) ) &
               !                                   + solar flux transmitted thru ice and the 1st ocean level (also not used by sea-ice)
               &             + SUM( a_i_b(ji,jj,:) * qtr_ice_bot(ji,jj,:) ) * ( 1._wp - frq_m(ji,jj) )
            !
         ELSE                                                       !-- cooling or no ice left
            qsr(ji,jj) = zqsr
         ENDIF
         !
         ! the non-solar is simply derived from the solar flux
         qns(ji,jj) = qt_oce_ai(ji,jj) - qsr(ji,jj)

         ! Mass flux at the atm. surface
         !-----------------------------------
         wfx_sub(ji,jj) = wfx_snw_sub(ji,jj) + wfx_ice_sub(ji,jj)

         ! Mass flux at the ocean surface
         !------------------------------------
         ! ice-ocean  mass flux
         wfx_ice(ji,jj) = wfx_bog(ji,jj) + wfx_bom(ji,jj) + wfx_sum(ji,jj) + wfx_sni(ji,jj)   &
            &           + wfx_opw(ji,jj) + wfx_dyn(ji,jj) + wfx_res(ji,jj) + wfx_lam(ji,jj)

         ! snw-ocean mass flux
         wfx_snw(ji,jj) = wfx_snw_sni(ji,jj) + wfx_snw_dyn(ji,jj) + wfx_snw_sum(ji,jj)

         ! total mass flux at the ocean/ice interface
         fmmflx(ji,jj) =                - wfx_ice(ji,jj) - wfx_snw(ji,jj) - wfx_pnd(ji,jj) - wfx_err_sub(ji,jj)   ! ice-ocean mass flux saved at least for biogeochemical model
         emp   (ji,jj) = emp_oce(ji,jj) - wfx_ice(ji,jj) - wfx_snw(ji,jj) - wfx_pnd(ji,jj) - wfx_err_sub(ji,jj)   ! atm-ocean + ice-ocean mass flux

         ! Salt flux at the ocean surface
         !------------------------------------------
         sfx(ji,jj) = sfx_bog(ji,jj) + sfx_bom(ji,jj) + sfx_sum(ji,jj) + sfx_sni(ji,jj) + sfx_opw(ji,jj)   &
            &       + sfx_res(ji,jj) + sfx_dyn(ji,jj) + sfx_bri(ji,jj) + sfx_sub(ji,jj) + sfx_lam(ji,jj)

         ! Mass of snow and ice per unit area
         !----------------------------------------
         snwice_mass_b(ji,jj) = snwice_mass(ji,jj)       ! save mass from the previous ice time step
         !                                               ! new mass per unit area
         snwice_mass  (ji,jj) = tmask(ji,jj,1) * ( rhos * vt_s(ji,jj) + rhoi * vt_i(ji,jj) + rhow * (vt_ip(ji,jj) + vt_il(ji,jj)) )
         !                                               ! time evolution of snow+ice mass
         snwice_fmass (ji,jj) = ( snwice_mass(ji,jj) - snwice_mass_b(ji,jj) ) * r1_Dt_ice

      END_2D

      ! Storing the transmitted variables
      !----------------------------------
      fr_i  (:,:)   = at_i(:,:)             ! Sea-ice fraction
      tn_ice(:,:,:) = t_su(:,:,:)           ! Ice surface temperature

      ! Snow/ice albedo (only if sent to coupler, useless in forced mode)
      !------------------------------------------------------------------
      CALL ice_alb( t_su, h_i, h_s, ln_pnd_alb, a_ip_eff, h_ip, cloud_fra, alb_ice ) ! ice albedo

      !
      IF( lrst_ice ) THEN                       !* write snwice_mass fields in the restart file
         CALL update_rst( 'WRITE', kt )
      ENDIF
      !
      ! output all fluxes
      !------------------
      !
      ! --- salt fluxes [kg/m2/s] --- !
      !                           ! sfxice =  sfxbog + sfxbom + sfxsum + sfxsni + sfxopw + sfxres + sfxdyn + sfxbri + sfxsub + sfxlam
      IF( iom_use('sfxice'  ) )   CALL iom_put( 'sfxice', sfx     * 1.e-03 )   ! salt flux from total ice growth/melt
      IF( iom_use('sfxbog'  ) )   CALL iom_put( 'sfxbog', sfx_bog * 1.e-03 )   ! salt flux from bottom growth
      IF( iom_use('sfxbom'  ) )   CALL iom_put( 'sfxbom', sfx_bom * 1.e-03 )   ! salt flux from bottom melting
      IF( iom_use('sfxsum'  ) )   CALL iom_put( 'sfxsum', sfx_sum * 1.e-03 )   ! salt flux from surface melting
      IF( iom_use('sfxlam'  ) )   CALL iom_put( 'sfxlam', sfx_lam * 1.e-03 )   ! salt flux from lateral melting
      IF( iom_use('sfxsni'  ) )   CALL iom_put( 'sfxsni', sfx_sni * 1.e-03 )   ! salt flux from snow ice formation
      IF( iom_use('sfxopw'  ) )   CALL iom_put( 'sfxopw', sfx_opw * 1.e-03 )   ! salt flux from open water formation
      IF( iom_use('sfxdyn'  ) )   CALL iom_put( 'sfxdyn', sfx_dyn * 1.e-03 )   ! salt flux from ridging rafting
      IF( iom_use('sfxbri'  ) )   CALL iom_put( 'sfxbri', sfx_bri * 1.e-03 )   ! salt flux from brines
      IF( iom_use('sfxres'  ) )   CALL iom_put( 'sfxres', sfx_res * 1.e-03 )   ! salt flux from undiagnosed processes
      IF( iom_use('sfxsub'  ) )   CALL iom_put( 'sfxsub', sfx_sub * 1.e-03 )   ! salt flux from sublimation

      ! --- mass fluxes [kg/m2/s] --- !
      CALL iom_put( 'emp_oce', emp_oce )   ! emp over ocean (taking into account the snow blown away from the ice)
      CALL iom_put( 'emp_ice', emp_ice )   ! emp over ice   (taking into account the snow blown away from the ice)

      !                           ! vfxice = vfxbog + vfxbom + vfxsum + vfxsni + vfxopw + vfxdyn + vfxres + vfxlam + vfxpnd
      CALL iom_put( 'vfxice'    , wfx_ice     )   ! mass flux from total ice growth/melt
      CALL iom_put( 'vfxbog'    , wfx_bog     )   ! mass flux from bottom growth
      CALL iom_put( 'vfxbom'    , wfx_bom     )   ! mass flux from bottom melt
      CALL iom_put( 'vfxsum'    , wfx_sum     )   ! mass flux from surface melt
      CALL iom_put( 'vfxlam'    , wfx_lam     )   ! mass flux from lateral melt
      CALL iom_put( 'vfxsni'    , wfx_sni     )   ! mass flux from snow-ice formation
      CALL iom_put( 'vfxopw'    , wfx_opw     )   ! mass flux from growth in open water
      CALL iom_put( 'vfxdyn'    , wfx_dyn     )   ! mass flux from dynamics (ridging)
      CALL iom_put( 'vfxres'    , wfx_res     )   ! mass flux from undiagnosed processes
      CALL iom_put( 'vfxpnd'    , wfx_pnd     )   ! mass flux from melt ponds
      CALL iom_put( 'vfxsub'    , wfx_ice_sub )   ! mass flux from ice sublimation (ice-atm.)
      CALL iom_put( 'vfxsub_err', wfx_err_sub )   ! "excess" of sublimation sent to ocean

      IF ( iom_use( 'vfxthin' ) ) THEN   ! mass flux from ice growth in open water + thin ice (<20cm) => comparable to observations
         ALLOCATE( z2d(jpi,jpj) )
         WHERE( hm_i(:,:) < 0.2 .AND. hm_i(:,:) > 0. ) ; z2d = wfx_bog
         ELSEWHERE                                     ; z2d = 0._wp
         END WHERE
         CALL iom_put( 'vfxthin', wfx_opw + z2d )
         DEALLOCATE( z2d )
      ENDIF

      !                            ! vfxsnw = vfxsnw_sni + vfxsnw_dyn + vfxsnw_sum
      CALL iom_put( 'vfxsnw'     , wfx_snw     )   ! mass flux from total snow growth/melt
      CALL iom_put( 'vfxsnw_sum' , wfx_snw_sum )   ! mass flux from snow melt at the surface
      CALL iom_put( 'vfxsnw_sni' , wfx_snw_sni )   ! mass flux from snow melt during snow-ice formation
      CALL iom_put( 'vfxsnw_dyn' , wfx_snw_dyn )   ! mass flux from dynamics (ridging)
      CALL iom_put( 'vfxsnw_sub' , wfx_snw_sub )   ! mass flux from snow sublimation (ice-atm.)
      CALL iom_put( 'vfxsnw_pre' , wfx_spr     )   ! snow precip

      ! --- heat fluxes [W/m2] --- !
      !                              ! qt_atm_oi - qt_oce_ai = hfxdhc - ( dihctrp + dshctrp )
      IF( iom_use('qsr_oce'    ) )   CALL iom_put( 'qsr_oce'    , qsr_oce * ( 1._wp - at_i_b )                               )   !     solar flux at ocean surface
      IF( iom_use('qns_oce'    ) )   CALL iom_put( 'qns_oce'    , qns_oce * ( 1._wp - at_i_b ) + qemp_oce                    )   ! non-solar flux at ocean surface
      IF( iom_use('qsr_ice'    ) )   CALL iom_put( 'qsr_ice'    , SUM( qsr_ice * a_i_b, dim=3 )                              )   !     solar flux at ice surface
      IF( iom_use('qns_ice'    ) )   CALL iom_put( 'qns_ice'    , SUM( qns_ice * a_i_b, dim=3 ) + qemp_ice                   )   ! non-solar flux at ice surface
      IF( iom_use('qtr_ice_bot') )   CALL iom_put( 'qtr_ice_bot', SUM( qtr_ice_bot * a_i_b, dim=3 )                          )   !     solar flux transmitted thru ice
      IF( iom_use('qtr_ice_top') )   CALL iom_put( 'qtr_ice_top', SUM( qtr_ice_top * a_i_b, dim=3 )                          )   !     solar flux transmitted thru ice surface
      IF( iom_use('qt_oce'     ) )   CALL iom_put( 'qt_oce'     ,      ( qsr_oce + qns_oce ) * ( 1._wp - at_i_b ) + qemp_oce )
      IF( iom_use('qt_ice'     ) )   CALL iom_put( 'qt_ice'     , SUM( ( qns_ice + qsr_ice ) * a_i_b, dim=3 )     + qemp_ice )
      IF( iom_use('qt_oce_ai'  ) )   CALL iom_put( 'qt_oce_ai'  , qt_oce_ai * tmask(:,:,1)                                   )   ! total heat flux at the ocean   surface: interface oce-(ice+atm)
      IF( iom_use('qt_atm_oi'  ) )   CALL iom_put( 'qt_atm_oi'  , qt_atm_oi * tmask(:,:,1)                                   )   ! total heat flux at the oce-ice surface: interface atm-(ice+oce)
      IF( iom_use('qemp_oce'   ) )   CALL iom_put( 'qemp_oce'   , qemp_oce                                                   )   ! Downward Heat Flux from E-P over ocean
      IF( iom_use('qemp_ice'   ) )   CALL iom_put( 'qemp_ice'   , qemp_ice                                                   )   ! Downward Heat Flux from E-P over ice

      ! heat fluxes from ice transformations
      !                            ! hfxdhc = hfxbog + hfxbom + hfxsum + hfxopw + hfxdif + hfxsnw - ( hfxthd + hfxdyn + hfxres + hfxsub + hfxspr )
      CALL iom_put ('hfxbog'     , hfx_bog     )   ! heat flux used for ice bottom growth
      CALL iom_put ('hfxbom'     , hfx_bom     )   ! heat flux used for ice bottom melt
      CALL iom_put ('hfxsum'     , hfx_sum     )   ! heat flux used for ice surface melt
      CALL iom_put ('hfxopw'     , hfx_opw     )   ! heat flux used for ice formation in open water
      CALL iom_put ('hfxdif'     , hfx_dif     )   ! heat flux used for ice temperature change
      CALL iom_put ('hfxsnw'     , hfx_snw     )   ! heat flux used for snow melt
      CALL iom_put ('hfxerr'     , hfx_err_dif )   ! heat flux error after heat diffusion

      ! heat fluxes associated with mass exchange (freeze/melt/precip...)
      CALL iom_put ('hfxthd'     , hfx_thd     )   !
      CALL iom_put ('hfxdyn'     , hfx_dyn     )   !
      CALL iom_put ('hfxres'     , hfx_res     )   !
      CALL iom_put ('hfxsub'     , hfx_sub     )   !
      CALL iom_put ('hfxspr'     , hfx_spr     )   ! Heat flux from snow precip heat content

      ! other heat fluxes
      IF( iom_use('hfxsensib'  ) )   CALL iom_put( 'hfxsensib'  ,      qsb_ice_bot * at_i_b         )   ! Sensible oceanic heat flux
      IF( iom_use('hfxcndbot'  ) )   CALL iom_put( 'hfxcndbot'  , SUM( qcn_ice_bot * a_i_b, dim=3 ) )   ! Bottom conduction flux
      IF( iom_use('hfxcndtop'  ) )   CALL iom_put( 'hfxcndtop'  , SUM( qcn_ice_top * a_i_b, dim=3 ) )   ! Surface conduction flux
      IF( iom_use('hfxmelt'    ) )   CALL iom_put( 'hfxmelt'    , SUM( qml_ice     * a_i_b, dim=3 ) )   ! Surface melt flux
      IF( iom_use('hfxldmelt'  ) )   CALL iom_put( 'hfxldmelt'  ,      fhld        * at_i_b         )   ! Heat in lead for ice melting 
      IF( iom_use('hfxldgrow'  ) )   CALL iom_put( 'hfxldgrow'  ,      qlead       * r1_Dt_ice      )   ! Heat in lead for ice growth

      ! controls
      !---------
#if ! defined key_agrif
      IF( ln_icediachk      )   CALL ice_cons_final('iceupdate')                                       ! conservation
#endif
      IF( ln_icectl         )   CALL ice_prt       (kt, iiceprt, jiceprt, 3, 'Final state ice_update') ! prints
      IF( sn_cfctl%l_prtctl )   CALL ice_prt3D     ('iceupdate')                                       ! prints
      IF( ln_timing         )   CALL timing_stop   ('iceupdate')                                       ! timing
      !
   END SUBROUTINE ice_update_flx


   SUBROUTINE ice_update_tau( kt, pu_oce, pv_oce )
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE ice_update_tau ***
      !!
      !! ** Purpose : Update the ocean surface stresses due to the ice
      !!
      !! ** Action  : * at each ice time step (every nn_fsbc time step):
      !!                - compute the modulus of ice-ocean relative velocity
      !!                  (*rho*Cd) at T-point (C-grid) or I-point (B-grid)
      !!                      tmod_io = rhoco * | U_ice-U_oce |
      !!                - update the modulus of stress at ocean surface
      !!                      taum = (1-a) * taum + a * tmod_io * | U_ice-U_oce |
      !!              * at each ocean time step (every kt):
      !!                  compute linearized ice-ocean stresses as
      !!                      Utau = tmod_io * | U_ice - pU_oce |
      !!                using instantaneous current ocean velocity (usually before)
      !!
      !!    NB: - ice-ocean rotation angle no more allowed
      !!        - here we make an approximation: taum is only computed every ice time step
      !!          This avoids mutiple average to pass from T -> U,V grids and next from U,V grids
      !!          to T grid. taum is used in TKE and GLS, which should not be too sensitive to this approximaton...
      !!
      !! ** Outputs : - utau, vtau   : surface ocean i- and j-stress (u- & v-pts) updated with ice-ocean fluxes
      !!              - taum         : modulus of the surface ocean stress (T-point) updated with ice-ocean fluxes
      !!---------------------------------------------------------------------
      INTEGER ,                     INTENT(in) ::   kt               ! ocean time-step index
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pu_oce, pv_oce   ! surface ocean currents
      !
      INTEGER  ::   ji, jj   ! dummy loop indices
      REAL(wp) ::   zat_u, zutau_ice, zu_t, zmodt   ! local scalar
      REAL(wp) ::   zat_v, zvtau_ice, zv_t, zrhoco  !   -      -
      REAL(wp) ::   zflagi                          !   -      -
      !!---------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('iceupdate')

      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*)'ice_update_tau: update stress at the ice-ocean interface'
         WRITE(numout,*)'~~~~~~~~~~~~~~'
      ENDIF

      zrhoco = rho0 * rn_cio
      !
      IF( MOD( kt-1, nn_fsbc ) == 0 ) THEN     !==  Ice time-step only  ==!   (i.e. surface module time-step)
         DO_2D( 0, 0, 0, 0 )                          !* update the modulus of stress at ocean surface (T-point)
            !                                               ! 2*(U_ice-U_oce) at T-point
            zu_t = u_ice(ji,jj) + u_ice(ji-1,jj) - u_oce(ji,jj) - u_oce(ji-1,jj)
            zv_t = v_ice(ji,jj) + v_ice(ji,jj-1) - v_oce(ji,jj) - v_oce(ji,jj-1)
            !                                              ! |U_ice-U_oce|^2
            zmodt =  0.25_wp * (  zu_t * zu_t + zv_t * zv_t  )
            !                                               ! update the ocean stress modulus
            taum(ji,jj) = ( 1._wp - at_i(ji,jj) ) * taum(ji,jj) + at_i(ji,jj) * zrhoco * zmodt
            tmod_io(ji,jj) = zrhoco * SQRT( zmodt )          ! rhoco * |U_ice-U_oce| at T-point
         END_2D
         CALL lbc_lnk( 'iceupdate', taum, 'T', 1.0_wp, tmod_io, 'T', 1.0_wp )
         !
         utau_oce(:,:) = utau(:,:)                    !* save the air-ocean stresses at ice time-step
         vtau_oce(:,:) = vtau(:,:)
         !
      ENDIF
      !
      !                                      !==  every ocean time-step  ==!
      IF ( ln_drgice_imp ) THEN
         ! Save drag with right sign to update top drag in the ocean implicit friction
         rCdU_ice(:,:) = -r1_rho0 * tmod_io(:,:) * at_i(:,:) * tmask(:,:,1)
         zflagi = 0._wp
      ELSE
         zflagi = 1._wp
      ENDIF
      !
      DO_2D( 0, 0, 0, 0 )                             !* update the stress WITHOUT an ice-ocean rotation angle
         ! ice area at u and v-points
         zat_u  = ( at_i(ji,jj) * tmask(ji,jj,1) + at_i (ji+1,jj    ) * tmask(ji+1,jj  ,1) )  &
            &     / MAX( 1.0_wp , tmask(ji,jj,1) + tmask(ji+1,jj  ,1) )
         zat_v  = ( at_i(ji,jj) * tmask(ji,jj,1) + at_i (ji  ,jj+1  ) * tmask(ji  ,jj+1,1) )  &
            &     / MAX( 1.0_wp , tmask(ji,jj,1) + tmask(ji  ,jj+1,1) )
         !                                                   ! linearized quadratic drag formulation
         zutau_ice   = 0.5_wp * ( tmod_io(ji,jj) + tmod_io(ji+1,jj) ) * ( u_ice(ji,jj) - zflagi * pu_oce(ji,jj) )
         zvtau_ice   = 0.5_wp * ( tmod_io(ji,jj) + tmod_io(ji,jj+1) ) * ( v_ice(ji,jj) - zflagi * pv_oce(ji,jj) )
         !                                                   ! stresses at the ocean surface
         utau(ji,jj) = ( 1._wp - zat_u ) * utau_oce(ji,jj) + zat_u * zutau_ice
         vtau(ji,jj) = ( 1._wp - zat_v ) * vtau_oce(ji,jj) + zat_v * zvtau_ice
      END_2D
      CALL lbc_lnk( 'iceupdate', utau, 'U', -1.0_wp, vtau, 'V', -1.0_wp )   ! lateral boundary condition
      !
      IF( ln_timing )   CALL timing_stop('iceupdate')
      !
   END SUBROUTINE ice_update_tau


   SUBROUTINE ice_update_init
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_update_init  ***
      !!
      !! ** Purpose :   allocate ice-ocean stress fields and read restarts
      !!                containing the snow & ice mass
      !!
      !!-------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk               ! dummy loop indices
      REAL(wp) ::   zcoefu, zcoefv, zcoeff   ! local scalar
      !!-------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'ice_update_init: ice-ocean stress init'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~'
      !
      !                                      ! allocate ice_update array
      IF( ice_update_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'ice_update_init : unable to allocate standard arrays' )
      !
      CALL update_rst( 'READ' )  !* read or initialize all required files
      !
   END SUBROUTINE ice_update_init


   SUBROUTINE update_rst( cdrw, kt )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE rhg_evp_rst  ***
      !!
      !! ** Purpose :   Read or write RHG file in restart file
      !!
      !! ** Method  :   use of IOM library
      !!----------------------------------------------------------------------
      CHARACTER(len=*) , INTENT(in) ::   cdrw   ! 'READ'/'WRITE' flag
      INTEGER, OPTIONAL, INTENT(in) ::   kt     ! ice time-step
      !
      INTEGER  ::   iter   ! local integer
      INTEGER  ::   id1    ! local integer
      !!----------------------------------------------------------------------
      !
      IF( TRIM(cdrw) == 'READ' ) THEN        ! Read/initialize
         !                                   ! ---------------
         IF( ln_rstart ) THEN                   !* Read the restart file
            !
            id1 = iom_varid( numrir, 'snwice_mass' , ldstop = .FALSE. )
            !
            IF( id1 > 0 ) THEN                       ! fields exist
               CALL iom_get( numrir, jpdom_auto, 'snwice_mass'  , snwice_mass   )
               CALL iom_get( numrir, jpdom_auto, 'snwice_mass_b', snwice_mass_b )
            ELSE                                     ! start from rest
               IF(lwp) WRITE(numout,*) '   ==>>   previous run without snow-ice mass output then set it'
               snwice_mass  (:,:) = tmask(:,:,1) * ( rhos * vt_s(:,:) + rhoi * vt_i(:,:) & 
                                                &  + rhow * (vt_ip(:,:) + vt_il(:,:))  )
               snwice_mass_b(:,:) = snwice_mass(:,:)
            ENDIF
         ELSE                                   !* Start from rest
!JC: I think this is useless with what is now done in ice_istate
            IF(lwp) WRITE(numout,*) '   ==>>   start from rest: set the snow-ice mass'
            snwice_mass  (:,:) = tmask(:,:,1) * ( rhos * vt_s(:,:) + rhoi * vt_i(:,:) & 
                                             &  + rhow * (vt_ip(:,:) + vt_il(:,:))  )
            snwice_mass_b(:,:) = snwice_mass(:,:)
         ENDIF
         !
      ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN   ! Create restart file
         !                                   ! -------------------
         IF(lwp) WRITE(numout,*) '---- update-rst ----'
         iter = kt + nn_fsbc - 1             ! ice restarts are written at kt == nitrst - nn_fsbc + 1
         !
         CALL iom_rstput( iter, nitrst, numriw, 'snwice_mass'  , snwice_mass   )
         CALL iom_rstput( iter, nitrst, numriw, 'snwice_mass_b', snwice_mass_b )
         !
      ENDIF
      !
   END SUBROUTINE update_rst

#else
   !!----------------------------------------------------------------------
   !!   Default option         Dummy module           NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE iceupdate
