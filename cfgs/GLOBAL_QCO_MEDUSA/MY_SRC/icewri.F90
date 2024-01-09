MODULE icewri
   !!======================================================================
   !!                     ***  MODULE  icewri  ***
   !!   sea-ice : output ice variables
   !!======================================================================
   !! History :  4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_wri       : write of the diagnostics variables in ouput file
   !!   ice_wri_state : write for initial state or/and abandon
   !!----------------------------------------------------------------------
   USE dianam         ! build name of file (routine)
   USE phycst         ! physical constant
   USE dom_oce        ! domain: ocean
   USE sbc_oce        ! surf. boundary cond.: ocean
   USE sbc_ice        ! Surface boundary condition: ice fields
   USE ice            ! sea-ice: variables
   USE icevar         ! sea-ice: operations
   USE icealb , ONLY : rn_alb_oce
   !
   USE ioipsl         !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)
   USE lbclnk         ! lateral boundary conditions (or mpp links)
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC ice_wri        ! called by ice_stp
   PUBLIC ice_wri_state  ! called by dia_wri_state

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/ICE 4.0 , NEMO Consortium (2018)
   !! $Id: icewri.F90 15388 2021-10-17 11:33:47Z clem $
   !! Software governed by the CeCILL licence     (./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_wri( kt )
      !!-------------------------------------------------------------------
      !!  This routine ouputs some (most?) of the sea ice fields
      !!-------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! time-step
      !
      INTEGER  ::   ji, jj, jk, jl  ! dummy loop indices
      REAL(wp) ::   z2da, z2db, zrho1, zrho2
      REAL(wp) ::   zmiss_val       ! missing value retrieved from xios
      REAL(wp), DIMENSION(jpi,jpj)     ::   z2d                            ! 2D workspace
      REAL(wp), DIMENSION(jpi,jpj)     ::   zmsk00, zmsk05, zmsk15, zmsksn ! O%, 5% and 15% concentration mask and snow mask
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   zmsk00l, zmsksnl               ! cat masks
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   zfast, zalb, zmskalb      ! 2D workspace
      !
      ! Global ice diagnostics (SIMIP)
      REAL(wp) ::   zdiag_area_nh, zdiag_extt_nh, zdiag_volu_nh   ! area, extent, volume
      REAL(wp) ::   zdiag_area_sh, zdiag_extt_sh, zdiag_volu_sh
      !!-------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('icewri')

      ! get missing value from xml
      CALL iom_miss_val( 'icetemp', zmiss_val )

      ! brine volume
      CALL ice_var_bv

      ! tresholds for outputs
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
         zmsk00(ji,jj) = MAX( 0._wp , SIGN( 1._wp , at_i(ji,jj) - epsi06  ) ) ! 1 if ice    , 0 if no ice
         zmsk05(ji,jj) = MAX( 0._wp , SIGN( 1._wp , at_i(ji,jj) - 0.05_wp ) ) ! 1 if 5% ice , 0 if less
         zmsk15(ji,jj) = MAX( 0._wp , SIGN( 1._wp , at_i(ji,jj) - 0.15_wp ) ) ! 1 if 15% ice, 0 if less
         zmsksn(ji,jj) = MAX( 0._wp , SIGN( 1._wp , vt_s(ji,jj) - epsi06  ) ) ! 1 if snow   , 0 if no snow
      END_2D
      DO jl = 1, jpl
         DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
            zmsk00l(ji,jj,jl)  = MAX( 0._wp , SIGN( 1._wp , a_i(ji,jj,jl) - epsi06 ) )
            zmsksnl(ji,jj,jl)  = MAX( 0._wp , SIGN( 1._wp , v_s(ji,jj,jl) - epsi06 ) )
         END_2D
      END DO

      !-----------------
      ! Standard outputs
      !-----------------
      zrho1 = ( rho0 - rhoi ) * r1_rho0 ; zrho2 = rhos * r1_rho0
      ! masks
      CALL iom_put( 'icemask'  , zmsk00 )   ! ice mask 0%
      CALL iom_put( 'icemask05', zmsk05 )   ! ice mask 5%
      CALL iom_put( 'icemask15', zmsk15 )   ! ice mask 15%
      CALL iom_put( 'icepres'  , zmsk00 )   ! Ice presence (1 or 0)
      !
      ! general fields
      IF( iom_use('icemass' ) )   CALL iom_put( 'icemass', vt_i * rhoi * zmsk00 )                                           ! Ice mass per cell area
      IF( iom_use('snwmass' ) )   CALL iom_put( 'snwmass', vt_s * rhos * zmsksn )                                           ! Snow mass per cell area
      IF( iom_use('iceconc' ) )   CALL iom_put( 'iceconc', at_i        * zmsk00 )                                           ! ice concentration
      IF( iom_use('icevolu' ) )   CALL iom_put( 'icevolu', vt_i        * zmsk00 )                                           ! ice volume = mean ice thickness over the cell
      IF( iom_use('icethic' ) )   CALL iom_put( 'icethic', hm_i        * zmsk00 )                                           ! ice thickness
      IF( iom_use('snwthic' ) )   CALL iom_put( 'snwthic', hm_s        * zmsk00 )                                           ! snw thickness
      IF( iom_use('icebrv'  ) )   CALL iom_put( 'icebrv' , bvm_i* 100. * zmsk00 )                                           ! brine volume
      IF( iom_use('iceage'  ) )   CALL iom_put( 'iceage' , om_i / rday * zmsk15 + zmiss_val * ( 1._wp - zmsk15 ) )          ! ice age
      IF( iom_use('icehnew' ) )   CALL iom_put( 'icehnew', ht_i_new             )                                           ! new ice thickness formed in the leads
      IF( iom_use('snwvolu' ) )   CALL iom_put( 'snwvolu', vt_s        * zmsksn )                                           ! snow volume
      IF( iom_use('icefrb'  ) ) THEN                                                                                        ! Ice freeboard
         z2d(:,:) = ( zrho1 * hm_i(:,:) - zrho2 * hm_s(:,:) )
         WHERE( z2d < 0._wp )   z2d = 0._wp
                                  CALL iom_put( 'icefrb' , z2d * zmsk00         )
      ENDIF
      ! melt ponds
      IF( iom_use('iceapnd' ) )   CALL iom_put( 'iceapnd', at_ip  * zmsk00      )                                           ! melt pond total fraction
      IF( iom_use('icehpnd' ) )   CALL iom_put( 'icehpnd', hm_ip  * zmsk00      )                                           ! melt pond depth
      IF( iom_use('icevpnd' ) )   CALL iom_put( 'icevpnd', vt_ip  * zmsk00      )                                           ! melt pond total volume per unit area
      IF( iom_use('icehlid' ) )   CALL iom_put( 'icehlid', hm_il  * zmsk00      )                                           ! melt pond lid depth
      IF( iom_use('icevlid' ) )   CALL iom_put( 'icevlid', vt_il  * zmsk00      )                                           ! melt pond lid total volume per unit area
      ! salt
      IF( iom_use('icesalt' ) )   CALL iom_put( 'icesalt', sm_i                 * zmsk00 + zmiss_val * ( 1._wp - zmsk00 ) ) ! mean ice salinity
      IF( iom_use('icesalm' ) )   CALL iom_put( 'icesalm', st_i * rhoi * 1.0e-3 * zmsk00 )                                  ! Mass of salt in sea ice per cell area
      IF( iom_use('iceepnd' ) )   CALL iom_put( 'iceepnd', SUM( a_ip_eff * a_i, dim=3 ) * zmsk00  )                         ! melt pond total effective fraction per cell area
      ! heat
      IF( iom_use('icetemp' ) )   CALL iom_put( 'icetemp', ( tm_i  - rt0 ) * zmsk00 + zmiss_val * ( 1._wp - zmsk00 ) )      ! ice mean temperature
      IF( iom_use('snwtemp' ) )   CALL iom_put( 'snwtemp', ( tm_s  - rt0 ) * zmsksn + zmiss_val * ( 1._wp - zmsksn ) )      ! snw mean temperature
      IF( iom_use('icettop' ) )   CALL iom_put( 'icettop', ( tm_su - rt0 ) * zmsk00 + zmiss_val * ( 1._wp - zmsk00 ) )      ! temperature at the ice surface
      IF( iom_use('icetbot' ) )   CALL iom_put( 'icetbot', ( t_bo  - rt0 ) * zmsk00 + zmiss_val * ( 1._wp - zmsk00 ) )      ! temperature at the ice bottom
      IF( iom_use('icetsni' ) )   CALL iom_put( 'icetsni', ( tm_si - rt0 ) * zmsk00 + zmiss_val * ( 1._wp - zmsk00 ) )      ! temperature at the snow-ice interface
      IF( iom_use('icehc'   ) )   CALL iom_put( 'icehc'  ,  -et_i          * zmsk00 )                                       ! ice heat content
      IF( iom_use('snwhc'   ) )   CALL iom_put( 'snwhc'  ,  -et_s          * zmsksn )                                       ! snow heat content
      ! momentum
      IF( iom_use('uice'    ) )   CALL iom_put( 'uice'   , u_ice    )                                                       ! ice velocity u
      IF( iom_use('vice'    ) )   CALL iom_put( 'vice'   , v_ice    )                                                       ! ice velocity v
      !
      IF( iom_use('icevel') .OR. iom_use('fasticepres') ) THEN                                                              ! module of ice velocity & fast ice
         ALLOCATE( zfast(jpi,jpj) )
         DO_2D( 0, 0, 0, 0 )
            z2da  = u_ice(ji,jj) + u_ice(ji-1,jj)
            z2db  = v_ice(ji,jj) + v_ice(ji,jj-1)
            z2d(ji,jj) = 0.5_wp * SQRT( z2da * z2da + z2db * z2db )
         END_2D
         CALL lbc_lnk( 'icewri', z2d, 'T', 1.0_wp )
         CALL iom_put( 'icevel', z2d )

         WHERE( z2d(:,:) < 5.e-04_wp .AND. zmsk15(:,:) == 1._wp ) ; zfast(:,:) = 1._wp                                      ! record presence of fast ice
         ELSEWHERE                                                ; zfast(:,:) = 0._wp
         END WHERE
         CALL iom_put( 'fasticepres', zfast )
         DEALLOCATE( zfast )
      ENDIF
      !
      IF( iom_use('icealb') .OR. iom_use('albedo') ) THEN                                                                   ! ice albedo and surface albedo
         ALLOCATE( zalb(jpi,jpj), zmskalb(jpi,jpj) )
         ! ice albedo
         WHERE( at_i_b < 1.e-03 )
            zmskalb(:,:) = 0._wp
            zalb   (:,:) = rn_alb_oce
         ELSEWHERE
            zmskalb(:,:) = 1._wp
            zalb   (:,:) = SUM( alb_ice * a_i_b, dim=3 ) / at_i_b
         END WHERE
         CALL iom_put( 'icealb' , zalb * zmskalb + zmiss_val * ( 1._wp - zmskalb ) )
         ! ice+ocean albedo
         zalb(:,:) = SUM( alb_ice * a_i_b, dim=3 ) + rn_alb_oce * ( 1._wp - at_i_b )
         CALL iom_put( 'albedo' , zalb )
         DEALLOCATE( zalb, zmskalb )
      ENDIF
      !
      ! --- category-dependent fields --- !
      IF( iom_use('icemask_cat' ) )   CALL iom_put( 'icemask_cat' ,                  zmsk00l                                   ) ! ice mask 0%
      IF( iom_use('iceconc_cat' ) )   CALL iom_put( 'iceconc_cat' , a_i            * zmsk00l                                   ) ! area for categories
      IF( iom_use('icethic_cat' ) )   CALL iom_put( 'icethic_cat' , h_i            * zmsk00l + zmiss_val * ( 1._wp - zmsk00l ) ) ! thickness for categories
      IF( iom_use('snwthic_cat' ) )   CALL iom_put( 'snwthic_cat' , h_s            * zmsksnl + zmiss_val * ( 1._wp - zmsksnl ) ) ! snow depth for categories
      IF( iom_use('icesalt_cat' ) )   CALL iom_put( 'icesalt_cat' , s_i            * zmsk00l + zmiss_val * ( 1._wp - zmsk00l ) ) ! salinity for categories
      IF( iom_use('iceage_cat'  ) )   CALL iom_put( 'iceage_cat'  , o_i / rday     * zmsk00l + zmiss_val * ( 1._wp - zmsk00l ) ) ! ice age
      IF( iom_use('icetemp_cat' ) )   CALL iom_put( 'icetemp_cat' , ( SUM( t_i, dim=3 ) * r1_nlay_i - rt0 ) &
         &                                                                         * zmsk00l + zmiss_val * ( 1._wp - zmsk00l ) ) ! ice temperature
      IF( iom_use('snwtemp_cat' ) )   CALL iom_put( 'snwtemp_cat' , ( SUM( t_s, dim=3 ) * r1_nlay_s - rt0 ) &
         &                                                                         * zmsksnl + zmiss_val * ( 1._wp - zmsksnl ) ) ! snow temperature
      IF( iom_use('icettop_cat' ) )   CALL iom_put( 'icettop_cat' , ( t_su - rt0 ) * zmsk00l + zmiss_val * ( 1._wp - zmsk00l ) ) ! surface temperature
      IF( iom_use('icebrv_cat'  ) )   CALL iom_put( 'icebrv_cat'  ,   bv_i * 100.  * zmsk00l + zmiss_val * ( 1._wp - zmsk00l ) ) ! brine volume
      IF( iom_use('iceapnd_cat' ) )   CALL iom_put( 'iceapnd_cat' ,   a_ip         * zmsk00l                                   ) ! melt pond frac for categories
      IF( iom_use('icevpnd_cat' ) )   CALL iom_put( 'icevpnd_cat' ,   v_ip         * zmsk00l                                   ) ! melt pond volume for categories
      IF( iom_use('icehpnd_cat' ) )   CALL iom_put( 'icehpnd_cat' ,   h_ip         * zmsk00l + zmiss_val * ( 1._wp - zmsk00l ) ) ! melt pond thickness for categories
      IF( iom_use('icehlid_cat' ) )   CALL iom_put( 'icehlid_cat' ,   h_il         * zmsk00l + zmiss_val * ( 1._wp - zmsk00l ) ) ! melt pond lid thickness for categories
      IF( iom_use('iceafpnd_cat') )   CALL iom_put( 'iceafpnd_cat',   a_ip_frac    * zmsk00l                                   ) ! melt pond frac per ice area for categories
      IF( iom_use('iceaepnd_cat') )   CALL iom_put( 'iceaepnd_cat',   a_ip_eff     * zmsk00l                                   ) ! melt pond effective frac for categories
      IF( iom_use('icealb_cat'  ) )   CALL iom_put( 'icealb_cat'  ,   alb_ice      * zmsk00l + zmiss_val * ( 1._wp - zmsk00l ) ) ! ice albedo for categories

      !------------------
      ! Add-ons for SIMIP
      !------------------
      ! trends
      IF( iom_use('dmithd') )   CALL iom_put( 'dmithd', - wfx_bog - wfx_bom - wfx_sum - wfx_sni - wfx_opw - wfx_lam - wfx_res ) ! Sea-ice mass change from thermodynamics
      IF( iom_use('dmidyn') )   CALL iom_put( 'dmidyn', - wfx_dyn + rhoi * diag_trp_vi                                        ) ! Sea-ice mass change from dynamics(kg/m2/s)
      IF( iom_use('dmiopw') )   CALL iom_put( 'dmiopw', - wfx_opw                                                             ) ! Sea-ice mass change through growth in open water
      IF( iom_use('dmibog') )   CALL iom_put( 'dmibog', - wfx_bog                                                             ) ! Sea-ice mass change through basal growth
      IF( iom_use('dmisni') )   CALL iom_put( 'dmisni', - wfx_sni                                                             ) ! Sea-ice mass change through snow-to-ice conversion
      IF( iom_use('dmisum') )   CALL iom_put( 'dmisum', - wfx_sum                                                             ) ! Sea-ice mass change through surface melting
      IF( iom_use('dmibom') )   CALL iom_put( 'dmibom', - wfx_bom                                                             ) ! Sea-ice mass change through bottom melting
      IF( iom_use('dmilam') )   CALL iom_put( 'dmilam', - wfx_lam                                                             ) ! Sea-ice mass change through lateral melting
      IF( iom_use('dmtsub') )   CALL iom_put( 'dmtsub', - wfx_sub                                                             ) ! Sea-ice mass change through evaporation and sublimation
      IF( iom_use('dmssub') )   CALL iom_put( 'dmssub', - wfx_snw_sub                                                         ) ! Snow mass change through sublimation
      IF( iom_use('dmisub') )   CALL iom_put( 'dmisub', - wfx_ice_sub                                                         ) ! Sea-ice mass change through sublimation
      IF( iom_use('dmsspr') )   CALL iom_put( 'dmsspr', - wfx_spr                                                             ) ! Snow mass change through snow fall
      IF( iom_use('dmsssi') )   CALL iom_put( 'dmsssi',   wfx_sni*rhos*r1_rhoi                                                ) ! Snow mass change through snow-to-ice conversion
      IF( iom_use('dmsmel') )   CALL iom_put( 'dmsmel', - wfx_snw_sum                                                         ) ! Snow mass change through melt
      IF( iom_use('dmsdyn') )   CALL iom_put( 'dmsdyn', - wfx_snw_dyn + rhos * diag_trp_vs                                    ) ! Snow mass change through dynamics(kg/m2/s)

      ! Global ice diagnostics
      IF(  iom_use('NH_icearea') .OR. iom_use('NH_icevolu') .OR. iom_use('NH_iceextt') .OR. &
         & iom_use('SH_icearea') .OR. iom_use('SH_icevolu') .OR. iom_use('SH_iceextt') ) THEN
         !
         WHERE( ff_t(:,:) > 0._wp )   ;   z2d(:,:) = 1._wp
         ELSEWHERE                    ;   z2d(:,:) = 0.
         END WHERE
         !
         IF( iom_use('NH_icearea') )   zdiag_area_nh = glob_sum( 'icewri', at_i *           z2d   * e1e2t * 1.e-12 )
         IF( iom_use('NH_icevolu') )   zdiag_volu_nh = glob_sum( 'icewri', vt_i *           z2d   * e1e2t * 1.e-12 )
         IF( iom_use('NH_iceextt') )   zdiag_extt_nh = glob_sum( 'icewri',                  z2d   * e1e2t * 1.e-12 * zmsk15 )
         !
         IF( iom_use('SH_icearea') )   zdiag_area_sh = glob_sum( 'icewri', at_i * ( 1._wp - z2d ) * e1e2t * 1.e-12 )
         IF( iom_use('SH_icevolu') )   zdiag_volu_sh = glob_sum( 'icewri', vt_i * ( 1._wp - z2d ) * e1e2t * 1.e-12 )
         IF( iom_use('SH_iceextt') )   zdiag_extt_sh = glob_sum( 'icewri',        ( 1._wp - z2d ) * e1e2t * 1.e-12 * zmsk15 )
         !
         CALL iom_put( 'NH_icearea' , zdiag_area_nh )
         CALL iom_put( 'NH_icevolu' , zdiag_volu_nh )
         CALL iom_put( 'NH_iceextt' , zdiag_extt_nh )
         CALL iom_put( 'SH_icearea' , zdiag_area_sh )
         CALL iom_put( 'SH_icevolu' , zdiag_volu_sh )
         CALL iom_put( 'SH_iceextt' , zdiag_extt_sh )
         !
      ENDIF
      !
!!CR      !     !  Create an output files (output.lim.abort.nc) if S < 0 or u > 20 m/s
!!CR      !     IF( kindic < 0 )   CALL ice_wri_state( 'output.abort' )
!!CR      !     not yet implemented
!!gm  idem for the ocean...  Ask Seb how to get rid of ioipsl....
      !
      IF( ln_timing )  CALL timing_stop('icewri')
      !
   END SUBROUTINE ice_wri


   SUBROUTINE ice_wri_state( kid )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE ice_wri_state  ***
      !!
      !! ** Purpose :   create a NetCDF file named cdfile_name which contains
      !!      the instantaneous ice state and forcing fields for ice model
      !!        Used to find errors in the initial state or save the last
      !!      ocean state in case of abnormal end of a simulation
      !!
      !! History :   4.0  !  2013-06  (C. Rousset)
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kid
      !!----------------------------------------------------------------------
      !
      !! The file is open in dia_wri_state (ocean routine)

      CALL iom_rstput( 0, 0, kid, 'sithic', hm_i         )   ! Ice thickness
      CALL iom_rstput( 0, 0, kid, 'siconc', at_i         )   ! Ice concentration
      CALL iom_rstput( 0, 0, kid, 'sitemp', tm_i - rt0   )   ! Ice temperature
      CALL iom_rstput( 0, 0, kid, 'sivelu', u_ice        )   ! i-Ice speed
      CALL iom_rstput( 0, 0, kid, 'sivelv', v_ice        )   ! j-Ice speed
      CALL iom_rstput( 0, 0, kid, 'sistru', utau_ice     )   ! i-Wind stress over ice
      CALL iom_rstput( 0, 0, kid, 'sistrv', vtau_ice     )   ! i-Wind stress over ice
      CALL iom_rstput( 0, 0, kid, 'sisflx', qsr          )   ! Solar flx over ocean
      CALL iom_rstput( 0, 0, kid, 'sinflx', qns          )   ! NonSolar flx over ocean
      CALL iom_rstput( 0, 0, kid, 'snwpre', sprecip      )   ! Snow precipitation
      CALL iom_rstput( 0, 0, kid, 'sisali', sm_i         )   ! Ice salinity
      CALL iom_rstput( 0, 0, kid, 'sivolu', vt_i         )   ! Ice volume
      CALL iom_rstput( 0, 0, kid, 'sidive', divu_i*1.0e8 )   ! Ice divergence
      CALL iom_rstput( 0, 0, kid, 'si_amp', at_ip        )   ! Melt pond fraction
      CALL iom_rstput( 0, 0, kid, 'si_vmp', vt_ip        )   ! Melt pond volume
      CALL iom_rstput( 0, 0, kid, 'sithicat', h_i        )   ! Ice thickness
      CALL iom_rstput( 0, 0, kid, 'siconcat', a_i        )   ! Ice concentration
      CALL iom_rstput( 0, 0, kid, 'sisalcat', s_i        )   ! Ice salinity
      CALL iom_rstput( 0, 0, kid, 'snthicat', h_s        )   ! Snw thickness

    END SUBROUTINE ice_wri_state

#else
   !!----------------------------------------------------------------------
   !!   Default option :         Empty module         NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE icewri
