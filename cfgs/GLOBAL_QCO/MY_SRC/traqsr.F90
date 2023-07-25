MODULE traqsr
   !!======================================================================
   !!                       ***  MODULE  traqsr  ***
   !! Ocean physics:   solar radiation penetration in the top ocean levels
   !!======================================================================
   !! History :  OPA  !  1990-10  (B. Blanke)  Original code
   !!            7.0  !  1991-11  (G. Madec)
   !!                 !  1996-01  (G. Madec)  s-coordinates
   !!   NEMO     1.0  !  2002-06  (G. Madec)  F90: Free form and module
   !!             -   !  2005-11  (G. Madec) zco, zps, sco coordinate
   !!            3.2  !  2009-04  (G. Madec & NEMO team)
   !!            3.6  !  2012-05  (C. Rousset) store attenuation coef for use in ice model
   !!            3.6  !  2015-12  (O. Aumont, J. Jouanno, C. Ethe) use vertical profile of chlorophyll
   !!            3.7  !  2015-11  (G. Madec, A. Coward)  remove optimisation for fix volume
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_qsr       : temperature trend due to the penetration of solar radiation
   !!   tra_qsr_init  : initialization of the qsr penetration
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and active tracers
   USE phycst         ! physical constants
   USE dom_oce        ! ocean space and time domain
   USE domtile
   USE sbc_oce        ! surface boundary condition: ocean
   USE trc_oce        ! share SMS/Ocean variables
   USE trd_oce        ! trends: ocean variables
   USE trdtra         ! trends manager: tracers
   !
   USE in_out_manager ! I/O manager
   USE prtctl         ! Print control
   USE iom            ! I/O library
   USE fldread        ! read input fields
   USE restart        ! ocean restart
   USE lib_mpp        ! MPP library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_qsr       ! routine called by step.F90 (ln_traqsr=T)
   PUBLIC   tra_qsr_init  ! routine called by nemogcm.F90

   !                                 !!* Namelist namtra_qsr: penetrative solar radiation
   LOGICAL , PUBLIC ::   ln_traqsr    !: light absorption (qsr) flag
   LOGICAL , PUBLIC ::   ln_qsr_rgb   !: Red-Green-Blue light absorption flag
   LOGICAL , PUBLIC ::   ln_qsr_2bd   !: 2 band         light absorption flag
   LOGICAL , PUBLIC ::   ln_qsr_bio   !: bio-model      light absorption flag
   INTEGER , PUBLIC ::   nn_chldta    !: use Chlorophyll data (=1) or not (=0)
   REAL(wp), PUBLIC ::   rn_chl_conc  !: Chlorophyll concentration (for nn_chldta=0)
   REAL(wp), PUBLIC ::   rn_abs       !: fraction absorbed in the very near surface (RGB & 2 bands)
   REAL(wp), PUBLIC ::   rn_si0       !: very near surface depth of extinction      (RGB & 2 bands)
   REAL(wp), PUBLIC ::   rn_si1       !: deepest depth of extinction (water type I)       (2 bands)
   !
   INTEGER , PUBLIC ::   nksr         !: levels below which the light cannot penetrate (depth larger than 391 m)

   INTEGER, PARAMETER ::   np_RGB  = 1   ! R-G-B     light penetration with constant Chlorophyll
   INTEGER, PARAMETER ::   np_RGBc = 2   ! R-G-B     light penetration with Chlorophyll data
   INTEGER, PARAMETER ::   np_2BD  = 3   ! 2 bands   light penetration
   INTEGER, PARAMETER ::   np_BIO  = 4   ! bio-model light penetration
   !
   INTEGER  ::   nqsr    ! user choice of the type of light penetration
   REAL(wp) ::   xsi0r   ! inverse of rn_si0
   REAL(wp) ::   xsi1r   ! inverse of rn_si1
   !
   REAL(wp) , PUBLIC, DIMENSION(3,61)   ::   rkrgb    ! tabulated attenuation coefficients for RGB absorption
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_chl   ! structure of input Chl (file informations, fields read)

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: traqsr.F90 14834 2021-05-11 09:24:44Z hadcv $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_qsr( kt, Kmm, pts, Krhs )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_qsr  ***
      !!
      !! ** Purpose :   Compute the temperature trend due to the solar radiation
      !!              penetration and add it to the general temperature trend.
      !!
      !! ** Method  : The profile of the solar radiation within the ocean is defined
      !!      through 2 wavebands (rn_si0,rn_si1) or 3 wavebands (RGB) and a ratio rn_abs
      !!      Considering the 2 wavebands case:
      !!         I(k) = Qsr*( rn_abs*EXP(z(k)/rn_si0) + (1.-rn_abs)*EXP(z(k)/rn_si1) )
      !!         The temperature trend associated with the solar radiation penetration
      !!         is given by : zta = 1/e3t dk[ I ] / (rho0*Cp)
      !!         At the bottom, boudary condition for the radiation is no flux :
      !!      all heat which has not been absorbed in the above levels is put
      !!      in the last ocean level.
      !!         The computation is only done down to the level where
      !!      I(k) < 1.e-15 W/m2 (i.e. over the top nksr levels) .
      !!
      !! ** Action  : - update ta with the penetrative solar radiation trend
      !!              - send  trend for further diagnostics (l_trdtra=T)
      !!
      !! Reference  : Jerlov, N. G., 1968 Optical Oceanography, Elsevier, 194pp.
      !!              Lengaigne et al. 2007, Clim. Dyn., V28, 5, 503-516.
      !!              Morel, A. et Berthon, JF, 1989, Limnol Oceanogr 34(8), 1545-1562
      !!----------------------------------------------------------------------
      INTEGER,                                   INTENT(in   ) :: kt            ! ocean time-step
      INTEGER,                                   INTENT(in   ) :: Kmm, Krhs     ! time level indices
      REAL(dp), DIMENSION(jpi,jpj,jpk,jpts,jpt), INTENT(inout) :: pts           ! active tracers and RHS of tracer equation
      !
      INTEGER  ::   ji, jj, jk               ! dummy loop indices
      INTEGER  ::   irgb                     ! local integers
      REAL(wp) ::   zchl, zcoef, z1_2        ! local scalars
      REAL(wp) ::   zc0 , zc1 , zc2 , zc3    !    -         -
      REAL(wp) ::   zz0 , zz1 , ze3t, zlui   !    -         -
      REAL(wp) ::   zCb, zCmax, zpsi, zpsimax, zrdpsi, zCze
      REAL(wp) ::   zlogc, zlogze, zlogCtot, zlogCze
      REAL(wp), ALLOCATABLE, DIMENSION(:,:)   :: ze0, ze1, ze2, ze3
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: ztrdt, zetot, ztmp3d
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('tra_qsr')
      !
      IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                       ! Do only on the first tile
         IF( kt == nit000 ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'tra_qsr : penetration of the surface solar radiation'
            IF(lwp) WRITE(numout,*) '~~~~~~~'
         ENDIF
      ENDIF
      !
      IF( l_trdtra ) THEN      ! trends diagnostic: save the input temperature trend
         ALLOCATE( ztrdt(jpi,jpj,jpk) )
         ztrdt(:,:,:) = pts(:,:,:,jp_tem,Krhs)
      ENDIF
      !
      !                         !-----------------------------------!
      !                         !  before qsr induced heat content  !
      !                         !-----------------------------------!
      IF( kt == nit000 ) THEN          !==  1st time step  ==!
         IF( ln_rstart .AND. .NOT.l_1st_euler ) THEN    ! read in restart
            z1_2 = 0.5_wp
            IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                        ! Do only on the first tile
               IF(lwp) WRITE(numout,*) '          nit000-1 qsr tracer content forcing field read in the restart file'
               CALL iom_get( numror, jpdom_auto, 'qsr_hc_b', qsr_hc_b )   ! before heat content trend due to Qsr flux
            ENDIF
         ELSE                                           ! No restart or Euler forward at 1st time step
            z1_2 = 1._wp
            DO_3D_OVR( nn_hls, nn_hls, nn_hls, nn_hls, 1, jpk )
               qsr_hc_b(ji,jj,jk) = 0._wp
            END_3D
         ENDIF
      ELSE                             !==  Swap of qsr heat content  ==!
         z1_2 = 0.5_wp
         DO_3D_OVR( nn_hls, nn_hls, nn_hls, nn_hls, 1, jpk )
            qsr_hc_b(ji,jj,jk) = qsr_hc(ji,jj,jk)
         END_3D
      ENDIF
      !
      !                         !--------------------------------!
      SELECT CASE( nqsr )       !  now qsr induced heat content  !
      !                         !--------------------------------!
      !
      CASE( np_BIO )                   !==  bio-model fluxes  ==!
         !
         DO_3D_OVR( nn_hls, nn_hls, nn_hls, nn_hls, 1, nksr )
            qsr_hc(ji,jj,jk) = r1_rho0_rcp * ( etot3(ji,jj,jk) - etot3(ji,jj,jk+1) )
         END_3D
         !
      CASE( np_RGB , np_RGBc )         !==  R-G-B fluxes  ==!
         !
         ALLOCATE( ze0 (A2D(nn_hls))           , ze1 (A2D(nn_hls)) ,   &
            &      ze2 (A2D(nn_hls))           , ze3 (A2D(nn_hls)) ,   &
            &      ztmp3d(A2D(nn_hls),nksr + 1)                     )
         !
         IF( nqsr == np_RGBc ) THEN          !*  Variable Chlorophyll
            IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                                         ! Do only for the full domain
               IF( ln_tile ) CALL dom_tile_stop( ldhold=.TRUE. )             ! Use full domain
               CALL fld_read( kt, 1, sf_chl )         ! Read Chl data and provides it at the current time step
               IF( ln_tile ) CALL dom_tile_start( ldhold=.TRUE. )            ! Revert to tile domain
            ENDIF
            !
            ! Separation in R-G-B depending on the surface Chl
            ! perform and store as many of the 2D calculations as possible
            ! before the 3D loop (use the temporary 2D arrays to replace the
            ! most expensive calculations)
            !
            DO_2D_OVR( nn_hls, nn_hls, nn_hls, nn_hls )
                       ! zlogc = log(zchl)
               zlogc = LOG ( MIN( 10. , MAX( 0.03, sf_chl(1)%fnow(ji,jj,1) ) ) )
                       ! zc1 : log(zCze)  = log (1.12  * zchl**0.803)
               zc1   = 0.113328685307 + 0.803 * zlogc
                       ! zc2 : log(zCtot) = log(40.6  * zchl**0.459)
               zc2   = 3.703768066608 + 0.459 * zlogc
                       ! zc3 : log(zze)   = log(568.2 * zCtot**(-0.746))
               zc3   = 6.34247346942  - 0.746 * zc2
                       ! IF( log(zze) > log(102.) ) log(zze) = log(200.0 * zCtot**(-0.293))
               IF( zc3 > 4.62497281328 ) zc3 = 5.298317366548 - 0.293 * zc2
               !
               ze0(ji,jj) = zlogc                                                 ! ze0 = log(zchl)
               ze1(ji,jj) = EXP( zc1 )                                            ! ze1 = zCze
               ze2(ji,jj) = 1._wp / ( 0.710 + zlogc * ( 0.159 + zlogc * 0.021 ) ) ! ze2 = 1/zdelpsi
               ze3(ji,jj) = EXP( - zc3 )                                          ! ze3 = 1/zze
            END_2D

!
            DO_3D_OVR( nn_hls, nn_hls, nn_hls, nn_hls, 1, nksr + 1 )
               ! zchl    = ALOG( ze0(ji,jj) )
               zlogc = ze0(ji,jj)
               !
               zCb       = 0.768 + zlogc * ( 0.087 - zlogc * ( 0.179 + zlogc * 0.025 ) )
               zCmax     = 0.299 - zlogc * ( 0.289 - zlogc * 0.579 )
               zpsimax   = 0.6   - zlogc * ( 0.640 - zlogc * ( 0.021 + zlogc * 0.115 ) )
               ! zdelpsi = 0.710 + zlogc * ( 0.159 + zlogc * 0.021 )
               !
               zCze   = ze1(ji,jj)
               zrdpsi = ze2(ji,jj)                                                 ! 1/zdelpsi
               zpsi   = ze3(ji,jj) * gdepw(ji,jj,jk,Kmm)                           ! gdepw/zze
               !
               ! NB. make sure zchl value is such that: zchl = MIN( 10. , MAX( 0.03, zchl ) )
               zchl = MIN( 10. , MAX( 0.03, zCze * ( zCb + zCmax * EXP( -( (zpsi - zpsimax) * zrdpsi )**2 ) ) ) )
               ! Convert chlorophyll value to attenuation coefficient look-up table index
               ztmp3d(ji,jj,jk) = 41 + 20.*LOG10(zchl) + 1.e-15
            END_3D
         ELSE                                !* constant chlorophyll
            zchl = rn_chl_conc
            ! NB. make sure constant value is such that:
            zchl = MIN( 10. , MAX( 0.03, zchl ) )
            ! Convert chlorophyll value to attenuation coefficient look-up table index
            zlui = 41 + 20.*LOG10(zchl) + 1.e-15
            DO jk = 1, nksr + 1
               ztmp3d(:,:,jk) = zlui
            END DO
         ENDIF
         !
         zcoef  = ( 1. - rn_abs ) / 3._wp    !* surface equi-partition in R-G-B
         DO_2D_OVR( nn_hls, nn_hls, nn_hls, nn_hls )
            ze0(ji,jj) = rn_abs * qsr(ji,jj)
            ze1(ji,jj) = zcoef  * qsr(ji,jj)
            ze2(ji,jj) = zcoef  * qsr(ji,jj)
            ze3(ji,jj) = zcoef  * qsr(ji,jj)
            ! store the surface SW radiation; re-use the surface ztmp3d array
            ! since the surface attenuation coefficient is not used
            ztmp3d(ji,jj,1) =       qsr(ji,jj)
         END_2D
         !
         !                                    !* interior equi-partition in R-G-B depending on vertical profile of Chl
         DO_3D_OVR( nn_hls, nn_hls, nn_hls, nn_hls, 2, nksr + 1 )
            ze3t = e3t(ji,jj,jk-1,Kmm)
            irgb = NINT( ztmp3d(ji,jj,jk) )
            zc0 = ze0(ji,jj) * EXP( - ze3t * xsi0r )
            zc1 = ze1(ji,jj) * EXP( - ze3t * rkrgb(1,irgb) )
            zc2 = ze2(ji,jj) * EXP( - ze3t * rkrgb(2,irgb) )
            zc3 = ze3(ji,jj) * EXP( - ze3t * rkrgb(3,irgb) )
            ze0(ji,jj) = zc0
            ze1(ji,jj) = zc1
            ze2(ji,jj) = zc2
            ze3(ji,jj) = zc3
            ztmp3d(ji,jj,jk) = ( zc0 + zc1 + zc2 + zc3 ) * wmask(ji,jj,jk)
         END_3D
         !
         DO_3D_OVR( nn_hls, nn_hls, nn_hls, nn_hls, 1, nksr )          !* now qsr induced heat content
            qsr_hc(ji,jj,jk) = r1_rho0_rcp * ( ztmp3d(ji,jj,jk) - ztmp3d(ji,jj,jk+1) )
         END_3D
         !
         DEALLOCATE( ze0 , ze1 , ze2 , ze3 , ztmp3d )
         !
      CASE( np_2BD  )            !==  2-bands fluxes  ==!
         !
         zz0 =        rn_abs   * r1_rho0_rcp      ! surface equi-partition in 2-bands
         zz1 = ( 1. - rn_abs ) * r1_rho0_rcp
         DO_3D_OVR( nn_hls, nn_hls, nn_hls, nn_hls, 1, nksr )          !* now qsr induced heat content
            zc0 = zz0 * EXP( -gdepw(ji,jj,jk  ,Kmm)*xsi0r ) + zz1 * EXP( -gdepw(ji,jj,jk  ,Kmm)*xsi1r )
            zc1 = zz0 * EXP( -gdepw(ji,jj,jk+1,Kmm)*xsi0r ) + zz1 * EXP( -gdepw(ji,jj,jk+1,Kmm)*xsi1r )
            qsr_hc(ji,jj,jk) = qsr(ji,jj) * ( zc0 * wmask(ji,jj,jk) - zc1 * wmask(ji,jj,jk+1) )
         END_3D
         !
      END SELECT
      !
      !                          !-----------------------------!
      !                          !  update to the temp. trend  !
      !                          !-----------------------------!
      DO_3D( 0, 0, 0, 0, 1, nksr )
         pts(ji,jj,jk,jp_tem,Krhs) = pts(ji,jj,jk,jp_tem,Krhs)   &
            &                      + z1_2 * ( qsr_hc_b(ji,jj,jk) + qsr_hc(ji,jj,jk) )   &
            &                             / e3t(ji,jj,jk,Kmm)
      END_3D
      !
      ! sea-ice: store the 1st ocean level attenuation coefficient
      DO_2D_OVR( nn_hls, nn_hls, nn_hls, nn_hls )
         zz0 = r1_rho0_rcp * qsr(ji,jj)   ! test zz0 and not qsr for rounding errors in single precision
         IF( zz0 /= 0._wp ) THEN   ;   fraqsr_1lev(ji,jj) = qsr_hc(ji,jj,1) / zz0
         ELSE                      ;   fraqsr_1lev(ji,jj) = 1._wp
         ENDIF
      END_2D
      !
      IF( iom_use('qsr3d') ) THEN      ! output the shortwave Radiation distribution
         ALLOCATE( zetot(A2D(0),jpk) )
         zetot(:,:,nksr+1:jpk) = 0._wp     ! below ~400m set to zero
         DO_3DS(0, 0, 0, 0, nksr, 1, -1)
            zetot(ji,jj,jk) = zetot(ji,jj,jk+1) + qsr_hc(ji,jj,jk) * rho0_rcp
         END_3D
         CALL iom_put( 'qsr3d', zetot )   ! 3D distribution of shortwave Radiation
         DEALLOCATE( zetot )
      ENDIF
      !
      IF( .NOT. l_istiled .OR. ntile == nijtile )  THEN                ! Do only on the last tile
         IF( lrst_oce ) THEN     ! write in the ocean restart file
            CALL iom_rstput( kt, nitrst, numrow, 'qsr_hc_b'   , qsr_hc      )
            CALL iom_rstput( kt, nitrst, numrow, 'fraqsr_1lev', fraqsr_1lev )
         ENDIF
      ENDIF
      !
      IF( l_trdtra ) THEN     ! qsr tracers trends saved for diagnostics
         ztrdt(:,:,:) = pts(:,:,:,jp_tem,Krhs) - ztrdt(:,:,:)
         CALL trd_tra( kt, Kmm, Krhs, 'TRA', jp_tem, jptra_qsr, ztrdt )
         DEALLOCATE( ztrdt )
      ENDIF
      !                       ! print mean trends (used for debugging)
      IF(sn_cfctl%l_prtctl)   CALL prt_ctl( tab3d_1=pts(:,:,:,jp_tem,Krhs), clinfo1=' qsr  - Ta: ', mask1=tmask, clinfo3='tra-ta' )
      
      IF( ln_timing )   CALL timing_stop('tra_qsr')
      !
   END SUBROUTINE tra_qsr


   SUBROUTINE tra_qsr_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_qsr_init  ***
      !!
      !! ** Purpose :   Initialization for the penetrative solar radiation
      !!
      !! ** Method  :   The profile of solar radiation within the ocean is set
      !!      from two length scale of penetration (rn_si0,rn_si1) and a ratio
      !!      (rn_abs). These parameters are read in the namtra_qsr namelist. The
      !!      default values correspond to clear water (type I in Jerlov'
      !!      (1968) classification.
      !!         called by tra_qsr at the first timestep (nit000)
      !!
      !! ** Action  : - initialize rn_si0, rn_si1 and rn_abs
      !!
      !! Reference : Jerlov, N. G., 1968 Optical Oceanography, Elsevier, 194pp.
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk                  ! dummy loop indices
      INTEGER  ::   ios, irgb, ierror, ioptio   ! local integer
      REAL(wp) ::   zz0, zc0 , zc1, zcoef       ! local scalars
      REAL(wp) ::   zz1, zc2 , zc3, zchl        !   -      -
      !
      CHARACTER(len=100) ::   cn_dir   ! Root directory for location of ssr files
      TYPE(FLD_N)        ::   sn_chl   ! informations about the chlorofyl field to be read
      !!
      NAMELIST/namtra_qsr/  sn_chl, cn_dir, ln_qsr_rgb, ln_qsr_2bd, ln_qsr_bio,  &
         &                  nn_chldta, rn_chl_conc, rn_abs, rn_si0, rn_si1
      !!----------------------------------------------------------------------
      !
      READ  ( numnam_ref, namtra_qsr, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namtra_qsr in reference namelist' )
      !
      READ  ( numnam_cfg, namtra_qsr, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namtra_qsr in configuration namelist' )
      IF(lwm) WRITE ( numond, namtra_qsr )
      !
      IF(lwp) THEN                ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'tra_qsr_init : penetration of the surface solar radiation'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namtra_qsr : set the parameter of penetration'
         WRITE(numout,*) '      RGB (Red-Green-Blue) light penetration       ln_qsr_rgb = ', ln_qsr_rgb
         WRITE(numout,*) '      2 band               light penetration       ln_qsr_2bd = ', ln_qsr_2bd
         WRITE(numout,*) '      bio-model            light penetration       ln_qsr_bio = ', ln_qsr_bio
         WRITE(numout,*) '      RGB : Chl data (=1) or cst value (=0)        nn_chldta  = ', nn_chldta
         WRITE(numout,*) '      Chlorophyll concentration (for nn_chldta=0)  rn_chl_conc = ', rn_chl_conc
         WRITE(numout,*) '      RGB & 2 bands: fraction of light (rn_si1)    rn_abs     = ', rn_abs
         WRITE(numout,*) '      RGB & 2 bands: shortess depth of extinction  rn_si0     = ', rn_si0
         WRITE(numout,*) '      2 bands: longest depth of extinction         rn_si1     = ', rn_si1
         WRITE(numout,*)
      ENDIF
      !
      ioptio = 0                    ! Parameter control
      IF( ln_qsr_rgb  )   ioptio = ioptio + 1
      IF( ln_qsr_2bd  )   ioptio = ioptio + 1
      IF( ln_qsr_bio  )   ioptio = ioptio + 1
      !
      IF( ioptio /= 1 )   CALL ctl_stop( 'Choose ONE type of light penetration in namelist namtra_qsr',  &
         &                               ' 2 bands, 3 RGB bands or bio-model light penetration' )
      !
      IF( ln_qsr_rgb .AND. nn_chldta == 0 )   nqsr = np_RGB
      IF( ln_qsr_rgb .AND. nn_chldta == 1 )   nqsr = np_RGBc
      IF( ln_qsr_2bd                      )   nqsr = np_2BD
      IF( ln_qsr_bio                      )   nqsr = np_BIO
      !
      !                             ! Initialisation
      xsi0r = 1._wp / rn_si0
      xsi1r = 1._wp / rn_si1
      !
      SELECT CASE( nqsr )
      !
      CASE( np_RGB , np_RGBc )         !==  Red-Green-Blue light penetration  ==!
         !
         IF(lwp)   WRITE(numout,*) '   ==>>>   R-G-B   light penetration '
         !
         CALL trc_oce_rgb( rkrgb )                 ! tabulated attenuation coef.
         !
         nksr = trc_oce_ext_lev( r_si2, 33._wp )   ! level of light extinction
         !
         IF(lwp) WRITE(numout,*) '        level of light extinction = ', nksr, ' ref depth = ', gdepw_1d(nksr+1), ' m'
         !
         IF( nqsr == np_RGBc ) THEN                ! Chl data : set sf_chl structure
            IF(lwp) WRITE(numout,*) '   ==>>>   Chlorophyll read in a file'
            ALLOCATE( sf_chl(1), STAT=ierror )
            IF( ierror > 0 ) THEN
               CALL ctl_stop( 'tra_qsr_init: unable to allocate sf_chl structure' )   ;   RETURN
            ENDIF
            ALLOCATE( sf_chl(1)%fnow(jpi,jpj,1)   )
            IF( sn_chl%ln_tint )   ALLOCATE( sf_chl(1)%fdta(jpi,jpj,1,2) )
            !                                        ! fill sf_chl with sn_chl and control print
            CALL fld_fill( sf_chl, (/ sn_chl /), cn_dir, 'tra_qsr_init',   &
               &           'Solar penetration function of read chlorophyll', 'namtra_qsr' , no_print )
         ENDIF
         IF( nqsr == np_RGB ) THEN                 ! constant Chl
            IF(lwp) WRITE(numout,*) '   ==>>>   Constant Chlorophyll concentration = ', rn_chl_conc
         ENDIF
         !
      CASE( np_2BD )                   !==  2 bands light penetration  ==!
         !
         IF(lwp)  WRITE(numout,*) '   ==>>>   2 bands light penetration'
         !
         nksr = trc_oce_ext_lev( rn_si1, 100._wp )    ! level of light extinction
         IF(lwp) WRITE(numout,*) '        level of light extinction = ', nksr, ' ref depth = ', gdepw_1d(nksr+1), ' m'
         !
      CASE( np_BIO )                   !==  BIO light penetration  ==!
         !
         IF(lwp) WRITE(numout,*) '   ==>>>   bio-model light penetration'
         IF( .NOT.lk_top )   CALL ctl_stop( 'No bio model : ln_qsr_bio = true impossible ' )
         !
         CALL trc_oce_rgb( rkrgb )                 ! tabulated attenuation coef.
         !
         nksr = trc_oce_ext_lev( r_si2, 33._wp )   ! level of light extinction
         !
         IF(lwp) WRITE(numout,*) '        level of light extinction = ', nksr, ' ref depth = ', gdepw_1d(nksr+1), ' m'
         !
      END SELECT
      !
      qsr_hc(:,:,:) = 0._wp     ! now qsr heat content set to zero where it will not be computed
      !
      ! 1st ocean level attenuation coefficient (used in sbcssm)
      IF( iom_varid( numror, 'fraqsr_1lev', ldstop = .FALSE. ) > 0 ) THEN
         CALL iom_get( numror, jpdom_auto, 'fraqsr_1lev'  , fraqsr_1lev  )
      ELSE
         fraqsr_1lev(:,:) = 1._wp   ! default : no penetration
      ENDIF
      !
   END SUBROUTINE tra_qsr_init

   !!======================================================================
END MODULE traqsr
