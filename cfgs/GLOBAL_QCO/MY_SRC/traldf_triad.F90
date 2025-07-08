MODULE traldf_triad
   !!======================================================================
   !!                   ***  MODULE  traldf_triad  ***
   !! Ocean  tracers:  horizontal component of the lateral tracer mixing trend
   !!======================================================================
   !! History :  3.3  ! 2010-10  (G. Nurser, C. Harris, G. Madec)  Griffies operator (original code)
   !!            3.7  ! 2013-12  (F. Lemarie, G. Madec)  triad operator (Griffies) + Method of Stabilizing Correction
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_ldf_triad : update the tracer trend with the iso-neutral laplacian triad-operator
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and active tracers
   USE dom_oce        ! ocean space and time domain
   USE domutl, ONLY : is_tile
   USE phycst         ! physical constants
   USE trc_oce        ! share passive tracers/Ocean variables
   USE zdf_oce        ! ocean vertical physics
   USE ldftra         ! lateral physics: eddy diffusivity
   USE ldfslp         ! lateral physics: iso-neutral slopes
   USE traldf_iso     ! lateral diffusion (Madec operator)         (tra_ldf_iso routine)
   USE diaptr         ! poleward transport diagnostics
   USE diaar5         ! AR5 diagnostics
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_ldf_triad   ! routine called by traldf.F90

   LOGICAL  ::   l_ptr   ! flag to compute poleward transport
   LOGICAL  ::   l_hst   ! flag to compute heat transport


   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: traldf_triad.F90 15062 2021-06-28 11:19:48Z jchanut $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_ldf_triad( kt, Kmm, kit000, cdtype, pahu, pahv,             &
      &                                               pgu , pgv , pgui, pgvi, &
      &                                               pt, pt2, pt_rhs, kjpt, kpass )
      !!
      INTEGER                     , INTENT(in   ) ::   kt         ! ocean time-step index
      INTEGER                     , INTENT(in   ) ::   kit000     ! first time step index
      CHARACTER(len=3)            , INTENT(in   ) ::   cdtype     ! =TRA or TRC (tracer indicator)
      INTEGER                     , INTENT(in   ) ::   kjpt       ! number of tracers
      INTEGER                     , INTENT(in   ) ::   kpass      ! =1/2 first or second passage
      INTEGER                     , INTENT(in   ) ::   Kmm        ! ocean time level indices
      REAL(wp), DIMENSION(:,:,:)  , INTENT(in   ) ::   pahu, pahv ! eddy diffusivity at u- and v-points  [m2/s]
      REAL(wp), DIMENSION(:,:,:)  , INTENT(in   ) ::   pgu , pgv  ! tracer gradient at pstep levels
      REAL(wp), DIMENSION(:,:,:)  , INTENT(in   ) ::   pgui, pgvi ! tracer gradient at top   levels
      REAL(dp), DIMENSION(:,:,:,:), INTENT(in   ) ::   pt         ! tracer (kpass=1) or laplacian of tracer (kpass=2)
      REAL(dp), DIMENSION(:,:,:,:), INTENT(in   ) ::   pt2        ! tracer (only used in kpass=2)
      REAL(dp), DIMENSION(:,:,:,:), INTENT(inout) ::   pt_rhs     ! tracer trend
      !!
      CALL tra_ldf_triad_t( kt, Kmm, kit000, cdtype, pahu, pahv, is_tile(pahu),                            &
      &                                              pgu , pgv , is_tile(pgu) , pgui, pgvi, is_tile(pgui), &
      &                                              pt, is_tile(pt), pt2, is_tile(pt2), pt_rhs, is_tile(pt_rhs), kjpt, kpass )
   END SUBROUTINE tra_ldf_triad


  SUBROUTINE tra_ldf_triad_t( kt, Kmm, kit000, cdtype, pahu, pahv, ktah,                   &
      &                                                pgu , pgv , ktg , pgui, pgvi, ktgi, &
      &                                                pt, ktt, pt2, ktt2, pt_rhs, ktt_rhs, kjpt, kpass )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_ldf_triad  ***
      !!
      !! ** Purpose :   Compute the before horizontal tracer (t & s) diffusive
      !!      trend for a laplacian tensor (ezxcept the dz[ dz[.] ] term) and
      !!      add it to the general trend of tracer equation.
      !!
      !! ** Method  :   The horizontal component of the lateral diffusive trends
      !!      is provided by a 2nd order operator rotated along neural or geopo-
      !!      tential surfaces to which an eddy induced advection can be added
      !!      It is computed using before fields (forward in time) and isopyc-
      !!      nal or geopotential slopes computed in routine ldfslp.
      !!
      !!      see documentation for the desciption
      !!
      !! ** Action :   pt_rhs   updated with the before rotated diffusion
      !!               ah_wslp2 ....
      !!               akz   stabilizing vertical diffusivity coefficient (used in trazdf_imp)
      !!----------------------------------------------------------------------
      INTEGER                              , INTENT(in   ) ::   kt         ! ocean time-step index
      INTEGER                              , INTENT(in   ) ::   kit000     ! first time step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype     ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt       ! number of tracers
      INTEGER                              , INTENT(in   ) ::   kpass      ! =1/2 first or second passage
      INTEGER                              , INTENT(in)    ::   Kmm        ! ocean time level indices
      INTEGER                              , INTENT(in   ) ::   ktah, ktg, ktgi, ktt, ktt2, ktt_rhs
      REAL(wp), DIMENSION(A2D_T(ktah),   JPK)     , INTENT(in   ) ::   pahu, pahv ! eddy diffusivity at u- and v-points  [m2/s]
      REAL(wp), DIMENSION(A2D_T(ktg),        KJPT), INTENT(in   ) ::   pgu , pgv  ! tracer gradient at pstep levels
      REAL(wp), DIMENSION(A2D_T(ktgi),       KJPT), INTENT(in   ) ::   pgui, pgvi ! tracer gradient at top   levels
      REAL(dp), DIMENSION(A2D_T(ktt),    JPK,KJPT), INTENT(in   ) ::   pt         ! tracer (kpass=1) or laplacian of tracer (kpass=2)
      REAL(dp), DIMENSION(A2D_T(ktt2),   JPK,KJPT), INTENT(in   ) ::   pt2        ! tracer (only used in kpass=2)
      REAL(dp), DIMENSION(A2D_T(ktt_rhs),JPK,KJPT), INTENT(inout) ::   pt_rhs     ! tracer trend
      !
      INTEGER  ::  ji, jj, jk, jn, kp, iij   ! dummy loop indices
      REAL(wp) ::  zcoef0, ze3w_2, zsign          !   -      -
      !
      REAL(wp) ::   zslope2, zbu, zbv, zbu1, zbv1, zslope21, zah, zah1, zah_ip1, zah_jp1, zbu_ip1, zbv_jp1
      REAL(wp) ::   ze1ur, ze2vr, ze3wr, zdxt, zdyt, zdzt, zdyt_jp1, ze3wr_jp1, zdzt_jp1, zah_slp1, zah_slp_jp1, zaei_slp_jp1
      REAL(wp) ::   zah_slp, zaei_slp, zdxt_ip1, ze3wr_ip1, zdzt_ip1, zah_slp_ip1, zaei_slp_ip1, zaei_slp1
      REAL(wp), DIMENSION(A2D(nn_hls),0:1) ::   zdkt3d                                           ! vertical tracer gradient at 2 levels
      REAL(wp), DIMENSION(A2D(nn_hls)    ) ::   z2d                                              ! 2D workspace
      REAL(wp), DIMENSION(A2D(nn_hls),jpk) ::   zdit, zdjt, zftu, zftv, ztfw, zpsi_uw, zpsi_vw   ! 3D     -
      !!----------------------------------------------------------------------
      !
      IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                       ! Do only on the first tile
         IF( kpass == 1 .AND. kt == kit000 )  THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'tra_ldf_triad : rotated laplacian diffusion operator on ', cdtype
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~'
         ENDIF
         !
         l_hst = .FALSE.
         l_ptr = .FALSE.
         IF( cdtype == 'TRA' ) THEN
            IF( iom_use( 'sophtldf' ) .OR. iom_use( 'sopstldf') )      l_ptr = .TRUE.
            IF( iom_use("uadv_heattr") .OR. iom_use("vadv_heattr") .OR.                   &
            &   iom_use("uadv_salttr") .OR. iom_use("vadv_salttr")  )   l_hst = .TRUE.
         ENDIF

         IF( ln_ldfeiv ) THEN
            CALL lbc_lnk( 'traldf_triad', triadi_g(:,:,:,1,:), 'U', -1.0_wp,  &
               &                          triadj_g(:,:,:,1,:), 'V', -1.0_wp )
            CALL lbc_lnk( 'traldf_triad', triadi_g(:,:,:,0,:), 'U', -1.0_wp,  &
               &                          triadj_g(:,:,:,0,:), 'V', -1.0_wp )
         ENDIF
         CALL lbc_lnk( 'traldf_triad', triadi(:,:,:,1,:), 'U', -1.0_wp,  &
            &                          triadj(:,:,:,1,:), 'V', -1.0_wp )
         CALL lbc_lnk( 'traldf_triad', triadi(:,:,:,0,:), 'U', -1.0_wp,  &
            &                          triadj(:,:,:,0,:), 'V', -1.0_wp )

      ENDIF
      !
      ! Define pt_rhs halo points for multi-point haloes in bilaplacian case
      IF( nldf_tra == np_blp_it .AND. kpass == 1 ) THEN ; iij = nn_hls
      ELSE                                              ; iij = 1
      ENDIF

      !
      IF( kpass == 1 ) THEN   ;   zsign =  1._wp      ! bilaplacian operator require a minus sign (eddy diffusivity >0)
      ELSE                    ;   zsign = -1._wp
      ENDIF
      !
      !!----------------------------------------------------------------------
      !!   0 - calculate  ah_wslp2, akz, and optionally zpsi_uw, zpsi_vw
      !!----------------------------------------------------------------------
      !
      IF( kpass == 1 ) THEN         !==  first pass only  and whatever the tracer is  ==!
         !
         DO_3D_OVR( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, 1, jpk )
            akz     (ji,jj,jk) = 0._wp
            ah_wslp2(ji,jj,jk) = 0._wp
         END_3D
         !
         DO kp = 0, 1                            ! i-k triads
            DO_3D_OVR( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, 1, jpkm1 )
               ze3wr = 1._wp / e3w(ji,jj,jk+kp,Kmm)
               zbu   = e1e2u(ji,jj) * e3u(ji,jj,jk,Kmm)
               zbu1  = e1e2u(ji-1,jj) * e3u(ji-1,jj,jk,Kmm)
               zah   = 0.25_wp * pahu(ji,jj,jk)
               zah1  = 0.25_wp * pahu(ji-1,jj,jk)
               ! Subtract s-coordinate slope at t-points to give slope rel to s-surfaces (do this by *adding* gradient of depth)
               zslope2 = triadi_g(ji,jj,jk,1,kp) + ( gdept(ji+1,jj,jk,Kmm) - gdept(ji,jj,jk,Kmm) ) * r1_e1u(ji,jj) * umask(ji,jj,jk+kp)
               zslope2 = zslope2 *zslope2
               zslope21 = triadi_g(ji,jj,jk,0,kp) + ( gdept(ji,jj,jk,Kmm) - gdept(ji-1,jj,jk,Kmm) ) * r1_e1u(ji-1,jj) * umask(ji-1,jj,jk+kp)
               zslope21 = zslope21 *zslope21
               ! round brackets added to fix the order of floating point operations
               ! needed to ensure halo 1 - halo 2 compatibility
               ah_wslp2(ji,jj,jk+kp) =  ah_wslp2(ji,jj,jk+kp) + ( zah * zbu * ze3wr * r1_e1e2t(ji,jj) * zslope2                    &
                        &                                       + zah1 * zbu1 * ze3wr * r1_e1e2t(ji,jj) * zslope21                 &
                        &                                       )                                                                  ! bracket for halo 1 - halo 2 compatibility
               akz     (ji,jj,jk+kp) =  akz     (ji,jj,jk+kp) + ( zah * r1_e1u(ji,jj) * r1_e1u(ji,jj) * umask(ji,jj,jk+kp)         &
                                                                + zah1 * r1_e1u(ji-1,jj) * r1_e1u(ji-1,jj) * umask(ji-1,jj,jk+kp)  &
                        &                                       )                                                                  ! bracket for halo 1 - halo 2 compatibility
            END_3D
         END DO
         !
         DO kp = 0, 1                            ! j-k triads
            DO_3D_OVR( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, 1, jpkm1 )
               ze3wr = 1.0_wp / e3w(ji,jj,jk+kp,Kmm)
               zbv   = e1e2v(ji,jj) * e3v(ji,jj,jk,Kmm)
               zbv1   = e1e2v(ji,jj-1) * e3v(ji,jj-1,jk,Kmm)
               zah   = 0.25_wp * pahv(ji,jj,jk)
               zah1   = 0.25_wp * pahv(ji,jj-1,jk)
               ! Subtract s-coordinate slope at t-points to give slope rel to s surfaces
               !    (do this by *adding* gradient of depth)
               zslope2 = triadj_g(ji,jj,jk,1,kp) + ( gdept(ji,jj+1,jk,Kmm) - gdept(ji,jj,jk,Kmm) ) * r1_e2v(ji,jj) * vmask(ji,jj,jk+kp)
               zslope2 = zslope2 * zslope2
               zslope21 = triadj_g(ji,jj,jk,0,kp) + ( gdept(ji,jj,jk,Kmm) - gdept(ji,jj-1,jk,Kmm) ) * r1_e2v(ji,jj-1) * vmask(ji,jj-1,jk+kp)
               zslope21 = zslope21 * zslope21
               ! round brackets added to fix the order of floating point operations
               ! needed to ensure halo 1 - halo 2 compatibility
               ah_wslp2(ji,jj,jk+kp) = ah_wslp2(ji,jj,jk+kp) + ( zah * zbv * ze3wr * r1_e1e2t(ji,jj) * zslope2                     &
                        &                                      + zah1 * zbv1 * ze3wr * r1_e1e2t(ji,jj) * zslope21                  &
                        &                                      )                                                                   ! bracket for halo 1 - halo 2 compatibility
               akz     (ji,jj,jk+kp) = akz     (ji,jj,jk+kp) + ( zah * r1_e2v(ji,jj) * r1_e2v(ji,jj) * vmask(ji,jj,jk+kp)          &
                        &                                      + zah1 * r1_e2v(ji,jj-1) * r1_e2v(ji,jj-1) * vmask(ji,jj-1,jk+kp)   &
                        &                                      )                                                                   ! bracket for halo 1 - halo 2 compatibility
            END_3D
         END DO
         !
         IF( ln_traldf_msc ) THEN                ! stabilizing vertical diffusivity coefficient
            !
            IF( ln_traldf_blp ) THEN                ! bilaplacian operator
               DO_3D_OVR( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, 2, jpkm1 )
                  akz(ji,jj,jk) = 16._wp           &
                     &   * ah_wslp2   (ji,jj,jk)   &
                     &   * (  akz     (ji,jj,jk)   &
                     &      + ah_wslp2(ji,jj,jk)   &
                     &        / ( e3w(ji,jj,jk,Kmm) * e3w(ji,jj,jk,Kmm) )  )
               END_3D
            ELSEIF( ln_traldf_lap ) THEN              ! laplacian operator
               DO_3D_OVR( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, 2, jpkm1 )
                  ze3w_2 = e3w(ji,jj,jk,Kmm) * e3w(ji,jj,jk,Kmm)
                  zcoef0 = rDt * (  akz(ji,jj,jk) + ah_wslp2(ji,jj,jk) / ze3w_2  )
                  akz(ji,jj,jk) = MAX( zcoef0 - 0.5_wp , 0._wp ) * ze3w_2 * r1_Dt
               END_3D
           ENDIF
           !
         ELSE                                    ! 33 flux set to zero with akz=ah_wslp2 ==>> computed in full implicit
            DO_3D_OVR( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, 1, jpk )
               akz(ji,jj,jk) = ah_wslp2(ji,jj,jk)
            END_3D
         ENDIF
         !
         IF( ln_ldfeiv_dia .AND. cdtype == 'TRA' ) THEN
            zpsi_uw(:,:,:) = 0._wp
            zpsi_vw(:,:,:) = 0._wp

            DO kp = 0, 1
               DO_3D( 1, 0, 1, 0, 1, jpkm1 )
                  ! round brackets added to fix the order of floating point operations
                  ! needed to ensure halo 1 - halo 2 compatibility
                  zpsi_uw(ji,jj,jk+kp) = zpsi_uw(ji,jj,jk+kp)                                     &
                     & + ( 0.25_wp * aeiu(ji,jj,jk) * e2u(ji,jj) * triadi_g(ji,jj,jk,1,kp)        &
                     &   + 0.25_wp * aeiu(ji,jj,jk) * e2u(ji,jj) * triadi_g(ji+1,jj,jk,0,kp)      &
                     &   )                                                                        ! bracket for halo 1 - halo 2 compatibility
                  zpsi_vw(ji,jj,jk+kp) = zpsi_vw(ji,jj,jk+kp)                                     &
                     & + ( 0.25_wp * aeiv(ji,jj,jk) * e1v(ji,jj) * triadj_g(ji,jj,jk,1,kp)        &
                     &   + 0.25_wp * aeiv(ji,jj,jk) * e1v(ji,jj) * triadj_g(ji,jj+1,jk,0,kp)      &
                     &   )                                                                        ! bracket for halo 1 - halo 2 compatibility
               END_3D
            END DO
            CALL ldf_eiv_dia( zpsi_uw, zpsi_vw, Kmm )
         ENDIF
         !
      ENDIF                                  !==  end 1st pass only  ==!
      !
      !                                                           ! ===========
      DO jn = 1, kjpt                                             ! tracer loop
         !                                                        ! ===========
         ! Zero fluxes for each tracer
!!gm  this should probably be done outside the jn loop
         ztfw(:,:,:) = 0._wp
         zftu(:,:,:) = 0._wp
         zftv(:,:,:) = 0._wp
         zdit(:,:,:) = 0._wp
         zdjt(:,:,:) = 0._wp
         !
         DO_3D( iij, iij-1, iij, iij-1, 1, jpkm1 )    !==  before lateral T & S gradients at T-level jk  ==!
            zdit(ji,jj,jk) = ( pt(ji+1,jj  ,jk,jn) - pt(ji,jj,jk,jn) ) * umask(ji,jj,jk)
            zdjt(ji,jj,jk) = ( pt(ji  ,jj+1,jk,jn) - pt(ji,jj,jk,jn) ) * vmask(ji,jj,jk)
         END_3D
         IF( ln_zps .AND. l_grad_zps ) THEN    ! partial steps: correction at top/bottom ocean level
            DO_2D( iij, iij-1, iij, iij-1 )                    ! bottom level
               zdit(ji,jj,mbku(ji,jj)) = pgu(ji,jj,jn)
               zdjt(ji,jj,mbkv(ji,jj)) = pgv(ji,jj,jn)
            END_2D
            IF( ln_isfcav ) THEN                   ! top level (ocean cavities only)
               DO_2D( iij, iij-1, iij, iij-1 )
                  IF( miku(ji,jj)  > 1 )   zdit(ji,jj,miku(ji,jj) ) = pgui(ji,jj,jn)
                  IF( mikv(ji,jj)  > 1 )   zdjt(ji,jj,mikv(ji,jj) ) = pgvi(ji,jj,jn)
               END_2D
            ENDIF
         ENDIF
         !
         !!----------------------------------------------------------------------
         !!   II - horizontal trend  (full)
         !!----------------------------------------------------------------------
         !
         DO jk = 1, jpkm1
            !                    !==  Vertical tracer gradient at level jk and jk+1
            DO_2D( iij, iij, iij, iij )
               zdkt3d(ji,jj,1) = ( pt(ji,jj,jk,jn) - pt(ji,jj,jk+1,jn) ) * tmask(ji,jj,jk+1)
            END_2D
            !
            !                    ! surface boundary condition: zdkt3d(jk=0)=zdkt3d(jk=1)
            IF( jk == 1 ) THEN   ;   zdkt3d(:,:,0) = zdkt3d(:,:,1)
            ELSE
               DO_2D( iij, iij, iij, iij )
                  zdkt3d(ji,jj,0) = ( pt(ji,jj,jk-1,jn) - pt(ji,jj,jk,jn) ) * tmask(ji,jj,jk)
               END_2D
            ENDIF
            !
            zaei_slp = 0._wp
            zaei_slp_ip1 = 0._wp
            zaei_slp_jp1 = 0._wp
            zaei_slp1 = 0._wp
            !
            IF( ln_botmix_triad ) THEN
               DO kp = 0, 1              !==  Horizontal & vertical fluxes
                  DO_2D( iij, iij-1, iij, iij-1 )
                     ze1ur = r1_e1u(ji,jj)
                     zdxt  = zdit(ji,jj,jk) * ze1ur
                     zdxt_ip1  = zdit(ji+1,jj,jk) * r1_e1u(ji+1,jj)
                     ze3wr = 1._wp / e3w(ji,jj,jk+kp,Kmm)
                     ze3wr_ip1 = 1._wp / e3w(ji+1,jj,jk+kp,Kmm)
                     zdzt  = zdkt3d(ji,jj,kp) * ze3wr
                     zdzt_ip1  = zdkt3d(ji+1,jj,kp) * ze3wr_ip1
                     !
                     zbu = 0.25_wp * e1e2u(ji,jj) * e3u(ji,jj,jk,Kmm)
                     zbu_ip1 = 0.25_wp * e1e2u(ji+1,jj) * e3u(ji+1,jj,jk,Kmm)
                     ! ln_botmix_triad is .T. don't mask zah for bottom half cells    !!gm ?????   ahu is masked....
                     zah = pahu(ji,jj,jk)
                     zah_ip1 = pahu(ji+1,jj,jk)
                     zah_slp  = zah * triadi(ji,jj,jk,1,kp)
                     zah_slp_ip1  = zah_ip1 * triadi(ji+1,jj,jk,1,kp)
                     zah_slp1  = zah * triadi(ji+1,jj,jk,0,kp)
                     IF( ln_ldfeiv )   THEN
                        zaei_slp = aeiu(ji,jj,jk) * triadi_g(ji,jj,jk,1,kp)
                        zaei_slp_ip1 = aeiu(ji+1,jj,jk) * triadi_g(ji+1,jj,jk,1,kp)
                        zaei_slp1 = aeiu(ji,jj,jk) * triadi_g(ji+1,jj,jk,0,kp)
                     ENDIF
                     ! round brackets added to fix the order of floating point operations
                     ! needed to ensure halo 1 - halo 2 compatibility
                     zftu(ji   ,jj,jk  ) =  zftu(ji   ,jj,jk )                                                               &
                                         &    - ( ( zah * zdxt + ( zah_slp - zaei_slp ) * zdzt ) * zbu * ze1ur               &
                                         &      + ( zah * zdxt + zah_slp1 * zdzt_ip1 - zaei_slp1 * zdzt_ip1 ) * zbu * ze1ur  &
                                         &      )                                                                            ! bracket for halo 1 - halo 2 compatibility
                     ztfw(ji+1,jj,jk+kp) =  ztfw(ji+1,jj,jk+kp)                                                              &
                                         &    - ( (zah_slp_ip1 + zaei_slp_ip1) * zdxt_ip1 * zbu_ip1 * ze3wr_ip1              &
                                         &      + ( zah_slp1 + zaei_slp1) * zdxt * zbu * ze3wr_ip1                           &
                                         &      )                                                                            ! bracket for halo 1 - halo 2 compatibility
                  END_2D
               END DO
               !
               DO kp = 0, 1
                  DO_2D( iij, iij-1, iij, iij-1 )
                     ze2vr = r1_e2v(ji,jj)
                     zdyt  = zdjt(ji,jj,jk) * ze2vr
                     zdyt_jp1  = zdjt(ji,jj+1,jk) * r1_e2v(ji,jj+1)
                     ze3wr = 1._wp / e3w(ji,jj,jk+kp,Kmm)
                     ze3wr_jp1 = 1._wp / e3w(ji,jj+1,jk+kp,Kmm)
                     zdzt  = zdkt3d(ji,jj,kp) * ze3wr
                     zdzt_jp1  = zdkt3d(ji,jj+1,kp) * ze3wr_jp1
                     zbv = 0.25_wp * e1e2v(ji,jj) * e3v(ji,jj,jk,Kmm)
                     zbv_jp1 = 0.25_wp * e1e2v(ji,jj+1) * e3v(ji,jj+1,jk,Kmm)
                     ! ln_botmix_triad is .T. don't mask zah for bottom half cells    !!gm ?????   ahu is masked....
                     zah = pahv(ji,jj,jk)          ! pahv(ji,jj+jp,jk)  ????
                     zah_jp1 = pahv(ji,jj+1,jk)
                     zah_slp = zah * triadj(ji,jj,jk,1,kp)
                     zah_slp1 = zah * triadj(ji,jj+1,jk,0,kp)
                     zah_slp_jp1 = zah_jp1 * triadj(ji,jj+1,jk,1,kp)
                     IF( ln_ldfeiv )   THEN
                        zaei_slp = aeiv(ji,jj,jk) * triadj_g(ji,jj,jk,1,kp)
                        zaei_slp_jp1 = aeiv(ji,jj+1,jk) * triadj_g(ji,jj+1,jk,1,kp)
                        zaei_slp1 = aeiv(ji,jj,jk) * triadj_g(ji,jj+1,jk,0,kp)
                     ENDIF
                     ! round brackets added to fix the order of floating point operations
                     ! needed to ensure halo 1 - halo 2 compatibility
                     zftv(ji,jj  ,jk   ) =  zftv(ji,jj  ,jk   )                                                              &
                                         &    - ( ( zah * zdyt + ( zah_slp - zaei_slp ) * zdzt ) * zbv * ze2vr               &
                                         &      + ( zah * zdyt + zah_slp1 * zdzt_jp1 - zaei_slp1 * zdzt_jp1 ) * zbv * ze2vr  &
                                         &      )                                                                            ! bracket for halo 1 - halo 2 compatibility
                     ztfw(ji,jj+1,jk+kp) =  ztfw(ji,jj+1,jk+kp)                                                              &
                                         &    - ( ( zah_slp_jp1 + zaei_slp_jp1) * zdyt_jp1 * zbv_jp1 * ze3wr_jp1             &
                                         &      + ( zah_slp1 + zaei_slp1) * zdyt * zbv * ze3wr_jp1                           &
                                         &      )                                                                            ! bracket for halo 1 - halo 2 compatibility
                  END_2D
               END DO
               !
            ELSE
               !
               DO kp = 0, 1               !==  Horizontal & vertical fluxes
                  DO_2D( iij, iij-1, iij, iij-1 )
                     ze1ur = r1_e1u(ji,jj)
                     zdxt  = zdit(ji,jj,jk) * ze1ur
                     zdxt_ip1  = zdit(ji+1,jj,jk) * r1_e1u(ji+1,jj)
                     ze3wr = 1._wp / e3w(ji,jj,jk+kp,Kmm)
                     ze3wr_ip1 = 1._wp / e3w(ji+1,jj,jk+kp,Kmm)
                     zdzt  = zdkt3d(ji,jj,kp) * ze3wr
                     zdzt_ip1  = zdkt3d(ji+1,jj,kp) * ze3wr_ip1
                     !
                     zbu = 0.25_wp * e1e2u(ji,jj) * e3u(ji,jj,jk,Kmm)
                     zbu_ip1 = 0.25_wp * e1e2u(ji+1,jj) * e3u(ji+1,jj,jk,Kmm)
                     ! ln_botmix_triad is .F. mask zah for bottom half cells
                     zah = pahu(ji,jj,jk) * umask(ji,jj,jk+kp)         ! pahu(ji+ip,jj,jk)   ===>>  ????
                     zah_ip1 = pahu(ji+1,jj,jk) * umask(ji+1,jj,jk+kp)
                     zah_slp  = zah * triadi(ji,jj,jk,1,kp)
                     zah_slp_ip1  = zah_ip1 * triadi(ji+1,jj,jk,1,kp)
                     zah_slp1  = zah * triadi(ji+1,jj,jk,0,kp)
                     IF( ln_ldfeiv )   THEN
                        zaei_slp = aeiu(ji,jj,jk) * triadi_g(ji,jj,jk,1,kp)
                        zaei_slp_ip1 = aeiu(ji+1,jj,jk) * triadi_g(ji+1,jj,jk,1,kp)
                        zaei_slp1 = aeiu(ji,jj,jk) * triadi_g(ji+1,jj,jk,0,kp)
                     ENDIF
                     ! round brackets added to fix the order of floating point operations
                     ! needed to ensure halo 1 - halo 2 compatibility
                     zftu(ji   ,jj,jk  ) =  zftu(ji   ,jj,jk )                                                               &
                                         &    - ( ( zah * zdxt + ( zah_slp - zaei_slp ) * zdzt ) * zbu * ze1ur               &
                                         &      + ( zah * zdxt + zah_slp1 * zdzt_ip1 - zaei_slp1 * zdzt_ip1 ) * zbu * ze1ur  &
                                         &      )                                                                            ! bracket for halo 1 - halo 2 compatibility
                     ztfw(ji+1,jj,jk+kp) =  ztfw(ji+1,jj,jk+kp)                                                              &
                                         &    - ( (zah_slp_ip1 + zaei_slp_ip1) * zdxt_ip1 * zbu_ip1 * ze3wr_ip1              &
                                         &      + ( zah_slp1 + zaei_slp1) * zdxt * zbu * ze3wr_ip1                           &
                                         &      )                                                                            ! bracket for halo 1 - halo 2 compatibility
                  END_2D
               END DO
               !
               DO kp = 0, 1
                  DO_2D( iij, iij-1, iij, iij-1 )
                     ze2vr = r1_e2v(ji,jj)
                     zdyt  = zdjt(ji,jj,jk) * ze2vr
                     zdyt_jp1  = zdjt(ji,jj+1,jk) * r1_e2v(ji,jj+1)
                     ze3wr = 1._wp / e3w(ji,jj,jk+kp,Kmm)
                     ze3wr_jp1 = 1._wp / e3w(ji,jj+1,jk+kp,Kmm)
                     zdzt  = zdkt3d(ji,jj,kp) * ze3wr
                     zdzt_jp1  = zdkt3d(ji,jj+1,kp) * ze3wr_jp1
                     zbv = 0.25_wp * e1e2v(ji,jj) * e3v(ji,jj,jk,Kmm)
                     zbv_jp1 = 0.25_wp * e1e2v(ji,jj+1) * e3v(ji,jj+1,jk,Kmm)
                     ! ln_botmix_triad is .F. mask zah for bottom half cells
                     zah = pahv(ji,jj,jk) * vmask(ji,jj,jk+kp)         ! pahv(ji,jj+jp,jk)  ????
                     zah_jp1 = pahv(ji,jj+1,jk) * vmask(ji,jj+1,jk+kp)
                     zah_slp = zah * triadj(ji,jj,jk,1,kp)
                     zah_slp1 = zah * triadj(ji,jj+1,jk,0,kp)
                     zah_slp_jp1 = zah_jp1 * triadj(ji,jj+1,jk,1,kp)
                     IF( ln_ldfeiv )   THEN
                        zaei_slp = aeiv(ji,jj,jk) * triadj_g(ji,jj,jk,1,kp)
                        zaei_slp_jp1 = aeiv(ji,jj+1,jk) * triadj_g(ji,jj+1,jk,1,kp)
                        zaei_slp1 = aeiv(ji,jj,jk) * triadj_g(ji,jj+1,jk,0,kp)
                     ENDIF
                     ! round brackets added to fix the order of floating point operations
                     ! needed to ensure halo 1 - halo 2 compatibility
                     zftv(ji,jj  ,jk   ) =  zftv(ji,jj  ,jk   )                                                              &
                                         &    - ( ( zah * zdyt + ( zah_slp - zaei_slp ) * zdzt ) * zbv * ze2vr               &
                                         &      + ( zah * zdyt + zah_slp1 * zdzt_jp1 - zaei_slp1 * zdzt_jp1 ) * zbv * ze2vr  &
                                         &      )                                                                            ! bracket for halo 1 - halo 2 compatibility
                     ztfw(ji,jj+1,jk+kp) =  ztfw(ji,jj+1,jk+kp)                                                              &
                                         &    - ( ( zah_slp_jp1 + zaei_slp_jp1) * zdyt_jp1 * zbv_jp1 * ze3wr_jp1             &
                                         &      + ( zah_slp1 + zaei_slp1) * zdyt * zbv * ze3wr_jp1                           &
                                         &      )                                                                            ! bracket for halo 1 - halo 2 compatibility
                  END_2D
               END DO
            ENDIF
            !                             !==  horizontal divergence and add to the general trend  ==!
            DO_2D( iij-1, iij-1, iij-1, iij-1 )
               ! round brackets added to fix the order of floating point operations
               ! needed to ensure halo 1 - halo 2 compatibility
               pt_rhs(ji,jj,jk,jn) = pt_rhs(ji,jj,jk,jn)                                                &
                  &                       + zsign * ( ( zftu(ji-1,jj  ,jk) - zftu(ji,jj,jk)             &
                  &                                   )                                                 & ! bracket for halo 1 - halo 2 compatibility
                  &                                 + ( zftv(ji,jj-1,jk) - zftv(ji,jj,jk)               &
                  &                                   )                                                 & ! bracket for halo 1 - halo 2 compatibility
                  &                                 ) / (  e1e2t(ji,jj) * e3t(ji,jj,jk,Kmm)  )
            END_2D
            !
         END DO
         !
         !                                !==  add the vertical 33 flux  ==!
         IF( ln_traldf_lap ) THEN               ! laplacian case: eddy coef = ah_wslp2 - akz
            DO_3D( iij-1, iij-1, iij-1, iij-1, 2, jpkm1 )
               ztfw(ji,jj,jk) = ztfw(ji,jj,jk) - e1e2t(ji,jj) / e3w(ji,jj,jk,Kmm) * tmask(ji,jj,jk)   &
                  &                            * ( ah_wslp2(ji,jj,jk) - akz(ji,jj,jk) )             &
                  &                            * (  pt(ji,jj,jk-1,jn) - pt(ji,jj,jk,jn) )
            END_3D
         ELSE                                   ! bilaplacian
            SELECT CASE( kpass )
            CASE(  1  )                            ! 1st pass : eddy coef = ah_wslp2
               DO_3D( iij-1, iij-1, iij-1, iij-1, 2, jpkm1 )
                  ztfw(ji,jj,jk) = ztfw(ji,jj,jk) - e1e2t(ji,jj) / e3w(ji,jj,jk,Kmm) * tmask(ji,jj,jk)             &
                     &                            * ah_wslp2(ji,jj,jk) * ( pt(ji,jj,jk-1,jn) - pt(ji,jj,jk,jn) )
               END_3D
            CASE(  2  )                            ! 2nd pass : eddy flux = ah_wslp2 and akz applied on pt  and pt2 gradients, resp.
               DO_3D( 0, 0, 0, 0, 2, jpkm1 )
                  ztfw(ji,jj,jk) = ztfw(ji,jj,jk) - e1e2t(ji,jj) / e3w(ji,jj,jk,Kmm) * tmask(ji,jj,jk)                      &
                     &                            * (  ah_wslp2(ji,jj,jk) * ( pt (ji,jj,jk-1,jn) - pt (ji,jj,jk,jn) )   &
                     &                               + akz     (ji,jj,jk) * ( pt2(ji,jj,jk-1,jn) - pt2(ji,jj,jk,jn) )   )
               END_3D
            END SELECT
         ENDIF
         !
         DO_3D( iij-1, iij-1, iij-1, iij-1, 1, jpkm1 )      !==  Divergence of vertical fluxes added to pta  ==!
            pt_rhs(ji,jj,jk,jn) = pt_rhs(ji,jj,jk,jn)    &
            &                                  + zsign * (  ztfw(ji,jj,jk+1) - ztfw(ji,jj,jk)  )   &
               &                                              / ( e1e2t(ji,jj) * e3t(ji,jj,jk,Kmm) )
         END_3D
         !
         IF( ( kpass == 1 .AND. ln_traldf_lap ) .OR.  &     !==  first pass only (  laplacian)  ==!
             ( kpass == 2 .AND. ln_traldf_blp ) ) THEN      !==  2nd   pass      (bilaplacian)  ==!
            !
            !                          ! "Poleward" diffusive heat or salt transports (T-S case only)
            IF( l_ptr )  CALL dia_ptr_hst( jn, 'ldf', zftv(:,:,:)  )
            !                          ! Diffusive heat transports
            IF( l_hst )  CALL dia_ar5_hst( jn, 'ldf', zftu(:,:,:), zftv(:,:,:) )
            !
         ENDIF                                                    !== end pass selection  ==!
         !
         !                                                        ! ===============
      END DO                                                      ! end tracer loop
      !                                                           ! ===============
   END SUBROUTINE tra_ldf_triad_t

   !!==============================================================================
END MODULE traldf_triad
