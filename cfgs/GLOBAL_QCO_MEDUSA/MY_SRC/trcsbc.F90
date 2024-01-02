MODULE trcsbc
   !!==============================================================================
   !!                       ***  MODULE  trcsbc  ***
   !! Ocean passive tracers:  surface boundary condition
   !!======================================================================
   !! History :  8.2  !  1998-10  (G. Madec, G. Roullet, M. Imbard)  Original code
   !!            8.2  !  2001-02  (D. Ludicone)  sea ice and free surface
   !!            8.5  !  2002-06  (G. Madec)  F90: Free form and module
   !!            9.0  !  2004-03  (C. Ethe)  adapted for passive tracers
   !!                 !  2006-08  (C. Deltel) Diagnose ML trends for passive tracers
   !!==============================================================================
#if defined key_top
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !!   trc_sbc      : update the tracer trend at ocean surface
   !!----------------------------------------------------------------------
   USE par_trc        ! need jptra, number of passive tracers
   USE oce_trc         ! ocean dynamics and active tracers variables
   USE trc             ! ocean  passive tracers variables
   USE prtctl          ! Print control for debbuging
   USE iom
   USE trd_oce
   USE trdtra

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sbc   ! routine called by step.F90

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcsbc.F90 15394 2021-10-18 10:55:29Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_sbc ( kt, Kmm, ptr, Krhs )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_sbc  ***
      !!                   
      !! ** Purpose :   Compute the tracer surface boundary condition trend of
      !!      (concentration/dilution effect) and add it to the general 
      !!       trend of tracer equations.
      !!
      !! ** Method :
      !!      * concentration/dilution effect:
      !!            The surface freshwater flux modify the ocean volume
      !!         and thus the concentration of a tracer as :
      !!            tr(Krhs) = tr(Krhs) + emp * tr(Kmm) / e3t_ + fmmflx * tri / e3t   for k=1
      !!          - tr(Kmm) , the concentration of tracer in the ocean
      !!          - tri, the concentration of tracer in the sea-ice
      !!          - emp, the surface freshwater budget (evaporation minus precipitation + fmmflx)
      !!            given in kg/m2/s is divided by 1035 kg/m3 (density of ocean water) to obtain m/s.
      !!          - fmmflx, the flux asscociated to freezing-melting of sea-ice 
      !!            In linear free surface case (ln_linssh=T), the volume of the
      !!            ocean does not change with the water exchanges at the (air+ice)-sea
      !!
      !! ** Action  : - Update the 1st level of tr(:,:,:,:,Krhs) with the trend associated
      !!                with the tracer surface boundary condition 
      !!
      !!----------------------------------------------------------------------
      INTEGER,                                    INTENT(in   ) :: kt        ! ocean time-step index
      INTEGER,                                    INTENT(in   ) :: Kmm, Krhs ! time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jptra,jpt), INTENT(inout) :: ptr       ! passive tracers and RHS of tracer equation
      !
      INTEGER  ::   ji, jj, jn                      ! dummy loop indices
      REAL(wp) ::   zse3t, zrtrn, zfact     ! local scalars
      REAL(wp) ::   zdtra          !   -      -
      CHARACTER (len=22) :: charout
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   ztrtrd
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('trc_sbc')
      !
      ! Allocate temporary workspace
      IF( l_trdtrc )  ALLOCATE( ztrtrd(jpi,jpj,jpk) )
      !
      zrtrn = 1.e-15_wp

      IF( kt == nittrc000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'trc_sbc : Passive tracers surface boundary condition'
         IF(lwp) WRITE(numout,*) '~~~~~~~ '
         !
         IF( ln_rsttr .AND. .NOT.ln_top_euler .AND.   &                     ! Restart: read in restart  file
            iom_varid( numrtr, 'sbc_'//TRIM(ctrcnm(1))//'_b', ldstop = .FALSE. ) > 0 ) THEN
            IF(lwp) WRITE(numout,*) '          nittrc000-1 surface tracer content forcing fields read in the restart file'
            zfact = 0.5_wp
            DO jn = 1, jptra
               IF( iom_varid( numrtr, 'sbc_'//TRIM(ctrcnm(jn))//'_b', ldstop = .FALSE. ) > 0 ) THEN
                  !! Jpalm -- ZMP -- need to read restart without oithona
                  !! ZMP is here :
                  CALL iom_get( numrtr, jpdom_auto, 'sbc_'//TRIM(ctrcnm(jn))//'_b', sbc_trc_b(:,:,jn) )   ! before tracer content sbc
               ELSE
                  !! Jpalm -- ZMP -- not sure what's a reasonable value... try with 0
                  IF(lwp) WRITE(numout,*) 'trc_sb ---> sbc_',TRIM(ctrcnm(jn)),'_b NOT in RESTART file '
                  IF(lwp) WRITE(numout,*) '         -- Initialize to 0.0'
                  sbc_trc_b(:,:,jn) = 0._wp
               ENDIF
            END DO
         ELSE                                         ! No restart or restart not found: Euler forward time stepping
           zfact = 1._wp
           sbc_trc_b(:,:,:) = 0._wp
         ENDIF
      ELSE                                         ! Swap of forcing fields
         IF( ln_top_euler ) THEN
            zfact = 1._wp
            sbc_trc_b(:,:,:) = 0._wp
         ELSE
            zfact = 0.5_wp
            sbc_trc_b(:,:,:) = sbc_trc(:,:,:)
         ENDIF
         !
      ENDIF

      ! 0. initialization
      SELECT CASE ( nn_ice_tr )

      CASE ( -1 ) ! ! No tracers in sea ice ( trc_i = 0 )
         !
         DO jn = 1, jptra
            DO_2D( 0, 0, 0, 1 )
               sbc_trc(ji,jj,jn) = 0._wp
            END_2D
         END DO
         !
         IF( ln_linssh ) THEN  !* linear free surface  
            DO jn = 1, jptra
               DO_2D( 0, 0, 0, 1 )
                  sbc_trc(ji,jj,jn) = sbc_trc(ji,jj,jn) + r1_rho0 * emp(ji,jj) * ptr(ji,jj,1,jn,Kmm) !==>> add concentration/dilution effect due to constant volume cell
               END_2D
            END DO
         ENDIF
         !
      CASE ( 0 )  ! Same concentration in sea ice and in the ocean ( trc_i = ptr(...,Kmm)  )
         !
         DO jn = 1, jptra
            DO_2D( 0, 0, 0, 1 )
               sbc_trc(ji,jj,jn) = - fmmflx(ji,jj) * r1_rho0 * ptr(ji,jj,1,jn,Kmm)
            END_2D
         END DO
         !
         IF( ln_linssh ) THEN  !* linear free surface  
            DO jn = 1, jptra
               DO_2D( 0, 0, 0, 1 )
                  sbc_trc(ji,jj,jn) = sbc_trc(ji,jj,jn) + r1_rho0 * emp(ji,jj) * ptr(ji,jj,1,jn,Kmm) !==>> add concentration/dilution effect due to constant volume cell
               END_2D
            END DO
         ENDIF
         !
      CASE ( 1 )  ! Specific treatment of sea ice fluxes with an imposed concentration in sea ice 
         !
         DO jn = 1, jptra
            DO_2D( 0, 0, 0, 1 )
               sbc_trc(ji,jj,jn) = - fmmflx(ji,jj) * r1_rho0 * trc_i(ji,jj,jn)
            END_2D
         END DO
         !
         IF( ln_linssh ) THEN  !* linear free surface  
            DO jn = 1, jptra
               DO_2D( 0, 0, 0, 1 )
                  sbc_trc(ji,jj,jn) = sbc_trc(ji,jj,jn) + r1_rho0 * emp(ji,jj) * ptr(ji,jj,1,jn,Kmm) !==>> add concentration/dilution effect due to constant volume cell
               END_2D
            END DO
         ENDIF
         !
         DO jn = 1, jptra
            DO_2D( 0, 0, 0, 1 )
               zse3t = rDt_trc / e3t(ji,jj,1,Kmm)
               zdtra = ptr(ji,jj,1,jn,Kmm) + sbc_trc(ji,jj,jn) * zse3t 
               IF( zdtra < 0. ) sbc_trc(ji,jj,jn) = MAX( zdtra, -ptr(ji,jj,1,jn,Kmm) / zse3t  ) ! avoid negative concentration that can occurs if trc_i > ptr 
            END_2D
         END DO
         !                             
      END SELECT
      !
      CALL lbc_lnk( 'trcsbc', sbc_trc(:,:,:), 'T', 1.0_wp )
      !                                       Concentration dilution effect on tracers due to evaporation & precipitation 
      DO jn = 1, jptra
         !
         IF( l_trdtrc )   ztrtrd(:,:,:) = ptr(:,:,:,jn,Krhs)  ! save trends
         !
         DO_2D( 0, 0, 0, 1 )
            zse3t = zfact / e3t(ji,jj,1,Kmm)
            ptr(ji,jj,1,jn,Krhs) = ptr(ji,jj,1,jn,Krhs) + ( sbc_trc_b(ji,jj,jn) + sbc_trc(ji,jj,jn) ) * zse3t
         END_2D
         !
         IF( l_trdtrc ) THEN
            ztrtrd(:,:,:) = ptr(:,:,:,jn,Krhs) - ztrtrd(:,:,:)
            CALL trd_tra( kt, Kmm, Krhs, 'TRC', jn, jptra_nsr, ztrtrd )
         END IF
         !                                                       ! ===========
      END DO                                                     ! tracer loop
      !                                                          ! ===========
      !
      !                                           Write in the tracer restar  file
      !                                          *******************************
      IF( lrst_trc .AND. .NOT.ln_top_euler ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'sbc : ocean surface tracer content forcing fields written in tracer restart file ',   &
            &                    'at it= ', kt,' date= ', ndastp
         IF(lwp) WRITE(numout,*) '~~~~'
         DO jn = 1, jptra
            CALL iom_rstput( kt, nitrst, numrtw, 'sbc_'//TRIM(ctrcnm(jn))//'_b', sbc_trc(:,:,jn) )
         END DO
      ENDIF
      !
      IF( sn_cfctl%l_prttrc )   THEN
         WRITE(charout, FMT="('sbc ')") ;  CALL prt_ctl_info( charout, cdcomp = 'top' )
                                           CALL prt_ctl( tab4d_1=ptr(:,:,:,:,Krhs), mask1=tmask, clinfo=ctrcnm, clinfo3='trd' )
      ENDIF
      IF( l_trdtrc )  DEALLOCATE( ztrtrd )
      !
      IF( ln_timing )   CALL timing_stop('trc_sbc')
      !
   END SUBROUTINE trc_sbc

#else
   !!----------------------------------------------------------------------
   !!   Dummy module :                      NO passive tracer
   !!----------------------------------------------------------------------
   USE par_oce
   USE par_trc
CONTAINS
   SUBROUTINE trc_sbc ( kt, Kmm, ptr, Krhs )      ! Empty routine
      INTEGER,                                    INTENT(in   ) :: kt        ! ocean time-step index
      INTEGER,                                    INTENT(in   ) :: Kmm, Krhs ! time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jptra,jpt), INTENT(inout) :: ptr       ! passive tracers and RHS of tracer equation
      WRITE(*,*) 'trc_sbc: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sbc
#endif
   
   !!======================================================================
END MODULE trcsbc
