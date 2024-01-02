MODULE diaprod
! Requires key_iom_put
# if defined key_xios
   !!======================================================================
   !!                     ***  MODULE  diaprod  ***
   !! Ocean diagnostics :  write ocean product diagnostics
   !!=====================================================================
   !! History :  3.4  ! 2012  (D. Storkey)  Original code
   !!            4.0  ! 2019  (D. Storkey)
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dia_prod      : calculate and write out product diagnostics
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers 
   USE dom_oce         ! ocean space and time domain
   USE domvvl          ! for thickness weighted diagnostics if key_vvl
   USE eosbn2          ! equation of state  (eos call)
   USE phycst          ! physical constants
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE in_out_manager  ! I/O manager
   USE iom
   USE ioipsl
   USE lib_mpp         ! MPP library
   USE timing          ! preformance summary

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dia_prod                 ! routines called by step.F90

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.4 , NEMO Consortium (2012)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dia_prod( kt, Kmm )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dia_prod  ***
      !!                   
      !! ** Purpose :   Write out product diagnostics (uT, vS etc.)
      !!
      !! ** Method  :  use iom_put
      !!               Product diagnostics are not thickness-weighted in this routine.
      !!               They should be thickness-weighted using XIOS if key_vvl is set. 
      !!----------------------------------------------------------------------
      !!
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      INTEGER, INTENT( in ) ::   Kmm     ! ocean time level index
      !!
      INTEGER                      ::   ji, jj, jk              ! dummy loop indices
      REAL(wp)                     ::   zztmp, zztmpx, zztmpy   ! 
      !!
      REAL(wp), POINTER, DIMENSION(:,:)             ::   z2d   ! 2D workspace
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   z3d   ! 3D workspace
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zrhop    ! potential density
      !!----------------------------------------------------------------------
      ! 
      IF( ln_timing )   CALL timing_start('dia_prod')
      ! 
      ALLOCATE( z2d(jpi,jpj), z3d(jpi,jpj,jpk), zrhop(jpi,jpj,jpk) )
      !

      IF( iom_use("urhop") .OR. iom_use("vrhop") .OR. iom_use("wrhop") &
#if ! defined key_diaar5
     &  .OR. iom_use("rhop") &
#endif
     & ) THEN
         CALL eos( ts(:,:,:,:,Kmm), z3d, zrhop )                 ! now in situ and potential density
         zrhop(:,:,:) = zrhop(:,:,:)-1000.e0         ! reference potential density to 1000 to avoid precision issues in rhop2 calculation
         zrhop(:,:,jpk) = 0._wp
#if ! defined key_diaar5
         CALL iom_put( 'rhop', zrhop )
#else
         ! If key_diaar5 set then there is already an iom_put call to output rhop.
         ! Really should be a standard diagnostics option?
#endif
      ENDIF

      IF( iom_use("ut") ) THEN
         z3d(:,:,:) = 0.e0 
         DO_3D( 0, 0, 0, 0, 1, jpk )
                  z3d(ji,jj,jk) = uu(ji,jj,jk,Kmm) * 0.5 * ( ts(ji,jj,jk,jp_tem,Kmm) + ts(ji+1,jj,jk,jp_tem,Kmm) )
         END_3D
         CALL iom_put( "ut", z3d )                  ! product of temperature and zonal velocity at U points
      ENDIF

      IF( iom_use("vt") ) THEN
         z3d(:,:,:) = 0.e0 
         DO_3D( 0, 0, 0, 0, 1, jpk )
                  z3d(ji,jj,jk) = vv(ji,jj,jk,Kmm) * 0.5 * ( ts(ji,jj,jk,jp_tem,Kmm) + ts(ji,jj+1,jk,jp_tem,Kmm) )
         END_3D
         CALL iom_put( "vt", z3d )                  ! product of temperature and meridional velocity at V points
      ENDIF

      IF( iom_use("wt") ) THEN
         z3d(:,:,:) = 0.e0 
         DO_2D( 0, 0, 0, 0 )
               z3d(ji,jj,1) = ww(ji,jj,1) * ts(ji,jj,1,jp_tem,Kmm)
         END_2D
         
         DO_3D( 0, 0, 0, 0, 1, jpk )
                  z3d(ji,jj,jk) = ww(ji,jj,jk) * 0.5 * ( ts(ji,jj,jk-1,jp_tem,Kmm) + ts(ji,jj,jk,jp_tem,Kmm) )
         END_3D
         CALL iom_put( "wt", z3d )                  ! product of temperature and vertical velocity at W points
      ENDIF

      IF( iom_use("us") ) THEN
         z3d(:,:,:) = 0.e0 
         DO_3D( 0, 0, 0, 0, 1, jpk )
                  z3d(ji,jj,jk) = uu(ji,jj,jk,Kmm) * 0.5 * ( ts(ji,jj,jk,jp_sal,Kmm) + ts(ji+1,jj,jk,jp_sal,Kmm) )
         END_3D
         CALL iom_put( "us", z3d )                  ! product of salinity and zonal velocity at U points
      ENDIF

      IF( iom_use("vs") ) THEN
         z3d(:,:,:) = 0.e0 
         DO_3D( 0, 0, 0, 0, 1, jpk )
                  z3d(ji,jj,jk) = vv(ji,jj,jk,Kmm) * 0.5 * ( ts(ji,jj,jk,jp_sal,Kmm) + ts(ji,jj+1,jk,jp_sal,Kmm) )
         END_3D
         CALL iom_put( "vs", z3d )                  ! product of salinity and meridional velocity at V points
      ENDIF

      IF( iom_use("ws") ) THEN
         z3d(:,:,:) = 0.e0 
         DO_2D( 0, 0, 0, 0 )
               z3d(ji,jj,1) = ww(ji,jj,1) * ts(ji,jj,1,jp_sal,Kmm)
         END_2D
         
         DO_3D( 0, 0, 0, 0, 1, jpk )
                  z3d(ji,jj,jk) = ww(ji,jj,jk) * 0.5 * ( ts(ji,jj,jk-1,jp_sal,Kmm) + ts(ji,jj,jk,jp_sal,Kmm) )
         END_3D
         
         CALL iom_put( "ws", z3d )                  ! product of salinity and vertical velocity at W points
      ENDIF

      IF( iom_use("uv") ) THEN
         z3d(:,:,:) = 0.e0 
         DO_3D( 0, 0, 0, 0, 1, jpk )
                  z3d(ji,jj,jk) = 0.25 * ( uu(ji-1,jj,jk,Kmm) + uu(ji,jj,jk,Kmm) ) * ( vv(ji,jj-1,jk,Kmm) + vv(ji,jj,jk,Kmm) ) 
         END_3D
         CALL iom_put( "uv", z3d )                  ! product of zonal velocity and meridional velocity at T points
      ENDIF

      IF( iom_use("uw") ) THEN
         z3d(:,:,:) = 0.e0 
         DO_2D( 0, 0, 0, 0 )
               z3d(ji,jj,1) = 0.5 * ( ww(ji,jj,1) + ww(ji+1,jj,1) ) * uu(ji,jj,1,Kmm) 
         END_2D
         
         DO_3D( 0, 0, 0, 0, 1, jpk )
                  z3d(ji,jj,jk) = 0.25 * ( ww(ji,jj,jk) + ww(ji+1,jj,jk) ) * ( uu(ji,jj,jk-1,Kmm) + uu(ji,jj,jk,Kmm) ) 
         END_3D
         CALL iom_put( "uw", z3d )                  ! product of zonal velocity and vertical velocity at UW points
      ENDIF

      IF( iom_use("vw") ) THEN
         z3d(:,:,:) = 0.e0 
         DO_2D( 0, 0, 0, 0 )
               z3d(ji,jj,1) = 0.5 * ( ww(ji,jj,1) + ww(ji,jj+1,1) ) * vv(ji,jj,1,Kmm) 
         END_2D
         
         DO_3D( 0, 0, 0, 0, 1, jpk )
                  z3d(ji,jj,jk) = 0.25 * ( ww(ji,jj,jk) + ww(ji,jj+1,jk) ) * ( vv(ji,jj,jk-1,Kmm) + vv(ji,jj,jk,Kmm) ) 
         END_3D
         CALL iom_put( "vw", z3d )                  ! product of meriodional velocity and vertical velocity at VW points
      ENDIF

      IF( iom_use("urhop") ) THEN
         z3d(:,:,:) = 0.e0 
         DO_3D( 0, 0, 0, 0, 1, jpk )
                  z3d(ji,jj,jk) = uu(ji,jj,jk,Kmm) * 0.5 * ( zrhop(ji,jj,jk) + zrhop(ji+1,jj,jk) )
         END_3D
         CALL iom_put( "urhop", z3d )                  ! product of density and zonal velocity at U points
      ENDIF

      IF( iom_use("vrhop") ) THEN
         z3d(:,:,:) = 0.e0 
         DO_3D( 0, 0, 0, 0, 1, jpk )
                  z3d(ji,jj,jk) = vv(ji,jj,jk,Kmm) * 0.5 * ( zrhop(ji,jj,jk) + zrhop(ji,jj+1,jk) )
         END_3D
         CALL iom_put( "vrhop", z3d )                  ! product of density and meridional velocity at V points
      ENDIF

      IF( iom_use("wrhop") ) THEN
         z3d(:,:,:) = 0.e0 
         DO_2D( 0, 0, 0, 0 )
               z3d(ji,jj,1) = ww(ji,jj,1) * zrhop(ji,jj,1)
         END_2D
         DO_3D( 0, 0, 0, 0, 1, jpk )
                  z3d(ji,jj,jk) = ww(ji,jj,jk) * 0.5 * ( zrhop(ji,jj,jk-1) + zrhop(ji,jj,jk) )
         END_3D
         CALL iom_put( "wrhop", z3d )                  ! product of density and vertical velocity at W points
      ENDIF

      !
      DEALLOCATE( z2d, z3d, zrhop )
      !
      IF( ln_timing )   CALL timing_stop('dia_prod')
      !
   END SUBROUTINE dia_prod
#else
   !!----------------------------------------------------------------------
   !!   Default option :                                         NO diaprod
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER :: lk_diaprod = .FALSE.   ! coupled flag
CONTAINS
   SUBROUTINE dia_prod( kt , Kmm)   ! Empty routine
      INTEGER ::   kt, Kmm
      !WRITE(*,*) 'dia_prod: You should not have seen this print! error?', kt
   END SUBROUTINE dia_prod
#endif
   !!======================================================================
END MODULE diaprod
