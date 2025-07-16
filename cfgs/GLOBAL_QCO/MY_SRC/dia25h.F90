MODULE dia25h 
   !!======================================================================
   !!                       ***  MODULE  diaharm  ***
   !! Harmonic analysis of tidal constituents 
   !!======================================================================
   !! History :  3.6  !  2014  (E O'Dea)  Original code
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain
   USE zdf_oce         ! ocean vertical physics    
   USE zdfgls   , ONLY : hmxl_n
   !
   USE in_out_manager  ! I/O units
   USE iom             ! I/0 library
   USE wet_dry

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dia_25h_init               ! routine called by nemogcm.F90
   PUBLIC   dia_25h                    ! routine called by diawri.F90

   LOGICAL, PUBLIC ::   ln_dia25h      !:  25h mean output

   ! variables for calculating 25-hourly means
   INTEGER , SAVE ::   cnt_25h           ! Counter for 25 hour means
   REAL(wp), SAVE ::   r1_25 = 0.04_wp   ! =1/25 
   REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::   tn_25h  , sn_25h
   REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:,:)   ::   sshn_25h 
   REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::   un_25h  , vn_25h  , wn_25h
   REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::   avt_25h , avm_25h
   REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::   en_25h  , rmxln_25h

!! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: dia25h.F90 15249 2021-09-13 09:59:09Z hadcv $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dia_25h_init( Kbb )
      !!---------------------------------------------------------------------------
      !!                  ***  ROUTINE dia_25h_init  ***
      !!     
      !! ** Purpose: Initialization of 25h mean namelist 
      !!        
      !! ** Method : Read namelist
      !!---------------------------------------------------------------------------
      INTEGER, INTENT(in) :: Kbb       ! Time level index
      !
      INTEGER ::   ios                 ! Local integer output status for namelist read
      INTEGER ::   ierror              ! Local integer for memory allocation
      INTEGER ::   ji, jj, jk
      !
      NAMELIST/nam_dia25h/ ln_dia25h
      !!----------------------------------------------------------------------
      !
      READ   ( numnam_ref, nam_dia25h, IOSTAT=ios, ERR= 901 )
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nam_dia25h in reference namelist' )
      READ  ( numnam_cfg, nam_dia25h, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'nam_dia25h in configuration namelist' )
      IF(lwm) WRITE ( numond, nam_dia25h )

      IF(lwp) THEN                   ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'dia_25h_init : Output 25 hour mean diagnostics'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist nam_dia25h : set 25h outputs '
         WRITE(numout,*) '      Switch for 25h diagnostics (T) or not (F)  ln_dia25h  = ', ln_dia25h
      ENDIF
      IF( .NOT. ln_dia25h )   RETURN
      ! ------------------- !
      ! 1 - Allocate memory !
      ! ------------------- !
      !                                ! ocean arrays
      ALLOCATE( tn_25h (A2D(0),jpk), sn_25h (A2D(0),jpk), sshn_25h(A2D(0))  ,     &
         &      un_25h (A2D(0),jpk), vn_25h (A2D(0),jpk), wn_25h(A2D(0),jpk),     &
         &      avt_25h(A2D(0),jpk), avm_25h(A2D(0),jpk),                      STAT=ierror )
      IF( ierror > 0 ) THEN
         CALL ctl_stop( 'dia_25h: unable to allocate ocean arrays' )   ;   RETURN
      ENDIF
      IF( ln_zdftke ) THEN             ! TKE physics
         ALLOCATE( en_25h(A2D(0),jpk), STAT=ierror )
         IF( ierror > 0 ) THEN
            CALL ctl_stop( 'dia_25h: unable to allocate en_25h' )   ;   RETURN
         ENDIF
      ENDIF
      IF( ln_zdfgls ) THEN             ! GLS physics
         ALLOCATE( en_25h(A2D(0),jpk), rmxln_25h(A2D(0),jpk), STAT=ierror )
         IF( ierror > 0 ) THEN
            CALL ctl_stop( 'dia_25h: unable to allocate en_25h and rmxln_25h' )   ;   RETURN
         ENDIF
      ENDIF
      ! ------------------------- !
      ! 2 - Assign Initial Values !
      ! ------------------------- !
      cnt_25h = 1  ! sets the first value of sum at timestep 1 (note - should strictly be at timestep zero so before values used where possible)
      DO_3D( 0, 0, 0, 0, 1, jpk )
         tn_25h (ji,jj,jk) = ts (ji,jj,jk,jp_tem,Kbb)
         sn_25h (ji,jj,jk) = ts (ji,jj,jk,jp_sal,Kbb)
         un_25h (ji,jj,jk) = uu (ji,jj,jk,Kbb)
         vn_25h (ji,jj,jk) = vv (ji,jj,jk,Kbb)
         avt_25h(ji,jj,jk) = avt(ji,jj,jk)
         avm_25h(ji,jj,jk) = avm(ji,jj,jk)
      END_3D
      DO_2D( 0, 0, 0, 0 )
         sshn_25h(ji,jj) = ssh(ji,jj,Kbb)
      END_2D
      IF( ln_zdftke ) THEN
         DO_3D( 0, 0, 0, 0, 1, jpk )
            en_25h(ji,jj,jk) = en(ji,jj,jk)
         END_3D
      ENDIF
      IF( ln_zdfgls ) THEN
         DO_3D( 0, 0, 0, 0, 1, jpk )
            en_25h   (ji,jj,jk) = en    (ji,jj,jk)
            rmxln_25h(ji,jj,jk) = hmxl_n(ji,jj,jk)
         END_3D
      ENDIF
!#if defined key_si3
!      CALL ctl_stop('STOP', 'dia_25h not setup yet to do tidemean ice')
!#endif 
      !
   END SUBROUTINE dia_25h_init


   SUBROUTINE dia_25h( kt, Kmm )  
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE dia_25h  ***
      !!         
      !! ** Purpose :   Write diagnostics with M2/S2 tide removed
      !!
      !! ** Method  :   25hr mean outputs for shelf seas
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      INTEGER, INTENT(in) ::   Kmm  ! ocean time level index
      !!
      INTEGER ::   ji, jj, jk
      INTEGER                          ::   iyear0, nimonth0,iday0            ! start year,imonth,day
      LOGICAL ::   ll_print = .FALSE.    ! =T print and flush numout
      REAL(dp)                         ::   zsto, zout, zmax, zjulian, zmdi   ! local scalars
      INTEGER                          ::   i_steps                           ! no of timesteps per hour
      REAL(wp), DIMENSION(A2D(0)    )  ::   zw2d, un_dm, vn_dm                ! workspace
      REAL(wp), DIMENSION(A2D(0),jpk)  ::   zw3d                              ! workspace
      REAL(wp), DIMENSION(A2D(0),3)    ::   zwtmb                             ! workspace
      !!----------------------------------------------------------------------

      ! 0. Initialisation
      ! -----------------
      ! Define frequency of summing to create 25 h mean
      IF( MOD( 3600,NINT(rn_Dt) ) == 0 ) THEN
         i_steps = 3600/NINT(rn_Dt)
      ELSE
         CALL ctl_stop('STOP', 'dia_wri_tide: timestep must give MOD(3600,rn_Dt) = 0 otherwise no hourly values are possible')
      ENDIF

      ! local variable for debugging
      ll_print = ll_print .AND. lwp

      ! wn_25h could not be initialised in dia_25h_init, so we do it here instead
      IF( kt == nn_it000 ) THEN
         DO_3D( 0, 0, 0, 0, 1, jpk )
            wn_25h(ji,jj,jk) = ww(ji,jj,jk)
         END_3D
      ENDIF

      ! Sum of 25 hourly instantaneous values to give a 25h mean from 24hours every day
      IF( MOD( kt, i_steps ) == 0  .AND. kt /= nn_it000 ) THEN

         IF (lwp) THEN
              WRITE(numout,*) 'dia_wri_tide : Summing instantaneous hourly diagnostics at timestep ',kt
              WRITE(numout,*) '~~~~~~~~~~~~ '
         ENDIF

         DO_3D( 0, 0, 0, 0, 1, jpk )
            tn_25h  (ji,jj,jk) = tn_25h  (ji,jj,jk) + ts (ji,jj,jk,jp_tem,Kmm)
            sn_25h  (ji,jj,jk) = sn_25h  (ji,jj,jk) + ts (ji,jj,jk,jp_sal,Kmm)
            un_25h  (ji,jj,jk) = un_25h  (ji,jj,jk) + uu (ji,jj,jk,Kmm)
            vn_25h  (ji,jj,jk) = vn_25h  (ji,jj,jk) + vv (ji,jj,jk,Kmm)
            wn_25h  (ji,jj,jk) = wn_25h  (ji,jj,jk) + ww (ji,jj,jk)
            avt_25h (ji,jj,jk) = avt_25h (ji,jj,jk) + avt(ji,jj,jk)
            avm_25h (ji,jj,jk) = avm_25h (ji,jj,jk) + avm(ji,jj,jk)
         END_3D
         DO_2D( 0, 0, 0, 0 )
            sshn_25h(ji,jj)    = sshn_25h(ji,jj)    + ssh(ji,jj,Kmm)
         END_2D
         IF( ln_zdftke ) THEN
            DO_3D( 0, 0, 0, 0, 1, jpk )
               en_25h(ji,jj,jk) = en_25h(ji,jj,jk) + en(ji,jj,jk)
            END_3D
         ENDIF
         IF( ln_zdfgls ) THEN
            DO_3D( 0, 0, 0, 0, 1, jpk )
               en_25h   (ji,jj,jk) = en_25h   (ji,jj,jk) + en    (ji,jj,jk)
               rmxln_25h(ji,jj,jk) = rmxln_25h(ji,jj,jk) + hmxl_n(ji,jj,jk)
            END_3D
         ENDIF
         cnt_25h = cnt_25h + 1
         !
         IF (lwp) THEN
            WRITE(numout,*) 'dia_tide : Summed the following number of hourly values so far',cnt_25h
         ENDIF
         !
      ENDIF ! MOD( kt, i_steps ) == 0

      ! Write data for 25 hour mean output streams
      IF( cnt_25h == 25 .AND.  MOD( kt, i_steps*24) == 0 .AND. kt /= nn_it000 ) THEN
         !
         IF(lwp) THEN
            WRITE(numout,*) 'dia_wri_tide : Writing 25 hour mean tide diagnostics at timestep', kt
            WRITE(numout,*) '~~~~~~~~~~~~ '
         ENDIF
         !
         tn_25h  (:,:,:) = tn_25h  (:,:,:) * r1_25
         sn_25h  (:,:,:) = sn_25h  (:,:,:) * r1_25
         sshn_25h(:,:)   = sshn_25h(:,:)   * r1_25
         un_25h  (:,:,:) = un_25h  (:,:,:) * r1_25
         vn_25h  (:,:,:) = vn_25h  (:,:,:) * r1_25
         wn_25h  (:,:,:) = wn_25h  (:,:,:) * r1_25
         avt_25h (:,:,:) = avt_25h (:,:,:) * r1_25
         avm_25h (:,:,:) = avm_25h (:,:,:) * r1_25
         IF( ln_zdftke ) THEN
            en_25h(:,:,:) = en_25h(:,:,:) * r1_25
         ENDIF
         IF( ln_zdfgls ) THEN
            en_25h   (:,:,:) = en_25h   (:,:,:) * r1_25
            rmxln_25h(:,:,:) = rmxln_25h(:,:,:) * r1_25
         ENDIF
         !
         IF(lwp)  WRITE(numout,*) 'dia_wri_tide : Mean calculated by dividing 25 hour sums and writing output'
         zmdi=1.e+20 !missing data indicator for masking
         ! write tracers (instantaneous)
         DO_3D( 0, 0, 0, 0, 1, jpk )
            zw3d(ji,jj,jk) = tn_25h(ji,jj,jk)*tmask(ji,jj,jk) + zmdi*(1.0-tmask(ji,jj,jk))
         END_3D
         CALL iom_put("temper25h", zw3d)   ! potential temperature
         DO_3D( 0, 0, 0, 0, 1, jpk )
            zw3d(ji,jj,jk) = sn_25h(ji,jj,jk)*tmask(ji,jj,jk) + zmdi*(1.0-tmask(ji,jj,jk))
         END_3D
         CALL iom_put( "salin25h", zw3d  )   ! salinity
         DO_2D( 0, 0, 0, 0 )
            zw2d(ji,jj) = sshn_25h(ji,jj)*tmask(ji,jj,1) + zmdi*(1.0-tmask(ji,jj,1))
         END_2D
         IF( ll_wd ) THEN
            CALL iom_put( "ssh25h", zw2d+ssh_ref )   ! sea surface 
         ELSE
            CALL iom_put( "ssh25h", zw2d )   ! sea surface
         ENDIF
         ! Write velocities (instantaneous)
         DO_3D( 0, 0, 0, 0, 1, jpk )
            zw3d(ji,jj,jk) = un_25h(ji,jj,jk)*umask(ji,jj,jk) + zmdi*(1.0-umask(ji,jj,jk))
         END_3D
         CALL iom_put("vozocrtx25h", zw3d)    ! i-current
         DO_3D( 0, 0, 0, 0, 1, jpk )
            zw3d(ji,jj,jk) = vn_25h(ji,jj,jk)*vmask(ji,jj,jk) + zmdi*(1.0-vmask(ji,jj,jk))
         END_3D
         CALL iom_put("vomecrty25h", zw3d  )   ! j-current
         DO_3D( 0, 0, 0, 0, 1, jpk )
            zw3d(ji,jj,jk) = wn_25h(ji,jj,jk)*wmask(ji,jj,jk) + zmdi*(1.0-tmask(ji,jj,jk))
         END_3D
         CALL iom_put("vovecrtz25h", zw3d )   ! k-current
         ! Write vertical physics
         DO_3D( 0, 0, 0, 0, 1, jpk )
            zw3d(ji,jj,jk) = avt_25h(ji,jj,jk)*wmask(ji,jj,jk) + zmdi*(1.0-tmask(ji,jj,jk))
         END_3D
         CALL iom_put("avt25h", zw3d )   ! diffusivity
         DO_3D( 0, 0, 0, 0, 1, jpk )
            zw3d(ji,jj,jk) = avm_25h(ji,jj,jk)*wmask(ji,jj,jk) + zmdi*(1.0-tmask(ji,jj,jk))
         END_3D
         CALL iom_put("avm25h", zw3d)   ! viscosity
         IF( ln_zdftke ) THEN
            DO_3D( 0, 0, 0, 0, 1, jpk )
               zw3d(ji,jj,jk) = en_25h(ji,jj,jk)*wmask(ji,jj,jk) + zmdi*(1.0-tmask(ji,jj,jk))
            END_3D
            CALL iom_put("tke25h", zw3d)   ! tke
         ENDIF
         IF( ln_zdfgls ) THEN
            DO_3D( 0, 0, 0, 0, 1, jpk )
               zw3d(ji,jj,jk) = en_25h(ji,jj,jk)*wmask(ji,jj,jk) + zmdi*(1.0-tmask(ji,jj,jk))
            END_3D
            CALL iom_put("tke25h", zw3d)   ! tke
            DO_3D( 0, 0, 0, 0, 1, jpk )
               zw3d(ji,jj,jk) = rmxln_25h(ji,jj,jk)*wmask(ji,jj,jk) + zmdi*(1.0-tmask(ji,jj,jk))
            END_3D
            CALL iom_put( "mxln25h",zw3d)
         ENDIF
         !
         ! After the write reset the values to cnt=1 and sum values equal current value 
         DO_3D( 0, 0, 0, 0, 1, jpk )
            tn_25h  (ji,jj,jk) = ts (ji,jj,jk,jp_tem,Kmm)
            sn_25h  (ji,jj,jk) = ts (ji,jj,jk,jp_sal,Kmm)
            un_25h  (ji,jj,jk) = uu (ji,jj,jk,Kmm)
            vn_25h  (ji,jj,jk) = vv (ji,jj,jk,Kmm)
            wn_25h  (ji,jj,jk) = ww (ji,jj,jk)
            avt_25h (ji,jj,jk) = avt(ji,jj,jk)
            avm_25h (ji,jj,jk) = avm(ji,jj,jk)
         END_3D
         DO_2D( 0, 0, 0, 0 )
            sshn_25h(ji,jj)    = ssh(ji,jj,Kmm)
         END_2D
         IF( ln_zdftke ) THEN
            DO_3D( 0, 0, 0, 0, 1, jpk )
               en_25h(ji,jj,jk) = en(ji,jj,jk)
            END_3D
         ENDIF
         IF( ln_zdfgls ) THEN
            DO_3D( 0, 0, 0, 0, 1, jpk )
               en_25h   (ji,jj,jk) = en    (ji,jj,jk)
               rmxln_25h(ji,jj,jk) = hmxl_n(ji,jj,jk)
            END_3D
         ENDIF
         cnt_25h = 1
         IF(lwp)  WRITE(numout,*) 'dia_wri_tide :   &
            &    After 25hr mean write, reset sum to current value and cnt_25h to one for overlapping average', cnt_25h
      ENDIF !  cnt_25h .EQ. 25 .AND.  MOD( kt, i_steps * 24) == 0 .AND. kt .NE. nn_it000
      !
   END SUBROUTINE dia_25h 

   !!======================================================================
END MODULE dia25h
