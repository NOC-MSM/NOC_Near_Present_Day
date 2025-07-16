MODULE diapea
   !!======================================================================
   !!                       ***  MODULE  diapea  ***
   !! Potential Energy Anomaly 
   !!======================================================================
   !! History :  3.6  !  12/2016  (J Tinker)  Original code
   !! 12/2024 update to calculate peas slwa
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain
   USE in_out_manager  ! I/O units
   USE iom             ! I/0 library
   USE eosbn2          ! Equation of state - in situ and potential density
   USE phycst          ! physical constant
   USE wet_dry

   IMPLICIT NONE
   PRIVATE 
   
   PUBLIC   dia_pea_init            ! routine called by nemogcm.F90
   PUBLIC   dia_pea                 ! routine called by diawri.F90
   REAL(wp), PUBLIC, SAVE, ALLOCATABLE,   DIMENSION(:,:)  ::   pea,peat,peas   
   REAL(wp), SAVE, ALLOCATABLE,   DIMENSION(:,:,:)  ::   wgt_co_mat   ! Weighting array for proportion of grid shallower than cut off depth
   REAL(wp), SAVE, ALLOCATABLE,   DIMENSION(:,:) ::   t_zmean, s_zmean  !Depth mean temperature and salinity: 2d fields
   REAL(wp), SAVE, ALLOCATABLE,   DIMENSION(:,:,:) ::   t_zmean_mat, s_zmean_mat  !Depth mean temperature and salinity: 3d fields
   REAL(wp) ::   zcutoff
   LOGICAL , PUBLIC ::   ln_pea  ! region mean calculation
   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.6 , NEMO Consortium (2014)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dia_pea_init(Kmm)
      ! Local variables
      INTEGER, INTENT(in) :: Kmm       ! Time level index
      INTEGER :: ji,jj,jk  ! Dummy loop indices
      REAL(wp) :: sumz,tmpsumz
      INTEGER  :: ierr                ! error integer for IOM_get
      INTEGER ::   ios                  ! Local integer output status for namelist read
      NAMELIST/nam_pea/ ln_pea
      
      zcutoff = 200.!200m
      
      
      
      !
      !REWIND ( numnam_ref )              ! Read Namelist nam_diatmb in referdiatmbence namelist : TMB diagnostics
      READ   ( numnam_ref, nam_pea, IOSTAT=ios, ERR= 901 )
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_pea in reference namelist' )

      !REWIND( numnam_cfg )              ! Namelist nam_diatmb in configuration namelist  TMB diagnostics
      READ  ( numnam_cfg, nam_pea, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_pea in configuration namelist' )
      IF(lwm) WRITE ( numond, nam_pea )

      IF(lwp) THEN                   ! Control print
          WRITE(numout,*)
          WRITE(numout,*) 'dia_pea_init : Output potential energy anomaly Diagnostics'
          WRITE(numout,*) '~~~~~~~~~~~~'
          WRITE(numout,*) 'Namelist nam_pea : set pea output '
          WRITE(numout,*) 'Switch for pea diagnostics (T) or not (F)  ln_dia  = ', ln_pea
      ENDIF
      
      
   ALLOCATE( pea(A2D(nn_hls)),  STAT= ierr )
   IF( ierr /= 0 )   CALL ctl_stop( 'dia_pea_init: failed to allocate pea array' )
   pea = 0.0_wp
   
   ALLOCATE( peat(A2D(nn_hls)),  STAT= ierr )
   IF( ierr /= 0 )   CALL ctl_stop( 'dia_pea_init: failed to allocate peat array' )
   peat = 0.0_wp
   
   ALLOCATE( peas(A2D(nn_hls)),  STAT= ierr )
   IF( ierr /= 0 )   CALL ctl_stop( 'dia_pea_init: failed to allocate peas array' )   
   peas = 0.0_wp
      
   ALLOCATE( t_zmean_mat(A2D(nn_hls),jpk),  STAT= ierr )
   IF( ierr /= 0 )   CALL ctl_stop( 'dia_pea_init: failed to allocate t_zmean_mat array' )
   t_zmean_mat = 0.0_wp
   
   ALLOCATE( s_zmean_mat(A2D(nn_hls),jpk),  STAT= ierr )
   IF( ierr /= 0 )   CALL ctl_stop( 'dia_pea_init: failed to allocate s_zmean_mat array' )
   s_zmean_mat = 0.0_wp
   
   ALLOCATE( t_zmean(A2D(nn_hls)),  STAT= ierr )
   IF( ierr /= 0 )   CALL ctl_stop( 'dia_pea_init: failed to allocate t_zmean array' )
   t_zmean = 0.0_wp

   ALLOCATE( s_zmean(A2D(nn_hls)),  STAT= ierr )
   IF( ierr /= 0 )   CALL ctl_stop( 'dia_pea_init: failed to allocate s_zmean array' )
   s_zmean = 0.0_wp
   
   ALLOCATE( wgt_co_mat(A2D(nn_hls),jpk),  STAT= ierr )
   IF( ierr /= 0 )   CALL ctl_stop( 'dia_pea_init: failed to allocate wgt_co_mat array' )
   wgt_co_mat = 0.0_wp
   
   if ( ln_pea ) THEN
          
      ! create wgt_co_mat mat, with the proportion of the grid (gdept_0) below cut off (200m)
      
      DO_2D( 0, 0, 0, 0 )
              IF ( tmask(ji,jj,1) == 1.0_wp ) THEN
                sumz = 0.
                DO jk = 1,jpk
                  IF ( tmask(ji,jj,jk) == 1.0_wp ) THEN
                    tmpsumz = sumz + e3t(ji,jj,jk,Kmm)
                    IF (sumz .ge. zcutoff) THEN
                      ! Already too deep
                      wgt_co_mat(ji,jj,jk) = 0.
                    ELSE IF (tmpsumz .le. zcutoff) THEN
                      ! Too shallow
                      wgt_co_mat(ji,jj,jk) = 1.
                    ELSE
                      !proprotion of grid box above cut off depth
                      wgt_co_mat(ji,jj,jk) = (zcutoff-Sumz)/e3t(ji,jj,jk,Kmm)
                    END IF
                    sumz = tmpsumz
                  endif
                END DO
                
              ELSE
                !if land, set to 0.
                DO jk = 1,jpk
                  wgt_co_mat(ji,jj,jk) = 0.
                END DO
                  
              ENDIF
       END_2D
   ENDIF
   

   END SUBROUTINE dia_pea_init
   

   SUBROUTINE dia_pea(kt, Kmm)
   INTEGER, INTENT(in) ::   kt      ! ocean time-step index
   INTEGER, INTENT(in) ::   Kmm  ! ocean time level index

    
   INTEGER :: ji,jj,jk,ierr  ! Dummy loop indices
   REAL(wp) :: tmpdenom, tmpnum, maxz
   !rau0
   
   REAL(wp), ALLOCATABLE,   DIMENSION(:,:,:,:) ::   ts_pea_mat,ts_pea_mat_TS_mean,ts_pea_mat_S_mean,ts_pea_mat_T_mean  !slwa
   REAL(wp), ALLOCATABLE,   DIMENSION(:,:,:) ::   tmp_pea_rho,tmp_pea_TS_mean_rho,tmp_pea_S_mean_rho,tmp_pea_T_mean_rho  !slwa
   REAL(wp), ALLOCATABLE,   DIMENSION(:,:,:) ::   zgdept
   REAL(wp) ::   int_y_pea,int_y_pea_t,int_y_pea_s  ! slwa
   
   
      
   ALLOCATE( ts_pea_mat(A2D(nn_hls),jpk,2),  STAT= ierr )
   IF( ierr /= 0 )   CALL ctl_stop( 'dia_pea_init: failed to allocate ts_pea_mat array' )
   ts_pea_mat = 0.0_wp    
   ALLOCATE( ts_pea_mat_TS_mean(A2D(nn_hls),jpk,2),  STAT= ierr )
   IF( ierr /= 0 )   CALL ctl_stop( 'dia_pea_init: failed to allocate ts_pea_mat_TS_mean array' )   
   ts_pea_mat_TS_mean = 0.0_wp    
   ALLOCATE( ts_pea_mat_S_mean(A2D(nn_hls),jpk,2),  STAT= ierr )
   IF( ierr /= 0 )   CALL ctl_stop( 'dia_pea_init: failed to allocate ts_pea_mat_S_mean array' )
   ts_pea_mat_S_mean = 0.0_wp    
   ALLOCATE( tmp_pea_rho(A2D(nn_hls),jpk),  STAT= ierr )
   IF( ierr /= 0 )   CALL ctl_stop( 'dia_pea_init: failed to allocate tmp_pea_rho array' )
   tmp_pea_rho = 0.0_wp
   ALLOCATE( tmp_pea_TS_mean_rho(A2D(nn_hls),jpk),  STAT= ierr )
   IF( ierr /= 0 )   CALL ctl_stop( 'dia_pea_init: failed to allocate tmp_pea_TS_mean_rho array' )
   tmp_pea_TS_mean_rho = 0.0_wp
   ALLOCATE( tmp_pea_S_mean_rho(A2D(nn_hls),jpk),  STAT= ierr )
   IF( ierr /= 0 )   CALL ctl_stop( 'dia_pea_init: failed to allocate tmp_pea_S_mean_rho array' )
   tmp_pea_S_mean_rho = 0.0

   !slwa
   ALLOCATE( ts_pea_mat_T_mean(A2D(nn_hls),jpk,2),  STAT= ierr )
   IF( ierr /= 0 )   CALL ctl_stop( 'dia_pea_init: failed to allocate ts_pea_mat_T_mean array' )
   ts_pea_mat_T_mean = 0.0_wp
   ALLOCATE( tmp_pea_T_mean_rho(A2D(nn_hls),jpk),  STAT= ierr )
   IF( ierr /= 0 )   CALL ctl_stop( 'dia_pea_init: failed to allocate tmp_pea_T_mean_rho array' )
   tmp_pea_T_mean_rho = 0.0
   !slwa

    pea(:,:)=0.0d0
    peaT(:,:)=0.0d0
    peaS(:,:)=0.0d0
    
    !calculate the depth mean temperature and salinity of the upper 200m. Save this into a 3d array. Set the value where tmask=0 to be tsn. 
    
       DO_2D( 0, 0, 0, 0 )
            IF ( tmask(ji,jj,1) == 1.0_wp ) THEN  ! if a sea point. 
            
             !Depth mean temperature
             tmpdenom = 0.
             tmpnum = 0.
             DO jk = 1,jpk
                IF ( tmask(ji,jj,jk) == 1.0_wp ) THEN
                  tmpnum = tmpnum + (wgt_co_mat(ji,jj,jk)*e3t(ji,jj,jk,Kmm)*ts (ji,jj,jk,jp_tem,Kmm))
                  tmpdenom = tmpdenom + (wgt_co_mat(ji,jj,jk)*e3t(ji,jj,jk,Kmm))
                endif
             END DO
             t_zmean(ji,jj) = tmpnum/tmpdenom
             
             !Depth mean salinity
             tmpdenom = 0.
             tmpnum = 0.
             DO jk = 1,jpk
                IF ( tmask(ji,jj,jk) == 1.0_wp ) THEN
                  tmpnum     = tmpnum + (wgt_co_mat(ji,jj,jk)*e3t(ji,jj,jk,Kmm)*ts (ji,jj,jk,jp_sal,Kmm))
                  tmpdenom = tmpdenom + (wgt_co_mat(ji,jj,jk)*e3t(ji,jj,jk,Kmm))
                endif
             END DO
             s_zmean(ji,jj) = tmpnum/tmpdenom
             
             !save into a 3d grid
             DO jk = 1,jpk
               IF ( tmask(ji,jj,jk) == 1.0_wp ) THEN
                 t_zmean_mat(ji,jj,jk) = t_zmean(ji,jj)
                 s_zmean_mat(ji,jj,jk) = s_zmean(ji,jj)
               else
                 t_zmean_mat(ji,jj,jk) = ts (ji,jj,jk,jp_tem,Kmm)
                 s_zmean_mat(ji,jj,jk) = ts (ji,jj,jk,jp_sal,Kmm)
               endif
             END DO
            else
              t_zmean(ji,jj) = ts (ji,jj,1,jp_tem,Kmm)
              s_zmean(ji,jj) = ts (ji,jj,1,jp_sal,Kmm)
              DO jk = 1,jpk
                t_zmean_mat(ji,jj,jk) = ts (ji,jj,jk,jp_tem,Kmm)
                s_zmean_mat(ji,jj,jk) = ts (ji,jj,jk,jp_sal,Kmm)
              END DO
            endif
            
      END_2D
      
   !Calculate the density from the depth varying, and depth average temperature and salinity
   !-----------------------------
   !-----------------------------
   
   ts_pea_mat(:,:,:,:) = ts (:,:,:,:,Kmm)
   
   ts_pea_mat_TS_mean(:,:,:,1) = t_zmean_mat(:,:,:)
   ts_pea_mat_TS_mean(:,:,:,2) = s_zmean_mat(:,:,:)
   
   ts_pea_mat_S_mean(:,:,:,1) = t_zmean_mat(:,:,:)
   ts_pea_mat_S_mean(:,:,:,2) = ts (:,:,:,jp_sal,Kmm)

   !slwa
   ts_pea_mat_T_mean(:,:,:,1) = ts (:,:,:,jp_tem,Kmm)
   ts_pea_mat_T_mean(:,:,:,2) = s_zmean_mat(:,:,:)
   !slwa

   ALLOCATE( zgdept(jpi,jpj,jpk) )

   DO jk = 1, jpk
      zgdept(:,:,jk) = gdept(:,:,jk,Kmm)
   END DO
   
!  CALL eos ( ts_pea_mat,         tmp_pea_rho,         gdept(:,:,:,Kmm) )
!  CALL eos ( ts_pea_mat_TS_mean, tmp_pea_TS_mean_rho, gdept(:,:,:,Kmm) )
!  CALL eos ( ts_pea_mat_S_mean,  tmp_pea_S_mean_rho,  gdept(:,:,:,Kmm) )
!  CALL eos ( ts_pea_mat_T_mean,  tmp_pea_T_mean_rho,  gdept(:,:,:,Kmm) )   ! slwa

   CALL eos ( ts_pea_mat,         tmp_pea_rho,         zgdept )
   CALL eos ( ts_pea_mat_TS_mean, tmp_pea_TS_mean_rho, zgdept )
   CALL eos ( ts_pea_mat_S_mean,  tmp_pea_S_mean_rho,  zgdept )
   CALL eos ( ts_pea_mat_T_mean,  tmp_pea_T_mean_rho,  zgdept )   ! slwa
   tmp_pea_rho = (tmp_pea_rho * rho0) + rho0
   tmp_pea_TS_mean_rho = (tmp_pea_TS_mean_rho * rho0) + rho0
   tmp_pea_S_mean_rho = (tmp_pea_S_mean_rho * rho0) + rho0
   tmp_pea_T_mean_rho = (tmp_pea_T_mean_rho * rho0) + rho0   !slwa
   
   
   ! to test the density calculation
   !CALL iom_put( "tmp_pea_rho" , tmp_pea_rho )                 ! pea
   !CALL iom_put( "tmp_pea_TS_mean_rho" , tmp_pea_TS_mean_rho )                 ! pea
   
   
   ! Caluclation of the PEA.
       DO_2D( 0, 0, 0, 0 )
        pea(ji,jj) = 0.
        peat(ji,jj) = 0.
        peas(ji,jj) = 0.
        maxz = 0.
        int_y_pea = 0.
        int_y_pea_t = 0.
        int_y_pea_s = 0.  ! slwa
        IF ( tmask(ji,jj,1) == 1.0_wp ) THEN ! for sea points
        
          ! the depth integrated calculation is summed up over the depths, and then divided by the depth
          DO jk = 1,jpk
            !for each level...
            
            IF ( tmask(ji,jj,jk) == 1.0_wp ) THEN ! if above the sea bed...
              int_y_pea =   -((tmp_pea_TS_mean_rho(ji,jj,jk)) - (tmp_pea_rho(ji,jj,jk)))*9.81*gdept(ji,jj,jk,Kmm)*wgt_co_mat(ji,jj,jk)*e3t(ji,jj,jk,Kmm)
              int_y_pea_t = -((tmp_pea_S_mean_rho(ji,jj,jk))  - (tmp_pea_rho(ji,jj,jk)))*9.81*gdept(ji,jj,jk,Kmm)*wgt_co_mat(ji,jj,jk)*e3t(ji,jj,jk,Kmm)
              int_y_pea_s = -((tmp_pea_T_mean_rho(ji,jj,jk))  - (tmp_pea_rho(ji,jj,jk)))*9.81*gdept(ji,jj,jk,Kmm)*wgt_co_mat(ji,jj,jk)*e3t(ji,jj,jk,Kmm) ! slwa
            else
              int_y_pea = 0.
              int_y_pea_t = 0.
              int_y_pea_s = 0.   !slwa
            endif
            
            ! check that the sum is not NaN. 
            if ( int_y_pea .ne.  int_y_pea    ) int_y_pea = 0.
            if ( int_y_pea_t .ne. int_y_pea_t )  int_y_pea_t = 0.
            if ( int_y_pea_s .ne. int_y_pea_s )  int_y_pea_s = 0.   ! slwa
            !if ( (int_y_pea*int_y_pea    ) .gt. 1.0e6 ) int_y_pea = 0.
            !if ( (int_y_pea_t*int_y_pea_t) .gt. 1.0e6 ) int_y_pea_t = 0.
            
            pea(ji,jj) =  pea(ji,jj) + int_y_pea
            peat(ji,jj) =  peat(ji,jj) + int_y_pea_t
            peas(ji,jj) =  peas(ji,jj) + int_y_pea_s   ! slwa
            maxz = maxz + (e3t(ji,jj,jk,Kmm)*wgt_co_mat(ji,jj,jk))
          enddo
            
            
          !divide by the depth
          pea(ji,jj) = pea(ji,jj)/maxz
          peat(ji,jj) = peat(ji,jj)/maxz
          peas(ji,jj) = peas(ji,jj)/maxz
!slwa     peas(ji,jj) = pea(ji,jj) - peat(ji,jj)
            
            
          else
            pea(ji,jj) = 0.
            peat(ji,jj) = 0.
            peas(ji,jj) = 0.
          endif
     END_2D
!    
     CALL iom_put( "pea" , pea )                 ! pea
     CALL iom_put( "peat" , peat )               ! pea
     CALL iom_put( "peas" , peas )               ! pea
           

    DEALLOCATE(ts_pea_mat,ts_pea_mat_TS_mean,ts_pea_mat_S_mean,tmp_pea_rho,tmp_pea_TS_mean_rho,tmp_pea_S_mean_rho)
    DEALLOCATE(ts_pea_mat_T_mean,tmp_pea_T_mean_rho)   ! slwa
    DEALLOCATE(zgdept)
   
   END SUBROUTINE dia_pea

END MODULE diapea
