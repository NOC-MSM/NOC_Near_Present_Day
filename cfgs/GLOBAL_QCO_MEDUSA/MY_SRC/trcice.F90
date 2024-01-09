MODULE trcice
   !!======================================================================
   !!                         ***  MODULE trcice  ***
   !! TOP :   Manage the communication between TOP and sea ice
   !!======================================================================
   !! History :  3.5  ! 2013    (M. Vancoppenolle, O. Aumont, G. Madec), original code
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !!   trc_ice       :  Call the appropriate sea ice tracer subroutine
   !!----------------------------------------------------------------------
   USE par_trc         ! need jptra, number of passive tracers
   USE oce_trc        ! shared variables between ocean and passive tracers
   USE trc            ! passive tracers common variables
   USE trcice_cfc     ! CFC      initialisation
   USE trcice_pisces  ! PISCES   initialisation
   USE trcice_medusa  ! MEDUSA   initialisation
   USE trcice_c14     ! C14 bomb initialisation
   USE trcice_age     ! AGE      initialisation
   USE trcice_my_trc  ! MY_TRC   initialisation
   
   IMPLICIT NONE
   PRIVATE
   
   PUBLIC   trc_ice_ini   ! called by trc_nam

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcice.F90 14086 2020-12-04 11:37:14Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS
   
   SUBROUTINE trc_ice_ini
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_ice_ini ***
      !!
      !! ** Purpose :   Initialization of the ice module for tracers
      !!
      !! ** Method  : - 
      !!---------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'trc_ice_ini : Initialize sea ice tracer boundary condition'
         WRITE(numout,*) '~~~~~~~~~~~~~'
      ENDIF
      !
      CALL trc_nam_ice
      !
      trc_i(:,:,:) = 0._wp   ! by default
      trc_o(:,:,:) = 0._wp   ! by default
      !
      IF ( nn_ice_tr == 1 ) THEN
         IF( ln_pisces  )    CALL trc_ice_ini_pisces       ! PISCES  bio-model
         IF( ln_medusa  )    CALL trc_ice_ini_medusa       ! MEDUSA  bio-model
         IF( ll_cfc     )    CALL trc_ice_ini_cfc          ! CFC     tracers
         IF( ln_c14     )    CALL trc_ice_ini_c14          ! C14     tracer
         IF( ln_age     )    CALL trc_ice_ini_age          ! AGE     tracer
         IF( ln_my_trc  )    CALL trc_ice_ini_my_trc       ! MY_TRC  tracers
      ENDIF
      !
   END SUBROUTINE trc_ice_ini


   SUBROUTINE trc_nam_ice
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_nam_ice ***
      !!
      !! ** Purpose :   Read the namelist for the ice effect on tracers
      !!
      !! ** Method  : -
      !!---------------------------------------------------------------------
      INTEGER :: jn      ! dummy loop indices
      INTEGER :: ios, ierr     ! Local integer output status for namelist read
      !
      TYPE(TRC_I_NML), DIMENSION(jpmaxtrc) ::   sn_tri_tracer
      !!
      NAMELIST/namtrc_ice/ nn_ice_tr, sn_tri_tracer
      !!---------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'trc_nam_ice : Read the namelist for trc_ice'
         WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      ENDIF
      !
      READ  ( numnat_ref, namtrc_ice, IOSTAT = ios, ERR = 901)
 901  IF( ios /= 0 )   CALL ctl_nam ( ios , ' namtrc_ice in reference namelist ' )
      READ  ( numnat_cfg, namtrc_ice, IOSTAT = ios, ERR = 902 )
 902  IF( ios >  0 )   CALL ctl_nam ( ios , 'namtrc_ice in configuration namelist' )

      IF( lwp ) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) '   Namelist : namtrc_ice'
         WRITE(numout,*) '      Sea ice tracers option (nn_ice_tr) : ', nn_ice_tr
      ENDIF
      !
      ! Assign namelist stuff
      DO jn = 1, jptra
         trc_ice_ratio (jn) = sn_tri_tracer(jn)%trc_ratio
         trc_ice_prescr(jn) = sn_tri_tracer(jn)%trc_prescr
         cn_trc_o      (jn) = sn_tri_tracer(jn)%ctrc_o
      END DO
      !
   END SUBROUTINE trc_nam_ice

#else
   !!----------------------------------------------------------------------
   !!  Empty module :                                     No passive tracer
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_ice_ini                   ! Dummy routine   
   END SUBROUTINE trc_ice_ini
   SUBROUTINE trc_nam_ice
   END SUBROUTINE trc_nam_ice
#endif

   !!======================================================================
END MODULE trcice
