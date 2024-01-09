MODULE trcsms
   !!======================================================================
   !!                         ***  MODULE trcsms  ***
   !! TOP :   Time loop of passive tracers sms
   !!======================================================================
   !! History :   1.0  !  2005-03 (O. Aumont, A. El Moussaoui) F90
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  revised architecture
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !!   trc_sms        :  Time loop of passive tracers sms
   !!----------------------------------------------------------------------
   USE oce_trc            !
   USE trc                !
   USE trcsms_pisces      ! PISCES biogeo-model
   USE trcsms_medusa      ! MEDUSA biogeo-model
   USE trcsms_cfc         ! CFC 11 &/or 12
   USE trcsms_c14         ! C14 
   USE trcsms_age         ! AGE
   USE trcsms_my_trc      ! MY_TRC  tracers
   USE prtctl             ! Print control for debbuging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sms    ! called in trcstp.F90

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcsms.F90 13286 2020-07-09 15:48:29Z smasson $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_sms( kt, Kbb, Kmm , Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_sms  ***
      !!
      !! ** Purpose :   Managment of the time loop of passive tracers sms 
      !!
      !! ** Method  : -  call the main routine of of each defined tracer model
      !! -------------------------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt        ! ocean time-step index      
      INTEGER, INTENT( in ) ::   Kbb, Kmm, Krhs ! time level indices
      !!
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('trc_sms')
      !
      IF( ln_pisces  )   CALL trc_sms_pisces ( kt, Kbb, Kmm, Krhs )    ! main program of PISCES 
      IF( ln_medusa  )   CALL trc_sms_medusa ( kt, Kbb, Kmm, Krhs )    ! main program of MEDUSA 
      IF( ll_cfc     )   CALL trc_sms_cfc    ( kt, Kbb, Kmm, Krhs )    ! surface fluxes of CFC
      IF( ln_c14     )   CALL trc_sms_c14    ( kt, Kbb, Kmm, Krhs )    ! surface fluxes of C14
      IF( ln_age     )   CALL trc_sms_age    ( kt, Kbb, Kmm, Krhs )    ! Age tracer
      IF( ln_my_trc  )   CALL trc_sms_my_trc ( kt, Kbb, Kmm, Krhs )    ! MY_TRC  tracers

      IF(sn_cfctl%l_prttrc) THEN                       ! print mean trends (used for debugging)
         WRITE(charout, FMT="('sms ')")
         CALL prt_ctl_info( charout, cdcomp = 'top' )
         CALL prt_ctl( tab4d_1=tr(:,:,:,:,Kmm), mask1=tmask, clinfo=ctrcnm )
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('trc_sms')
      !
   END SUBROUTINE trc_sms

#else
   !!======================================================================
   !!  Dummy module :                                     No passive tracer
   !!======================================================================
CONTAINS
   SUBROUTINE trc_sms( kt )                   ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'trc_sms: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sms
#endif 

   !!======================================================================
END MODULE trcsms
