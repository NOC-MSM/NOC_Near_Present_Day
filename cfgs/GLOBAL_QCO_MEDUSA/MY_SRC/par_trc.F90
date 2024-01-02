MODULE par_trc
   !!======================================================================
   !!                        ***  par_trc  ***
   !! TOP :   set the passive tracers parameters
   !!======================================================================
   !! History :    -   !  1996-01  (M. Levy)  original code
   !!              -   !  2000-04  (O. Aumont, M.A. Foujols)  HAMOCC3 and P3ZD
   !!             1.0  !  2004-03  (C. Ethe) Free form and module
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  revised architecture
   !!----------------------------------------------------------------------
   USE par_kind          ! kind parameters
   USE par_pisces        ! PISCES model  parameters
   USE par_medusa        ! MEDUSA model  parameters
   USE par_cfc           ! CFCs  tracers parameters
   USE par_c14           ! C14 tracer    parameters
   USE par_age           ! AGE tracer    parameters
   USE par_my_trc        ! MY_TRC model  parameters
   !

   IMPLICIT NONE

   INTEGER, PUBLIC,  PARAMETER :: jpmaxtrc = 100  ! Maximum number of tracers

   INTEGER, PUBLIC             :: jptra           !: Total number of passive tracers
   INTEGER, PUBLIC             :: jp_pisces       !: number of passive tracers in PISCES model
   INTEGER, PUBLIC             :: jp_cfc          !: number of CFC passive tracers 
   INTEGER, PUBLIC             :: jp_my_trc       !: number of passive tracers in MY_TRC model
   INTEGER, PUBLIC             :: jp_bgc          !: number of passive tracers for the BGC model

   INTEGER, PUBLIC             :: jp_dia3d        !: number of 3D diagnostic variables
   INTEGER, PUBLIC             :: jp_dia2d        !: number of 2D diagnostic variables

   LOGICAL, PUBLIC             :: ln_pisces       !: PISCES flag 
   LOGICAL, PUBLIC             :: ln_medusa       !: MEDUSA flag 
   LOGICAL, PUBLIC             :: ln_age          !: AGE flag 
   LOGICAL, PUBLIC             :: ln_cfc11        !: CFC11 flag 
   LOGICAL, PUBLIC             :: ln_cfc12        !: CFC12 flag 
   LOGICAL, PUBLIC             :: ln_sf6          !: SF6 flag 
   LOGICAL, PUBLIC             :: ll_cfc          !: CFC flag 
   LOGICAL, PUBLIC             :: ln_c14          !: C14 flag 
   LOGICAL, PUBLIC             :: ln_my_trc       !: MY_TRC flag 

   REAL(wp), PUBLIC            :: rtrn  = 0.5 * EPSILON( 1.e0 )    !: truncation value

END MODULE par_trc
