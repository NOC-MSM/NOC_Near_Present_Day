MODULE trcrst
   !!======================================================================
   !!                         ***  MODULE trcrst  ***
   !! TOP :   Manage the passive tracer restart
   !!======================================================================
   !! History :    -   !  1991-03  ()  original code
   !!             1.0  !  2005-03 (O. Aumont, A. El Moussaoui) F90
   !!              -   !  2005-10 (C. Ethe) print control
   !!             2.0  !  2005-10 (C. Ethe, G. Madec) revised architecture
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   trc_rst        : Restart for passive tracer
   !!   trc_rst_opn    : open  restart file
   !!   trc_rst_read   : read  restart file
   !!   trc_rst_wri    : write restart file
   !!----------------------------------------------------------------------
   USE par_trc        ! need jptra, number of passive tracers
   USE oce_trc
   USE trc
   USE iom
   USE daymod
   USE lib_mpp
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_rst_opn       ! called by ???
   PUBLIC   trc_rst_read      ! called by ???
   PUBLIC   trc_rst_wri       ! called by ???
   PUBLIC   trc_rst_cal
   PUBLIC   trc_rst_stat

#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcrst.F90 14239 2020-12-23 08:57:16Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS
   
   SUBROUTINE trc_rst_opn( kt )
      !!----------------------------------------------------------------------
      !!                    ***  trc_rst_opn  ***
      !!
      !! ** purpose  :   output of sea-trc variable in a netcdf file
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt       ! number of iteration
      !
      CHARACTER(LEN=20)   ::   clkt     ! ocean time-step define as a character
      CHARACTER(LEN=50)   ::   clname   ! trc output restart file name
      CHARACTER(LEN=256)  ::   clpath   ! full path to ocean output restart file
      CHARACTER(LEN=50)   ::   clpname  ! trc output restart file name including AGRIF
      !!----------------------------------------------------------------------
      !
      IF( l_offline ) THEN
         IF( kt == nittrc000 ) THEN
            lrst_trc = .FALSE.
            IF( ln_rst_list ) THEN
               nrst_lst = 1
               nitrst = nn_stocklist( nrst_lst )
            ELSE
               nitrst = nitend
            ENDIF
         ENDIF

         IF( .NOT. ln_rst_list .AND. MOD( kt - 1, nn_stock ) == 0 ) THEN
            ! we use kt - 1 and not kt - nittrc000 to keep the same periodicity from the beginning of the experiment
            nitrst = kt + nn_stock - 1                  ! define the next value of nitrst for restart writing
            IF( nitrst > nitend )   nitrst = nitend   ! make sure we write a restart at the end of the run
         ENDIF
      ELSE
         IF( kt == nittrc000 ) lrst_trc = .FALSE.
      ENDIF

      IF( .NOT. ln_rst_list .AND. nn_stock == -1 )   RETURN   ! we will never do any restart

      ! to get better performances with NetCDF format:
      ! we open and define the tracer restart file one tracer time step before writing the data (-> at nitrst - 1)
      ! except if we write tracer restart files every tracer time step or if a tracer restart file was writen at nitend - 1
      IF( kt == nitrst - 1 .OR. nn_stock == 1 .OR. ( kt == nitend - 1 .AND. .NOT. lrst_trc ) ) THEN
         ! beware of the format used to write kt (default is i8.8, that should be large enough)
         IF( nitrst > 1.0e9 ) THEN   ;   WRITE(clkt,*       ) nitrst
         ELSE                        ;   WRITE(clkt,'(i8.8)') nitrst
         ENDIF
         ! create the file
         IF(lwp) WRITE(numout,*)
         clname = TRIM(cexper)//"_"//TRIM(ADJUSTL(clkt))//"_"//TRIM(cn_trcrst_out)
         clpath = TRIM(cn_trcrst_outdir)
         IF( clpath(LEN_TRIM(clpath):) /= '/' ) clpath = TRIM(clpath) // '/'
         IF(lwp) WRITE(numout,*) &
             '             open trc restart.output NetCDF file: ',TRIM(clpath)//clname
         IF(.NOT.lwxios) THEN
            CALL iom_open( TRIM(clpath)//TRIM(clname), numrtw, ldwrt = .TRUE. )
         ELSE
#if defined key_xios
            cw_toprst_cxt = "rstwt_"//TRIM(ADJUSTL(clkt))
            IF( TRIM(Agrif_CFixed()) == '0' ) THEN
               clpname = clname
            ELSE
               clpname = TRIM(Agrif_CFixed())//"_"//clname
            ENDIF
            numrtw = iom_xios_setid(TRIM(clpath)//TRIM(clpname))
            CALL iom_init( cw_toprst_cxt, kdid = numrtw, ld_closedef = .FALSE. )
#else
            CALL ctl_stop( 'Can not use XIOS in trc_rst_opn' )
#endif
         ENDIF
         lrst_trc = .TRUE.
      ENDIF
      !
   END SUBROUTINE trc_rst_opn

   SUBROUTINE trc_rst_read( Kbb, Kmm )
      !!----------------------------------------------------------------------
      !!                    ***  trc_rst_opn  ***
      !!
      !! ** purpose  :   read passive tracer fields in restart files
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   Kbb, Kmm  ! time level indices
      INTEGER  ::  jn     

      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'trc_rst_read : read data in the TOP restart file'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'

      ! READ prognostic variables and computes diagnostic variable
      DO jn = 1, jptra      ! MLF only : Now time step
         IF( iom_varid( numrtr, 'TRN'//ctrcnm(jn), ldstop = .FALSE. ) > 0 ) THEN 
            !! Jpalm -- ZMP -- need to read restart without oithona
            !! ZMP is here :
            CALL iom_get( numrtr, jpdom_auto, 'TRN'//ctrcnm(jn), tr(:,:,:,jn,Kmm) )
            !! if not there, we first try to leave the default init values
         ELSE
            IF(lwp) WRITE(numout,*) 'trc_rst_read ---> TRN_',TRIM(ctrcnm(jn)),' NOT in RESTART file '
           IF(lwp) WRITE(numout,*) '                -- Keep default value. '
         ENDIF
      END DO

      IF( l_1st_euler .OR. ln_top_euler ) THEN
         IF(lwp) WRITE(numout,*) '               + adjustment for forward Euler time stepping'
         tr(:,:,:,1:jptra,Kbb) = tr(:,:,:,1:jptra,Kmm)
      ELSE
         DO jn = 1, jptra      ! RK3 and MLF : Before time step
            IF( iom_varid( numrtr, 'TRB'//ctrcnm(jn), ldstop = .FALSE. ) > 0 ) THEN 
               !! Jpalm -- ZMP -- need to read restart without oithona
               !! ZMP is here :
               CALL iom_get( numrtr, jpdom_auto, 'TRB'//ctrcnm(jn), tr(:,:,:,jn,Kbb) )
               !! if not there, we first try to leave the default init values
            ELSE
               IF(lwp) WRITE(numout,*) 'trc_rst_read ---> TRB_',TRIM(ctrcnm(jn)),' NOT in RESTART file '
               IF(lwp) WRITE(numout,*) '                -- Keep default value. '
            ENDIF
         END DO
      END IF
      !
      IF(.NOT.lrxios) CALL iom_delay_rst( 'READ', 'TOP', numrtr )   ! read only TOP delayed global communication variables
   END SUBROUTINE trc_rst_read

   SUBROUTINE trc_rst_wri( kt, Kbb, Kmm, Krhs )
      !!----------------------------------------------------------------------
      !!                    ***  trc_rst_wri  ***
      !!
      !! ** purpose  :   write passive tracer fields in restart files
      !!----------------------------------------------------------------------
      USE trcrst_medusa, ONLY: trc_rst_conserv
      INTEGER, INTENT( in ) ::   kt              ! ocean time-step index
      INTEGER, INTENT( in ) ::   Kbb, Kmm, Krhs  ! time level indices
      !!
      INTEGER  :: jn
      !!----------------------------------------------------------------------
      !
      CALL iom_rstput( kt, nitrst, numrtw, 'rdttrc1', rn_Dt )   ! passive tracer time step (= ocean time step)
      ! prognostic variables 
      ! -------------------- 
      DO jn = 1, jptra
         CALL iom_rstput( kt, nitrst, numrtw, 'TRN'//ctrcnm(jn), tr(:,:,:,jn,Kmm) )
      END DO

      DO jn = 1, jptra
         CALL iom_rstput( kt, nitrst, numrtw, 'TRB'//ctrcnm(jn), tr(:,:,:,jn,Kbb) )
      END DO

      IF( .NOT. lwxios ) CALL iom_delay_rst( 'WRITE', 'TOP', numrtw )   ! save only TOP delayed global communication variables
    
      IF( kt == nitrst ) THEN
          CALL trc_rst_stat( Kmm, Krhs )             ! statistics
          CALL trc_rst_conserv( Kbb, Kmm, Krhs )     ! Conservation Checks
          IF(lwxios) THEN
             CALL iom_context_finalize(      cw_toprst_cxt          )
             iom_file(numrtw)%nfid       = 0
             numrtw = 0
          ELSE
             CALL iom_close( numrtw )     ! close the restart file (only at last time step)
          ENDIF
#if ! defined key_trdmxl_trc
          lrst_trc = .FALSE.
#endif
          IF( l_offline .AND. ln_rst_list ) THEN
             nrst_lst = nrst_lst + 1
             nitrst = nn_stocklist( nrst_lst )
          ENDIF
      ENDIF
      !
   END SUBROUTINE trc_rst_wri 


   SUBROUTINE trc_rst_cal( kt, cdrw )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE trc_rst_cal  ***
      !!
      !!  ** Purpose : Read or write calendar in restart file:
      !!
      !!  WRITE(READ) mode:
      !!       kt        : number of time step since the begining of the experiment at the
      !!                   end of the current(previous) run
      !!       adatrj(0) : number of elapsed days since the begining of the experiment at the
      !!                   end of the current(previous) run (REAL -> keep fractions of day)
      !!       ndastp    : date at the end of the current(previous) run (coded as yyyymmdd integer)
      !!
      !!   According to namelist parameter nrstdt,
      !!       nn_rsttr = 0  no control on the date (nittrc000 is  arbitrary).
      !!       nn_rsttr = 1  we verify that nittrc000 is equal to the last
      !!                   time step of previous run + 1.
      !!       In both those options, the  exact duration of the experiment
      !!       since the beginning (cumulated duration of all previous restart runs)
      !!       is not stored in the restart and is assumed to be (nittrc000-1)*rn_Dt.
      !!       This is valid is the time step has remained constant.
      !!
      !!       nn_rsttr = 2  the duration of the experiment in days (adatrj)
      !!                    has been stored in the restart file.
      !!----------------------------------------------------------------------
      INTEGER         , INTENT(in) ::   kt         ! ocean time-step
      CHARACTER(len=*), INTENT(in) ::   cdrw       ! "READ"/"WRITE" flag
      !
      LOGICAL  ::  llok
      REAL(wp) ::  zrdttrc1, zkt, zndastp, zdayfrac, ksecs, ktime
      INTEGER  ::   ihour, iminute
      CHARACTER(len=82) :: clpname

      ! Time domain : restart
      ! ---------------------

      IF( TRIM(cdrw) == 'READ' ) THEN

         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'trc_rst_cal : read the TOP restart file for calendar'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'

         IF( ln_rsttr ) THEN
            lxios_sini = .FALSE.
            CALL iom_open( TRIM(cn_trcrst_indir)//'/'//cn_trcrst_in, numrtr )
            IF( lrxios) THEN
                cr_toprst_cxt = 'top_rst'
                IF(lwp) WRITE(numout,*) 'Enable restart reading by XIOS for TOP'
!               IF( TRIM(Agrif_CFixed()) == '0' ) THEN
!                  clpname = cn_trcrst_in
!               ELSE
!                  clpname = TRIM(Agrif_CFixed())//"_"//cn_trcrst_in   
!               ENDIF
                CALL iom_init( cr_toprst_cxt, kdid = numrtr, ld_closedef = .TRUE. )
            ENDIF

            CALL iom_get ( numrtr, 'kt', zkt )   ! last time-step of previous run

            IF(lwp) THEN
               WRITE(numout,*) ' *** Info read in restart : '
               WRITE(numout,*) '   previous time-step                               : ', NINT( zkt )
               WRITE(numout,*) ' *** restart option'
               SELECT CASE ( nn_rsttr )
               CASE ( 0 )   ;   WRITE(numout,*) ' nn_rsttr = 0 : no control of nittrc000'
               CASE ( 1 )   ;   WRITE(numout,*) ' nn_rsttr = 1 : no control the date at nittrc000 (use ndate0 read in the namelist)'
               CASE ( 2 )   ;   WRITE(numout,*) ' nn_rsttr = 2 : calendar parameters read in restart'
               END SELECT
               WRITE(numout,*)
            ENDIF
            ! Control of date 
            IF( nittrc000  - NINT( zkt ) /= 1 .AND.  nn_rsttr /= 0 )                                  &
               &   CALL ctl_stop( ' ===>>>> : problem with nittrc000 for the restart',                 &
               &                  ' verify the restart file or rerun with nn_rsttr = 0 (namelist)' )
         ENDIF
         !
         IF( l_offline ) THEN    
            !                                          ! set the date in offline mode
            IF( ln_rsttr .AND. nn_rsttr == 2 ) THEN
               CALL iom_get( numrtr, 'ndastp', zndastp )
               ndastp = NINT( zndastp )
               CALL iom_get( numrtr, 'adatrj', adatrj  )
               CALL iom_get( numrtr, 'ntime' , ktime   )
               nn_time0=INT(ktime)
               ! calculate start time in hours and minutes
               zdayfrac=adatrj-INT(adatrj)
               ksecs = NINT(zdayfrac*86400)            ! Nearest second to catch rounding errors in adatrj
               ihour = INT(ksecs/3600)
               iminute = ksecs/60-ihour*60
                
               ! Add to nn_time0
               nhour   =   nn_time0 / 100
               nminute = ( nn_time0 - nhour * 100 )
               nminute=nminute+iminute
               
               IF( nminute >= 60 ) THEN
                  nminute=nminute-60
                  nhour=nhour+1
               ENDIF
               nhour=nhour+ihour
               IF( nhour >= 24 ) THEN
                  nhour=nhour-24
                  adatrj=adatrj+1
               ENDIF           
               nn_time0 = nhour * 100 + nminute
               adatrj = INT(adatrj)                    ! adatrj set to integer as nn_time0 updated            
             ELSE
               ndt05 = NINT( 0.5 * rn_Dt  )   !  --- WARNING --- not defined yet are we did not go through day_init
               ! parameters corresponding to nit000 - 1 (as we start the step
               ! loop with a call to day)
               ndastp = ndate0        ! ndate0 read in the namelist in dom_nam
               nhour   =   nn_time0 / 100
               nminute = ( nn_time0 - nhour * 100 )
               IF( nhour*3600+nminute*60-ndt05 .lt. 0 )  ndastp=ndastp-1      ! Start hour is specified in the namelist (default 0)
               adatrj = ( REAL( nit000-1, wp ) * rn_Dt ) / rday
               ! note this is wrong if time step has changed during run
            ENDIF
            IF( ABS(adatrj  - REAL(NINT(adatrj),wp)) < 0.1 / rday )   adatrj = REAL(NINT(adatrj),wp)   ! avoid truncation error
            !
            IF(lwp) THEN
              WRITE(numout,*) ' *** Info used values : '
              WRITE(numout,*) '   date ndastp                                      : ', ndastp
              WRITE(numout,*) '   number of elapsed days since the begining of run : ', adatrj
              WRITE(numout,*) '   nn_time0                                         : ', nn_time0
              WRITE(numout,*)
            ENDIF
            !
            IF( ln_rsttr )  THEN   ;    l_1st_euler = .false.
            ELSE                   ;    l_1st_euler = .true.
            ENDIF
            !
            CALL day_init          ! compute calendar
            !
         ENDIF
         !
      ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN
         !
         IF(  kt == nitrst ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'trc_wri : write the TOP restart file (NetCDF) at it= ', kt, ' date= ', ndastp
            IF(lwp) WRITE(numout,*) '~~~~~~~'
         ENDIF
         CALL iom_rstput( kt, nitrst, numrtw, 'kt'     , REAL( kt    , wp)   )   ! time-step
         CALL iom_rstput( kt, nitrst, numrtw, 'ndastp' , REAL( ndastp, wp)   )   ! date
         CALL iom_rstput( kt, nitrst, numrtw, 'adatrj' , adatrj              )   ! number of elapsed days since
         !                                                                     ! the begining of the run [s]
         CALL iom_rstput( kt, nitrst, numrtw, 'ntime'  , REAL( nn_time0, wp) ) ! time
      ENDIF

   END SUBROUTINE trc_rst_cal


   SUBROUTINE trc_rst_stat( Kmm, Krhs )
      !!----------------------------------------------------------------------
      !!                    ***  trc_rst_stat  ***
      !!
      !! ** purpose  :   Compute tracers statistics
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   Kmm, Krhs  ! time level indices
      INTEGER  :: jk, jn
      REAL(wp) :: ztraf, zmin, zmax, zmean, zdrift
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zvol
      !!----------------------------------------------------------------------

      IF( lwp ) THEN
         WRITE(numout,*) 
         WRITE(numout,*) '           ----TRACER STAT----             '
         WRITE(numout,*) 
      ENDIF
      !
      DO jk = 1, jpk
         zvol(:,:,jk) = e1e2t(:,:) * e3t(:,:,jk,Kmm) * tmask(:,:,jk)
      END DO
      !
      DO jn = 1, jptra
         ztraf = glob_sum( 'trcrst', tr(:,:,:,jn,Kmm) * zvol(:,:,:) )
         zmin  = MINVAL( tr(:,:,:,jn,Kmm), mask= ((tmask*SPREAD(tmask_i,DIM=3,NCOPIES=jpk).NE.0.)) )
         zmax  = MAXVAL( tr(:,:,:,jn,Kmm), mask= ((tmask*SPREAD(tmask_i,DIM=3,NCOPIES=jpk).NE.0.)) )
         IF( lk_mpp ) THEN
            CALL mpp_min( 'trcrst', zmin )      ! min over the global domain
            CALL mpp_max( 'trcrst', zmax )      ! max over the global domain
         END IF
         zmean  = ztraf / areatot
         zdrift = ( ( ztraf - trai(jn) ) / ( trai(jn) + 1.e-12 )  ) * 100._wp
         IF(lwp) WRITE(numout,9000) jn, TRIM( ctrcnm(jn) ), zmean, zmin, zmax, zdrift
      END DO
      IF(lwp) WRITE(numout,*) 
9000  FORMAT(' tracer nb :',i2,'    name :',a10,'    mean :',e18.10,'    min :',e18.10, &
      &      '    max :',e18.10,'    drift :',e18.10, ' %')
      !
   END SUBROUTINE trc_rst_stat

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                     No passive tracer
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_rst_read( Kbb, Kmm)                      ! Empty routines
      INTEGER, INTENT( in ) :: Kbb, Kmm  ! time level indices
   END SUBROUTINE trc_rst_read
   SUBROUTINE trc_rst_wri( kt, Kbb, Kmm, Krhs )
      INTEGER, INTENT( in ) :: kt
      INTEGER, INTENT( in ) :: Kbb, Kmm, Krhs  ! time level indices
      WRITE(*,*) 'trc_rst_wri: You should not have seen this print! error?', kt
   END SUBROUTINE trc_rst_wri   
#endif

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcrst.F90 14239 2020-12-23 08:57:16Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!======================================================================
END MODULE trcrst
