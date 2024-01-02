MODULE icerst
   !!======================================================================
   !!                     ***  MODULE  icerst  ***
   !!   sea-ice :  write/read the ice restart file
   !!======================================================================
   !! History:   4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_rst_opn   : open  restart file
   !!   ice_rst_write : write restart file
   !!   ice_rst_read  : read  restart file
   !!----------------------------------------------------------------------
   USE ice            ! sea-ice: variables
   USE dom_oce        ! ocean domain
   USE phycst  , ONLY : rt0
   USE sbc_oce , ONLY : nn_fsbc, ln_cpl
   USE sbc_oce , ONLY : nn_components, jp_iam_sas   ! SAS ss[st]_m init
   USE sbc_oce , ONLY : sst_m, sss_m                ! SAS ss[st]_m init
   USE oce     , ONLY : ts                          ! SAS ss[st]_m init
   USE eosbn2  , ONLY : l_useCT, eos_pt_from_ct     ! SAS ss[st]_m init
   USE iceistate      ! sea-ice: initial state
   USE icectl         ! sea-ice: control
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE ioipsl  , ONLY : ju2ymds                     ! for calendar
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_rst_opn     ! called by icestp
   PUBLIC   ice_rst_write   ! called by icestp
   PUBLIC   ice_rst_read    ! called by ice_init

   !!----------------------------------------------------------------------
   !! NEMO/ICE 4.0 , NEMO Consortium (2018)
   !! $Id: icerst.F90 14239 2020-12-23 08:57:16Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_rst_opn( kt )
      !!----------------------------------------------------------------------
      !!                    ***  ice_rst_opn  ***
      !!
      !! ** purpose  :   open restart file
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt       ! number of iteration
      !
      INTEGER             ::   iyear, imonth, iday
      REAL (wp)           ::   zsec
      REAL (wp)           ::   zfjulday
      CHARACTER(len=20)   ::   clkt     ! ocean time-step define as a character
      CHARACTER(len=50)   ::   clname   ! ice output restart file name
      CHARACTER(len=256)  ::   clpath   ! full path to ice output restart file
      CHARACTER(LEN=52)   ::   clpname   ! ocean output restart file name including prefix for AGRIF
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 )   lrst_ice = .FALSE.   ! default definition

      IF( ln_rst_list .OR. nn_stock /= -1 ) THEN
      ! in order to get better performances with NetCDF format, we open and define the ice restart file
      ! one ice time step before writing the data (-> at nitrst - 2*nn_fsbc + 1), except if we write ice
      ! restart files every ice time step or if an ice restart file was writen at nitend - 2*nn_fsbc + 1
      IF( kt == nitrst - 2*nn_fsbc + 1 .OR. nn_stock == nn_fsbc    &
         &                             .OR. ( kt == nitend - nn_fsbc + 1 .AND. .NOT. lrst_ice ) ) THEN
         IF( nitrst <= nitend .AND. nitrst > 0 ) THEN
            ! beware of the format used to write kt (default is i8.8, that should be large enough...)
            IF ( ln_rstdate ) THEN
               zfjulday = fjulday + (2*nn_fsbc+1)*rdt / rday
               IF( ABS(zfjulday - REAL(NINT(zfjulday),wp)) < 0.1 / rday )   zfjulday = REAL(NINT(zfjulday),wp)   ! avoid truncation error
               CALL ju2ymds( zfjulday, iyear, imonth, iday, zsec )           
               WRITE(clkt, '(i4.4,2i2.2)') iyear, imonth, iday
            ELSE
               IF( nitrst > 99999999 ) THEN   ;   WRITE(clkt, *       ) nitrst
               ELSE                           ;   WRITE(clkt, '(i8.8)') nitrst
               ENDIF
            ENDIF
            ! create the file
            clname = TRIM(cexper)//"_"//TRIM(ADJUSTL(clkt))//"_"//TRIM(cn_icerst_out)
            clpath = TRIM(cn_icerst_outdir)
            IF( clpath(LEN_TRIM(clpath):) /= '/' ) clpath = TRIM(clpath)//'/'
            IF(lwp) THEN
               WRITE(numout,*)
               WRITE(numout,*) '             open ice restart NetCDF file: ',TRIM(clpath)//clname
               IF( kt == nitrst - 2*nn_fsbc + 1 ) THEN
                  WRITE(numout,*) '             kt = nitrst - 2*nn_fsbc + 1 = ', kt,' date= ', ndastp
               ELSE
                  WRITE(numout,*) '             kt = '                         , kt,' date= ', ndastp
               ENDIF
            ENDIF
            !
            IF(.NOT.lwxios) THEN
               CALL iom_open( TRIM(clpath)//TRIM(clname), numriw, ldwrt = .TRUE., kdlev = jpl, cdcomp = 'ICE' )
            ELSE
#if defined key_xios
               cw_icerst_cxt = "rstwi_"//TRIM(ADJUSTL(clkt))
               IF( TRIM(Agrif_CFixed()) == '0' ) THEN
                  clpname = clname
               ELSE
                  clpname = TRIM(Agrif_CFixed())//"_"//clname
               ENDIF
               numriw = iom_xios_setid(TRIM(clpath)//TRIM(clpname))
               CALL iom_init( cw_icerst_cxt, kdid = numriw, ld_closedef = .FALSE. )
               CALL iom_swap( cxios_context )
#else
               CALL ctl_stop( 'Can not use XIOS in rst_opn' )
#endif
            ENDIF
            lrst_ice = .TRUE.
         ENDIF
      ENDIF
      ENDIF
      !
      IF( ln_icectl )   CALL ice_prt( kt, iiceprt, jiceprt, 1, ' - Beginning the time step - ' )   ! control print
      !
   END SUBROUTINE ice_rst_opn


   SUBROUTINE ice_rst_write( kt )
      !!----------------------------------------------------------------------
      !!                    ***  ice_rst_write  ***
      !!
      !! ** purpose  :   write restart file
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt     ! number of iteration
      !!
      INTEGER ::   jk    ! dummy loop indices
      INTEGER ::   iter
      CHARACTER(len=25) ::   znam
      CHARACTER(len=2)  ::   zchar, zchar1
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   z3d   ! 3D workspace
      !!----------------------------------------------------------------------

      iter = kt + nn_fsbc - 1   ! ice restarts are written at kt == nitrst - nn_fsbc + 1

      IF( iter == nitrst ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'ice_rst_write : write ice restart file  kt =', kt
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~'
      ENDIF

      ! Write in numriw (if iter == nitrst)
      ! ------------------
      !                                                                        ! calendar control
      CALL iom_rstput( iter, nitrst, numriw, 'nn_fsbc', REAL( nn_fsbc, wp ) )      ! time-step
      CALL iom_rstput( iter, nitrst, numriw, 'kt_ice' , REAL( iter   , wp ) )      ! date

      IF(.NOT.lwxios) CALL iom_delay_rst( 'WRITE', 'ICE', numriw )   ! save only ice delayed global communication variables

      ! Prognostic variables
      CALL iom_rstput( iter, nitrst, numriw, 'v_i'  , v_i   )
      CALL iom_rstput( iter, nitrst, numriw, 'v_s'  , v_s   )
      CALL iom_rstput( iter, nitrst, numriw, 'sv_i' , sv_i  )
      CALL iom_rstput( iter, nitrst, numriw, 'a_i'  , a_i   )
      CALL iom_rstput( iter, nitrst, numriw, 't_su' , t_su  )
      CALL iom_rstput( iter, nitrst, numriw, 'u_ice', u_ice )
      CALL iom_rstput( iter, nitrst, numriw, 'v_ice', v_ice )
      CALL iom_rstput( iter, nitrst, numriw, 'oa_i' , oa_i  )
      CALL iom_rstput( iter, nitrst, numriw, 'a_ip' , a_ip  )
      CALL iom_rstput( iter, nitrst, numriw, 'v_ip' , v_ip  )
      CALL iom_rstput( iter, nitrst, numriw, 'v_il' , v_il  )
      ! Snow enthalpy
      DO jk = 1, nlay_s
         WRITE(zchar1,'(I2.2)') jk
         znam = 'e_s'//'_l'//zchar1
         z3d(:,:,:) = e_s(:,:,jk,:)
         CALL iom_rstput( iter, nitrst, numriw, znam , z3d )
      END DO
      ! Ice enthalpy
      DO jk = 1, nlay_i
         WRITE(zchar1,'(I2.2)') jk
         znam = 'e_i'//'_l'//zchar1
         z3d(:,:,:) = e_i(:,:,jk,:)
         CALL iom_rstput( iter, nitrst, numriw, znam , z3d )
      END DO
      ! fields needed for Met Office (Jules) coupling
      IF( ln_cpl ) THEN
         CALL iom_rstput( iter, nitrst, numriw, 'cnd_ice', cnd_ice )
         CALL iom_rstput( iter, nitrst, numriw, 't1_ice' , t1_ice )
      ENDIF
      !

      ! close restart file
      ! ------------------
      IF( iter == nitrst ) THEN
         IF(.NOT.lwxios) THEN
            CALL iom_close( numriw )
         ELSE
            CALL iom_context_finalize(      cw_icerst_cxt          )
            iom_file(numriw)%nfid       = 0
            numriw = 0
         ENDIF
         lrst_ice = .FALSE.
      ENDIF
      !
   END SUBROUTINE ice_rst_write


   SUBROUTINE ice_rst_read( Kbb, Kmm, Kaa )
      !!----------------------------------------------------------------------
      !!                    ***  ice_rst_read  ***
      !!
      !! ** purpose  :   read restart file
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) :: Kbb, Kmm, Kaa ! ocean time level indices
      INTEGER           ::   jk
      LOGICAL           ::   llok
      INTEGER           ::   id0, id1, id2, id3, id4, id5   ! local integer
      CHARACTER(len=25) ::   znam
      CHARACTER(len=2)  ::   zchar, zchar1
      REAL(wp)          ::   zfice, ziter
      CHARACTER(lc)     ::   clpname
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   z3d   ! 3D workspace
      !!----------------------------------------------------------------------

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'ice_rst_read: read ice NetCDF restart file'
         WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF

      lxios_sini = .FALSE.
      CALL iom_open ( TRIM(cn_icerst_indir)//'/'//cn_icerst_in, numrir )

      IF( lrxios) THEN
          cr_icerst_cxt = 'si3_rst'
          IF(lwp) WRITE(numout,*) 'Enable restart reading by XIOS for SI3'
!         IF( TRIM(Agrif_CFixed()) == '0' ) THEN
!            clpname = cn_icerst_in
!         ELSE
!            clpname = TRIM(Agrif_CFixed())//"_"//cn_icerst_in
!         ENDIF
          CALL iom_init( cr_icerst_cxt, kdid = numrir, ld_closedef = .TRUE. )
      ENDIF

      ! test if v_i exists
      id0 = iom_varid( numrir, 'v_i' , ldstop = .FALSE. )

      !                    ! ------------------------------ !
      IF( id0 > 0 ) THEN   ! == case of a normal restart == !
         !                 ! ------------------------------ !
         ! Time info
         CALL iom_get( numrir, 'nn_fsbc', zfice )
         CALL iom_get( numrir, 'kt_ice' , ziter )
         IF(lwp) WRITE(numout,*) '   read ice restart file at time step    : ', ziter
         IF(lwp) WRITE(numout,*) '   in any case we force it to nit000 - 1 : ', nit000 - 1

         ! Control of date
         IF( ( nit000 - NINT(ziter) ) /= 1 .AND. ABS( nrstdt ) == 1 )   &
            &     CALL ctl_stop( 'ice_rst_read ===>>>> : problem with nit000 in ice restart',  &
            &                   '   verify the file or rerun with the value 0 for the',        &
            &                   '   control of time parameter  nrstdt' )
         IF( NINT(zfice) /= nn_fsbc          .AND. ABS( nrstdt ) == 1 )   &
            &     CALL ctl_stop( 'ice_rst_read ===>>>> : problem with nn_fsbc in ice restart',  &
            &                   '   verify the file or rerun with the value 0 for the',         &
            &                   '   control of time parameter  nrstdt' )

         ! --- mandatory fields --- !
         CALL iom_get( numrir, jpdom_auto, 'v_i'  , v_i   )
         CALL iom_get( numrir, jpdom_auto, 'v_s'  , v_s   )
         CALL iom_get( numrir, jpdom_auto, 'sv_i' , sv_i  )
         CALL iom_get( numrir, jpdom_auto, 'a_i'  , a_i   )
         CALL iom_get( numrir, jpdom_auto, 't_su' , t_su  )
         CALL iom_get( numrir, jpdom_auto, 'u_ice', u_ice, cd_type = 'U', psgn = -1._wp )
         CALL iom_get( numrir, jpdom_auto, 'v_ice', v_ice, cd_type = 'V', psgn = -1._wp )
         ! Snow enthalpy
         DO jk = 1, nlay_s
            WRITE(zchar1,'(I2.2)') jk
            znam = 'e_s'//'_l'//zchar1
            CALL iom_get( numrir, jpdom_auto, znam , z3d )
            e_s(:,:,jk,:) = z3d(:,:,:)
         END DO
         ! Ice enthalpy
         DO jk = 1, nlay_i
            WRITE(zchar1,'(I2.2)') jk
            znam = 'e_i'//'_l'//zchar1
            CALL iom_get( numrir, jpdom_auto, znam , z3d )
            e_i(:,:,jk,:) = z3d(:,:,:)
         END DO
         ! -- optional fields -- !
         ! ice age
         id1 = iom_varid( numrir, 'oa_i' , ldstop = .FALSE. )
         IF( id1 > 0 ) THEN                       ! fields exist
            CALL iom_get( numrir, jpdom_auto, 'oa_i', oa_i )
         ELSE                                     ! start from rest
            IF(lwp) WRITE(numout,*) '   ==>>   previous run without ice age output then set it to zero'
            oa_i(:,:,:) = 0._wp
         ENDIF
         ! melt ponds
         id2 = iom_varid( numrir, 'a_ip' , ldstop = .FALSE. )
         IF( id2 > 0 ) THEN                       ! fields exist
            CALL iom_get( numrir, jpdom_auto, 'a_ip' , a_ip )
            CALL iom_get( numrir, jpdom_auto, 'v_ip' , v_ip )
         ELSE                                     ! start from rest
            IF(lwp) WRITE(numout,*) '   ==>>   previous run without melt ponds output then set it to zero'
            a_ip(:,:,:) = 0._wp
            v_ip(:,:,:) = 0._wp
         ENDIF
         ! melt pond lids
         id3 = iom_varid( numrir, 'v_il' , ldstop = .FALSE. )
         IF( id3 > 0 ) THEN
            CALL iom_get( numrir, jpdom_auto, 'v_il', v_il)
         ELSE
            IF(lwp) WRITE(numout,*) '   ==>>   previous run without melt ponds lids output then set it to zero'
            v_il(:,:,:) = 0._wp
         ENDIF
         ! fields needed for Met Office (Jules) coupling
         IF( ln_cpl ) THEN
            id4 = iom_varid( numrir, 'cnd_ice' , ldstop = .FALSE. )
            id5 = iom_varid( numrir, 't1_ice'  , ldstop = .FALSE. )
            IF( id4 > 0 .AND. id5 > 0 ) THEN         ! fields exist
               CALL iom_get( numrir, jpdom_auto, 'cnd_ice', cnd_ice )
               CALL iom_get( numrir, jpdom_auto, 't1_ice' , t1_ice  )
            ELSE                                     ! start from rest
               IF(lwp) WRITE(numout,*) '   ==>>   previous run without conductivity output then set it to zero'
               cnd_ice(:,:,:) = 0._wp
               t1_ice (:,:,:) = rt0
            ENDIF
         ENDIF

         IF(.NOT.lrxios) CALL iom_delay_rst( 'READ', 'ICE', numrir )   ! read only ice delayed global communication variables
         !                 ! ---------------------------------- !
      ELSE                 ! == case of a simplified restart == !
         !                 ! ---------------------------------- !
         CALL ctl_warn('ice_rst_read: you are attempting to use an unsuitable ice restart')
         !
         IF( .NOT. ln_iceini .OR. nn_iceini_file == 2 ) THEN
            CALL ctl_stop('STOP', 'ice_rst_read: you need ln_ice_ini=T and nn_iceini_file=0 or 1')
         ELSE
            CALL ctl_warn('ice_rst_read: using ice_istate to set initial conditions instead')
         ENDIF
         !
         IF( nn_components == jp_iam_sas ) THEN   ! SAS case: ss[st]_m were not initialized by sbc_ssm_init
            !
            IF(lwp) WRITE(numout,*) '  SAS: default initialisation of ss[st]_m arrays used in ice_istate'
            IF( l_useCT )  THEN    ;   sst_m(:,:) = eos_pt_from_ct( ts(:,:,1,jp_tem, Kmm), ts(:,:,1,jp_sal, Kmm) )
            ELSE                   ;   sst_m(:,:) = ts(:,:,1,jp_tem, Kmm)
            ENDIF
            sss_m(:,:) = ts(:,:,1,jp_sal, Kmm)
         ENDIF
         !
         CALL ice_istate( nit000, Kbb, Kmm, Kaa )
         !
      ENDIF

   END SUBROUTINE ice_rst_read

#else
   !!----------------------------------------------------------------------
   !!   Default option :       Empty module           NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE icerst
