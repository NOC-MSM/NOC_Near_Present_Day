MODULE sbcrnf
   !!======================================================================
   !!                       ***  MODULE  sbcrnf  ***
   !! Ocean forcing:  river runoff
   !!=====================================================================
   !! History :  OPA  ! 2000-11  (R. Hordoir, E. Durand)  NetCDF FORMAT
   !!   NEMO     1.0  ! 2002-09  (G. Madec)  F90: Free form and module
   !!            3.0  ! 2006-07  (G. Madec)  Surface module
   !!            3.2  ! 2009-04  (B. Lemaire)  Introduce iom_put
   !!            3.3  ! 2010-10  (R. Furner, G. Madec) runoff distributed over ocean levels
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_rnf       : monthly runoffs read in a NetCDF file
   !!   sbc_rnf_init  : runoffs initialisation
   !!   rnf_mouth     : set river mouth mask
   !!----------------------------------------------------------------------
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE sbc_oce        ! surface boundary condition variables
   USE eosbn2         ! Equation Of State
   USE closea, ONLY: l_clo_rnf, clo_rnf ! closed seas
   !
   USE in_out_manager ! I/O manager
   USE fldread        ! read input field at current time step
   USE iom            ! I/O module
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_rnf       ! called in sbcmod module
   PUBLIC   sbc_rnf_div   ! called in divhor module
   PUBLIC   sbc_rnf_alloc ! called in sbcmod module
   PUBLIC   sbc_rnf_init  ! called in sbcmod module

   !                                                !!* namsbc_rnf namelist *
   CHARACTER(len=100)         ::   cn_dir            !: Root directory for location of rnf files
   LOGICAL           , PUBLIC ::   ln_rnf_depth      !: depth       river runoffs attribute specified in a file
   LOGICAL                    ::      ln_rnf_depth_ini  !: depth       river runoffs  computed at the initialisation
   REAL(wp)                   ::      rn_rnf_max        !: maximum value of the runoff climatologie (ln_rnf_depth_ini =T)
   REAL(wp)                   ::      rn_dep_max        !: depth over which runoffs is spread       (ln_rnf_depth_ini =T)
   INTEGER                    ::      nn_rnf_depth_file !: create (=1) a runoff depth file or not (=0)
   LOGICAL           , PUBLIC ::   ln_rnf_icb        !: iceberg flux is specified in a file
   LOGICAL                    ::   ln_rnf_tem        !: temperature river runoffs attribute specified in a file
   LOGICAL           , PUBLIC ::   ln_rnf_sal        !: salinity    river runoffs attribute specified in a file
   TYPE(FLD_N)       , PUBLIC ::   sn_rnf            !: information about the runoff file to be read
   TYPE(FLD_N)                ::   sn_cnf            !: information about the runoff mouth file to be read
   TYPE(FLD_N)                ::   sn_i_rnf        !: information about the iceberg flux file to be read
   TYPE(FLD_N)                ::   sn_s_rnf          !: information about the salinities of runoff file to be read
   TYPE(FLD_N)                ::   sn_t_rnf          !: information about the temperatures of runoff file to be read
   TYPE(FLD_N)                ::   sn_dep_rnf        !: information about the depth which river inflow affects
   LOGICAL           , PUBLIC ::   ln_rnf_mouth      !: specific treatment in mouths vicinity
   REAL(wp)                   ::   rn_hrnf           !: runoffs, depth over which enhanced vertical mixing is used
   REAL(wp)          , PUBLIC ::   rn_avt_rnf        !: runoffs, value of the additional vertical mixing coef. [m2/s]
   REAL(wp)          , PUBLIC ::   rn_rfact          !: multiplicative factor for runoff

   LOGICAL , PUBLIC ::   l_rnfcpl = .false.   !: runoffs recieved from oasis
   INTEGER , PUBLIC ::   nkrnf = 0            !: nb of levels over which Kz is increased at river mouths

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   rnfmsk              !: river mouth mask (hori.)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)     ::   rnfmsk_z            !: river mouth mask (vert.)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   h_rnf               !: depth of runoff in m
   INTEGER,  PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   nk_rnf              !: depth of runoff in model levels
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   rnf_tsc_b, rnf_tsc  !: before and now T & S runoff contents   [K.m/s & PSU.m/s]

   TYPE(FLD), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   sf_rnf       ! structure: river runoff (file information, fields read)
   TYPE(FLD),        ALLOCATABLE, DIMENSION(:) ::   sf_i_rnf     ! structure: iceberg flux (file information, fields read)
   TYPE(FLD),        ALLOCATABLE, DIMENSION(:) ::   sf_s_rnf     ! structure: river runoff salinity (file information, fields read)
   TYPE(FLD),        ALLOCATABLE, DIMENSION(:) ::   sf_t_rnf     ! structure: river runoff temperature (file information, fields read)

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "read_nml_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: sbcrnf.F90 14993 2021-06-14 22:35:18Z techene $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION sbc_rnf_alloc()
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE sbc_rnf_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( rnfmsk(jpi,jpj)        , rnfmsk_z(jpk)         ,     & ! needed over the whole domain by muscl (traadv_muscl)
         &      h_rnf (A2D(1))         , nk_rnf  (A2D(1))      ,     &
         &      rnf_tsc_b(A2D(0),jpts) , rnf_tsc (A2D(0),jpts) , STAT=sbc_rnf_alloc )
         !
      CALL mpp_sum ( 'sbcrnf', sbc_rnf_alloc )
      IF( sbc_rnf_alloc > 0 )   CALL ctl_warn('sbc_rnf_alloc: allocation of arrays failed')
   END FUNCTION sbc_rnf_alloc


   SUBROUTINE sbc_rnf( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_rnf  ***
      !!
      !! ** Purpose :   Introduce a climatological run off forcing
      !!
      !! ** Method  :   Set each river mouth with a monthly climatology
      !!                provided from different data.
      !!                CAUTION : upward water flux, runoff forced to be < 0
      !!
      !! ** Action  :   runoff updated runoff field at time-step kt
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt          ! ocean time step
      !
      INTEGER  ::   ji, jj    ! dummy loop indices
      INTEGER  ::   z_err = 0 ! dummy integer for error handling
      !!----------------------------------------------------------------------
      !
      !                                            !-------------------!
      !                                            !   Update runoff   !
      !                                            !-------------------!
      !
      !
      IF( .NOT. l_rnfcpl )  THEN
                             CALL fld_read ( kt, nn_fsbc, sf_rnf   )    ! Read Runoffs data and provide it at kt ( runoffs + iceberg )
         IF( ln_rnf_icb )    CALL fld_read ( kt, nn_fsbc, sf_i_rnf )    ! idem for iceberg flux if required
      ENDIF
      IF(   ln_rnf_tem   )   CALL fld_read ( kt, nn_fsbc, sf_t_rnf )    ! idem for runoffs temperature if required
      IF(   ln_rnf_sal   )   CALL fld_read ( kt, nn_fsbc, sf_s_rnf )    ! idem for runoffs salinity    if required
      !
      IF( MOD( kt - 1, nn_fsbc ) == 0 ) THEN
         !
         IF( .NOT. l_rnfcpl ) THEN
            rnf(A2D(0)) = rn_rfact * ( sf_rnf(1)%fnow(:,:,1) ) * smask0(:,:)  ! updated runoff value at time step kt
            IF( ln_rnf_icb ) THEN
               fwficb(:,:) = rn_rfact * ( sf_i_rnf(1)%fnow(:,:,1) ) * smask0(:,:)  ! updated runoff value at time step kt
               rnf(A2D(0)) = rnf(A2D(0)) + fwficb(:,:)
               qns(:,:) = qns(:,:) - fwficb(:,:) * rLfus
               !!qns_tot(:,:) = qns_tot(:,:) - fwficb(:,:) * rLfus
               !!qns_oce(:,:) = qns_oce(:,:) - fwficb(:,:) * rLfus
               CALL iom_put( 'iceberg_cea'  ,  fwficb(:,:)  )          ! output iceberg flux
               CALL iom_put( 'hflx_icb_cea' , -fwficb(:,:) * rLfus )   ! output Heat Flux into Sea Water due to Iceberg Thermodynamics -->
            ENDIF
         ENDIF
         !
         !                                                           ! set temperature & salinity content of runoffs
         IF( ln_rnf_tem ) THEN                                       ! use runoffs temperature data
            rnf_tsc(:,:,jp_tem) = ( sf_t_rnf(1)%fnow(:,:,1) ) * rnf(A2D(0)) * r1_rho0
            WHERE( sf_t_rnf(1)%fnow(:,:,1) == -999._wp )             ! if missing data value use SST as runoffs temperature
               rnf_tsc(:,:,jp_tem) = sst_m(A2D(0)) * rnf(A2D(0)) * r1_rho0
            END WHERE
         ELSE                                                        ! use SST as runoffs temperature
            !CEOD River is fresh water so must at least be 0 unless we consider ice
            rnf_tsc(:,:,jp_tem) = MAX( sst_m(A2D(0)), 0.0_wp ) * rnf(A2D(0)) * r1_rho0
         ENDIF
         !                                                           ! use runoffs salinity data
         IF( ln_rnf_sal )   rnf_tsc(:,:,jp_sal) = ( sf_s_rnf(1)%fnow(:,:,1) ) * rnf(A2D(0)) * r1_rho0
         !                                                           ! else use S=0 for runoffs (done one for all in the init)
                                         CALL iom_put( 'runoffs'     , rnf(:,:)                         )   ! output runoff mass flux
         IF( iom_use('hflx_rnf_cea') )   CALL iom_put( 'hflx_rnf_cea', rnf_tsc(:,:,jp_tem) * rho0 * rcp )   ! output runoff sensible heat (W/m2)
         IF( iom_use('sflx_rnf_cea') )   CALL iom_put( 'sflx_rnf_cea', rnf_tsc(:,:,jp_sal) * rho0       )   ! output runoff salt flux (g/m2/s)
      ENDIF
      !
#if ! defined key_RK3
      !                                                ! ---------------------------------------- !
      IF( kt == nit000 ) THEN                          !   set the forcing field at nit000 - 1    !
         !                                             ! ---------------------------------------- !
         IF( ln_rstart .AND. .NOT.l_1st_euler ) THEN         !* Restart: read in restart file
            IF(lwp) WRITE(numout,*) '          nit000-1 runoff forcing fields red in the restart file', lrxios
            CALL iom_get( numror, jpdom_auto, 'rnf_b'   , rnf_b                 )   ! before runoff
            CALL iom_get( numror, jpdom_auto, 'rnf_hc_b', rnf_tsc_b(:,:,jp_tem) )   ! before heat content of runoff
            CALL iom_get( numror, jpdom_auto, 'rnf_sc_b', rnf_tsc_b(:,:,jp_sal) )   ! before salinity content of runoff
         ELSE                                                !* no restart: set from nit000 values
            IF(lwp) WRITE(numout,*) '          nit000-1 runoff forcing fields set to nit000'
            CALL lbc_lnk( 'sbcrnf', rnf, 'T', 1.0_wp )
            rnf_b    (:,:)    = rnf    (:,:)
            rnf_tsc_b(:,:,:)  = rnf_tsc(:,:,:)
         ENDIF
      ENDIF
      !                                                ! ---------------------------------------- !
      IF( lrst_oce ) THEN                              !      Write in the ocean restart file     !
         !                                             ! ---------------------------------------- !
         !
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'sbcrnf : runoff forcing fields written in ocean restart file ',   &
            &                    'at it= ', kt,' date= ', ndastp
         IF(lwp) WRITE(numout,*) '~~~~'
         CALL iom_rstput( kt, nitrst, numrow, 'rnf_b'   , rnf                 )
         CALL iom_rstput( kt, nitrst, numrow, 'rnf_hc_b', rnf_tsc(:,:,jp_tem) )
         CALL iom_rstput( kt, nitrst, numrow, 'rnf_sc_b', rnf_tsc(:,:,jp_sal) )
      ENDIF
      !
#endif
   END SUBROUTINE sbc_rnf


   SUBROUTINE sbc_rnf_div( phdivn, Kmm )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_rnf  ***
      !!
      !! ** Purpose :   update the horizontal divergence with the runoff inflow
      !!
      !! ** Method  :
      !!                CAUTION : rnf is positive (inflow) decreasing the
      !!                          divergence and expressed in m/s
      !!
      !! ** Action  :   phdivn   decreased by the runoff inflow
      !!----------------------------------------------------------------------
      INTEGER                        , INTENT(in   ) ::   Kmm      ! ocean time level index
      REAL(wp), DIMENSION(A2D(1),jpk), INTENT(inout) ::   phdivn   ! horizontal divergence
      !!
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zfact     ! local scalar
      !!----------------------------------------------------------------------
      !
      zfact = 0.5_wp * r1_rho0
      !
      IF( ln_rnf_depth .OR. ln_rnf_depth_ini ) THEN      !==   runoff distributed over several levels   ==!
         IF( lk_linssh ) THEN    !* constant volume case : just apply the runoff input flow
            DO_2D( 1, 1, 1, 1 )
               DO jk = 1, nk_rnf(ji,jj)
#if defined key_RK3
                  phdivn(ji,jj,jk) = phdivn(ji,jj,jk) - rnf(ji,jj) * r1_rho0 / h_rnf(ji,jj)                    ! RK3: rnf forcing at n+1/2
#else
                  phdivn(ji,jj,jk) = phdivn(ji,jj,jk) - ( rnf(ji,jj) + rnf_b(ji,jj) ) * zfact / h_rnf(ji,jj)   ! MLF: rnf forcing at Kmm (n)
#endif
               END DO
            END_2D
         ELSE                    !* variable volume case
            DO_2D( 1, 1, 1, 1 )         ! compute the depth over which runoffs are distributed
               h_rnf(ji,jj) = 0._wp
               DO jk = 1, nk_rnf(ji,jj)
                  h_rnf(ji,jj)     = h_rnf(ji,jj) + e3t(ji,jj,jk,Kmm) ! recalculates h_rnf to be the depth in metres to the bottom of the relevant grid box
               END DO
            END_2D
            DO_2D( 1, 1, 1, 1 )         ! update phdivn
               DO jk = 1, nk_rnf(ji,jj)
#if defined key_RK3
                  phdivn(ji,jj,jk) = phdivn(ji,jj,jk) - rnf(ji,jj) * r1_rho0 / h_rnf(ji,jj)                    ! RK3: rnf forcing at n+1/2
#else
                  phdivn(ji,jj,jk) = phdivn(ji,jj,jk) - ( rnf(ji,jj) + rnf_b(ji,jj) ) * zfact / h_rnf(ji,jj)   ! MLF: rnf forcing at Kmm (n)
#endif
               END DO
            END_2D
         ENDIF
      ELSE                       !==   runoff put only at the surface   ==!
         DO_2D( 1, 1, 1, 1 )
            h_rnf (ji,jj)   = e3t(ji,jj,1,Kmm)        ! update h_rnf to be depth of top box
#if defined key_RK3
            phdivn(ji,jj,1) = phdivn(ji,jj,1) - rnf(ji,jj) * r1_rho0 / e3t(ji,jj,1,Kmm)                    ! RK3: rnf forcing at n+1/2
#else
            phdivn(ji,jj,1) = phdivn(ji,jj,1) - ( rnf(ji,jj) + rnf_b(ji,jj) ) * zfact / e3t(ji,jj,1,Kmm)   ! MLF: rnf forcing at Kmm (n)
#endif
         END_2D
      ENDIF
      !
   END SUBROUTINE sbc_rnf_div


   SUBROUTINE sbc_rnf_init( Kmm )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_rnf_init  ***
      !!
      !! ** Purpose :   Initialisation of the runoffs if (ln_rnf=T)
      !!
      !! ** Method  : - read the runoff namsbc_rnf namelist
      !!
      !! ** Action  : - read parameters
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) :: Kmm           ! ocean time level index
      CHARACTER(len=32) ::   rn_dep_file   ! runoff file name
      INTEGER           ::   ji, jj, jk, jm    ! dummy loop indices
      INTEGER           ::   ierror, inum  ! temporary integer
      INTEGER           ::   ios           ! Local integer output status for namelist read
      INTEGER           ::   nbrec         ! temporary integer
      REAL(wp)          ::   zacoef
      REAL(wp), DIMENSION(A2D(1),2) :: zrnfcl
      !!
      NAMELIST/namsbc_rnf/ cn_dir            , ln_rnf_depth, ln_rnf_tem, ln_rnf_sal, ln_rnf_icb,   &
         &                 sn_rnf, sn_cnf    , sn_i_rnf, sn_s_rnf    , sn_t_rnf  , sn_dep_rnf,   &
         &                 ln_rnf_mouth      , rn_hrnf     , rn_avt_rnf, rn_rfact,     &
         &                 ln_rnf_depth_ini  , rn_dep_max  , rn_rnf_max, nn_rnf_depth_file
      !!----------------------------------------------------------------------
      !
      !                                         !==  allocate runoff arrays
      IF( sbc_rnf_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'sbc_rnf_alloc : unable to allocate arrays' )
      !
      !                                   ! ============
      !                                   !   Namelist
      !                                   ! ============
      !
      READ_NML_REF(numnam,namsbc_rnf)
      READ_NML_CFG(numnam,namsbc_rnf)
      IF(lwm) WRITE ( numond, namsbc_rnf )
      !
      IF( .NOT. ln_rnf ) THEN                      ! no specific treatment in vicinity of river mouths
         ln_rnf_mouth  = .FALSE.                   ! default definition needed for example by sbc_ssr or by tra_adv_muscl
         ln_rnf_tem    = .FALSE.
         ln_rnf_sal    = .FALSE.
         ln_rnf_icb    = .FALSE.
         nkrnf         = 0
         rnf     (:,:) = 0.0_wp
         rnf_b   (:,:) = 0.0_wp
         rnfmsk  (:,:) = 0.0_wp
         rnfmsk_z(:)   = 0.0_wp
         RETURN
      ENDIF
      !
      !                                         ! Control print
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'sbc_rnf_init : runoff '
         WRITE(numout,*) '~~~~~~~~~~~~ '
         WRITE(numout,*) '   Namelist namsbc_rnf'
         WRITE(numout,*) '      specific river mouths treatment            ln_rnf_mouth = ', ln_rnf_mouth
         WRITE(numout,*) '      river mouth additional Kz                  rn_avt_rnf   = ', rn_avt_rnf
         WRITE(numout,*) '      depth of river mouth additional mixing     rn_hrnf      = ', rn_hrnf
         WRITE(numout,*) '      multiplicative factor for runoff           rn_rfact     = ', rn_rfact
      ENDIF
      !                                   ! ==================
      !                                   !   Type of runoff
      !                                   ! ==================
      !
      IF( .NOT. l_rnfcpl ) THEN
         ALLOCATE( sf_rnf(1), STAT=ierror )         ! Create sf_rnf structure (runoff inflow)
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   ==>>>   runoffs inflow read in a file'
         IF( ierror > 0 ) THEN
            CALL ctl_stop( 'sbc_rnf_init: unable to allocate sf_rnf structure' )   ;   RETURN
         ENDIF
         ALLOCATE( sf_rnf(1)%fnow(A2D(0),1)   )
         IF( sn_rnf%ln_tint ) ALLOCATE( sf_rnf(1)%fdta(A2D(0),1,2) )
         CALL fld_fill( sf_rnf, (/ sn_rnf /), cn_dir, 'sbc_rnf_init', 'read runoffs data', 'namsbc_rnf', no_print )
         !
         IF( ln_rnf_icb ) THEN                      ! Create (if required) sf_i_rnf structure
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) '          iceberg flux read in a file'
            ALLOCATE( sf_i_rnf(1), STAT=ierror  )
            IF( ierror > 0 ) THEN
               CALL ctl_stop( 'sbc_rnf_init: unable to allocate sf_i_rnf structure' )   ;   RETURN
            ENDIF
            ALLOCATE( sf_i_rnf(1)%fnow(A2D(0),1)   )
            IF( sn_i_rnf%ln_tint ) ALLOCATE( sf_i_rnf(1)%fdta(A2D(0),1,2) )
            CALL fld_fill (sf_i_rnf, (/ sn_i_rnf /), cn_dir, 'sbc_rnf_init', 'read iceberg flux data', 'namsbc_rnf' )
         ELSE
            DO_2D( 0, 0, 0, 0 )
               fwficb(ji,jj) = 0._wp
            END_2D
         ENDIF

      ENDIF
      !
      IF( ln_rnf_tem ) THEN                      ! Create (if required) sf_t_rnf structure
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   ==>>>   runoffs temperatures read in a file'
         ALLOCATE( sf_t_rnf(1), STAT=ierror  )
         IF( ierror > 0 ) THEN
            CALL ctl_stop( 'sbc_rnf_init: unable to allocate sf_t_rnf structure' )   ;   RETURN
         ENDIF
         ALLOCATE( sf_t_rnf(1)%fnow(A2D(0),1)   )
         IF( sn_t_rnf%ln_tint ) ALLOCATE( sf_t_rnf(1)%fdta(A2D(0),1,2) )
         CALL fld_fill (sf_t_rnf, (/ sn_t_rnf /), cn_dir, 'sbc_rnf_init', 'read runoff temperature data', 'namsbc_rnf', no_print )
      ENDIF
      !
      IF( ln_rnf_sal  ) THEN                     ! Create (if required) sf_s_rnf and sf_t_rnf structures
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   ==>>>   runoffs salinities read in a file'
         ALLOCATE( sf_s_rnf(1), STAT=ierror  )
         IF( ierror > 0 ) THEN
            CALL ctl_stop( 'sbc_rnf_init: unable to allocate sf_s_rnf structure' )   ;   RETURN
         ENDIF
         ALLOCATE( sf_s_rnf(1)%fnow(A2D(0),1)   )
         IF( sn_s_rnf%ln_tint ) ALLOCATE( sf_s_rnf(1)%fdta(A2D(0),1,2) )
         CALL fld_fill (sf_s_rnf, (/ sn_s_rnf /), cn_dir, 'sbc_rnf_init', 'read runoff salinity data', 'namsbc_rnf', no_print )
      ENDIF
      !
      IF( ln_rnf_depth ) THEN                    ! depth of runoffs set from a file
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   ==>>>   runoffs depth read in a file'
         rn_dep_file = TRIM( cn_dir )//TRIM( sn_dep_rnf%clname )
         IF( .NOT. sn_dep_rnf%ln_clim ) THEN   ;   WRITE(rn_dep_file, '(a,"_y",i4)' ) TRIM( rn_dep_file ), nyear    ! add year
            IF( sn_dep_rnf%clftyp == 'monthly' )   WRITE(rn_dep_file, '(a,"m",i2)'  ) TRIM( rn_dep_file ), nmonth   ! add month
         ENDIF
         CALL iom_open ( rn_dep_file, inum )                                                 ! open file
         CALL iom_get  ( inum, jpdom_global, sn_dep_rnf%clvar, h_rnf, kfill = jpfillcopy )   ! read the river mouth. no 0 on halos!
         CALL iom_close( inum )                                                              ! close file
         !
         DO_2D( 1, 1, 1, 1 )      ! set the number of level over which river runoffs are applied
            IF( h_rnf(ji,jj) > 0._wp ) THEN
               jk = 2
               DO WHILE ( jk < mbkt(ji,jj) .AND. gdept_0(ji,jj,jk) < h_rnf(ji,jj) ) ;  jk = jk + 1
               END DO
               nk_rnf(ji,jj) = jk
            ELSEIF( h_rnf(ji,jj) == -1._wp   ) THEN   ;  nk_rnf(ji,jj) = 1
            ELSEIF( h_rnf(ji,jj) == -999._wp ) THEN   ;  nk_rnf(ji,jj) = mbkt(ji,jj)
            ELSE
               CALL ctl_stop( 'sbc_rnf_init: runoff depth not positive, and not -999 or -1, rnf value in file fort.999'  )
               WRITE(999,*) 'ji, jj, h_rnf(ji,jj) :', ji, jj, h_rnf(ji,jj)
            ENDIF
         END_2D
         !
         DO_2D( 1, 1, 1, 1 )      ! set the associated depth
            h_rnf(ji,jj) = 0._wp
            DO jk = 1, nk_rnf(ji,jj)
               h_rnf(ji,jj) = h_rnf(ji,jj) + e3t(ji,jj,jk,Kmm)
            END DO
         END_2D
         !
      ELSE IF( ln_rnf_depth_ini ) THEN           ! runoffs applied at the surface
         !
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   ==>>>   depth of runoff computed once from max value of runoff'
         IF(lwp) WRITE(numout,*) '        max value of the runoff climatologie (over global domain) rn_rnf_max = ', rn_rnf_max
         IF(lwp) WRITE(numout,*) '        depth over which runoffs is spread                        rn_dep_max = ', rn_dep_max
         IF(lwp) WRITE(numout,*) '        create (=1) a runoff depth file or not (=0)      nn_rnf_depth_file  = ', nn_rnf_depth_file

         CALL iom_open( TRIM( sn_rnf%clname ), inum ) !  open runoff file
         nbrec = iom_getszuld( inum )
         zrnfcl(:,:,1) = 0._wp                                                            ! init the max to 0. in 1
         DO jm = 1, nbrec
            CALL iom_get( inum, jpdom_global, TRIM( sn_rnf%clvar ), zrnfcl(:,:,2), jm )   ! read the value in 2
            zrnfcl(:,:,1) = MAXVAL( zrnfcl(:,:,:), DIM=3 )                                ! store the maximum value in time in 1
         END DO
         CALL iom_close( inum )
         !
         zacoef = rn_dep_max / rn_rnf_max             ! coef of linear relation between runoff and its depth (150m for max of runoff)
         !
         WHERE( zrnfcl(:,:,1) > 0._wp ) ; h_rnf(:,:) = zacoef * zrnfcl(:,:,1)   ! compute depth for all runoffs
         ELSEWHERE                      ; h_rnf(:,:) = 1._wp
         ENDWHERE
         !
         DO_2D( 1, 1, 1, 1 )      ! take in account min depth of ocean rn_hmin
            IF( zrnfcl(ji,jj,1) > 0._wp ) THEN
               jk = mbkt(ji,jj)
               h_rnf(ji,jj) = MIN( h_rnf(ji,jj), gdept_0(ji,jj,jk ) )
            ENDIF
         END_2D
         !
         DO_2D( 1, 1, 1, 1 )      ! number of levels on which runoffs are distributed
            IF( zrnfcl(ji,jj,1) > 0._wp ) THEN
               jk = 2
               DO WHILE ( jk < mbkt(ji,jj) .AND. gdept_0(ji,jj,jk) < h_rnf(ji,jj) ) ;  jk = jk + 1
               END DO
               nk_rnf(ji,jj) = jk
            ELSE
               nk_rnf(ji,jj) = 1
            ENDIF
         END_2D
         !
         DO_2D( 1, 1, 1, 1 )      ! set the associated depth
            h_rnf(ji,jj) = 0._wp
            DO jk = 1, nk_rnf(ji,jj)
               h_rnf(ji,jj) = h_rnf(ji,jj) + e3t(ji,jj,jk,Kmm)
            END DO
         END_2D
         !
         IF( nn_rnf_depth_file == 1 ) THEN            !  save  output nb levels for runoff
            IF(lwp) WRITE(numout,*) '   ==>>>   create runoff depht file'
            CALL iom_open  ( TRIM( sn_dep_rnf%clname ), inum, ldwrt = .TRUE. )
            CALL iom_rstput( 0, 0, inum, 'rodepth', h_rnf )
            CALL iom_close ( inum )
         ENDIF
      ELSE                                            ! runoffs applied at the surface
         DO_2D( 1, 1, 1, 1 )
            nk_rnf(ji,jj) = 1
            h_rnf (ji,jj) = e3t(ji,jj,1,Kmm)
         END_2D
      ENDIF
      !
      DO_2D( 0, 0, 0, 0 )
         rnf(ji,jj) =  0._wp                          ! runoff initialisation
      END_2D
      DO_3D( 0, 0, 0, 0, 1, jpts )
         rnf_tsc(ji,jj,jk) = 0._wp                    ! runoffs temperature & salinty contents initilisation
      END_3D
         !
      !                                   ! ========================
      !                                   !   River mouth vicinity
      !                                   ! ========================
      !
      IF( ln_rnf_mouth ) THEN                   ! Specific treatment in vicinity of river mouths :
         !                                      !    - Increase Kz in surface layers ( rn_hrnf > 0 )
         !                                      !    - set to zero SSS damping (ln_ssr=T)
         !                                      !    - mixed upstream-centered (ln_traadv_cen2=T)
         !
         IF( ln_rnf_depth )   CALL ctl_warn( 'sbc_rnf_init: increased mixing turned on but effects may already',   &
            &                                              'be spread through depth by ln_rnf_depth'               )
         !
         nkrnf = 0                                  ! Number of level over which Kz increase
         IF( rn_hrnf > 0._wp ) THEN
            nkrnf = 2
            DO WHILE( nkrnf /= jpkm1 .AND. gdepw_1d(nkrnf+1) < rn_hrnf )   ;   nkrnf = nkrnf + 1
            END DO
            IF( l_sco )   CALL ctl_warn( 'sbc_rnf_init: number of levels over which Kz is increased is computed for zco...' )
         ENDIF
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   ==>>>   Specific treatment used in vicinity of river mouths :'
         IF(lwp) WRITE(numout,*) '             - Increase Kz in surface layers (if rn_hrnf > 0 )'
         IF(lwp) WRITE(numout,*) '               by ', rn_avt_rnf,' m2/s  over ', nkrnf, ' w-levels'
         IF(lwp) WRITE(numout,*) '             - set to zero SSS damping       (if ln_ssr=T)'
         IF(lwp) WRITE(numout,*) '             - mixed upstream-centered       (if ln_traadv_cen2=T)'
         !
         CALL rnf_mouth                             ! set river mouth mask
         !
      ELSE                                      ! No treatment at river mouths
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   ==>>>   No specific treatment at river mouths'
         DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
            rnfmsk(ji,jj) = 0._wp
         END_2D
         rnfmsk_z(:)   = 0._wp
         nkrnf = 0
      ENDIF
     !
   END SUBROUTINE sbc_rnf_init


   SUBROUTINE rnf_mouth
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE rnf_mouth  ***
      !!
      !! ** Purpose :   define the river mouths mask
      !!
      !! ** Method  :   read the river mouth mask (=0/1) in the river runoff
      !!                climatological file. Defined a given vertical structure.
      !!                CAUTION, the vertical structure is hard coded on the
      !!                first 5 levels.
      !!                This fields can be used to:
      !!                 - set an upstream advection scheme
      !!                   (ln_rnf_mouth=T and ln_traadv_cen2=T)
      !!                 - increase vertical on the top nn_krnf vertical levels
      !!                   at river runoff input grid point (nn_krnf>=2, see step.F90)
      !!                 - set to zero SSS restoring flux at river mouth grid points
      !!
      !! ** Action  :   rnfmsk   set to 1 at river runoff input, 0 elsewhere
      !!                rnfmsk_z vertical structure
      !!----------------------------------------------------------------------
      INTEGER            ::   inum        ! temporary integers
      CHARACTER(len=140) ::   cl_rnfile   ! runoff file name
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '   rnf_mouth : river mouth mask'
      IF(lwp) WRITE(numout,*) '   ~~~~~~~~~ '
      !
      cl_rnfile = TRIM( cn_dir )//TRIM( sn_cnf%clname )
      IF( .NOT. sn_cnf%ln_clim ) THEN   ;   WRITE(cl_rnfile, '(a,"_y",i4.4)' ) TRIM( cl_rnfile ), nyear    ! add year
         IF( sn_cnf%clftyp == 'monthly' )   WRITE(cl_rnfile, '(a,"m" ,i2.2)' ) TRIM( cl_rnfile ), nmonth   ! add month
      ENDIF
      !
      ! horizontal mask (read in NetCDF file)
      CALL iom_open ( cl_rnfile, inum )                                     ! open file
      CALL iom_get  ( inum, jpdom_global, sn_cnf%clvar, rnfmsk(A2D(0)) )    ! read the river mouth array
      CALL iom_close( inum )                                                ! close file
      !
      IF( l_clo_rnf )   CALL clo_rnf( rnfmsk )   ! closed sea inflow set as river mouth
      !
      rnfmsk_z(:)   = 0._wp                                       ! vertical structure
      rnfmsk_z(1)   = 1.0
      rnfmsk_z(2)   = 1.0                                         ! **********
      rnfmsk_z(3)   = 0.5                                         ! HARD CODED on the 5 first levels
      rnfmsk_z(4)   = 0.25                                        ! **********
      rnfmsk_z(5)   = 0.125
      !
      CALL lbc_lnk( 'rnf_mouth', rnfmsk, 'T', 1._wp )   ! needed by muscl scheme (traadv_muscl)
      !
   END SUBROUTINE rnf_mouth

   !!======================================================================
END MODULE sbcrnf
