MODULE fldread
   !!======================================================================
   !!                       ***  MODULE  fldread  ***
   !! Ocean forcing:  read input field for surface boundary condition
   !!=====================================================================
   !! History :  2.0  !  2006-06  (S. Masson, G. Madec)  Original code
   !!            3.0  !  2008-05  (S. Alderson)  Modified for Interpolation in memory from input grid to model grid
   !!            3.4  !  2013-10  (D. Delrosso, P. Oddo)  suppression of land point prior to interpolation
   !!                 !  12-2015  (J. Harle) Adding BDY on-the-fly interpolation
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   fld_read      : read input fields used for the computation of the surface boundary condition
   !!   fld_init      : initialization of field read
   !!   fld_def       : define the record(s) of the file and its name
   !!   fld_get       : read the data
   !!   fld_map       : read global data from file and map onto local data using a general mapping (use for open boundaries)
   !!   fld_rot       : rotate the vector fields onto the local grid direction
   !!   fld_clopn     : close/open the files
   !!   fld_fill      : fill the data structure with the associated information read in namelist
   !!   wgt_list      : manage the weights used for interpolation
   !!   wgt_print     : print the list of known weights
   !!   fld_weight    : create a WGT structure and fill in data from file, restructuring as required
   !!   apply_seaoverland : fill land with ocean values
   !!   seaoverland   : create shifted matrices for seaoverland application
   !!   fld_interp    : apply weights to input gridded data to create data on model grid
   !!   fld_filename  : define the filename according to a given date
   !!   ksec_week     : function returning seconds between 00h of the beginning of the week and half of the current time step
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constant
   USE sbc_oce        ! surface boundary conditions : fields
   USE geo2ocean      ! for vector rotation on to model grid
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE ioipsl  , ONLY : ymds2ju, ju2ymds   ! for calendar
   USE lib_mpp        ! MPP library
   USE lbclnk         ! ocean lateral boundary conditions (online interpolation case)
   
   IMPLICIT NONE
   PRIVATE   
 
   PUBLIC   fld_map    ! routine called by tides_init
   PUBLIC   fld_read, fld_fill   ! called by sbc... modules
   PUBLIC   fld_def, fld_clopn

   TYPE, PUBLIC ::   FLD_N      !: Namelist field informations
      CHARACTER(len = 256) ::   clname      ! generic name of the NetCDF flux file
      REAL(wp)             ::   freqh       ! frequency of each flux file
      CHARACTER(len = 34)  ::   clvar       ! generic name of the variable in the NetCDF flux file
      LOGICAL              ::   ln_tint     ! time interpolation or not (T/F)
      LOGICAL              ::   ln_clim     ! climatology or not (T/F)
      CHARACTER(len = 8)   ::   clftyp      ! type of data file 'daily', 'monthly' or yearly'
      CHARACTER(len = 256) ::   wname       ! generic name of a NetCDF weights file to be used, blank if not
      CHARACTER(len = 34)  ::   vcomp       ! symbolic component name if a vector that needs rotation
      !                                     ! a string starting with "U" or "V" for each component   
      !                                     ! chars 2 onwards identify which components go together  
      CHARACTER(len = 34)  ::   lname       ! generic name of a NetCDF land/sea mask file to be used, blank if not 
      !                                     ! 0=sea 1=land
   END TYPE FLD_N

   TYPE, PUBLIC ::   FLD        !: Input field related variables
      CHARACTER(len = 256)            ::   clrootname   ! generic name of the NetCDF file
      CHARACTER(len = 256)            ::   clname       ! current name of the NetCDF file
      REAL(wp)                        ::   freqh        ! frequency of each flux file
      CHARACTER(len = 34)             ::   clvar        ! generic name of the variable in the NetCDF flux file
      LOGICAL                         ::   ln_tint      ! time interpolation or not (T/F)
      REAL(wp)                        ::   rec_shft     ! record shift (from -0.5 to +0.5) when ln_tint=T
      LOGICAL                         ::   ln_clim      ! climatology or not (T/F)
      CHARACTER(len = 8)              ::   clftyp       ! type of data file 'daily', 'monthly' or yearly'
      CHARACTER(len = 1)              ::   cltype       ! nature of grid-points: T, U, V...
      REAL(wp)                        ::   zsgn         ! -1. the sign change across the north fold, =  1. otherwise
      INTEGER                         ::   num          ! iom id of the jpfld files to be read
      INTEGER , DIMENSION(2,2)        ::   nrec         ! before/after record (1: index, 2: second since Jan. 1st 00h of yr nit000)
      INTEGER                         ::   nbb          ! index of before values
      INTEGER                         ::   naa          ! index of after  values
      INTEGER , ALLOCATABLE, DIMENSION(:) ::   nrecsec  ! records time in seconds
      REAL(wp), POINTER, DIMENSION(:,:,:  ) ::   fnow   ! input fields interpolated to now time step
      REAL(wp), POINTER, DIMENSION(:,:,:,:) ::   fdta   ! 2 consecutive record of input fields
      CHARACTER(len = 256)            ::   wgtname      ! current name of the NetCDF weight file acting as a key
      !                                                 ! into the WGTLIST structure
      CHARACTER(len = 34)             ::   vcomp        ! symbolic name for a vector component that needs rotation
      LOGICAL, DIMENSION(2)           ::   rotn         ! flag to indicate whether before/after field has been rotated
      INTEGER                         ::   nreclast     ! last record to be read in the current file
      CHARACTER(len = 256)            ::   lsmname      ! current name of the NetCDF mask file acting as a key
      !                                                 ! 
      !                                                 ! Variables related to BDY
      INTEGER                         ::   igrd         !   grid type for bdy data
      INTEGER                         ::   ibdy         !   bdy set id number
      INTEGER, POINTER, DIMENSION(:)  ::   imap         !   Array of integer pointers to 1D arrays
      LOGICAL                         ::   ltotvel      !   total velocity or not (T/F)
      LOGICAL                         ::   lzint        !   T if it requires a vertical interpolation
   END TYPE FLD

!$AGRIF_DO_NOT_TREAT

   !! keep list of all weights variables so they're only read in once
   !! need to add AGRIF directives not to process this structure
   !! also need to force wgtname to include AGRIF nest number
   TYPE         ::   WGT        !: Input weights related variables
      CHARACTER(len = 256)                    ::   wgtname      ! current name of the NetCDF weight file
      INTEGER , DIMENSION(2)                  ::   ddims        ! shape of input grid
      INTEGER , DIMENSION(2)                  ::   botleft      ! top left corner of box in input grid containing 
      !                                                         ! current processor grid
      INTEGER , DIMENSION(2)                  ::   topright     ! top right corner of box 
      INTEGER                                 ::   jpiwgt       ! width of box on input grid
      INTEGER                                 ::   jpjwgt       ! height of box on input grid
      INTEGER                                 ::   numwgt       ! number of weights (4=bilinear, 16=bicubic)
      INTEGER                                 ::   nestid       ! for agrif, keep track of nest we're in
      INTEGER                                 ::   overlap      ! =0 when cyclic grid has no overlapping EW columns
      !                                                         ! =>1 when they have one or more overlapping columns      
      !                                                         ! =-1 not cyclic
      LOGICAL                                 ::   cyclic       ! east-west cyclic or not
      INTEGER,  DIMENSION(:,:,:), POINTER     ::   data_jpi     ! array of source integers
      INTEGER,  DIMENSION(:,:,:), POINTER     ::   data_jpj     ! array of source integers
      REAL(wp), DIMENSION(:,:,:), POINTER     ::   data_wgt     ! array of weights on model grid
      REAL(wp), DIMENSION(:,:,:), POINTER     ::   fly_dta      ! array of values on input grid
      REAL(wp), DIMENSION(:,:,:), POINTER     ::   col          ! temporary array for reading in columns
   END TYPE WGT

   INTEGER,     PARAMETER             ::   tot_wgts = 20
   TYPE( WGT ), DIMENSION(tot_wgts)   ::   ref_wgts     ! array of wgts
   INTEGER                            ::   nxt_wgt = 1  ! point to next available space in ref_wgts array
   INTEGER                            ::   nflag = 0
   REAL(wp), PARAMETER                ::   undeff_lsm = -999.00_wp

!$AGRIF_END_DO_NOT_TREAT

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "single_precision_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: fldread.F90 15023 2021-06-18 14:35:25Z gsamson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE fld_read( kt, kn_fsbc, sd, kit, pt_offset, Kmm )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE fld_read  ***
      !!                   
      !! ** Purpose :   provide at each time step the surface ocean fluxes
      !!                (momentum, heat, freshwater and runoff) 
      !!
      !! ** Method  :   READ each input fields in NetCDF files using IOM
      !!      and intepolate it to the model time-step.
      !!         Several assumptions are made on the input file:
      !!      blahblahblah....
      !!----------------------------------------------------------------------
      INTEGER  , INTENT(in   )               ::   kt        ! ocean time step
      INTEGER  , INTENT(in   )               ::   kn_fsbc   ! sbc computation period (in time step) 
      TYPE(FLD), INTENT(inout), DIMENSION(:) ::   sd        ! input field related variables
      INTEGER  , INTENT(in   ), OPTIONAL     ::   kit       ! subcycle timestep for timesplitting option
      REAL(wp) , INTENT(in   ), OPTIONAL     ::   pt_offset ! provide fields at time other than "now"
      INTEGER  , INTENT(in   ), OPTIONAL     ::   Kmm       ! ocean time level index
      !!
      INTEGER  ::   imf          ! size of the structure sd
      INTEGER  ::   jf           ! dummy indices
      INTEGER  ::   isecsbc      ! number of seconds between Jan. 1st 00h of nit000 year and the middle of sbc time step
      INTEGER  ::   ibb, iaa     ! shorter name for sd(jf)%nbb and sd(jf)%naa
      LOGICAL  ::   ll_firstcall ! true if this is the first call to fld_read for this set of fields
      REAL(wp) ::   zt_offset    ! local time offset variable
      REAL(wp) ::   ztinta       ! ratio applied to after  records when doing time interpolation
      REAL(wp) ::   ztintb       ! ratio applied to before records when doing time interpolation
      CHARACTER(LEN=1000) ::   clfmt  ! write format
      !!---------------------------------------------------------------------
      ll_firstcall = kt == nit000
      IF( PRESENT(kit) )   ll_firstcall = ll_firstcall .and. kit == 1

      IF( nn_components == jp_iam_sas ) THEN   ;   zt_offset = REAL( nn_fsbc, wp )
      ELSE                                     ;   zt_offset = 0.
      ENDIF
      IF( PRESENT(pt_offset) )   zt_offset = pt_offset

      ! Note that all varibles starting by nsec_* are shifted time by +1/2 time step to be centrered
      IF( PRESENT(kit) ) THEN   ! ignore kn_fsbc in this case
         isecsbc = nsec_year + nsec1jan000 + NINT( (     REAL(      kit,wp) + zt_offset ) * rn_Dt / REAL(nn_e,wp) )
      ELSE                      ! middle of sbc time step
         ! note: we use kn_fsbc-1 because nsec_year is defined at the middle of the current time step
         isecsbc = nsec_year + nsec1jan000 + NINT( ( 0.5*REAL(kn_fsbc-1,wp) + zt_offset ) * rn_Dt )
      ENDIF
      imf = SIZE( sd )
      !
      IF( ll_firstcall ) THEN                      ! initialization
         DO jf = 1, imf 
            IF( TRIM(sd(jf)%clrootname) == 'NOT USED' )   CYCLE
            CALL fld_init( isecsbc, sd(jf) )       ! read each before field (put them in after as they will be swapped)
         END DO
         IF( lwp ) CALL wgt_print()                ! control print
      ENDIF
      !                                            ! ====================================== !
      IF( MOD( kt-1, kn_fsbc ) == 0 ) THEN         ! update field at each kn_fsbc time-step !
         !                                         ! ====================================== !
         !
         DO jf = 1, imf                            ! ---   loop over field   --- !
            !
            IF( TRIM(sd(jf)%clrootname) == 'NOT USED' )   CYCLE
            CALL fld_update( isecsbc, sd(jf), Kmm )
            !
         END DO                                    ! --- end loop over field --- !

         CALL fld_rot( kt, sd )                    ! rotate vector before/now/after fields if needed

         DO jf = 1, imf                            ! ---   loop over field   --- !
            !
            IF( TRIM(sd(jf)%clrootname) == 'NOT USED' )   CYCLE
            !
            ibb = sd(jf)%nbb   ;   iaa = sd(jf)%naa
            !
            IF( sd(jf)%ln_tint ) THEN              ! temporal interpolation
               IF(lwp .AND. ( kt - nit000 <= 30 .OR. nitend - kt <= 30 ) ) THEN 
                  clfmt = "('   fld_read: var ', a, ' kt = ', i8, ' (', f9.4,' days), Y/M/D = ', i4.4,'/', i2.2,'/', i2.2," //   &
                     &    "', records b/a: ', i6.4, '/', i6.4, ' (days ', f9.4,'/', f9.4, ')')"
                  WRITE(numout, clfmt)  TRIM( sd(jf)%clvar ), kt, REAL(isecsbc,wp)/rday, nyear, nmonth, nday,   &            
                     & sd(jf)%nrec(1,ibb), sd(jf)%nrec(1,iaa), REAL(sd(jf)%nrec(2,ibb),wp)/rday, REAL(sd(jf)%nrec(2,iaa),wp)/rday
                  IF( zt_offset /= 0._wp )   WRITE(numout, *) '      zt_offset is : ', zt_offset
               ENDIF
               ! temporal interpolation weights
               ztinta =  REAL( isecsbc - sd(jf)%nrec(2,ibb), wp ) / REAL( sd(jf)%nrec(2,iaa) - sd(jf)%nrec(2,ibb), wp )
               ztintb =  1. - ztinta
               sd(jf)%fnow(:,:,:) = ztintb * sd(jf)%fdta(:,:,:,ibb) + ztinta * sd(jf)%fdta(:,:,:,iaa)
            ELSE   ! nothing to do...
               IF(lwp .AND. ( kt - nit000 <= 30 .OR. nitend - kt <= 30 ) ) THEN
                  clfmt = "('   fld_read: var ', a, ' kt = ', i8,' (', f9.4,' days), Y/M/D = ', i4.4,'/', i2.2,'/', i2.2," //   &
                     &    "', record: ', i6.4, ' (days ', f9.4, ' <-> ', f9.4, ')')"
                  WRITE(numout, clfmt) TRIM(sd(jf)%clvar), kt, REAL(isecsbc,wp)/rday, nyear, nmonth, nday,    &
                     &                 sd(jf)%nrec(1,iaa), REAL(sd(jf)%nrec(2,ibb),wp)/rday, REAL(sd(jf)%nrec(2,iaa),wp)/rday
               ENDIF
            ENDIF
            !
            IF( kt == nitend - kn_fsbc + 1 )   CALL iom_close( sd(jf)%num )   ! Close the input files

         END DO                                    ! --- end loop over field --- !
         !
      ENDIF
      !
   END SUBROUTINE fld_read


   SUBROUTINE fld_init( ksecsbc, sdjf )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE fld_init  ***
      !!
      !! ** Purpose :  - first call(s) to fld_def to define before values
      !!               - open file
      !!----------------------------------------------------------------------
      INTEGER  , INTENT(in   ) ::   ksecsbc   ! 
      TYPE(FLD), INTENT(inout) ::   sdjf         ! input field related variables
      !!---------------------------------------------------------------------
      !
      IF( nflag == 0 )   nflag = -HUGE(0)
      !
      CALL fld_def( sdjf )
      IF( sdjf%ln_tint .AND. ksecsbc < sdjf%nrecsec(1) )   CALL fld_def( sdjf, ldprev = .TRUE. )
      !
      CALL fld_clopn( sdjf )
      sdjf%nrec(:,sdjf%naa) = (/ 1, nflag /)  ! default definition to force flp_update to read the file.
      !
   END SUBROUTINE fld_init


   SUBROUTINE fld_update( ksecsbc, sdjf, Kmm )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE fld_update  ***
      !!
      !! ** Purpose : Compute
      !!              if sdjf%ln_tint = .TRUE.
      !!                  nrec(:,iaa): record number and its time (nrec(:,ibb) is obtained from nrec(:,iaa) when swapping)
      !!              if sdjf%ln_tint = .FALSE.
      !!                  nrec(1,iaa): record number
      !!                  nrec(2,ibb) and nrec(2,iaa): time of the beginning and end of the record
      !!----------------------------------------------------------------------
      INTEGER  ,           INTENT(in   ) ::   ksecsbc   ! 
      TYPE(FLD),           INTENT(inout) ::   sdjf      ! input field related variables
      INTEGER  , OPTIONAL, INTENT(in   ) ::   Kmm    ! ocean time level index
      !
      INTEGER  ::   ja           ! end of this record (in seconds)
      INTEGER  ::   ibb, iaa     ! shorter name for sdjf%nbb and sdjf%naa
      !!----------------------------------------------------------------------
      ibb = sdjf%nbb   ;   iaa = sdjf%naa
      !
      IF( ksecsbc > sdjf%nrec(2,iaa) ) THEN     ! --> we need to update after data
        
         ! find where is the new after record... (it is not necessary sdjf%nrec(1,iaa)+1 )
         ja = sdjf%nrec(1,iaa)
         DO WHILE ( ksecsbc >= sdjf%nrecsec(ja) .AND. ja < sdjf%nreclast )   ! Warning: make sure ja <= sdjf%nreclast in this test
            ja = ja + 1
         END DO
         IF( ksecsbc > sdjf%nrecsec(ja) )   ja = ja + 1   ! in case ksecsbc > sdjf%nrecsec(sdjf%nreclast)

         ! if ln_tint and if the new after is not ja+1, we need also to update after data before the swap
         ! so, after the swap, sdjf%nrec(2,ibb) will still be the closest value located just before ksecsbc
         IF( sdjf%ln_tint .AND. ( ja > sdjf%nrec(1,iaa) + 1 .OR. sdjf%nrec(2,iaa) == nflag ) ) THEN
            sdjf%nrec(:,iaa) = (/ ja-1, sdjf%nrecsec(ja-1) /)   ! update nrec(:,iaa) with before information
            CALL fld_get( sdjf, Kmm )                           ! read after data that will be used as before data
         ENDIF
            
         ! if after is in the next file...
         IF( ja > sdjf%nreclast ) THEN
            
            CALL fld_def( sdjf )
            IF( ksecsbc > sdjf%nrecsec(sdjf%nreclast) )   CALL fld_def( sdjf, ldnext = .TRUE. )
            CALL fld_clopn( sdjf )           ! open next file
            
            ! find where is after in this new file
            ja = 1
            DO WHILE ( ksecsbc > sdjf%nrecsec(ja) .AND. ja < sdjf%nreclast )
               ja = ja + 1
            END DO
            IF( ksecsbc > sdjf%nrecsec(ja) )   ja = ja + 1   ! in case ksecsbc > sdjf%nrecsec(sdjf%nreclast)
            
            IF( ja > sdjf%nreclast ) THEN
               CALL ctl_stop( "STOP", "fld_def: need next-next file? we should not be there... file: "//TRIM(sdjf%clrootname) )
            ENDIF
            
            ! if ln_tint and if after is not the first record, we must (potentially again) update after data before the swap
            IF( sdjf%ln_tint .AND. ja > 1 ) THEN
               IF( sdjf%nrecsec(0) /= nflag ) THEN                    ! no trick used: after file is not the current file
                  sdjf%nrec(:,iaa) = (/ ja-1, sdjf%nrecsec(ja-1) /)   ! update nrec(:,iaa) with before information
                  CALL fld_get( sdjf, Kmm )                           ! read after data that will be used as before data
               ENDIF
            ENDIF
            
         ENDIF

         IF( sdjf%ln_tint ) THEN                                ! Swap data
            sdjf%nbb = sdjf%naa                                 !    swap indices
            sdjf%naa = 3 - sdjf%naa                             !    = 2(1) if naa == 1(2)
         ELSE                                                   ! No swap
            sdjf%nrec(:,ibb) = (/ ja-1, sdjf%nrecsec(ja-1) /)   !    only for print 
         ENDIF
            
         ! read new after data
         sdjf%nrec(:,sdjf%naa) = (/ ja, sdjf%nrecsec(ja) /)     ! update nrec(:,naa) as it is used by fld_get
         CALL fld_get( sdjf, Kmm )                              ! read after data (with nrec(:,naa) informations)
        
      ENDIF
      !
   END SUBROUTINE fld_update


   SUBROUTINE fld_get( sdjf, Kmm )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE fld_get  ***
      !!
      !! ** Purpose :   read the data
      !!----------------------------------------------------------------------
      TYPE(FLD),           INTENT(inout) ::   sdjf   ! input field related variables
      INTEGER  , OPTIONAL, INTENT(in   ) ::   Kmm    ! ocean time level index
      !
      INTEGER ::   ipk      ! number of vertical levels of sdjf%fdta ( 2D: ipk=1 ; 3D: ipk=jpk )
      INTEGER ::   iaa      ! shorter name for sdjf%naa
      INTEGER ::   iw       ! index into wgts array
      INTEGER ::   idvar    ! variable ID
      INTEGER ::   idmspc   ! number of spatial dimensions
      REAL(wp)                            ::   zsgn        ! sign used in the call to lbc_lbk called by iom_get
      REAL(wp), DIMENSION(:,:,:), POINTER ::   dta_alias   ! short cut
      !!---------------------------------------------------------------------
      iaa = sdjf%naa
      !
      IF( sdjf%ln_tint ) THEN   ;   dta_alias => sdjf%fdta(:,:,:,iaa)
      ELSE                      ;   dta_alias => sdjf%fnow(:,:,:    )
      ENDIF
      ipk = SIZE( dta_alias, 3 )
      !
      IF( LEN_TRIM(sdjf%vcomp) > 0 ) THEN   ;   zsgn = 1._wp        ! geographical vectors -> sign change done later when rotating
      ELSE                                  ;   zsgn = sdjf%zsgn
      ENDIF
      !
      IF( ASSOCIATED(sdjf%imap) ) THEN              ! BDY case 
         CALL fld_map( sdjf%num, sdjf%clvar, dta_alias(:,:,:), sdjf%nrec(1,iaa),   &
            &          sdjf%imap, sdjf%igrd, sdjf%ibdy, sdjf%ltotvel, sdjf%lzint, Kmm )
      ELSE IF( LEN(TRIM(sdjf%wgtname)) > 0 ) THEN   ! On-the-fly interpolation
         CALL wgt_list( sdjf, iw )
         CALL fld_interp( sdjf%num, sdjf%clvar, iw, ipk, dta_alias(:,:,:), sdjf%nrec(1,iaa), sdjf%lsmname )
         CALL lbc_lnk( 'fldread', dta_alias(:,:,:), sdjf%cltype, zsgn, kfillmode = jpfillcopy )
      ELSE                                          ! default case
         idvar  = iom_varid( sdjf%num, sdjf%clvar )
         idmspc = iom_file ( sdjf%num )%ndims( idvar )
         IF( iom_file( sdjf%num )%luld( idvar ) )   idmspc = idmspc - 1   ! id of the last spatial dimension
         CALL iom_get( sdjf%num,  jpdom_global, sdjf%clvar, dta_alias(:,:,:), sdjf%nrec(1,iaa),   &
            &          sdjf%cltype, zsgn, kfill = jpfillcopy )
      ENDIF
      !
      sdjf%rotn(iaa) = .false.   ! vector not yet rotated
      !
   END SUBROUTINE fld_get

   
   SUBROUTINE fld_map( knum, cdvar, pdta, krec, kmap, kgrd, kbdy, ldtotvel, ldzint, Kmm )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE fld_map  ***
      !!
      !! ** Purpose :   read global data from file and map onto local data
      !!                using a general mapping (for open boundaries)
      !!----------------------------------------------------------------------
      INTEGER                   , INTENT(in   ) ::   knum         ! stream number
      CHARACTER(LEN=*)          , INTENT(in   ) ::   cdvar        ! variable name
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pdta         ! bdy output field on model grid
      INTEGER                   , INTENT(in   ) ::   krec         ! record number to read (ie time slice)
      INTEGER , DIMENSION(:)    , INTENT(in   ) ::   kmap         ! global-to-local bdy mapping indices
      ! optional variables used for vertical interpolation:
      INTEGER, OPTIONAL         , INTENT(in   ) ::   kgrd         ! grid type (t, u, v)
      INTEGER, OPTIONAL         , INTENT(in   ) ::   kbdy         ! bdy number
      LOGICAL, OPTIONAL         , INTENT(in   ) ::   ldtotvel     ! true if total ( = barotrop + barocline) velocity
      LOGICAL, OPTIONAL         , INTENT(in   ) ::   ldzint       ! true if 3D variable requires a vertical interpolation
      INTEGER, OPTIONAL         , INTENT(in   ) ::   Kmm          ! ocean time level index 
      !!
      INTEGER                                   ::   ipi          ! length of boundary data on local process
      INTEGER                                   ::   ipj          ! length of dummy dimension ( = 1 )
      INTEGER                                   ::   ipk          ! number of vertical levels of pdta ( 2D: ipk=1 ; 3D: ipk=jpk )
      INTEGER                                   ::   ipkb         ! number of vertical levels in boundary data file
      INTEGER                                   ::   idvar        ! variable ID
      INTEGER                                   ::   indims       ! number of dimensions of the variable
      INTEGER, DIMENSION(4)                     ::   idimsz       ! size of variable dimensions 
      REAL(wp)                                  ::   zfv          ! fillvalue 
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   ::   zz_read      ! work space for global boundary data
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   ::   zdta_read    ! work space local data requiring vertical interpolation
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   ::   zdta_read_z  ! work space local data requiring vertical interpolation
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   ::   zdta_read_dz ! work space local data requiring vertical interpolation
      CHARACTER(LEN=1),DIMENSION(3)             ::   cltype
      LOGICAL                                   ::   lluld        ! is the variable using the unlimited dimension
      LOGICAL                                   ::   llzint       ! local value of ldzint
      !!---------------------------------------------------------------------
      !
      cltype = (/'t','u','v'/)
      !
      ipi = SIZE( pdta, 1 )
      ipj = SIZE( pdta, 2 )   ! must be equal to 1
      ipk = SIZE( pdta, 3 )
      !
      llzint = .FALSE.
      IF( PRESENT(ldzint) )   llzint = ldzint
      !
      idvar = iom_varid( knum, cdvar, kndims = indims, kdimsz = idimsz, lduld = lluld  )
      IF( indims == 4 .OR. ( indims == 3 .AND. .NOT. lluld ) ) THEN   ;   ipkb = idimsz(3)   ! xy(zl)t or xy(zl)
      ELSE                                                            ;   ipkb = 1           ! xy or xyt
      ENDIF
      !
      ALLOCATE( zz_read( idimsz(1), idimsz(2), ipkb ) )  ! ++++++++ !!! this can be very big...         
      !
      IF( ipk == 1 ) THEN

         IF( ipkb /= 1 ) CALL ctl_stop( 'fld_map : we must have ipkb = 1 to read surface data' )
         CALL iom_get ( knum, jpdom_unknown, cdvar, zz_read(:,:,1), krec )   ! call iom_get with a 2D file
         CALL fld_map_core( zz_read, kmap, pdta )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Do we include something here to adjust barotropic velocities !
      ! in case of a depth difference between bdy files and          !
      ! bathymetry in the case ln_totvel = .false. and ipkb>0?       !
      ! [as the enveloping and parital cells could change H]         !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ELSE
         !
         CALL iom_get ( knum, jpdom_unknown, cdvar, zz_read(:,:,:), krec )   ! call iom_get with a 3D file
         !
         IF( ipkb /= ipk .OR. llzint ) THEN   ! boundary data not on model vertical grid : vertical interpolation
            !
            IF( ipk == jpk .AND. iom_varid(knum,'gdep'//cltype(kgrd)) /= -1 .AND. iom_varid(knum,'e3'//cltype(kgrd)) /= -1 ) THEN
               
               ALLOCATE( zdta_read(ipi,ipj,ipkb), zdta_read_z(ipi,ipj,ipkb), zdta_read_dz(ipi,ipj,ipkb) )
                
               CALL fld_map_core( zz_read, kmap, zdta_read )
               CALL iom_get ( knum, jpdom_unknown, 'gdep'//cltype(kgrd), zz_read )   ! read only once? Potential temporal evolution?
               CALL fld_map_core( zz_read, kmap, zdta_read_z )
               CALL iom_get ( knum, jpdom_unknown,   'e3'//cltype(kgrd), zz_read )   ! read only once? Potential temporal evolution?
               CALL fld_map_core( zz_read, kmap, zdta_read_dz )
               
               CALL iom_getatt(knum, '_FillValue', zfv, cdvar=cdvar )
               CALL fld_bdy_interp(zdta_read, zdta_read_z, zdta_read_dz, pdta, kgrd, kbdy, zfv, ldtotvel, Kmm)
               DEALLOCATE( zdta_read, zdta_read_z, zdta_read_dz )
               
            ELSE
               IF( ipk /= jpk ) CALL ctl_stop( 'fld_map : this should be an impossible case...' )
               WRITE(ctmp1,*) 'fld_map : vertical interpolation for bdy variable '//TRIM(cdvar)//' requires ' 
               IF( iom_varid(knum, 'gdep'//cltype(kgrd)) == -1 ) CALL ctl_stop( ctmp1//'gdep'//cltype(kgrd)//' variable' )
               IF( iom_varid(knum,   'e3'//cltype(kgrd)) == -1 ) CALL ctl_stop( ctmp1//  'e3'//cltype(kgrd)//' variable' )

            ENDIF
            !
         ELSE                            ! bdy data assumed to be the same levels as bdy variables
            !
            CALL fld_map_core( zz_read, kmap, pdta )
            !
         ENDIF   ! ipkb /= ipk
      ENDIF   ! ipk == 1
      
      DEALLOCATE( zz_read )

   END SUBROUTINE fld_map

     
   SUBROUTINE fld_map_core( pdta_read, kmap, pdta_bdy )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE fld_map_core  ***
      !!
      !! ** Purpose :  inner core of fld_map
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pdta_read    ! global boundary data
      INTEGER,  DIMENSION(:    ), INTENT(in   ) ::   kmap         ! global-to-local bdy mapping indices
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pdta_bdy     ! bdy output field on model grid
      !!
      INTEGER,  DIMENSION(3) ::   idim_read,  idim_bdy            ! arrays dimensions
      INTEGER                ::   ji, jj, jk, jb                  ! loop counters
      INTEGER                ::   im1
      !!---------------------------------------------------------------------
      !
      idim_read = SHAPE( pdta_read )
      idim_bdy  = SHAPE( pdta_bdy  )
      !
      ! in all cases: idim_bdy(2) == 1 .AND. idim_read(1) * idim_read(2) == idim_bdy(1)
      ! structured BDY with rimwidth > 1                     : idim_read(2) == rimwidth /= 1
      ! structured BDY with rimwidth == 1 or unstructured BDY: idim_read(2) == 1
      !
      IF( idim_read(2) > 1 ) THEN    ! structured BDY with rimwidth > 1  
         DO jk = 1, idim_bdy(3)
            DO jb = 1, idim_bdy(1)
               im1 = kmap(jb) - 1
               jj = im1 / idim_read(1) + 1
               ji = MOD( im1, idim_read(1) ) + 1
               pdta_bdy(jb,1,jk) =  pdta_read(ji,jj,jk)
            END DO
         END DO
      ELSE
         DO jk = 1, idim_bdy(3)
            DO jb = 1, idim_bdy(1)   ! horizontal remap of bdy data on the local bdy 
               pdta_bdy(jb,1,jk) = pdta_read(kmap(jb),1,jk)
            END DO
         END DO
      ENDIF
      
   END SUBROUTINE fld_map_core
   
   SUBROUTINE fld_bdy_interp(pdta_read, pdta_read_z, pdta_read_dz, pdta, kgrd, kbdy, pfv, ldtotvel, Kmm )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE fld_bdy_interp  ***
      !!
      !! ** Purpose :   on the fly vertical interpolation to allow the use of
      !!                boundary data from non-native vertical grid
      !!----------------------------------------------------------------------
      USE bdy_oce, ONLY:  idx_bdy         ! indexing for map <-> ij transformation

      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pdta_read       ! data read in bdy file
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pdta_read_z     ! depth of the data read in bdy file
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pdta_read_dz    ! thickness of the levels in bdy file
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pdta            ! output field on model grid (2 dimensional)
      REAL(wp)                  , INTENT(in   ) ::   pfv             ! fillvalue of the data read in bdy file
      LOGICAL                   , INTENT(in   ) ::   ldtotvel        ! true if toal ( = barotrop + barocline) velocity
      INTEGER                   , INTENT(in   ) ::   kgrd            ! grid type (t, u, v)
      INTEGER                   , INTENT(in   ) ::   kbdy            ! bdy number
      INTEGER, OPTIONAL         , INTENT(in   ) ::   Kmm             ! ocean time level index
      !!
      INTEGER                  ::   ipi                 ! length of boundary data on local process
      INTEGER                  ::   ipkb                ! number of vertical levels in boundary data file
      INTEGER                  ::   ipkmax              ! number of vertical levels in boundary data file where no mask
      INTEGER                  ::   jb, ji, jj, jk, jkb ! loop counters
      REAL(wp)                 ::   zcoef, zi           ! 
      REAL(wp)                 ::   ztrans, ztrans_new  ! transports
      REAL(wp), DIMENSION(jpk) ::   zdepth, zdhalf      ! level and half-level depth
      !!---------------------------------------------------------------------
     
      ipi  = SIZE( pdta, 1 )
      ipkb = SIZE( pdta_read, 3 )
      
      DO jb = 1, ipi
         ji = idx_bdy(kbdy)%nbi(jb,kgrd)
         jj = idx_bdy(kbdy)%nbj(jb,kgrd)
         !
         ! --- max jk where input data /= FillValue --- !
         ipkmax = 1
         DO jkb = 2, ipkb
            IF( pdta_read(jb,1,jkb) /= pfv )   ipkmax = MAX( ipkmax, jkb )
         END DO
         !
         ! --- calculate depth at t,u,v points --- !
         SELECT CASE( kgrd )                         
         CASE(1)            ! depth of T points:
            zdepth(:) = gdept(ji,jj,:,Kmm)
         CASE(2)            ! depth of U points: we must not use gdept_n as we don't want to do a communication
            !                 --> copy what is done for gdept_n in domvvl...
            zdhalf(1) = 0.0_wp
            zdepth(1) = 0.5_wp * e3uw(ji,jj,1,Kmm)
            DO jk = 2, jpk                               ! vertical sum
               !    zcoef = umask - wumask    ! 0 everywhere tmask = wmask, ie everywhere expect at jk = mikt
               !                              ! 1 everywhere from mbkt to mikt + 1 or 1 (if no isf)
               !                              ! 0.5 where jk = mikt     
               !!gm ???????   BUG ?  gdept_n as well as gde3w_n  does not include the thickness of ISF ??
               zcoef = ( umask(ji,jj,jk) - wumask(ji,jj,jk) )
               zdhalf(jk) = zdhalf(jk-1) + e3u(ji,jj,jk-1,Kmm)
               zdepth(jk) =          zcoef  * ( zdhalf(jk  ) + 0.5_wp * e3uw(ji,jj,jk,Kmm))  &
                  &         + (1._wp-zcoef) * ( zdepth(jk-1) +          e3uw(ji,jj,jk,Kmm))
            END DO
         CASE(3)            ! depth of V points: we must not use gdept_n as we don't want to do a communication
            !                 --> copy what is done for gdept_n in domvvl...
            zdhalf(1) = 0.0_wp
            zdepth(1) = 0.5_wp * e3vw(ji,jj,1,Kmm)
            DO jk = 2, jpk                               ! vertical sum
               !    zcoef = vmask - wvmask    ! 0 everywhere tmask = wmask, ie everywhere expect at jk = mikt
               !                              ! 1 everywhere from mbkt to mikt + 1 or 1 (if no isf)
               !                              ! 0.5 where jk = mikt     
               !!gm ???????   BUG ?  gdept_n as well as gde3w_n  does not include the thickness of ISF ??
               zcoef = ( vmask(ji,jj,jk) - wvmask(ji,jj,jk) )
               zdhalf(jk) = zdhalf(jk-1) + e3v(ji,jj,jk-1,Kmm)
               zdepth(jk) =          zcoef  * ( zdhalf(jk  ) + 0.5_wp * e3vw(ji,jj,jk,Kmm))  &
                     + (1._wp-zcoef) * ( zdepth(jk-1) +          e3vw(ji,jj,jk,Kmm))
            END DO
         END SELECT
         !         
         ! --- interpolate bdy data on the model grid --- !
         DO jk = 1, jpk
            IF(     zdepth(jk) <= pdta_read_z(jb,1,1)      ) THEN   ! above the first level of external data
               pdta(jb,1,jk) = pdta_read(jb,1,1)
            ELSEIF( zdepth(jk) >  pdta_read_z(jb,1,ipkmax) ) THEN   ! below the last level of external data /= FillValue
               pdta(jb,1,jk) = pdta_read(jb,1,ipkmax)
            ELSE                                                    ! inbetween: vertical interpolation between jkb & jkb+1
               DO jkb = 1, ipkmax-1
                  IF( ( ( zdepth(jk) - pdta_read_z(jb,1,jkb) ) * ( zdepth(jk) - pdta_read_z(jb,1,jkb+1) ) ) <= 0._wp ) THEN ! linear interpolation between 2 levels
                     zi = ( zdepth(jk) - pdta_read_z(jb,1,jkb) ) / ( pdta_read_z(jb,1,jkb+1) - pdta_read_z(jb,1,jkb) )
                     pdta(jb,1,jk) = pdta_read(jb,1,jkb) + zi * ( pdta_read(jb,1,jkb+1) - pdta_read(jb,1,jkb) )
                  ENDIF
               END DO
            ENDIF
         END DO   ! jpk
         !
      END DO   ! ipi

      ! --- mask data and adjust transport --- !
      SELECT CASE( kgrd )                         

      CASE(1)                                 ! mask data (probably unecessary)
         DO jb = 1, ipi
            ji = idx_bdy(kbdy)%nbi(jb,kgrd)
            jj = idx_bdy(kbdy)%nbj(jb,kgrd)
            DO jk = 1, jpk                      
               pdta(jb,1,jk) = pdta(jb,1,jk) * tmask(ji,jj,jk)
            END DO
         END DO
         
      CASE(2)                                 ! adjust the U-transport term
         DO jb = 1, ipi
            ji = idx_bdy(kbdy)%nbi(jb,kgrd)
            jj = idx_bdy(kbdy)%nbj(jb,kgrd)
            ztrans = 0._wp
            DO jkb = 1, ipkb                              ! calculate transport on input grid
               IF( pdta_read(jb,1,jkb) /= pfv )   ztrans = ztrans + pdta_read(jb,1,jkb) * pdta_read_dz(jb,1,jkb)
            ENDDO
            ztrans_new = 0._wp
            DO jk = 1, jpk                                ! calculate transport on model grid
               ztrans_new = ztrans_new +      pdta(jb,1,jk ) * e3u(ji,jj,jk,Kmm ) * umask(ji,jj,jk)
            ENDDO
            DO jk = 1, jpk                                ! make transport correction
               IF(ldtotvel) THEN ! bdy data are total velocity so adjust bt transport term to match input data
                  pdta(jb,1,jk) = ( pdta(jb,1,jk) + ( ztrans - ztrans_new ) * r1_hu(ji,jj,Kmm) ) * umask(ji,jj,jk)
               ELSE              ! we're just dealing with bc velocity so bt transport term should sum to zero
                  pdta(jb,1,jk) =   pdta(jb,1,jk) + (  0._wp - ztrans_new ) * r1_hu(ji,jj,Kmm)   * umask(ji,jj,jk)
               ENDIF
            ENDDO
         ENDDO

      CASE(3)                                 ! adjust the V-transport term
         DO jb = 1, ipi
            ji = idx_bdy(kbdy)%nbi(jb,kgrd)
            jj = idx_bdy(kbdy)%nbj(jb,kgrd)
            ztrans = 0._wp
            DO jkb = 1, ipkb                              ! calculate transport on input grid
               IF( pdta_read(jb,1,jkb) /= pfv )   ztrans = ztrans + pdta_read(jb,1,jkb) * pdta_read_dz(jb,1,jkb)
            ENDDO
            ztrans_new = 0._wp
            DO jk = 1, jpk                                ! calculate transport on model grid
               ztrans_new = ztrans_new +      pdta(jb,1,jk ) * e3v(ji,jj,jk,Kmm ) * vmask(ji,jj,jk)
            ENDDO
            DO jk = 1, jpk                                ! make transport correction
               IF(ldtotvel) THEN ! bdy data are total velocity so adjust bt transport term to match input data
                  pdta(jb,1,jk) = ( pdta(jb,1,jk) + ( ztrans - ztrans_new ) * r1_hv(ji,jj,Kmm) ) * vmask(ji,jj,jk)
               ELSE              ! we're just dealing with bc velocity so bt transport term should sum to zero
                  pdta(jb,1,jk) =   pdta(jb,1,jk) + (  0._wp - ztrans_new ) * r1_hv(ji,jj,Kmm)   * vmask(ji,jj,jk)
               ENDIF
            ENDDO
         ENDDO
      END SELECT
      
   END SUBROUTINE fld_bdy_interp


   SUBROUTINE fld_rot( kt, sd )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE fld_rot  ***
      !!
      !! ** Purpose :   Vector fields may need to be rotated onto the local grid direction
      !!----------------------------------------------------------------------
      INTEGER                , INTENT(in   ) ::   kt   ! ocean time step
      TYPE(FLD), DIMENSION(:), INTENT(inout) ::   sd   ! input field related variables
      !
      INTEGER ::   ju, jv, jk, jn  ! loop indices
      INTEGER ::   imf             ! size of the structure sd
      INTEGER ::   ill             ! character length
      INTEGER ::   iv              ! indice of V component
      CHARACTER (LEN=100)          ::   clcomp       ! dummy weight name
      REAL(wp), DIMENSION(jpi,jpj) ::   utmp, vtmp   ! temporary arrays for vector rotation
      REAL(wp), DIMENSION(:,:,:), POINTER ::   dta_u, dta_v    ! short cut
      !!---------------------------------------------------------------------
      !
      !! (sga: following code should be modified so that pairs arent searched for each time
      !
      imf = SIZE( sd )
      DO ju = 1, imf
         IF( TRIM(sd(ju)%clrootname) == 'NOT USED' )   CYCLE
         ill = LEN_TRIM( sd(ju)%vcomp )
         DO jn = 2-COUNT((/sd(ju)%ln_tint/)), 2
            IF( ill > 0 .AND. .NOT. sd(ju)%rotn(jn) ) THEN   ! find vector rotations required             
               IF( sd(ju)%vcomp(1:1) == 'U' ) THEN      ! east-west component has symbolic name starting with 'U'
                  ! look for the north-south component which has same symbolic name but with 'U' replaced with 'V'
                  clcomp = 'V' // sd(ju)%vcomp(2:ill)   ! works even if ill == 1
                  iv = -1
                  DO jv = 1, imf
                     IF( TRIM(sd(jv)%clrootname) == 'NOT USED' )   CYCLE
                     IF( TRIM(sd(jv)%vcomp) == TRIM(clcomp) )   iv = jv
                  END DO
                  IF( iv > 0 ) THEN   ! fields ju and iv are two components which need to be rotated together
                     IF( sd(ju)%ln_tint ) THEN   ;   dta_u => sd(ju)%fdta(:,:,:,jn)   ;   dta_v => sd(iv)%fdta(:,:,:,jn) 
                     ELSE                        ;   dta_u => sd(ju)%fnow(:,:,:   )   ;   dta_v => sd(iv)%fnow(:,:,:   )
                     ENDIF
                     DO jk = 1, SIZE( sd(ju)%fnow, 3 )
                        CALL rot_rep( dta_u(:,:,jk), dta_v(:,:,jk), 'T', 'en->i', utmp(:,:) )
                        CALL rot_rep( dta_u(:,:,jk), dta_v(:,:,jk), 'T', 'en->j', vtmp(:,:) )
                        dta_u(:,:,jk) = utmp(:,:)   ;   dta_v(:,:,jk) = vtmp(:,:)
                     END DO
                     sd(ju)%rotn(jn) = .TRUE.               ! vector was rotated 
                     IF( lwp .AND. kt == nit000 )   WRITE(numout,*)   &
                        &   'fld_read: vector pair ('//TRIM(sd(ju)%clvar)//', '//TRIM(sd(iv)%clvar)//') rotated on to model grid'
                  ENDIF
               ENDIF
            ENDIF
         END DO
       END DO
      !
   END SUBROUTINE fld_rot


   SUBROUTINE fld_def( sdjf, ldprev, ldnext )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE fld_def  ***
      !!
      !! ** Purpose :   define the record(s) of the file and its name
      !!----------------------------------------------------------------------
      TYPE(FLD)        , INTENT(inout) ::   sdjf       ! input field related variables
      LOGICAL, OPTIONAL, INTENT(in   ) ::   ldprev     ! 
      LOGICAL, OPTIONAL, INTENT(in   ) ::   ldnext     ! 
      !
      INTEGER  :: jt
      INTEGER  :: idaysec               ! number of seconds in 1 day = NINT(rday)
      INTEGER  :: iyr, imt, idy, isecwk
      INTEGER  :: indexyr, indexmt
      INTEGER  :: ireclast, irecshft
      INTEGER  :: ishift, istart
      INTEGER, DIMENSION(2)  :: isave
      REAL(wp) :: zfreqs
      LOGICAL  :: llprev, llnext, llstop
      LOGICAL  :: llprevmt, llprevyr
      LOGICAL  :: llnextmt, llnextyr
      !!----------------------------------------------------------------------
      idaysec = NINT(rday)
      !
      IF( PRESENT(ldprev) ) THEN   ;   llprev = ldprev
      ELSE                         ;   llprev = .FALSE.
      ENDIF
      IF( PRESENT(ldnext) ) THEN   ;   llnext = ldnext
      ELSE                         ;   llnext = .FALSE.
      ENDIF

      ! current file parameters
      IF( sdjf%clftyp(1:4) == 'week' ) THEN         ! find the day of the beginning of the current week
         isecwk = ksec_week( sdjf%clftyp(6:8) )     ! seconds between the beginning of the week and half of current time step
         llprevmt = isecwk > nsec_month             ! longer time since beginning of the current week than the current month
         llprevyr = llprevmt .AND. nmonth == 1
         iyr = nyear  - COUNT((/llprevyr/))
         imt = nmonth - COUNT((/llprevmt/)) + 12 * COUNT((/llprevyr/))
         idy = nday + nmonth_len(nmonth-1) * COUNT((/llprevmt/)) - isecwk / idaysec
         isecwk = nsec_year - isecwk                ! seconds between 00h jan 1st of current year and current week beginning
      ELSE
         iyr = nyear
         imt = nmonth
         idy = nday
         isecwk  = 0
      ENDIF

      ! previous file parameters
      IF( llprev ) THEN
         IF( sdjf%clftyp(1:4) == 'week'    ) THEN   ! find the day of the beginning of previous week
            isecwk = ksec_week( sdjf%clftyp(6:8) )  ! seconds between the beginning of the week and half of current time step
            isecwk = isecwk + 7 * idaysec           ! seconds between the beginning of previous week and half of the time step
            llprevmt = isecwk > nsec_month          ! longer time since beginning of the previous week than the current month
            llprevyr = llprevmt .AND. nmonth == 1
            iyr = nyear  - COUNT((/llprevyr/))
            imt = nmonth - COUNT((/llprevmt/)) + 12 * COUNT((/llprevyr/))
            idy = nday + nmonth_len(nmonth-1) * COUNT((/llprevmt/)) - isecwk / idaysec
            isecwk = nsec_year - isecwk             ! seconds between 00h jan 1st of current year and previous week beginning
         ELSE
            idy = nday   - COUNT((/ sdjf%clftyp == 'daily'                 /))
            imt = nmonth - COUNT((/ sdjf%clftyp == 'monthly' .OR. idy == 0 /))
            iyr = nyear  - COUNT((/ sdjf%clftyp == 'yearly'  .OR. imt == 0 /))
            IF( idy == 0 ) idy = nmonth_len(imt)
            IF( imt == 0 ) imt = 12
            isecwk = 0
         ENDIF
      ENDIF

      ! next file parameters
      IF( llnext ) THEN
         IF( sdjf%clftyp(1:4) == 'week'    ) THEN   ! find the day of the beginning of next week
            isecwk = ksec_week( sdjf%clftyp(6:8) )  ! seconds between the beginning of the week and half of current time step
            isecwk = 7 * idaysec - isecwk           ! seconds between half of the time step and the beginning of next week
            llnextmt = isecwk >= ( nmonth_len(nmonth)*idaysec - nsec_month )   ! larger than the seconds to the end of the month
            llnextyr = llnextmt .AND. nmonth == 12
            iyr = nyear  + COUNT((/llnextyr/))
            imt = nmonth + COUNT((/llnextmt/)) - 12 * COUNT((/llnextyr/))
            idy = nday - nmonth_len(nmonth) * COUNT((/llnextmt/)) + isecwk / idaysec + 1
            isecwk = nsec_year + isecwk             ! seconds between 00h jan 1st of current year and next week beginning
         ELSE
            idy = nday   + COUNT((/ sdjf%clftyp == 'daily'                                 /))
            imt = nmonth + COUNT((/ sdjf%clftyp == 'monthly' .OR. idy > nmonth_len(nmonth) /))
            iyr = nyear  + COUNT((/ sdjf%clftyp == 'yearly'  .OR. imt == 13                /))
            IF( idy > nmonth_len(nmonth) )   idy = 1
            IF( imt == 13                )   imt = 1
            isecwk = 0
         ENDIF
      ENDIF
      !
      ! find the last record to be read -> update sdjf%nreclast
      indexyr = iyr - nyear + 1                 ! which  year are we looking for? previous(0), current(1) or next(2)?
      indexmt = imt + 12 * ( indexyr - 1 )      ! which month are we looking for (relatively to current year)? 
      !
      ! Last record to be read in the current file
      ! Predefine the number of record in the file according of its type.
      ! We could compare this number with the number of records in the file and make a stop if the 2 numbers do not match...
      ! However this would be much less fexible (e.g. for tests) and will force to rewite input files according to nleapy...
      IF    ( NINT(sdjf%freqh) == -12 ) THEN            ;   ireclast = 1    ! yearly mean: consider only 1 record
      ELSEIF( NINT(sdjf%freqh) ==  -1 ) THEN                                ! monthly mean:
         IF(     sdjf%clftyp      == 'monthly' ) THEN   ;   ireclast = 1    !  consider that the file has  1 record
         ELSE                                           ;   ireclast = 12   !  consider that the file has 12 record
         ENDIF
      ELSE                                                                  ! higher frequency mean (in hours)
         IF(     sdjf%clftyp      == 'monthly' ) THEN   ;   ireclast = NINT( 24. * REAL(nmonth_len(indexmt), wp) / sdjf%freqh )
         ELSEIF( sdjf%clftyp(1:4) == 'week'    ) THEN   ;   ireclast = NINT( 24. * 7.                            / sdjf%freqh )
         ELSEIF( sdjf%clftyp      == 'daily'   ) THEN   ;   ireclast = NINT( 24.                                 / sdjf%freqh )
         ELSE                                           ;   ireclast = NINT( 24. * REAL( nyear_len(indexyr), wp) / sdjf%freqh )
         ENDIF
      ENDIF

      sdjf%nreclast = ireclast
      ! Allocate arrays for beginning/middle/end of each record (seconds since Jan. 1st 00h of nit000 year)
      IF( ALLOCATED(sdjf%nrecsec) )   DEALLOCATE( sdjf%nrecsec )
      ALLOCATE( sdjf%nrecsec( 0:ireclast ) )
      !
      IF    ( NINT(sdjf%freqh) == -12 ) THEN                                     ! yearly mean and yearly file
         SELECT CASE( indexyr )
         CASE(0)   ;   sdjf%nrecsec(0) = nsec1jan000 - nyear_len( 0 ) * idaysec
         CASE(1)   ;   sdjf%nrecsec(0) = nsec1jan000
         CASE(2)   ;   sdjf%nrecsec(0) = nsec1jan000 + nyear_len( 1 ) * idaysec
         ENDSELECT
         sdjf%nrecsec(1) = sdjf%nrecsec(0) + nyear_len( indexyr ) * idaysec
      ELSEIF( NINT(sdjf%freqh) ==  -1 ) THEN                                     ! monthly mean:
         IF(     sdjf%clftyp      == 'monthly' ) THEN                            !    monthly file
            sdjf%nrecsec(0   ) = nsec1jan000 + nmonth_beg(indexmt  )
            sdjf%nrecsec(1   ) = nsec1jan000 + nmonth_beg(indexmt+1)
         ELSE                                                                    !    yearly  file
            ishift = 12 * ( indexyr - 1 )
            sdjf%nrecsec(0:12) = nsec1jan000 + nmonth_beg(1+ishift:13+ishift)
         ENDIF
      ELSE                                                                       ! higher frequency mean (in hours)
         IF(     sdjf%clftyp      == 'monthly' ) THEN   ;   istart = nsec1jan000 + nmonth_beg(indexmt)
         ELSEIF( sdjf%clftyp(1:4) == 'week'    ) THEN   ;   istart = nsec1jan000 + isecwk
         ELSEIF( sdjf%clftyp      == 'daily'   ) THEN   ;   istart = nsec1jan000 + nmonth_beg(indexmt) + ( idy - 1 ) * idaysec
         ELSEIF( indexyr          == 0         ) THEN   ;   istart = nsec1jan000 - nyear_len( 0 ) * idaysec
         ELSEIF( indexyr          == 2         ) THEN   ;   istart = nsec1jan000 + nyear_len( 1 ) * idaysec
         ELSE                                           ;   istart = nsec1jan000
         ENDIF
         zfreqs = sdjf%freqh * rhhmm * rmmss
         DO jt = 0, sdjf%nreclast
            sdjf%nrecsec(jt) = istart + NINT( zfreqs * REAL(jt,wp) )
         END DO
      ENDIF
      !
      IF( sdjf%ln_tint ) THEN   ! record time defined in the middle of the record, computed using an implementation
                                ! of the rounded average that is valid over the full integer range
         irecshft = NINT( sdjf%rec_shft * (sdjf%nrecsec(1) - sdjf%nrecsec(0)) )
         sdjf%nrecsec(1:sdjf%nreclast) = sdjf%nrecsec(0:sdjf%nreclast-1) / 2 + sdjf%nrecsec(1:sdjf%nreclast) / 2 + &
            & MAX( MOD( sdjf%nrecsec(0:sdjf%nreclast-1), 2 ), MOD( sdjf%nrecsec(1:sdjf%nreclast), 2 ) ) + irecshft
      END IF
      !
      sdjf%clname = fld_filename( sdjf, idy, imt, iyr )
      !
   END SUBROUTINE fld_def

   
   SUBROUTINE fld_clopn( sdjf )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE fld_clopn  ***
      !!
      !! ** Purpose :   close/open the files
      !!----------------------------------------------------------------------
      TYPE(FLD)        , INTENT(inout) ::   sdjf       ! input field related variables
      !
      INTEGER  :: isave
      LOGICAL  :: llprev, llnext, llstop
      !!----------------------------------------------------------------------
      !
      llprev = sdjf%nrecsec(sdjf%nreclast) < nsec000_1jan000   ! file ends before the beginning of the job -> file may not exist
      llnext = sdjf%nrecsec(      1      ) > nsecend_1jan000   ! file begins after the end of the job -> file may not exist 

      llstop = sdjf%ln_clim .OR. .NOT. ( llprev .OR. llnext )

      IF( sdjf%num <= 0 .OR. .NOT. sdjf%ln_clim  ) THEN
         IF( sdjf%num > 0 )   CALL iom_close( sdjf%num )   ! close file if already open
         CALL iom_open( sdjf%clname, sdjf%num, ldstop = llstop, ldiof = LEN_TRIM(sdjf%wgtname) > 0 )
      ENDIF
      !
      IF( sdjf%num <= 0 .AND. .NOT. llstop ) THEN   ! file not found but we do accept this...
         !
         IF( llprev ) THEN   ! previous file does not exist : go back to current and accept to read only the first record
            CALL ctl_warn('previous file: '//TRIM(sdjf%clname)//' not found -> go back to current year/month/week/day file')
            isave = sdjf%nrecsec(sdjf%nreclast)   ! save previous file info
            CALL fld_def( sdjf )                  ! go back to current file
            sdjf%nreclast = 1                     ! force to use only the first record (do as if other were not existing...)
         ENDIF
         !
         IF( llnext ) THEN   ! next     file does not exist : go back to current and accept to read only the last  record 
            CALL ctl_warn('next file: '//TRIM(sdjf%clname)//' not found -> go back to current year/month/week/day file')
            isave = sdjf%nrecsec(1)               ! save next file info
            CALL fld_def( sdjf )                  ! go back to current file
         ENDIF
         ! -> read "last" record but keep record info from the first record of next file
         sdjf%nrecsec(  sdjf%nreclast  ) = isave
         sdjf%nrecsec(0:sdjf%nreclast-1) = nflag
         !
         CALL iom_open( sdjf%clname, sdjf%num, ldiof = LEN_TRIM(sdjf%wgtname) > 0 )   
         !
      ENDIF
      !
   END SUBROUTINE fld_clopn


   SUBROUTINE fld_fill( sdf, sdf_n, cdir, cdcaller, cdtitle, cdnam, knoprint )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE fld_fill  ***
      !!
      !! ** Purpose :   fill the data structure (sdf) with the associated information 
      !!              read in namelist (sdf_n) and control print
      !!----------------------------------------------------------------------
      TYPE(FLD)  , DIMENSION(:)          , INTENT(inout) ::   sdf        ! structure of input fields (file informations, fields read)
      TYPE(FLD_N), DIMENSION(:)          , INTENT(in   ) ::   sdf_n      ! array of namelist information structures
      CHARACTER(len=*)                   , INTENT(in   ) ::   cdir       ! Root directory for location of flx files
      CHARACTER(len=*)                   , INTENT(in   ) ::   cdcaller   ! name of the calling routine
      CHARACTER(len=*)                   , INTENT(in   ) ::   cdtitle    ! description of the calling routine 
      CHARACTER(len=*)                   , INTENT(in   ) ::   cdnam      ! name of the namelist from which sdf_n comes
      INTEGER                  , OPTIONAL, INTENT(in   ) ::   knoprint   ! no calling routine information printed
      !
      INTEGER  ::   jf   ! dummy indices
      !!---------------------------------------------------------------------
      !
      DO jf = 1, SIZE(sdf)
         sdf(jf)%clrootname = sdf_n(jf)%clname
         IF( TRIM(sdf_n(jf)%clname) /= 'NOT USED' )   sdf(jf)%clrootname = TRIM( cdir )//sdf(jf)%clrootname
         sdf(jf)%clname     = "not yet defined"
         sdf(jf)%freqh      = sdf_n(jf)%freqh
         sdf(jf)%clvar      = sdf_n(jf)%clvar
         sdf(jf)%ln_tint    = sdf_n(jf)%ln_tint
         sdf(jf)%rec_shft   = 0._wp
         IF( sdf(jf)%ln_tint ) THEN
            IF( sdf(jf)%clvar(1:1) == '+' ) THEN
               sdf(jf)%rec_shft = +0.5_wp
               sdf(jf)%clvar = sdf(jf)%clvar(2:)
            ELSE IF( sdf(jf)%clvar(1:1) == '-' ) THEN
               sdf(jf)%rec_shft = -0.5_wp
               sdf(jf)%clvar = sdf(jf)%clvar(2:)
            ENDIF
         ENDIF
         sdf(jf)%ln_clim    = sdf_n(jf)%ln_clim
         sdf(jf)%clftyp     = sdf_n(jf)%clftyp
         sdf(jf)%cltype     = 'T'   ! by default don't do any call to lbc_lnk in iom_get
         sdf(jf)%zsgn       = 1.    ! by default don't do change signe across the north fold
         sdf(jf)%num        = -1
         sdf(jf)%nbb        = 1  ! start with before data in 1
         sdf(jf)%naa        = 2  ! start with after  data in 2
         sdf(jf)%wgtname    = " "
         IF( LEN( TRIM(sdf_n(jf)%wname) ) > 0 )   sdf(jf)%wgtname = TRIM( cdir )//sdf_n(jf)%wname
         sdf(jf)%lsmname = " "
         IF( LEN( TRIM(sdf_n(jf)%lname) ) > 0 )   sdf(jf)%lsmname = TRIM( cdir )//sdf_n(jf)%lname
         sdf(jf)%vcomp      = sdf_n(jf)%vcomp
         sdf(jf)%rotn(:)    = .TRUE.   ! pretend to be rotated -> won't try to rotate data before the first call to fld_get
         IF( sdf(jf)%clftyp(1:4) == 'week' .AND. nn_leapy == 0  )   &
            &   CALL ctl_stop('fld_clopn: weekly file ('//TRIM(sdf(jf)%clrootname)//') needs nn_leapy = 1')
         IF( sdf(jf)%clftyp(1:4) == 'week' .AND. sdf(jf)%ln_clim )   &
            &   CALL ctl_stop('fld_clopn: weekly file ('//TRIM(sdf(jf)%clrootname)//') needs ln_clim = .FALSE.')
         sdf(jf)%nreclast   = -1 ! Set to non zero default value to avoid errors, is updated to meaningful value during fld_clopn
         sdf(jf)%igrd       = 0
         sdf(jf)%ibdy       = 0
         sdf(jf)%imap       => NULL()
         sdf(jf)%ltotvel    = .FALSE.
         sdf(jf)%lzint      = .FALSE.
      END DO
      !
      IF(lwp) THEN      ! control print
         WRITE(numout,*)
         IF( .NOT.PRESENT( knoprint) ) THEN
            WRITE(numout,*) TRIM( cdcaller )//' : '//TRIM( cdtitle )
            WRITE(numout,*) (/ ('~', jf = 1, LEN_TRIM( cdcaller ) ) /)
         ENDIF
         WRITE(numout,*) '   fld_fill : fill data structure with information from namelist '//TRIM( cdnam )
         WRITE(numout,*) '   ~~~~~~~~'
         WRITE(numout,*) '      list of files and frequency (>0: in hours ; <0 in months)'
         DO jf = 1, SIZE(sdf)
            WRITE(numout,*) '      root filename: '  , TRIM( sdf(jf)%clrootname ), '   variable name: ', TRIM( sdf(jf)%clvar )
            WRITE(numout,*) '         frequency: '      ,       sdf(jf)%freqh       ,   &
               &                  '   time interp: '    ,       sdf(jf)%ln_tint     ,   &
               &                  '   climatology: '    ,       sdf(jf)%ln_clim
            WRITE(numout,*) '         weights: '        , TRIM( sdf(jf)%wgtname    ),   &
               &                  '   pairing: '        , TRIM( sdf(jf)%vcomp      ),   &
               &                  '   data type: '      ,       sdf(jf)%clftyp      ,   &
               &                  '   land/sea mask:'   , TRIM( sdf(jf)%lsmname    )
            call flush(numout)
         END DO
      ENDIF
      !
   END SUBROUTINE fld_fill


   SUBROUTINE wgt_list( sd, kwgt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE wgt_list  ***
      !!
      !! ** Purpose :   search array of WGTs and find a weights file entry,
      !!                or return a new one adding it to the end if new entry.
      !!                the weights data is read in and restructured (fld_weight)
      !!----------------------------------------------------------------------
      TYPE( FLD ), INTENT(in   ) ::   sd        ! field with name of weights file
      INTEGER    , INTENT(  out) ::   kwgt      ! index of weights
      !
      INTEGER ::   kw, nestid   ! local integer
      !!----------------------------------------------------------------------
      !
      !! search down linked list 
      !! weights filename is either present or we hit the end of the list
      !
      !! because agrif nest part of filenames are now added in iom_open
      !! to distinguish between weights files on the different grids, need to track
      !! nest number explicitly
      nestid = 0
#if defined key_agrif
      nestid = Agrif_Fixed()
#endif
      DO kw = 1, nxt_wgt-1
         IF( ref_wgts(kw)%wgtname == sd%wgtname .AND. &
             ref_wgts(kw)%nestid  == nestid) THEN
            kwgt = kw
            RETURN
         ENDIF
      END DO
      kwgt = nxt_wgt
      CALL fld_weight( sd )
      !
   END SUBROUTINE wgt_list


   SUBROUTINE wgt_print( )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE wgt_print  ***
      !!
      !! ** Purpose :   print the list of known weights
      !!----------------------------------------------------------------------
      INTEGER ::   kw   !
      !!----------------------------------------------------------------------
      !
      DO kw = 1, nxt_wgt-1
         WRITE(numout,*) 'weight file:  ',TRIM(ref_wgts(kw)%wgtname)
         WRITE(numout,*) '      ddims:  ',ref_wgts(kw)%ddims(1),ref_wgts(kw)%ddims(2)
         WRITE(numout,*) '     numwgt:  ',ref_wgts(kw)%numwgt
         WRITE(numout,*) '     jpiwgt:  ',ref_wgts(kw)%jpiwgt
         WRITE(numout,*) '     jpjwgt:  ',ref_wgts(kw)%jpjwgt
         WRITE(numout,*) '    botleft:  ',ref_wgts(kw)%botleft
         WRITE(numout,*) '   topright:  ',ref_wgts(kw)%topright
         IF( ref_wgts(kw)%cyclic ) THEN
            WRITE(numout,*) '       cyclical'
            IF( ref_wgts(kw)%overlap > 0 ) WRITE(numout,*) '              with overlap of ', ref_wgts(kw)%overlap
         ELSE
            WRITE(numout,*) '       not cyclical'
         ENDIF
         IF( ASSOCIATED(ref_wgts(kw)%data_wgt) )  WRITE(numout,*) '       allocated'
      END DO
      !
   END SUBROUTINE wgt_print


   SUBROUTINE fld_weight( sd )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE fld_weight  ***
      !!
      !! ** Purpose :   create a new WGT structure and fill in data from file,
      !!              restructuring as required
      !!----------------------------------------------------------------------
      TYPE( FLD ), INTENT(in) ::   sd   ! field with name of weights file
      !!
      INTEGER ::   ji,jj,jn   ! dummy loop indices
      INTEGER ::   inum       ! local logical unit
      INTEGER ::   id         ! local variable id
      INTEGER ::   ipk        ! local vertical dimension
      INTEGER ::   zwrap      ! local integer
      LOGICAL ::   cyclical   ! 
      CHARACTER (len=5) ::   clname   !
      INTEGER , DIMENSION(4) ::   ddims
      INTEGER                ::   isrc
      REAL(dp), DIMENSION(jpi,jpj) ::   data_tmp
      !!----------------------------------------------------------------------
      !
      IF( nxt_wgt > tot_wgts ) THEN
        CALL ctl_stop("fld_weight: weights array size exceeded, increase tot_wgts")
      ENDIF
      !
      !! new weights file entry, add in extra information
      !! a weights file represents a 2D grid of a certain shape, so we assume that the current
      !! input data file is representative of all other files to be opened and processed with the
      !! current weights file

      !! get data grid dimensions
      id = iom_varid( sd%num, sd%clvar, ddims )

      !! now open the weights file
      CALL iom_open ( sd%wgtname, inum )   ! interpolation weights
      IF( inum > 0 ) THEN

         !! determine whether we have an east-west cyclic grid
         !! from global attribute called "ew_wrap" in the weights file
         !! note that if not found, iom_getatt returns -999 and cyclic with no overlap is assumed
         !! since this is the most common forcing configuration

         CALL iom_getatt(inum, 'ew_wrap', zwrap)
         IF( zwrap >= 0 ) THEN
            cyclical = .TRUE.
         ELSE IF( zwrap == -999 ) THEN
            cyclical = .TRUE.
            zwrap = 0
         ELSE
            cyclical = .FALSE.
         ENDIF

         ref_wgts(nxt_wgt)%ddims(1) = ddims(1)
         ref_wgts(nxt_wgt)%ddims(2) = ddims(2)
         ref_wgts(nxt_wgt)%wgtname = sd%wgtname
         ref_wgts(nxt_wgt)%overlap = zwrap
         ref_wgts(nxt_wgt)%cyclic = cyclical
         ref_wgts(nxt_wgt)%nestid = 0
#if defined key_agrif
         ref_wgts(nxt_wgt)%nestid = Agrif_Fixed()
#endif
         !! weights file is stored as a set of weights (wgt01->wgt04 or wgt01->wgt16)
         !! for each weight wgtNN there is an integer array srcNN which gives the point in
         !! the input data grid which is to be multiplied by the weight
         !! they are both arrays on the model grid so the result of the multiplication is
         !! added into an output array on the model grid as a running sum

         !! two possible cases: bilinear (4 weights) or bicubic (16 weights)
         id = iom_varid(inum, 'src05', ldstop=.FALSE.)
         IF( id <= 0 ) THEN   ;   ref_wgts(nxt_wgt)%numwgt = 4
         ELSE                 ;   ref_wgts(nxt_wgt)%numwgt = 16
         ENDIF

         ALLOCATE( ref_wgts(nxt_wgt)%data_jpi(Nis0:Nie0,Njs0:Nje0,4) )
         ALLOCATE( ref_wgts(nxt_wgt)%data_jpj(Nis0:Nie0,Njs0:Nje0,4) )
         ALLOCATE( ref_wgts(nxt_wgt)%data_wgt(Nis0:Nie0,Njs0:Nje0,ref_wgts(nxt_wgt)%numwgt) )

         DO jn = 1,4
            WRITE(clname,'(a3,i2.2)') 'src',jn
            CALL iom_get ( inum, jpdom_global, clname, data_tmp(:,:), cd_type = 'Z' )   !  no call to lbc_lnk
            DO_2D( 0, 0, 0, 0 )
               isrc = NINT(data_tmp(ji,jj)) - 1
               ref_wgts(nxt_wgt)%data_jpi(ji,jj,jn) = 1 + MOD(isrc,  ref_wgts(nxt_wgt)%ddims(1))
               ref_wgts(nxt_wgt)%data_jpj(ji,jj,jn) = 1 +     isrc / ref_wgts(nxt_wgt)%ddims(1)
            END_2D
         END DO

         DO jn = 1, ref_wgts(nxt_wgt)%numwgt
            WRITE(clname,'(a3,i2.2)') 'wgt',jn
            CALL iom_get ( inum, jpdom_global, clname, data_tmp(:,:), cd_type = 'Z' )   !  no call to lbc_lnk
            DO_2D( 0, 0, 0, 0 )
               ref_wgts(nxt_wgt)%data_wgt(ji,jj,jn) = data_tmp(ji,jj)
            END_2D
         END DO
         CALL iom_close (inum)
 
         ! find min and max indices in grid
         ref_wgts(nxt_wgt)%botleft( 1) = MINVAL(ref_wgts(nxt_wgt)%data_jpi(:,:,:))
         ref_wgts(nxt_wgt)%botleft( 2) = MINVAL(ref_wgts(nxt_wgt)%data_jpj(:,:,:))
         ref_wgts(nxt_wgt)%topright(1) = MAXVAL(ref_wgts(nxt_wgt)%data_jpi(:,:,:))
         ref_wgts(nxt_wgt)%topright(2) = MAXVAL(ref_wgts(nxt_wgt)%data_jpj(:,:,:))

         ! and therefore dimensions of the input box
         ref_wgts(nxt_wgt)%jpiwgt = ref_wgts(nxt_wgt)%topright(1) - ref_wgts(nxt_wgt)%botleft(1) + 1
         ref_wgts(nxt_wgt)%jpjwgt = ref_wgts(nxt_wgt)%topright(2) - ref_wgts(nxt_wgt)%botleft(2) + 1

         ! shift indexing of source grid
         ref_wgts(nxt_wgt)%data_jpi(:,:,:) = ref_wgts(nxt_wgt)%data_jpi(:,:,:) - ref_wgts(nxt_wgt)%botleft(1) + 1
         ref_wgts(nxt_wgt)%data_jpj(:,:,:) = ref_wgts(nxt_wgt)%data_jpj(:,:,:) - ref_wgts(nxt_wgt)%botleft(2) + 1

         ! create input grid, give it a halo to allow gradient calculations
         ! SA: +3 stencil is a patch to avoid out-of-bound computation in some configuration. 
         ! a more robust solution will be given in next release
         ipk =  SIZE(sd%fnow, 3)
         ALLOCATE( ref_wgts(nxt_wgt)%fly_dta(ref_wgts(nxt_wgt)%jpiwgt+3, ref_wgts(nxt_wgt)%jpjwgt+3 ,ipk) )
         IF( ref_wgts(nxt_wgt)%cyclic ) ALLOCATE( ref_wgts(nxt_wgt)%col(1,ref_wgts(nxt_wgt)%jpjwgt+3,ipk) )
         !
         nxt_wgt = nxt_wgt + 1
         !
      ELSE 
         CALL ctl_stop( '    fld_weight : unable to read the file ' )
      ENDIF
      !
   END SUBROUTINE fld_weight


   SUBROUTINE apply_seaoverland( clmaskfile, zfieldo, jpi1_lsm, jpi2_lsm, jpj1_lsm,   &
      &                          jpj2_lsm, itmpi, itmpj, itmpz, rec1_lsm, recn_lsm )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE apply_seaoverland  ***
      !!
      !! ** Purpose :   avoid spurious fluxes in coastal or near-coastal areas
      !!                due to the wrong usage of "land" values from the coarse
      !!                atmospheric model when spatial interpolation is required
      !!      D. Delrosso INGV          
      !!---------------------------------------------------------------------- 
      INTEGER,                   INTENT(in   ) :: itmpi,itmpj,itmpz                    ! lengths
      INTEGER,                   INTENT(in   ) :: jpi1_lsm,jpi2_lsm,jpj1_lsm,jpj2_lsm  ! temporary indices
      INTEGER, DIMENSION(3),     INTENT(in   ) :: rec1_lsm,recn_lsm                    ! temporary arrays for start and length
      REAL(wp),DIMENSION (:,:,:),INTENT(inout) :: zfieldo                              ! input/output array for seaoverland application
      CHARACTER (len=100),       INTENT(in   ) :: clmaskfile                           ! land/sea mask file name
      !
      INTEGER :: inum,jni,jnj,jnz,jc   ! local indices
      REAL(wp),DIMENSION (:,:,:),ALLOCATABLE :: zslmec1             ! local array for land point detection
      REAL(wp),DIMENSION (:,:),  ALLOCATABLE :: zfieldn   ! array of forcing field with undeff for land points
      REAL(wp),DIMENSION (:,:),  ALLOCATABLE :: zfield    ! array of forcing field
      !!---------------------------------------------------------------------
      !
      ALLOCATE ( zslmec1(itmpi,itmpj,itmpz), zfieldn(itmpi,itmpj), zfield(itmpi,itmpj) )
      !
      ! Retrieve the land sea mask data
      CALL iom_open( clmaskfile, inum )
      SELECT CASE( SIZE(zfieldo(jpi1_lsm:jpi2_lsm,jpj1_lsm:jpj2_lsm,:),3) )
      CASE(1)
         CALL iom_get( inum, jpdom_unknown, 'LSM', zslmec1(jpi1_lsm:jpi2_lsm,jpj1_lsm:jpj2_lsm,1),   &
            &          1, kstart = rec1_lsm, kcount = recn_lsm)
      CASE DEFAULT
         CALL iom_get( inum, jpdom_unknown, 'LSM', zslmec1(jpi1_lsm:jpi2_lsm,jpj1_lsm:jpj2_lsm,:),   &
            &          1, kstart = rec1_lsm, kcount = recn_lsm)
      END SELECT
      CALL iom_close( inum )
      !
      DO jnz=1,rec1_lsm(3)             !! Loop over k dimension
         !
         DO jni = 1, itmpi                               !! copy the original field into a tmp array
            DO jnj = 1, itmpj                            !! substituting undeff over land points
               zfieldn(jni,jnj) = zfieldo(jni,jnj,jnz)
               IF( zslmec1(jni,jnj,jnz) == 1. )   zfieldn(jni,jnj) = undeff_lsm
            END DO
         END DO
         !
         CALL seaoverland( zfieldn, itmpi, itmpj, zfield )
         DO jc = 1, nn_lsm
            CALL seaoverland( zfield, itmpi, itmpj, zfield )
         END DO
         !
         !   Check for Undeff and substitute original values
         IF( ANY(zfield==undeff_lsm) ) THEN
            DO jni = 1, itmpi
               DO jnj = 1, itmpj
                  IF( zfield(jni,jnj)==undeff_lsm )   zfield(jni,jnj) = zfieldo(jni,jnj,jnz)
               END DO
            END DO
         ENDIF
         !
         zfieldo(:,:,jnz) = zfield(:,:)
         !
      END DO                           !! End Loop over k dimension
      !
      DEALLOCATE ( zslmec1, zfieldn, zfield )
      !
   END SUBROUTINE apply_seaoverland 


   SUBROUTINE seaoverland( zfieldn, ileni, ilenj, zfield )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE seaoverland  ***
      !!
      !! ** Purpose :   create shifted matrices for seaoverland application  
      !!      D. Delrosso INGV
      !!---------------------------------------------------------------------- 
      INTEGER                          , INTENT(in   ) :: ileni,ilenj   ! lengths 
      REAL(wp), DIMENSION (ileni,ilenj), INTENT(in   ) :: zfieldn       ! array of forcing field with undeff for land points
      REAL(wp), DIMENSION (ileni,ilenj), INTENT(  out) :: zfield        ! array of forcing field
      !
      REAL(wp) , DIMENSION (ileni,ilenj)   :: zmat1, zmat2, zmat3, zmat4  ! local arrays 
      REAL(wp) , DIMENSION (ileni,ilenj)   :: zmat5, zmat6, zmat7, zmat8  !   -     - 
      REAL(wp) , DIMENSION (ileni,ilenj)   :: zlsm2d                      !   -     - 
      REAL(wp) , DIMENSION (ileni,ilenj,8) :: zlsm3d                      !   -     -
      LOGICAL  , DIMENSION (ileni,ilenj,8) :: ll_msknan3d                 ! logical mask for undeff detection
      LOGICAL  , DIMENSION (ileni,ilenj)   :: ll_msknan2d                 ! logical mask for undeff detection
      !!---------------------------------------------------------------------- 
      zmat8 = eoshift( zfieldn , SHIFT=-1 , BOUNDARY = (/zfieldn(:,1)/)     , DIM=2 )
      zmat1 = eoshift( zmat8   , SHIFT=-1 , BOUNDARY = (/zmat8(1,:)/)       , DIM=1 )
      zmat2 = eoshift( zfieldn , SHIFT=-1 , BOUNDARY = (/zfieldn(1,:)/)     , DIM=1 )
      zmat4 = eoshift( zfieldn , SHIFT= 1 , BOUNDARY = (/zfieldn(:,ilenj)/) , DIM=2 )
      zmat3 = eoshift( zmat4   , SHIFT=-1 , BOUNDARY = (/zmat4(1,:)/)       , DIM=1 )
      zmat5 = eoshift( zmat4   , SHIFT= 1 , BOUNDARY = (/zmat4(ileni,:)/)   , DIM=1 )
      zmat6 = eoshift( zfieldn , SHIFT= 1 , BOUNDARY = (/zfieldn(ileni,:)/) , DIM=1 )
      zmat7 = eoshift( zmat8   , SHIFT= 1 , BOUNDARY = (/zmat8(ileni,:)/)   , DIM=1 )
      !
      zlsm3d  = RESHAPE( (/ zmat1, zmat2, zmat3, zmat4, zmat5, zmat6, zmat7, zmat8 /), (/ ileni, ilenj, 8 /))
      ll_msknan3d = .NOT.( zlsm3d  == undeff_lsm )
      ll_msknan2d = .NOT.( zfieldn == undeff_lsm )  ! FALSE where is Undeff (land)
      zlsm2d = SUM( zlsm3d, 3 , ll_msknan3d ) / MAX( 1 , COUNT( ll_msknan3d , 3 ) )
      WHERE( COUNT( ll_msknan3d , 3 ) == 0._wp )   zlsm2d = undeff_lsm
      zfield = MERGE( zfieldn, zlsm2d, ll_msknan2d )
      !
   END SUBROUTINE seaoverland


   SUBROUTINE fld_interp( num, clvar, kw, kk, dta, nrec, lsmfile)      
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE fld_interp  ***
      !!
      !! ** Purpose :   apply weights to input gridded data to create data
      !!                on model grid
      !!----------------------------------------------------------------------
      INTEGER                   , INTENT(in   ) ::   num     ! stream number
      CHARACTER(LEN=*)          , INTENT(in   ) ::   clvar   ! variable name
      INTEGER                   , INTENT(in   ) ::   kw      ! weights number
      INTEGER                   , INTENT(in   ) ::   kk      ! vertical dimension of kk
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   dta     ! output field on model grid
      INTEGER                   , INTENT(in   ) ::   nrec    ! record number to read (ie time slice)
      CHARACTER(LEN=*)          , INTENT(in   ) ::   lsmfile ! land sea mask file name
      !
      INTEGER, DIMENSION(3) ::   rec1, recn           ! temporary arrays for start and length
      INTEGER, DIMENSION(3) ::   rec1_lsm, recn_lsm   ! temporary arrays for start and length in case of seaoverland
      INTEGER ::   ii_lsm1,ii_lsm2,ij_lsm1,ij_lsm2    ! temporary indices
      INTEGER ::   ji, jj, jk, jn, jir, jjr           ! loop counters
      INTEGER ::   ipk
      INTEGER ::   ni, nj                             ! lengths
      INTEGER ::   jpimin,jpiwid                      ! temporary indices
      INTEGER ::   jpimin_lsm,jpiwid_lsm              ! temporary indices
      INTEGER ::   jpjmin,jpjwid                      ! temporary indices
      INTEGER ::   jpjmin_lsm,jpjwid_lsm              ! temporary indices
      INTEGER ::   jpi1,jpi2,jpj1,jpj2                ! temporary indices
      INTEGER ::   jpi1_lsm,jpi2_lsm,jpj1_lsm,jpj2_lsm   ! temporary indices
      INTEGER ::   itmpi,itmpj,itmpz                     ! lengths
      REAL(wp),DIMENSION(:,:,:), ALLOCATABLE ::   ztmp_fly_dta                 ! local array of values on input grid     
      !!----------------------------------------------------------------------
      ipk = SIZE(dta, 3)
      !
      !! for weighted interpolation we have weights at four corners of a box surrounding 
      !! a model grid point, each weight is multiplied by a grid value (bilinear case)
      !! or by a grid value and gradients at the corner point (bicubic case) 
      !! so we need to have a 4 by 4 subgrid surrounding each model point to cover both cases

      !! sub grid from non-model input grid which encloses all grid points in this nemo process
      jpimin = ref_wgts(kw)%botleft(1)
      jpjmin = ref_wgts(kw)%botleft(2)
      jpiwid = ref_wgts(kw)%jpiwgt
      jpjwid = ref_wgts(kw)%jpjwgt

      !! when reading in, expand this sub-grid by one halo point all the way round for calculating gradients
      rec1(1) = MAX( jpimin-1, 1 )
      rec1(2) = MAX( jpjmin-1, 1 )
      rec1(3) = 1
      recn(1) = MIN( jpiwid+2, ref_wgts(kw)%ddims(1)-rec1(1)+1 )
      recn(2) = MIN( jpjwid+2, ref_wgts(kw)%ddims(2)-rec1(2)+1 )
      recn(3) = kk

      !! where we need to put it in the non-nemo grid fly_dta
      !! note that jpi1 and jpj1 only differ from 1 when jpimin and jpjmin are 1
      !! (ie at the extreme west or south of the whole input grid) and similarly for jpi2 and jpj2
      jpi1 = 2 + rec1(1) - jpimin
      jpj1 = 2 + rec1(2) - jpjmin
      jpi2 = jpi1 + recn(1) - 1
      jpj2 = jpj1 + recn(2) - 1


      IF( LEN_TRIM(lsmfile) > 0 ) THEN
      !! indeces for ztmp_fly_dta
      ! --------------------------
         rec1_lsm(1)=MAX(rec1(1)-nn_lsm,1)  ! starting index for enlarged external data, x direction
         rec1_lsm(2)=MAX(rec1(2)-nn_lsm,1)  ! starting index for enlarged external data, y direction
         rec1_lsm(3) = 1                    ! vertical dimension
         recn_lsm(1)=MIN(rec1(1)-rec1_lsm(1)+recn(1)+nn_lsm,ref_wgts(kw)%ddims(1)-rec1_lsm(1)) ! n points in x direction
         recn_lsm(2)=MIN(rec1(2)-rec1_lsm(2)+recn(2)+nn_lsm,ref_wgts(kw)%ddims(2)-rec1_lsm(2)) ! n points in y direction
         recn_lsm(3) = kk                   ! number of vertical levels in the input file

      !  Avoid out of bound
         jpimin_lsm = MAX( rec1_lsm(1)+1, 1 )
         jpjmin_lsm = MAX( rec1_lsm(2)+1, 1 )
         jpiwid_lsm = MIN( recn_lsm(1)-2,ref_wgts(kw)%ddims(1)-rec1(1)+1)
         jpjwid_lsm = MIN( recn_lsm(2)-2,ref_wgts(kw)%ddims(2)-rec1(2)+1)

         jpi1_lsm = 2+rec1_lsm(1)-jpimin_lsm
         jpj1_lsm = 2+rec1_lsm(2)-jpjmin_lsm
         jpi2_lsm = jpi1_lsm + recn_lsm(1) - 1
         jpj2_lsm = jpj1_lsm + recn_lsm(2) - 1


         itmpi=jpi2_lsm-jpi1_lsm+1
         itmpj=jpj2_lsm-jpj1_lsm+1
         itmpz=kk
         ALLOCATE(ztmp_fly_dta(itmpi,itmpj,itmpz))
         ztmp_fly_dta(:,:,:) = 0.0
         SELECT CASE( SIZE(ztmp_fly_dta(jpi1_lsm:jpi2_lsm,jpj1_lsm:jpj2_lsm,:),3) )
         CASE(1)
              CALL iom_get( num, jpdom_unknown, clvar, ztmp_fly_dta(jpi1_lsm:jpi2_lsm,jpj1_lsm:jpj2_lsm,1),   &
                 &          nrec, kstart = rec1_lsm, kcount = recn_lsm)
         CASE DEFAULT
              CALL iom_get( num, jpdom_unknown, clvar, ztmp_fly_dta(jpi1_lsm:jpi2_lsm,jpj1_lsm:jpj2_lsm,:),   &
                 &          nrec, kstart = rec1_lsm, kcount = recn_lsm)
         END SELECT
         CALL apply_seaoverland(lsmfile,ztmp_fly_dta(jpi1_lsm:jpi2_lsm,jpj1_lsm:jpj2_lsm,:),                  &
                 &                                      jpi1_lsm,jpi2_lsm,jpj1_lsm,jpj2_lsm,                  &
                 &                                      itmpi,itmpj,itmpz,rec1_lsm,recn_lsm)


         ! Relative indeces for remapping
         ii_lsm1 = (rec1(1)-rec1_lsm(1))+1
         ii_lsm2 = (ii_lsm1+recn(1))-1
         ij_lsm1 = (rec1(2)-rec1_lsm(2))+1
         ij_lsm2 = (ij_lsm1+recn(2))-1

         ref_wgts(kw)%fly_dta(:,:,:) = 0.0
         ref_wgts(kw)%fly_dta(jpi1:jpi2,jpj1:jpj2,:) = ztmp_fly_dta(ii_lsm1:ii_lsm2,ij_lsm1:ij_lsm2,:)
         DEALLOCATE(ztmp_fly_dta)

      ELSE
         
         ref_wgts(kw)%fly_dta(:,:,:) = 0.0
         CALL iom_get( num, jpdom_unknown, clvar, ref_wgts(kw)%fly_dta(jpi1:jpi2,jpj1:jpj2,:), nrec, kstart = rec1, kcount = recn)
      ENDIF
      

      !! first four weights common to both bilinear and bicubic
      !! data_jpi, data_jpj have already been shifted to (1,1) corresponding to botleft
      !! note that we have to offset by 1 into fly_dta array because of halo added to fly_dta (rec1 definition)
      dta(:,:,:) = 0._wp
      DO jn = 1,4
         DO_3D( 0, 0, 0, 0, 1,ipk )
            ni = ref_wgts(kw)%data_jpi(ji,jj,jn) + 1
            nj = ref_wgts(kw)%data_jpj(ji,jj,jn) + 1
            dta(ji,jj,jk) = dta(ji,jj,jk) + ref_wgts(kw)%data_wgt(ji,jj,jn) * ref_wgts(kw)%fly_dta(ni,nj,jk)
         END_3D
      END DO

      IF(ref_wgts(kw)%numwgt .EQ. 16) THEN

         !! fix up halo points that we couldnt read from file
         IF( jpi1 == 2 ) THEN
            ref_wgts(kw)%fly_dta(jpi1-1,:,:) = ref_wgts(kw)%fly_dta(jpi1,:,:)
         ENDIF
         IF( jpi2 + jpimin - 1 == ref_wgts(kw)%ddims(1)+1 ) THEN
            ref_wgts(kw)%fly_dta(jpi2+1,:,:) = ref_wgts(kw)%fly_dta(jpi2,:,:)
         ENDIF
         IF( jpj1 == 2 ) THEN
            ref_wgts(kw)%fly_dta(:,jpj1-1,:) = ref_wgts(kw)%fly_dta(:,jpj1,:)
         ENDIF
         IF( jpj2 + jpjmin - 1 == ref_wgts(kw)%ddims(2)+1 .AND. jpj2 .LT. jpjwid+2 ) THEN
            ref_wgts(kw)%fly_dta(:,jpj2+1,:) = 2.0*ref_wgts(kw)%fly_dta(:,jpj2,:) - ref_wgts(kw)%fly_dta(:,jpj2-1,:)
         ENDIF
         
         !! if data grid is cyclic we can do better on east-west edges
         !! but have to allow for whether first and last columns are coincident
         IF( ref_wgts(kw)%cyclic ) THEN
            rec1(2) = MAX( jpjmin-1, 1 )
            recn(1) = 1
            recn(2) = MIN( jpjwid+2, ref_wgts(kw)%ddims(2)-rec1(2)+1 )
            jpj1 = 2 + rec1(2) - jpjmin
            jpj2 = jpj1 + recn(2) - 1
            IF( jpi1 == 2 ) THEN
               rec1(1) = ref_wgts(kw)%ddims(1) - ref_wgts(kw)%overlap
               CALL iom_get( num, jpdom_unknown, clvar, ref_wgts(kw)%col(:,jpj1:jpj2,:), nrec, kstart = rec1, kcount = recn)
               ref_wgts(kw)%fly_dta(jpi1-1,jpj1:jpj2,:) = ref_wgts(kw)%col(1,jpj1:jpj2,:)
            ENDIF
            IF( jpi2 + jpimin - 1 == ref_wgts(kw)%ddims(1)+1 ) THEN
               rec1(1) = 1 + ref_wgts(kw)%overlap
               CALL iom_get( num, jpdom_unknown, clvar, ref_wgts(kw)%col(:,jpj1:jpj2,:), nrec, kstart = rec1, kcount = recn)
               ref_wgts(kw)%fly_dta(jpi2+1,jpj1:jpj2,:) = ref_wgts(kw)%col(1,jpj1:jpj2,:)
            ENDIF
         ENDIF
         !
!!$         DO jn = 1,4
!!$            DO_3D( 0, 0, 0, 0, 1,ipk )
!!$               ni = ref_wgts(kw)%data_jpi(ji,jj,jn) + 1
!!$               nj = ref_wgts(kw)%data_jpj(ji,jj,jn) + 1
!!$               dta(ji,jj,jk) = dta(ji,jj,jk)   &
!!$                  ! gradient in the i direction
!!$                  &            + ref_wgts(kw)%data_wgt(ji,jj,jn+4) * 0.5_wp *                                    &
!!$                  &                (ref_wgts(kw)%fly_dta(ni+1,nj  ,jk) - ref_wgts(kw)%fly_dta(ni-1,nj  ,jk))     &
!!$                  ! gradient in the j direction
!!$                  &            + ref_wgts(kw)%data_wgt(ji,jj,jn+8) * 0.5_wp *                                    &
!!$                  &                (ref_wgts(kw)%fly_dta(ni  ,nj+1,jk) - ref_wgts(kw)%fly_dta(ni  ,nj-1,jk))     &
!!$                  ! gradient in the ij direction
!!$                  &            + ref_wgts(kw)%data_wgt(ji,jj,jn+12) * 0.25_wp *                                  &
!!$                  &               ((ref_wgts(kw)%fly_dta(ni+1,nj+1,jk) - ref_wgts(kw)%fly_dta(ni-1,nj+1,jk)) -   &
!!$                  &                (ref_wgts(kw)%fly_dta(ni+1,nj-1,jk) - ref_wgts(kw)%fly_dta(ni-1,nj-1,jk)))
!!$            END_3D
!!$         END DO
         !
         DO jn = 1,4
            DO_3D( 0, 0, 0, 0, 1,ipk )
               ni = ref_wgts(kw)%data_jpi(ji,jj,jn)
               nj = ref_wgts(kw)%data_jpj(ji,jj,jn)
               ! gradient in the i direction
               dta(ji,jj,jk) = dta(ji,jj,jk) + ref_wgts(kw)%data_wgt(ji,jj,jn+4) * 0.5_wp *         &
                  &                (ref_wgts(kw)%fly_dta(ni+2,nj+1,jk) - ref_wgts(kw)%fly_dta(ni  ,nj+1,jk))
            END_3D
         END DO
         DO jn = 1,4
            DO_3D( 0, 0, 0, 0, 1,ipk )
               ni = ref_wgts(kw)%data_jpi(ji,jj,jn)
               nj = ref_wgts(kw)%data_jpj(ji,jj,jn)
               ! gradient in the j direction
               dta(ji,jj,jk) = dta(ji,jj,jk) + ref_wgts(kw)%data_wgt(ji,jj,jn+8) * 0.5_wp *         &
                  &                (ref_wgts(kw)%fly_dta(ni+1,nj+2,jk) - ref_wgts(kw)%fly_dta(ni+1,nj  ,jk))
            END_3D
         END DO
         DO jn = 1,4
            DO_3D( 0, 0, 0, 0, 1,ipk )
               ni = ref_wgts(kw)%data_jpi(ji,jj,jn)
               nj = ref_wgts(kw)%data_jpj(ji,jj,jn)
               ! gradient in the ij direction
               dta(ji,jj,jk) = dta(ji,jj,jk) + ref_wgts(kw)%data_wgt(ji,jj,jn+12) * 0.25_wp * (     &
                  &                (ref_wgts(kw)%fly_dta(ni+2,nj+2,jk) - ref_wgts(kw)%fly_dta(ni  ,nj+2,jk)) -   &
                  &                (ref_wgts(kw)%fly_dta(ni+2,nj  ,jk) - ref_wgts(kw)%fly_dta(ni  ,nj  ,jk)))
            END_3D
         END DO
         !
      ENDIF
      !
   END SUBROUTINE fld_interp


   FUNCTION fld_filename( sdjf, kday, kmonth, kyear )
      !!---------------------------------------------------------------------
      !!                    ***  FUNCTION fld_filename *** 
      !!
      !! ** Purpose :   define the filename according to a given date
      !!---------------------------------------------------------------------
      TYPE(FLD), INTENT(in) ::   sdjf         ! input field related variables
      INTEGER  , INTENT(in) ::   kday, kmonth, kyear
      !
      CHARACTER(len = 256) ::   clname, fld_filename
      !!---------------------------------------------------------------------

      
      ! build the new filename if not climatological data
      clname=TRIM(sdjf%clrootname)
      !
      ! note that sdjf%ln_clim is is only acting on the presence of the year in the file name
      IF( .NOT. sdjf%ln_clim ) THEN   
                                         WRITE(clname, '(a,"_y",i4.4)' ) TRIM( sdjf%clrootname ), kyear    ! add year
         IF( sdjf%clftyp /= 'yearly' )   WRITE(clname, '(a, "m",i2.2)' ) TRIM( clname          ), kmonth   ! add month
      ELSE
         ! build the new filename if climatological data
         IF( sdjf%clftyp /= 'yearly' )   WRITE(clname, '(a,"_m",i2.2)' ) TRIM( sdjf%clrootname ), kmonth   ! add month
      ENDIF
      IF(    sdjf%clftyp == 'daily' .OR. sdjf%clftyp(1:4) == 'week' ) &
         &                               WRITE(clname, '(a,"d" ,i2.2)' ) TRIM( clname          ), kday     ! add day

      fld_filename = clname
      
   END FUNCTION fld_filename


   FUNCTION ksec_week( cdday )
      !!---------------------------------------------------------------------
      !!                    ***  FUNCTION ksec_week *** 
      !!
      !! ** Purpose :   seconds between 00h of the beginning of the week and half of the current time step
      !!---------------------------------------------------------------------
      CHARACTER(len=*), INTENT(in)   ::   cdday   ! first 3 letters of the first day of the weekly file
      !!
      INTEGER                        ::   ksec_week      ! output variable
      INTEGER                        ::   ijul, ishift   ! local integer
      CHARACTER(len=3),DIMENSION(7)  ::   cl_week 
      !!----------------------------------------------------------------------
      cl_week = (/"sun","sat","fri","thu","wed","tue","mon"/)
      DO ijul = 1, 7
         IF( cl_week(ijul) == TRIM(cdday) ) EXIT
      END DO
      IF( ijul .GT. 7 )   CALL ctl_stop( 'ksec_week: wrong day for sdjf%clftyp(6:8): '//TRIM(cdday) )
      !
      ishift = ijul * NINT(rday)
      ! 
      ksec_week = nsec_monday + ishift
      ksec_week = MOD( ksec_week, 7*NINT(rday) )
      ! 
   END FUNCTION ksec_week

   !!======================================================================
END MODULE fldread
