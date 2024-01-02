MODULE iom
   !!======================================================================
   !!                    ***  MODULE  iom ***
   !! Input/Output manager :  Library to read input files
   !!======================================================================
   !! History :  2.0  ! 2005-12  (J. Belier) Original code
   !!            2.0  ! 2006-02  (S. Masson) Adaptation to NEMO
   !!            3.0  ! 2007-07  (D. Storkey) Changes to iom_gettime
   !!            3.4  ! 2012-12  (R. Bourdalle-Badie and G. Reffray)  add C1D case
   !!            3.6  ! 2014-15  DIMG format removed
   !!            3.6  ! 2015-15  (J. Harle) Added procedure to read REAL attributes
   !!            4.0  ! 2017-11  (M. Andrejczuk) Extend IOM interface to write any 3D fields
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   iom_open       : open a file read only
   !!   iom_close      : close a file or all files opened by iom
   !!   iom_get        : read a field (interfaced to several routines)
   !!   iom_varid      : get the id of a variable in a file
   !!   iom_rstput     : write a field in a restart file (interfaced to several routines)
   !!----------------------------------------------------------------------
   USE dom_oce         ! ocean space and time domain
   USE domutl          !
   USE flo_oce         ! floats module declarations
   USE lbclnk          ! lateal boundary condition / mpp exchanges
   USE iom_def         ! iom variables definitions
   USE iom_nf90        ! NetCDF format with native NetCDF library
   USE in_out_manager  ! I/O manager
   USE lib_mpp           ! MPP library
   USE sbc_oce  , ONLY :   nn_fsbc, ght_abl, ghw_abl, e3t_abl, e3w_abl, jpka, jpkam1
   USE icb_oce  , ONLY :   nclasses, class_num       !  !: iceberg classes
#if defined key_si3
   USE ice      , ONLY :   jpl
#endif
   USE phycst          ! physical constants
   USE dianam          ! build name of file
#if defined key_xios
   USE xios
# endif
   USE ioipsl, ONLY :  ju2ymds    ! for calendar
   USE crs             ! Grid coarsening
#if defined key_top
   USE trc, ONLY    :  profsed
#endif
   USE lib_fortran
   USE diu_bulk, ONLY : ln_diurnal_only, ln_diurnal
   USE iom_nf90
   USE netcdf

   IMPLICIT NONE
   PUBLIC   !   must be public to be able to access iom_def through iom

#if defined key_xios
   LOGICAL, PUBLIC, PARAMETER ::   lk_iomput = .TRUE.        !: iom_put flag
#else
   LOGICAL, PUBLIC, PARAMETER ::   lk_iomput = .FALSE.       !: iom_put flag
#endif
   PUBLIC iom_init, iom_init_closedef, iom_swap, iom_open, iom_close, iom_setkt, iom_varid, iom_get, iom_get_var
   PUBLIC iom_chkatt, iom_getatt, iom_putatt, iom_getszuld, iom_rstput, iom_delay_rst, iom_put
   PUBLIC iom_use, iom_context_finalize, iom_update_file_name, iom_miss_val
   PUBLIC iom_xios_setid

   PRIVATE iom_rp0d_sp, iom_rp1d_sp, iom_rp2d_sp, iom_rp3d_sp
   PRIVATE iom_rp0d_dp, iom_rp1d_dp, iom_rp2d_dp, iom_rp3d_dp
   PRIVATE iom_get_123d
   PRIVATE iom_g0d_sp, iom_g1d_sp, iom_g2d_sp, iom_g3d_sp
   PRIVATE iom_g0d_dp, iom_g1d_dp, iom_g2d_dp, iom_g3d_dp
   PRIVATE iom_p1d_sp, iom_p2d_sp, iom_p3d_sp, iom_p4d_sp
   PRIVATE iom_p1d_dp, iom_p2d_dp, iom_p3d_dp, iom_p4d_dp
#if defined key_xios
   PRIVATE iom_set_domain_attr, iom_set_axis_attr, iom_set_field_attr, iom_set_file_attr, iom_get_file_attr, iom_set_grid_attr
   PRIVATE set_grid, set_grid_bounds, set_scalar, set_xmlatt, set_mooring, iom_sdate
   PRIVATE iom_set_rst_context, iom_set_vars_active
# endif
   PRIVATE set_xios_context
   PRIVATE iom_set_rstw_active

   INTERFACE iom_get
      MODULE PROCEDURE iom_g0d_sp, iom_g1d_sp, iom_g2d_sp, iom_g3d_sp
      MODULE PROCEDURE iom_g0d_dp, iom_g1d_dp, iom_g2d_dp, iom_g3d_dp
   END INTERFACE
   INTERFACE iom_getatt
      MODULE PROCEDURE iom_g0d_iatt, iom_g1d_iatt, iom_g0d_ratt, iom_g1d_ratt, iom_g0d_catt
   END INTERFACE
   INTERFACE iom_putatt
      MODULE PROCEDURE iom_p0d_iatt, iom_p1d_iatt, iom_p0d_ratt, iom_p1d_ratt, iom_p0d_catt
   END INTERFACE
   INTERFACE iom_rstput
      MODULE PROCEDURE iom_rp0d_sp, iom_rp1d_sp, iom_rp2d_sp, iom_rp3d_sp
      MODULE PROCEDURE iom_rp0d_dp, iom_rp1d_dp, iom_rp2d_dp, iom_rp3d_dp
   END INTERFACE
   INTERFACE iom_put
      MODULE PROCEDURE iom_p0d_sp, iom_p1d_sp, iom_p2d_sp, iom_p3d_sp, iom_p4d_sp
      MODULE PROCEDURE iom_p0d_dp, iom_p1d_dp, iom_p2d_dp, iom_p3d_dp, iom_p4d_dp
   END INTERFACE iom_put

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: iom.F90 15033 2021-06-21 10:24:45Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE iom_init( cdname, kdid, ld_closedef )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE   ***
      !!
      !! ** Purpose :
      !!
      !!----------------------------------------------------------------------
      CHARACTER(len=*),           INTENT(in)  :: cdname
      INTEGER         , OPTIONAL, INTENT(in)  :: kdid
      LOGICAL         , OPTIONAL, INTENT(in)  :: ld_closedef
#if defined key_xios
      !
      TYPE(xios_duration) :: dtime    = xios_duration(0, 0, 0, 0, 0, 0)
      TYPE(xios_date)     :: start_date
      CHARACTER(len=lc) :: clname, cltmpn
      INTEGER             :: irefyear, irefmonth, irefday
      INTEGER           :: ji
      LOGICAL           :: llrst_context              ! is context related to restart
      LOGICAL           :: llrstr, llrstw
      INTEGER           :: inum
      INTEGER           :: iln
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zt_bnds, zw_bnds
      REAL(wp), DIMENSION(2,jpkam1)         :: za_bnds   ! ABL vertical boundaries
      LOGICAL ::   ll_closedef
      LOGICAL ::   ll_exist
      !!----------------------------------------------------------------------
      !
      ll_closedef = .TRUE.
      IF ( PRESENT(ld_closedef) ) ll_closedef = ld_closedef
      !
      ALLOCATE( zt_bnds(2,jpk), zw_bnds(2,jpk) )
      !
      clname = TRIM(cdname)
      IF ( .NOT. Agrif_Root() ) THEN
         iln    = INDEX(clname,'/', BACK=.TRUE.)
         cltmpn = clname(1:iln)
         clname = clname(iln+1:LEN_TRIM(clname))
         clname = TRIM(cltmpn)//TRIM(Agrif_CFixed())//'_'//TRIM(clname)
      ENDIF

      CALL xios_context_initialize(TRIM(clname), mpi_comm_oce)
      CALL iom_swap( cdname )

      llrstr = (cdname == cr_ocerst_cxt) .OR. (cdname == cr_icerst_cxt)
      llrstr = llrstr .OR. (cdname == cr_ablrst_cxt)
      llrstr = llrstr .OR. (cdname == cr_toprst_cxt)
      llrstr = llrstr .OR. (cdname == cr_sedrst_cxt)

      llrstw = (cdname == cw_ocerst_cxt) .OR. (cdname == cw_icerst_cxt)
      llrstw = llrstw .OR. (cdname == cw_ablrst_cxt)
      llrstw = llrstw .OR. (cdname == cw_toprst_cxt)
      llrstw = llrstw .OR. (cdname == cw_sedrst_cxt)

      llrst_context = llrstr .OR. llrstw

      ! Calendar type is now defined in xml file
      IF (.NOT.(xios_getvar('ref_year' ,irefyear ))) irefyear  = 1900
      IF (.NOT.(xios_getvar('ref_month',irefmonth))) irefmonth = 01
      IF (.NOT.(xios_getvar('ref_day'  ,irefday  ))) irefday   = 01

      SELECT CASE ( nleapy )        ! Choose calendar for IOIPSL
      CASE ( 1)   ;   CALL xios_define_calendar( TYPE = "Gregorian", time_origin = xios_date(irefyear,irefmonth,irefday,0,0,0),   &
          &                                                          start_date  = xios_date(   nyear,   nmonth,   nday,0,0,0) )
      CASE ( 0)   ;   CALL xios_define_calendar( TYPE = "NoLeap"   , time_origin = xios_date(irefyear,irefmonth,irefday,0,0,0),   &
          &                                                          start_date  = xios_date(   nyear,   nmonth,   nday,0,0,0) )
      CASE (30)   ;   CALL xios_define_calendar( TYPE = "D360"     , time_origin = xios_date(irefyear,irefmonth,irefday,0,0,0),   &
          &                                                          start_date  = xios_date(   nyear,   nmonth,   nday,0,0,0) )
      END SELECT

      ! horizontal grid definition
      IF(.NOT.llrst_context) CALL set_scalar
      !
      IF( cdname == cxios_context ) THEN
         CALL set_grid( "T", glamt, gphit, .FALSE., .FALSE. )
         CALL set_grid( "U", glamu, gphiu, .FALSE., .FALSE. )
         CALL set_grid( "V", glamv, gphiv, .FALSE., .FALSE. )
         CALL set_grid( "W", glamt, gphit, .FALSE., .FALSE. )
         CALL set_grid( "F", glamf, gphif, .FALSE., .FALSE. )
         CALL set_grid_znl( gphit )
         !
         IF( ln_cfmeta ) THEN   ! Add additional grid metadata
            CALL iom_set_domain_attr("grid_T", area = real( e1e2t(Nis0:Nie0, Njs0:Nje0), dp))
            CALL iom_set_domain_attr("grid_U", area = real( e1e2u(Nis0:Nie0, Njs0:Nje0), dp))
            CALL iom_set_domain_attr("grid_V", area = real( e1e2v(Nis0:Nie0, Njs0:Nje0), dp))
            CALL iom_set_domain_attr("grid_W", area = REAL( e1e2t(Nis0:Nie0, Njs0:Nje0), dp))
            CALL iom_set_domain_attr("grid_F", area = real( e1e2f(Nis0:Nie0, Njs0:Nje0), dp))
            CALL set_grid_bounds( "T", glamf, gphif, glamt, gphit )
            CALL set_grid_bounds( "U", glamv, gphiv, glamu, gphiu )
            CALL set_grid_bounds( "V", glamu, gphiu, glamv, gphiv )
            CALL set_grid_bounds( "W", glamf, gphif, glamt, gphit )
            CALL set_grid_bounds( "F", glamt, gphit, glamf, gphif )
         ENDIF
      ENDIF
      !
      IF( TRIM(cdname) == TRIM(cxios_context)//"_crs" ) THEN
         CALL dom_grid_crs   ! Save the parent grid information  & Switch to coarse grid domain
         !
         CALL set_grid( "T", glamt_crs, gphit_crs, .FALSE., .FALSE. )
         CALL set_grid( "U", glamu_crs, gphiu_crs, .FALSE., .FALSE. )
         CALL set_grid( "V", glamv_crs, gphiv_crs, .FALSE., .FALSE. )
         CALL set_grid( "W", glamt_crs, gphit_crs, .FALSE., .FALSE. )
         CALL set_grid_znl( gphit_crs )
          !
         CALL dom_grid_glo   ! Return to parent grid domain
         !
         IF( ln_cfmeta .AND. .NOT. llrst_context) THEN   ! Add additional grid metadata
            CALL iom_set_domain_attr("grid_T", area = real(e1e2t_crs(Nis0:Nie0, Njs0:Nje0), dp))
            CALL iom_set_domain_attr("grid_U", area = real(e1u_crs(Nis0:Nie0, Njs0:Nje0) * e2u_crs(Nis0:Nie0, Njs0:Nje0), dp))
            CALL iom_set_domain_attr("grid_V", area = real(e1v_crs(Nis0:Nie0, Njs0:Nje0) * e2v_crs(Nis0:Nie0, Njs0:Nje0), dp))
            CALL iom_set_domain_attr("grid_W", area = real(e1e2t_crs(Nis0:Nie0, Njs0:Nje0), dp))
            CALL set_grid_bounds( "T", glamf_crs, gphif_crs, glamt_crs, gphit_crs )
            CALL set_grid_bounds( "U", glamv_crs, gphiv_crs, glamu_crs, gphiu_crs )
            CALL set_grid_bounds( "V", glamu_crs, gphiu_crs, glamv_crs, gphiv_crs )
            CALL set_grid_bounds( "W", glamf_crs, gphif_crs, glamt_crs, gphit_crs )
         ENDIF
      ENDIF
      !
      ! vertical grid definition
      IF(.NOT.llrst_context) THEN
         CALL iom_set_axis_attr(  "deptht", paxis = gdept_1d )
         CALL iom_set_axis_attr(  "depthu", paxis = gdept_1d )
         CALL iom_set_axis_attr(  "depthv", paxis = gdept_1d )
         CALL iom_set_axis_attr(  "depthw", paxis = gdepw_1d )
          CALL iom_set_axis_attr(  "depthf", paxis = gdept_1d )

          ! ABL
         IF( .NOT. ALLOCATED(ght_abl) ) THEN   ! force definition for xml files (xios)
            ALLOCATE( ght_abl(jpka), ghw_abl(jpka), e3t_abl(jpka), e3w_abl(jpka) )   ! default allocation needed by iom
            ght_abl(:) = -1._wp   ;   ghw_abl(:) = -1._wp
            e3t_abl(:) = -1._wp   ;   e3w_abl(:) = -1._wp
         ENDIF
         CALL iom_set_axis_attr( "ght_abl", ght_abl(2:jpka) )
         CALL iom_set_axis_attr( "ghw_abl", ghw_abl(2:jpka) )

         ! Add vertical grid bounds
         zt_bnds(2,:      ) = gdept_1d(:)
         zt_bnds(1,2:jpk  ) = gdept_1d(1:jpkm1)
         zt_bnds(1,1      ) = gdept_1d(1) - e3w_1d(1)
         zw_bnds(1,:      ) = gdepw_1d(:)
         zw_bnds(2,1:jpkm1) = gdepw_1d(2:jpk)
         zw_bnds(2,jpk:   ) = gdepw_1d(jpk) + e3t_1d(jpk)
         CALL iom_set_axis_attr(  "deptht", bounds=zw_bnds )
         CALL iom_set_axis_attr(  "depthu", bounds=zw_bnds )
         CALL iom_set_axis_attr(  "depthv", bounds=zw_bnds )
         CALL iom_set_axis_attr(  "depthw", bounds=zt_bnds )
          CALL iom_set_axis_attr(  "depthf", bounds=zw_bnds )

         ! ABL
         za_bnds(1,:) = ghw_abl(1:jpkam1)
         za_bnds(2,:) = ghw_abl(2:jpka  )
         CALL iom_set_axis_attr( "ght_abl", bounds=za_bnds )
         za_bnds(1,:) = ght_abl(2:jpka  )
         za_bnds(2,:) = ght_abl(2:jpka  ) + e3w_abl(2:jpka)
         CALL iom_set_axis_attr( "ghw_abl", bounds=za_bnds )

         CALL iom_set_axis_attr(  "nfloat", (/ (REAL(ji,wp), ji=1,jpnfl) /) )
# if defined key_si3
         CALL iom_set_axis_attr( "ncatice", (/ (REAL(ji,wp), ji=1,jpl) /) )
         ! SIMIP diagnostics (4 main arctic straits)
!         remove from si3 case to general case to use dct diagnostics over
!         oce+ice
!         CALL iom_set_axis_attr( "nstrait", (/ (REAL(ji,wp), ji=1,4) /) )
# endif
#if defined key_top
         IF( ALLOCATED(profsed) ) CALL iom_set_axis_attr( "profsed", paxis = profsed )
#endif
         CALL iom_set_axis_attr( "icbcla", class_num )
         CALL iom_set_axis_attr( "iax_20C", (/ REAL(20,wp) /) )   ! strange syntaxe and idea...
         CALL iom_set_axis_attr( "iax_26C", (/ REAL(26,wp) /) )   ! strange syntaxe and idea...
         CALL iom_set_axis_attr( "iax_28C", (/ REAL(28,wp) /) )   ! strange syntaxe and idea...
         ! for diaprt, we need to define an axis which size can be 1 (default) or 5 (if the file subbasins.nc exists)
         INQUIRE( FILE = 'subbasins.nc', EXIST = ll_exist )
         nbasin = 1 + 4 * COUNT( (/ll_exist/) )
         CALL iom_set_axis_attr( "basin"  , (/ (REAL(ji,wp), ji=1,nbasin) /) )
         ! Transport diagnostics diadct - max number of section 150 in diadct -> ! pb to ensure consistency here
         IF( ln_diadct  ) THEN
            CALL iom_set_axis_attr( "nstrait"  , (/ (REAL(ji,wp), ji=1,nb_sec) /) )
         ELSE
            CALL iom_set_axis_attr( "nstrait", (/ (REAL(ji,wp), ji=1,4) /) )
         ENDIF
      ENDIF
      !
      ! automatic definitions of some of the xml attributs
      IF(llrstr) THEN
         IF(PRESENT(kdid)) THEN
            CALL iom_set_rst_context(.TRUE.)
!set which fields will be read from restart file
            CALL iom_set_vars_active(kdid)
         ELSE
            CALL ctl_stop( 'iom_init:', 'restart read with XIOS: missing pointer to NETCDF file' )
         ENDIF
      ELSE IF(llrstw) THEN
         CALL iom_set_rstw_file(iom_file(kdid)%name)
      ELSE
         CALL set_xmlatt
      ENDIF
      !
      ! set time step length
      dtime%second = rn_Dt
      CALL xios_set_timestep( dtime )
      !
      ! conditional closure of context definition
      IF ( ll_closedef ) CALL iom_init_closedef
      !
      DEALLOCATE( zt_bnds, zw_bnds )
      !
#endif
      !
   END SUBROUTINE iom_init

   SUBROUTINE iom_init_closedef(cdname)
      !!----------------------------------------------------------------------
      !!            ***  SUBROUTINE iom_init_closedef  ***
      !!----------------------------------------------------------------------
      !!
      !! ** Purpose : Closure of context definition
      !!
      !!----------------------------------------------------------------------
      CHARACTER(len=*), OPTIONAL, INTENT(IN) :: cdname
#if defined key_xios
      LOGICAL :: llrstw

      llrstw = .FALSE.
      IF(PRESENT(cdname)) THEN
         llrstw = (cdname == cw_ocerst_cxt)
         llrstw = llrstw .OR. (cdname == cw_icerst_cxt)
         llrstw = llrstw .OR. (cdname == cw_ablrst_cxt)
         llrstw = llrstw .OR. (cdname == cw_toprst_cxt)
         llrstw = llrstw .OR. (cdname == cw_sedrst_cxt)
      ENDIF

      IF( llrstw ) THEN
!set names of the fields in restart file IF using XIOS to write data
         CALL iom_set_rst_context(.FALSE.)
         CALL xios_close_context_definition()
      ELSE
         CALL xios_close_context_definition()
         CALL xios_update_calendar( 0 )
      ENDIF
#else
      IF( .FALSE. )   WRITE(numout,*) 'iom_init_closedef: should not see this'   ! useless statement to avoid compilation warnings
#endif

   END SUBROUTINE iom_init_closedef

   SUBROUTINE iom_set_vars_active(idnum)
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE  iom_set_vars_active  ***
      !!
      !! ** Purpose :  define filename in XIOS context for reading file,
      !!               enable variables present in a file for reading with XIOS
      !!               id of the file is assumed to be rrestart.
      !!---------------------------------------------------------------------
      INTEGER, INTENT(IN) :: idnum

#if defined key_xios
      INTEGER                                    :: ndims, nvars, natts, unlimitedDimId, dimlen, xtype,mdims
      TYPE(xios_field)                           :: field_hdl
      TYPE(xios_file)                            :: file_hdl
      TYPE(xios_filegroup)                       :: filegroup_hdl
      INTEGER                                    :: dimids(4), jv,i, idim
      CHARACTER(LEN=256)                         :: clinfo               ! info character
      INTEGER, ALLOCATABLE                       :: indimlens(:)
      CHARACTER(LEN=nf90_max_name), ALLOCATABLE  :: indimnames(:)
      CHARACTER(LEN=nf90_max_name)               :: dimname, varname
      INTEGER                                    :: iln
      CHARACTER(LEN=lc)                          :: fname
      LOGICAL                                    :: lmeta
!metadata in restart file for restart read with XIOS
      INTEGER, PARAMETER                         :: NMETA = 11
      CHARACTER(LEN=lc)                          :: meta(NMETA)


      meta(1) = "nav_lat"
      meta(2) = "nav_lon"
      meta(3) = "nav_lev"
      meta(4) = "time_instant"
      meta(5) = "time_instant_bounds"
      meta(6) = "time_counter"
      meta(7) = "time_counter_bounds"
      meta(8) = "x"
      meta(9) = "y"
      meta(10) = "numcat"
      meta(11) = "nav_hgt"

      clinfo = '          iom_set_vars_active, file: '//TRIM(iom_file(idnum)%name)

      iln = INDEX( iom_file(idnum)%name, '.nc' )
!XIOS doee not need .nc
      IF(iln > 0) THEN
        fname =  iom_file(idnum)%name(1:iln-1)
      ELSE
        fname =  iom_file(idnum)%name
      ENDIF

!set name of the restart file and enable available fields
      CALL xios_get_handle("file_definition", filegroup_hdl )
      CALL xios_add_child(filegroup_hdl, file_hdl, 'rrestart')
      CALL xios_set_file_attr( "rrestart", name=fname, type="one_file",      &
           par_access="collective", enabled=.TRUE., mode="read",              &
                                                    output_freq=xios_timestep )

      CALL iom_nf90_check( nf90_inquire(iom_file(idnum)%nfid, ndims, nvars, natts ), clinfo )
      ALLOCATE(indimlens(ndims), indimnames(ndims))
      CALL iom_nf90_check( nf90_inquire(iom_file(idnum)%nfid, unlimitedDimId = unlimitedDimId ), clinfo )

      DO idim = 1, ndims
         CALL iom_nf90_check( nf90_inquire_dimension(iom_file(idnum)%nfid, idim, dimname, dimlen ), clinfo )
         indimlens(idim) = dimlen
         indimnames(idim) = dimname
      ENDDO

      DO jv =1, nvars
         lmeta = .FALSE.
         CALL iom_nf90_check( nf90_inquire_variable(iom_file(idnum)%nfid, jv, varname, xtype, ndims, dimids, natts ), clinfo )
         DO i = 1, NMETA
           IF(varname == meta(i)) THEN
             lmeta = .TRUE.
           ENDIF
         ENDDO
         IF(.NOT.lmeta) THEN
            CALL xios_add_child(file_hdl, field_hdl, varname)
            mdims = ndims

            IF(ANY(dimids(1:ndims) == unlimitedDimId)) THEN
               mdims = mdims - 1
            ENDIF

            IF(mdims == 3) THEN
               CALL xios_set_attr (field_hdl, enabled = .TRUE., name = varname,   &
                                   domain_ref="grid_N",                           &
                                   axis_ref=iom_axis(indimlens(dimids(mdims))),   &
                                   prec = 8, operation = "instant"                )
            ELSEIF(mdims == 2) THEN
               CALL xios_set_attr (field_hdl, enabled = .TRUE., name = varname,  &
                                   domain_ref="grid_N", prec = 8,                &
                                   operation = "instant"                         )
            ELSEIF(mdims == 1) THEN
               CALL xios_set_attr (field_hdl, enabled = .TRUE., name = varname, &
                                   axis_ref=iom_axis(indimlens(dimids(mdims))), &
                                   prec = 8, operation = "instant"              )
            ELSEIF(mdims == 0) THEN
               CALL xios_set_attr (field_hdl, enabled = .TRUE., name = varname, &
                                   scalar_ref = "grid_scalar", prec = 8,        &
                                   operation = "instant"                        )
            ELSE
               WRITE(ctmp1,*) 'iom_set_vars_active: variable ', TRIM(varname) ,' incorrect number of dimensions'
               CALL ctl_stop( 'iom_set_vars_active:', ctmp1 )
            ENDIF
         ENDIF
      ENDDO
      DEALLOCATE(indimlens, indimnames)
#endif
   END SUBROUTINE iom_set_vars_active

   SUBROUTINE iom_set_rstw_file(cdrst_file)
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE iom_set_rstw_file   ***
      !!
      !! ** Purpose :  define file name in XIOS context for writing restart
      !!---------------------------------------------------------------------
      CHARACTER(len=*) :: cdrst_file
#if defined key_xios
      TYPE(xios_file) :: file_hdl
      TYPE(xios_filegroup) :: filegroup_hdl

!set name of the restart file and enable available fields
      IF(lwp) WRITE(numout,*) 'Setting restart filename (for XIOS write) to: ', TRIM(cdrst_file)
      CALL xios_get_handle("file_definition", filegroup_hdl )
      CALL xios_add_child(filegroup_hdl, file_hdl, 'wrestart')
      IF(nxioso.eq.1) THEN
         CALL xios_set_file_attr( "wrestart", type="one_file", enabled=.TRUE.,&
                                       mode="write", output_freq=xios_timestep)
         IF(lwp) write(numout,*) 'OPEN ', TRIM(cdrst_file), ' in one_file mode'
      ELSE
         CALL xios_set_file_attr( "wrestart", type="multiple_file", enabled=.TRUE.,&
                                            mode="write", output_freq=xios_timestep)
         IF(lwp) write(numout,*) 'OPEN ', TRIM(cdrst_file), ' in multiple_file mode'
      ENDIF
      CALL xios_set_file_attr( "wrestart", name=TRIM(cdrst_file))
#endif
   END SUBROUTINE iom_set_rstw_file


   SUBROUTINE iom_set_rstw_active(sdfield, rd0, rs0, rd1, rs1, rd2, rs2, rd3, rs3)
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE iom_set_rstw_active   ***
      !!
      !! ** Purpose :  define file name in XIOS context for writing restart
      !!               enable variables present in restart file for writing
      !!---------------------------------------------------------------------
!sets enabled = .TRUE. for each field in restart file
      CHARACTER(len = *), INTENT(IN)                     :: sdfield
      REAL(dp), OPTIONAL, INTENT(IN)                     :: rd0
      REAL(sp), OPTIONAL, INTENT(IN)                     :: rs0
      REAL(dp), OPTIONAL, INTENT(IN), DIMENSION(:)       :: rd1
      REAL(sp), OPTIONAL, INTENT(IN), DIMENSION(:)       :: rs1
      REAL(dp), OPTIONAL, INTENT(IN), DIMENSION(:, :)    :: rd2
      REAL(sp), OPTIONAL, INTENT(IN), DIMENSION(:, :)    :: rs2
      REAL(dp), OPTIONAL, INTENT(IN), DIMENSION(:, :, :) :: rd3
      REAL(sp), OPTIONAL, INTENT(IN), DIMENSION(:, :, :) :: rs3
#if defined key_xios
      TYPE(xios_field) :: field_hdl
      TYPE(xios_file) :: file_hdl

      CALL xios_get_handle("wrestart", file_hdl)
!define fields for restart context
      CALL xios_add_child(file_hdl, field_hdl, sdfield)

      IF(PRESENT(rd3)) THEN
         CALL xios_set_attr (field_hdl, enabled = .TRUE., name = sdfield, &
                             domain_ref = "grid_N",                       &
                             axis_ref = iom_axis(size(rd3, 3)),           &
                             prec = 8, operation = "instant"              )
      ELSEIF(PRESENT(rs3)) THEN
         CALL xios_set_attr (field_hdl, enabled = .TRUE., name = sdfield, &
                             domain_ref = "grid_N",                       &
                             axis_ref = iom_axis(size(rd3, 3)),           &
                             prec = 4, operation = "instant"              )
      ELSEIF(PRESENT(rd2)) THEN
         CALL xios_set_attr (field_hdl, enabled = .TRUE., name = sdfield, &
                             domain_ref = "grid_N", prec = 8,             &
                             operation = "instant"                        )
      ELSEIF(PRESENT(rs2)) THEN
         CALL xios_set_attr (field_hdl, enabled = .TRUE., name = sdfield, &
                             domain_ref = "grid_N", prec = 4,             &
                             operation = "instant"                        )
      ELSEIF(PRESENT(rd1)) THEN
         CALL xios_set_attr (field_hdl, enabled = .TRUE., name = sdfield, &
                             axis_ref = iom_axis(size(rd1, 1)),           &
                             prec = 8, operation = "instant"              )
      ELSEIF(PRESENT(rs1)) THEN
         CALL xios_set_attr (field_hdl, enabled = .TRUE., name = sdfield, &
                             axis_ref = iom_axis(size(rd1, 1)),           &
                             prec = 4, operation = "instant"              )
      ELSEIF(PRESENT(rd0)) THEN
         CALL xios_set_attr (field_hdl, enabled = .TRUE., name = sdfield, &
                             scalar_ref = "grid_scalar", prec = 8,        &
                             operation = "instant"                        )
      ELSEIF(PRESENT(rs0)) THEN
         CALL xios_set_attr (field_hdl, enabled = .TRUE., name = sdfield, &
                             scalar_ref = "grid_scalar", prec = 4,        &
                             operation = "instant"                        )
      ENDIF
#endif
   END SUBROUTINE iom_set_rstw_active

   FUNCTION iom_axis(idlev) result(axis_ref)
      !!---------------------------------------------------------------------
      !!                   ***  FUNCTION  iom_axis  ***
      !!
      !! ** Purpose : Used for grid definition when XIOS is used to read/write
      !!              restart. Returns axis corresponding to the number of levels
      !!              given as an input variable. Axes are defined in routine
      !!              iom_set_rst_context
      !!---------------------------------------------------------------------
      INTEGER, INTENT(IN) :: idlev
      CHARACTER(len=lc)   :: axis_ref
      CHARACTER(len=12)   :: str
      IF(idlev == jpk) THEN
         axis_ref="nav_lev"
      ELSEIF(idlev == jpka) THEN
         axis_ref="nav_hgt"
#if defined key_si3
      ELSEIF(idlev == jpl) THEN
         axis_ref="numcat"
#endif
      ELSE
         write(str, *) idlev
         CALL ctl_stop( 'iom_axis', 'Definition for axis with '//TRIM(ADJUSTL(str))//' levels missing')
      ENDIF
   END FUNCTION iom_axis

   FUNCTION iom_xios_setid(cdname) result(kid)
     !!---------------------------------------------------------------------
      !!                   ***  FUNCTION    ***
      !!
      !! ** Purpose : this function returns first available id to keep information about file
      !!              sets filename in iom_file structure and sets name
      !!              of XIOS context depending on cdcomp
      !!              corresponds to iom_nf90_open
      !!---------------------------------------------------------------------
      CHARACTER(len=*), INTENT(in   ) :: cdname      ! File name
      INTEGER                         :: kid      ! identifier of the opened file
      INTEGER                         :: jl

      kid = 0
      DO jl = jpmax_files, 1, -1
         IF( iom_file(jl)%nfid == 0 )   kid = jl
      ENDDO

      iom_file(kid)%name   = TRIM(cdname)
      iom_file(kid)%nfid   = 1
      iom_file(kid)%nvars  = 0
      iom_file(kid)%irec   = -1

   END FUNCTION iom_xios_setid

   SUBROUTINE iom_set_rst_context(ld_rstr)
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE  iom_set_rst_context  ***
      !!
      !! ** Purpose : Define domain, axis and grid for restart (read/write)
      !!              context
      !!
      !!---------------------------------------------------------------------
      LOGICAL, INTENT(IN)               :: ld_rstr
      INTEGER :: ji
#if defined key_xios
      TYPE(xios_domaingroup)            :: domaingroup_hdl
      TYPE(xios_domain)                 :: domain_hdl
      TYPE(xios_axisgroup)              :: axisgroup_hdl
      TYPE(xios_axis)                   :: axis_hdl
      TYPE(xios_scalar)                 :: scalar_hdl
      TYPE(xios_scalargroup)            :: scalargroup_hdl

      CALL xios_get_handle("domain_definition",domaingroup_hdl)
      CALL xios_add_child(domaingroup_hdl, domain_hdl, "grid_N")
      CALL set_grid("N", glamt, gphit, .TRUE., ld_rstr)

      CALL xios_get_handle("axis_definition",axisgroup_hdl)
      CALL xios_add_child(axisgroup_hdl, axis_hdl, "nav_lev")
!AGRIF fails to compile when unit= is in call to xios_set_axis_attr
!     CALL xios_set_axis_attr( "nav_lev", long_name="Vertical levels",  unit="m", positive="down")
      CALL xios_set_axis_attr( "nav_lev", long_name = "Vertical levels in meters", positive = "down")
      CALL iom_set_axis_attr( "nav_lev", paxis = gdept_1d )
#if defined key_si3
      CALL xios_add_child(axisgroup_hdl, axis_hdl, "numcat")
      CALL iom_set_axis_attr( "numcat", (/ (REAL(ji,wp), ji=1,jpl) /) )
#endif
      CALL xios_add_child(axisgroup_hdl, axis_hdl, "nav_hgt")
      CALL iom_set_axis_attr( "nav_hgt", (/ (REAL(ji,wp), ji=1,jpka) /) )
      CALL xios_get_handle("scalar_definition", scalargroup_hdl)
      CALL xios_add_child(scalargroup_hdl, scalar_hdl, "grid_scalar")
#endif
   END SUBROUTINE iom_set_rst_context


   SUBROUTINE set_xios_context(kdid, cdcont)
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE  iom_set_rst_context  ***
      !!
      !! ** Purpose : set correct XIOS context based on kdid
      !!
      !!---------------------------------------------------------------------
      INTEGER,           INTENT(IN)     :: kdid           ! Identifier of the file
      CHARACTER(LEN=lc), INTENT(OUT)    :: cdcont         ! name of the context for XIOS read/write

      cdcont = "NONE"

      IF(lrxios) THEN
         IF(kdid == numror) THEN
            cdcont = cr_ocerst_cxt
         ELSEIF(kdid == numrir) THEN
            cdcont = cr_icerst_cxt
         ELSEIF(kdid == numrar) THEN
            cdcont = cr_ablrst_cxt
         ELSEIF(kdid == numrtr) THEN
            cdcont = cr_toprst_cxt
         ELSEIF(kdid == numrsr) THEN
            cdcont = cr_sedrst_cxt
         ENDIF
      ENDIF

      IF(lwxios) THEN
         IF(kdid == numrow) THEN
            cdcont = cw_ocerst_cxt
         ELSEIF(kdid == numriw) THEN
            cdcont = cw_icerst_cxt
         ELSEIF(kdid == numraw) THEN
            cdcont = cw_ablrst_cxt
         ELSEIF(kdid == numrtw) THEN
            cdcont = cw_toprst_cxt
         ELSEIF(kdid == numrsw) THEN
            cdcont = cw_sedrst_cxt
         ENDIF
      ENDIF
   END SUBROUTINE set_xios_context


   SUBROUTINE iom_swap( cdname )
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE  iom_swap  ***
      !!
      !! ** Purpose :  swap context between different agrif grid for xmlio_server
      !!---------------------------------------------------------------------
      CHARACTER(len=*), INTENT(in) :: cdname
      CHARACTER(len=256)           :: clname, cltmpn
      INTEGER                      :: iln
#if defined key_xios
      TYPE(xios_context) :: nemo_hdl

      clname = TRIM(cdname)
      IF ( .NOT. Agrif_Root() ) THEN
         iln    = INDEX(clname,'/', BACK=.TRUE.)
         cltmpn = clname(1:iln)
         clname = clname(iln+1:LEN_TRIM(clname))
         clname = TRIM(cltmpn)//TRIM(Agrif_CFixed())//'_'//TRIM(clname)
      ENDIF
      !
      CALL xios_get_handle(clname,nemo_hdl)
      !
      CALL xios_set_current_context(nemo_hdl)
#endif
      !
   END SUBROUTINE iom_swap


   SUBROUTINE iom_open( cdname, kiomid, ldwrt, ldstop, ldiof, kdlev, cdcomp )
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE  iom_open  ***
      !!
      !! ** Purpose :  open an input file (return 0 if not found)
      !!---------------------------------------------------------------------
      CHARACTER(len=*), INTENT(in   )           ::   cdname   ! File name
      INTEGER         , INTENT(  out)           ::   kiomid   ! iom identifier of the opened file
      LOGICAL         , INTENT(in   ), OPTIONAL ::   ldwrt    ! open in write modeb          (default = .FALSE.)
      LOGICAL         , INTENT(in   ), OPTIONAL ::   ldstop   ! stop if open to read a non-existing file (default = .TRUE.)
      LOGICAL         , INTENT(in   ), OPTIONAL ::   ldiof    ! Interp On the Fly, needed for AGRIF (default = .FALSE.)
      INTEGER         , INTENT(in   ), OPTIONAL ::   kdlev    ! number of vertical levels
      CHARACTER(len=3), INTENT(in   ), OPTIONAL ::   cdcomp   ! name of component calling iom_nf90_open
      !
      CHARACTER(LEN=256)    ::   clname    ! the name of the file based on cdname [[+clcpu]+clcpu]
      CHARACTER(LEN=256)    ::   cltmpn    ! tempory name to store clname (in writting mode)
      CHARACTER(LEN=10)     ::   clsuffix  ! ".nc"
      CHARACTER(LEN=15)     ::   clcpu     ! the cpu number (max jpmax_digits digits)
      CHARACTER(LEN=256)    ::   clinfo    ! info character
      LOGICAL               ::   llok      ! check the existence
      LOGICAL               ::   llwrt     ! local definition of ldwrt
      LOGICAL               ::   llstop    ! local definition of ldstop
      LOGICAL               ::   lliof     ! local definition of ldiof
      INTEGER               ::   icnt      ! counter for digits in clcpu (max = jpmax_digits)
      INTEGER               ::   iln, ils  ! lengths of character
      INTEGER               ::   istop     !
      ! local number of points for x,y dimensions
      ! position of first local point for x,y dimensions
      ! position of last local point for x,y dimensions
      ! start halo size for x,y dimensions
      ! end halo size for x,y dimensions
      !---------------------------------------------------------------------
      ! Initializations and control
      ! =============
      kiomid = -1
      clinfo = '                    iom_open ~~~  '
      istop = nstop
      ! if iom_open is called for the first time: initialize iom_file(:)%nfid to 0
      ! (could be done when defining iom_file in f95 but not in f90)
      IF( Agrif_Root() ) THEN
         IF( iom_open_init == 0 ) THEN
            iom_file(:)%nfid = 0
            iom_open_init = 1
         ENDIF
      ENDIF
      ! do we read or write the file?
      IF( PRESENT(ldwrt) ) THEN   ;   llwrt = ldwrt
      ELSE                        ;   llwrt = .FALSE.
      ENDIF
      ! do we call ctl_stop if we try to open a non-existing file in read mode?
      IF( PRESENT(ldstop) ) THEN   ;   llstop = ldstop
      ELSE                         ;   llstop = .TRUE.
      ENDIF
      ! are we using interpolation on the fly?
      IF( PRESENT(ldiof) ) THEN   ;   lliof = ldiof
      ELSE                        ;   lliof = .FALSE.
      ENDIF
      ! create the file name by added, if needed, TRIM(Agrif_CFixed()) and TRIM(clsuffix)
      ! =============
      clname   = TRIM(cdname)
      IF ( .NOT. Agrif_Root() .AND. .NOT. lliof ) THEN
         iln    = INDEX(clname,'/', BACK=.TRUE.)
         cltmpn = clname(1:iln)
         clname = clname(iln+1:LEN_TRIM(clname))
         clname = TRIM(cltmpn)//TRIM(Agrif_CFixed())//'_'//TRIM(clname)
      ENDIF
      ! which suffix should we use?
      clsuffix = '.nc'
      ! Add the suffix if needed
      iln = LEN_TRIM(clname)
      ils = LEN_TRIM(clsuffix)
      IF( iln <= ils .OR. INDEX( TRIM(clname), TRIM(clsuffix), back = .TRUE. ) /= iln - ils + 1 )   &
         &   clname = TRIM(clname)//TRIM(clsuffix)
      cltmpn = clname   ! store this name
      ! try to find if the file to be opened already exist
      ! =============
      INQUIRE( FILE = clname, EXIST = llok )
      IF( .NOT.llok ) THEN
         ! we try to add the cpu number to the name
         WRITE(clcpu,*) narea-1

         clcpu  = TRIM(ADJUSTL(clcpu))
         iln = INDEX(clname,TRIM(clsuffix), back = .TRUE.)
         clname = clname(1:iln-1)//'_'//TRIM(clcpu)//TRIM(clsuffix)
         icnt = 0
         INQUIRE( FILE = clname, EXIST = llok )
         ! we try different formats for the cpu number by adding 0
         DO WHILE( .NOT.llok .AND. icnt < jpmax_digits )
            clcpu  = "0"//TRIM(clcpu)
            clname = clname(1:iln-1)//'_'//TRIM(clcpu)//TRIM(clsuffix)
            INQUIRE( FILE = clname, EXIST = llok )
            icnt = icnt + 1
         END DO
      ELSE
         lxios_sini = .TRUE.
      ENDIF
      ! Open the NetCDF file
      ! =============
      ! do we have some free file identifier?
      IF( MINVAL(iom_file(:)%nfid) /= 0 )   &
         &   CALL ctl_stop( TRIM(clinfo), 'No more free file identifier', 'increase jpmax_files in iom_def' )
      ! if no file was found...
      IF( .NOT. llok ) THEN
         IF( .NOT. llwrt ) THEN   ! we are in read mode
            IF( llstop ) THEN   ;   CALL ctl_stop( TRIM(clinfo), 'File '//TRIM(cltmpn)//'* not found' )
            ELSE                ;   istop = nstop + 1   ! make sure that istop /= nstop so we don't open the file
            ENDIF
         ELSE                     ! we are in write mode so we
            clname = cltmpn       ! get back the file name without the cpu number
         ENDIF
      ELSE
         IF( llwrt .AND. .NOT. ln_clobber ) THEN   ! we stop as we want to write in a new file
            CALL ctl_stop( TRIM(clinfo), 'We want to write in a new file but '//TRIM(clname)//' already exists...' )
            istop = nstop + 1                      ! make sure that istop /= nstop so we don't open the file
         ELSEIF( llwrt ) THEN     ! the file exists and we are in write mode with permission to
            clname = cltmpn       ! overwrite so get back the file name without the cpu number
         ENDIF
      ENDIF
      IF( istop == nstop ) THEN   ! no error within this routine
         CALL iom_nf90_open( clname, kiomid, llwrt, llok, kdlev = kdlev, cdcomp = cdcomp )
      ENDIF
      !
   END SUBROUTINE iom_open


   SUBROUTINE iom_close( kiomid )
      !!--------------------------------------------------------------------
      !!                   ***  SUBROUTINE  iom_close  ***
      !!
      !! ** Purpose : close an input file, or all files opened by iom
      !!--------------------------------------------------------------------
      INTEGER, INTENT(inout), OPTIONAL ::   kiomid   ! iom identifier of the file to be closed
      !                                              ! return 0 when file is properly closed
      !                                              ! No argument: all files opened by iom are closed

      INTEGER ::   jf         ! dummy loop indices
      INTEGER ::   i_s, i_e   ! temporary integer
      CHARACTER(LEN=100)    ::   clinfo    ! info character
      !---------------------------------------------------------------------
      !
      IF( iom_open_init == 0 )   RETURN   ! avoid to use iom_file(jf)%nfid that us not yet initialized
      !
      clinfo = '                    iom_close ~~~  '
      IF( PRESENT(kiomid) ) THEN
         i_s = kiomid
         i_e = kiomid
      ELSE
         i_s = 1
         i_e = jpmax_files
      ENDIF

      IF( i_s > 0 ) THEN
         DO jf = i_s, i_e
            IF( iom_file(jf)%nfid > 0 ) THEN
               CALL iom_nf90_close( jf )
               iom_file(jf)%nfid       = 0          ! free the id
               IF( PRESENT(kiomid) )   kiomid = 0   ! return 0 as id to specify that the file was closed
               IF(lwp) WRITE(numout,*) TRIM(clinfo)//' close file: '//TRIM(iom_file(jf)%name)//' ok'
            ELSEIF( PRESENT(kiomid) ) THEN
               WRITE(ctmp1,*) '--->',  kiomid
               CALL ctl_stop( TRIM(clinfo)//' Invalid file identifier', ctmp1 )
            ENDIF
         END DO
      ENDIF
      !
   END SUBROUTINE iom_close


   FUNCTION iom_varid ( kiomid, cdvar, kdimsz, kndims, lduld, ldstop )
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION  iom_varid  ***
      !!
      !! ** Purpose : get the id of a variable in a file (return 0 if not found)
      !!-----------------------------------------------------------------------
      INTEGER              , INTENT(in   )           ::   kiomid   ! file Identifier
      CHARACTER(len=*)     , INTENT(in   )           ::   cdvar    ! name of the variable
      INTEGER, DIMENSION(:), INTENT(  out), OPTIONAL ::   kdimsz   ! size of each dimension
      INTEGER              , INTENT(  out), OPTIONAL ::   kndims   ! number of dimensions
      LOGICAL              , INTENT(  out), OPTIONAL ::   lduld    ! true if the last dimension is unlimited (time)
      LOGICAL              , INTENT(in   ), OPTIONAL ::   ldstop   ! stop if looking for non-existing variable (default = .TRUE.)
      !
      INTEGER                        ::   iom_varid, iiv, i_nvd
      LOGICAL                        ::   ll_fnd
      CHARACTER(LEN=100)             ::   clinfo                   ! info character
      LOGICAL                        ::   llstop                   ! local definition of ldstop
      !!-----------------------------------------------------------------------
      iom_varid = 0                         ! default definition
      ! do we call ctl_stop if we look for non-existing variable?
      IF( PRESENT(ldstop) ) THEN   ;   llstop = ldstop
      ELSE                         ;   llstop = .TRUE.
      ENDIF
      !
      IF( kiomid > 0 ) THEN
         clinfo = 'iom_varid, file: '//TRIM(iom_file(kiomid)%name)//', var: '//TRIM(cdvar)
         IF( iom_file(kiomid)%nfid == 0 ) THEN
            CALL ctl_stop( TRIM(clinfo), 'the file is not open' )
         ELSE
            ll_fnd  = .FALSE.
            iiv = 0
            !
            DO WHILE ( .NOT.ll_fnd .AND. iiv < iom_file(kiomid)%nvars )
               iiv = iiv + 1
               ll_fnd  = ( TRIM(cdvar) == TRIM(iom_file(kiomid)%cn_var(iiv)) )
            END DO
            !
            IF( .NOT.ll_fnd ) THEN
               iiv = iiv + 1
               IF( iiv <= jpmax_vars ) THEN
                  iom_varid = iom_nf90_varid( kiomid, cdvar, iiv, kdimsz, kndims, lduld )
               ELSE
                  CALL ctl_stop( TRIM(clinfo), 'Too many variables in the file '//iom_file(kiomid)%name,   &
                        &                      'increase the parameter jpmax_vars')
               ENDIF
               IF( llstop .AND. iom_varid == -1 )   CALL ctl_stop( TRIM(clinfo)//' not found' )
            ELSE
               iom_varid = iiv
               IF( PRESENT(kdimsz) ) THEN
                  i_nvd = iom_file(kiomid)%ndims(iiv)
                  IF( i_nvd <= size(kdimsz) ) THEN
                     kdimsz(1:i_nvd) = iom_file(kiomid)%dimsz(1:i_nvd,iiv)
                  ELSE
                     WRITE(ctmp1,*) i_nvd, size(kdimsz)
                     CALL ctl_stop( TRIM(clinfo), 'error in kdimsz size'//TRIM(ctmp1) )
                  ENDIF
               ENDIF
               IF( PRESENT(kndims) )  kndims = iom_file(kiomid)%ndims(iiv)
               IF( PRESENT( lduld) )  lduld  = iom_file(kiomid)%luld( iiv)
            ENDIF
         ENDIF
      ENDIF
      !
   END FUNCTION iom_varid


   !!----------------------------------------------------------------------
   !!                   INTERFACE iom_get
   !!----------------------------------------------------------------------
   SUBROUTINE iom_g0d_sp( kiomid, cdvar, pvar, ktime )
      INTEGER         , INTENT(in   )                 ::   kiomid    ! Identifier of the file
      CHARACTER(len=*), INTENT(in   )                 ::   cdvar     ! Name of the variable
      REAL(sp)        , INTENT(  out)                 ::   pvar      ! read field
      REAL(dp)                                        ::   ztmp_pvar ! tmp var to read field
      INTEGER         , INTENT(in   ),     OPTIONAL   ::   ktime     ! record number
      !
      INTEGER                                         ::   idvar     ! variable id
      INTEGER                                         ::   idmspc    ! number of spatial dimensions
      INTEGER         , DIMENSION(1)                  ::   itime     ! record number
      CHARACTER(LEN=100)                              ::   clinfo    ! info character
      CHARACTER(LEN=100)                              ::   clname    ! file name
      CHARACTER(LEN=1)                                ::   cldmspc   !
      CHARACTER(LEN=lc)                               ::   context
      !
      CALL set_xios_context(kiomid, context)

      IF(context == "NONE") THEN  ! read data using default library
         itime = 1
         IF( PRESENT(ktime) ) itime = ktime
         !
         clname = iom_file(kiomid)%name
         clinfo = '          iom_g0d, file: '//TRIM(clname)//', var: '//TRIM(cdvar)
         !
         IF( kiomid > 0 ) THEN
            idvar = iom_varid( kiomid, cdvar )
            IF( iom_file(kiomid)%nfid > 0 .AND. idvar > 0 ) THEN
               idmspc = iom_file ( kiomid )%ndims( idvar )
               IF( iom_file(kiomid)%luld(idvar) )  idmspc = idmspc - 1
               WRITE(cldmspc , fmt='(i1)') idmspc
               IF( idmspc > 0 )  CALL ctl_stop( TRIM(clinfo), 'When reading to a 0D array, we do not accept data', &
                                    &                         'with 1 or more spatial dimensions: '//cldmspc//' were found.' , &
                                    &                         'Use ncwa -a to suppress the unnecessary dimensions' )
               CALL iom_nf90_get( kiomid, idvar, ztmp_pvar, itime )
               pvar = ztmp_pvar
            ENDIF
         ENDIF
      ELSE
#if defined key_xios
         IF(lwp) WRITE(numout,*) 'XIOS RST READ (0D): ', TRIM(cdvar)
         CALL iom_swap(context)
         CALL xios_recv_field( TRIM(cdvar), pvar)
         CALL iom_swap(cxios_context)
#else
         WRITE(ctmp1,*) 'Can not use XIOS in iom_g0d, file: '//TRIM(clname)//', var:'//TRIM(cdvar)
         CALL ctl_stop( 'iom_g0d', ctmp1 )
#endif
      ENDIF
   END SUBROUTINE iom_g0d_sp

   SUBROUTINE iom_g0d_dp( kiomid, cdvar, pvar, ktime )
      INTEGER         , INTENT(in   )                 ::   kiomid    ! Identifier of the file
      CHARACTER(len=*), INTENT(in   )                 ::   cdvar     ! Name of the variable
      REAL(dp)        , INTENT(  out)                 ::   pvar      ! read field
      INTEGER         , INTENT(in   ),     OPTIONAL   ::   ktime     ! record number
      !
      INTEGER                                         ::   idvar     ! variable id
      INTEGER                                         ::   idmspc    ! number of spatial dimensions
      INTEGER         , DIMENSION(1)                  ::   itime     ! record number
      CHARACTER(LEN=100)                              ::   clinfo    ! info character
      CHARACTER(LEN=100)                              ::   clname    ! file name
      CHARACTER(LEN=1)                                ::   cldmspc   !
      CHARACTER(LEN=lc)                               ::   context
      !
      CALL set_xios_context(kiomid, context)

      IF(context == "NONE") THEN  ! read data using default library
         itime = 1
         IF( PRESENT(ktime) ) itime = ktime
         !
         clname = iom_file(kiomid)%name
         clinfo = '          iom_g0d, file: '//TRIM(clname)//', var: '//TRIM(cdvar)
         !
         IF( kiomid > 0 ) THEN
            idvar = iom_varid( kiomid, cdvar )
            IF( iom_file(kiomid)%nfid > 0 .AND. idvar > 0 ) THEN
               idmspc = iom_file ( kiomid )%ndims( idvar )
               IF( iom_file(kiomid)%luld(idvar) )  idmspc = idmspc - 1
               WRITE(cldmspc , fmt='(i1)') idmspc
               IF( idmspc > 0 )  CALL ctl_stop( TRIM(clinfo), 'When reading to a 0D array, we do not accept data', &
                                    &                         'with 1 or more spatial dimensions: '//cldmspc//' were found.' , &
                                    &                         'Use ncwa -a to suppress the unnecessary dimensions' )
               CALL iom_nf90_get( kiomid, idvar, pvar, itime )
            ENDIF
         ENDIF
      ELSE
#if defined key_xios
         IF(lwp) WRITE(numout,*) 'XIOS RST READ (0D): ', TRIM(cdvar)
         CALL iom_swap(context)
         CALL xios_recv_field( TRIM(cdvar), pvar)
         CALL iom_swap(cxios_context)
#else
         WRITE(ctmp1,*) 'Can not use XIOS in iom_g0d, file: '//TRIM(clname)//', var:'//TRIM(cdvar)
         CALL ctl_stop( 'iom_g0d', ctmp1 )
#endif
      ENDIF
   END SUBROUTINE iom_g0d_dp

   SUBROUTINE iom_g1d_sp( kiomid, kdom, cdvar, pvar, ktime, kstart, kcount )
      INTEGER         , INTENT(in   )                         ::   kiomid    ! Identifier of the file
      INTEGER         , INTENT(in   )                         ::   kdom      ! Type of domain to be read
      CHARACTER(len=*), INTENT(in   )                         ::   cdvar     ! Name of the variable
      REAL(sp)        , INTENT(  out), DIMENSION(:)           ::   pvar      ! read field
      INTEGER         , INTENT(in   )              , OPTIONAL ::   ktime     ! record number
      INTEGER         , INTENT(in   ), DIMENSION(1), OPTIONAL ::   kstart    ! start axis position of the reading
      INTEGER         , INTENT(in   ), DIMENSION(1), OPTIONAL ::   kcount    ! number of points in each axis
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 ) CALL iom_get_123d( kiomid, kdom, cdvar, pvsp1d = pvar,   &
               &                                            ktime = ktime, kstart = kstart, kcount = kcount )
      ENDIF
   END SUBROUTINE iom_g1d_sp


   SUBROUTINE iom_g1d_dp( kiomid, kdom, cdvar, pvar, ktime, kstart, kcount )
      INTEGER         , INTENT(in   )                         ::   kiomid    ! Identifier of the file
      INTEGER         , INTENT(in   )                         ::   kdom      ! Type of domain to be read
      CHARACTER(len=*), INTENT(in   )                         ::   cdvar     ! Name of the variable
      REAL(dp)        , INTENT(  out), DIMENSION(:)           ::   pvar      ! read field
      INTEGER         , INTENT(in   )              , OPTIONAL ::   ktime     ! record number
      INTEGER         , INTENT(in   ), DIMENSION(1), OPTIONAL ::   kstart    ! start axis position of the reading
      INTEGER         , INTENT(in   ), DIMENSION(1), OPTIONAL ::   kcount    ! number of points in each axis
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 ) CALL iom_get_123d( kiomid, kdom, cdvar, pvdp1d = pvar,   &
            &                                               ktime = ktime, kstart = kstart, kcount = kcount)
      ENDIF
   END SUBROUTINE iom_g1d_dp

   SUBROUTINE iom_g2d_sp( kiomid, kdom, cdvar, pvar, ktime, cd_type, psgn, kfill, kstart, kcount)
      INTEGER         , INTENT(in   )                         ::   kiomid    ! Identifier of the file
      INTEGER         , INTENT(in   )                         ::   kdom      ! Type of domain to be read
      CHARACTER(len=*), INTENT(in   )                         ::   cdvar     ! Name of the variable
      REAL(sp)        , INTENT(  out), DIMENSION(:,:)         ::   pvar      ! read field
      INTEGER         , INTENT(in   )              , OPTIONAL ::   ktime     ! record number
      CHARACTER(len=1), INTENT(in   )              , OPTIONAL ::   cd_type   ! nature of grid-points (T, U, V, F, W)
      REAL(sp)        , INTENT(in   )              , OPTIONAL ::   psgn      ! -1.(1.): (not) change sign across the north fold
      INTEGER         , INTENT(in   )              , OPTIONAL ::   kfill     ! value of kfillmode in lbc_lbk
      INTEGER         , INTENT(in   ), DIMENSION(2), OPTIONAL ::   kstart    ! start axis position of the reading
      INTEGER         , INTENT(in   ), DIMENSION(2), OPTIONAL ::   kcount    ! number of points in each axis
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 ) CALL iom_get_123d( kiomid, kdom, cdvar, pvsp2d = pvar,                 &
            &                                               cd_type = cd_type, psgn_sp = psgn, kfill = kfill,   &
            &                                               ktime = ktime, kstart = kstart, kcount = kcount )
      ENDIF
   END SUBROUTINE iom_g2d_sp

   SUBROUTINE iom_g2d_dp( kiomid, kdom, cdvar, pvar, ktime, cd_type, psgn, kfill, kstart, kcount)
      INTEGER         , INTENT(in   )                         ::   kiomid    ! Identifier of the file
      INTEGER         , INTENT(in   )                         ::   kdom      ! Type of domain to be read
      CHARACTER(len=*), INTENT(in   )                         ::   cdvar     ! Name of the variable
      REAL(dp)        , INTENT(  out), DIMENSION(:,:)         ::   pvar      ! read field
      INTEGER         , INTENT(in   )              , OPTIONAL ::   ktime     ! record number
      CHARACTER(len=1), INTENT(in   )              , OPTIONAL ::   cd_type   ! nature of grid-points (T, U, V, F, W)
      REAL(dp)        , INTENT(in   )              , OPTIONAL ::   psgn      ! -1.(1.): (not) change sign across the north fold
      INTEGER         , INTENT(in   )              , OPTIONAL ::   kfill     ! value of kfillmode in lbc_lbk
      INTEGER         , INTENT(in   ), DIMENSION(2), OPTIONAL ::   kstart    ! start axis position of the reading
      INTEGER         , INTENT(in   ), DIMENSION(2), OPTIONAL ::   kcount    ! number of points in each axis
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 ) CALL iom_get_123d( kiomid, kdom, cdvar, pvdp2d = pvar,                 &
            &                                               cd_type = cd_type, psgn_dp = psgn, kfill = kfill,   &
            &                                               ktime = ktime, kstart = kstart, kcount = kcount )
      ENDIF
   END SUBROUTINE iom_g2d_dp

   SUBROUTINE iom_g3d_sp( kiomid, kdom, cdvar, pvar, ktime, cd_type, psgn, kfill, kstart, kcount )
      INTEGER         , INTENT(in   )                         ::   kiomid    ! Identifier of the file
      INTEGER         , INTENT(in   )                         ::   kdom      ! Type of domain to be read
      CHARACTER(len=*), INTENT(in   )                         ::   cdvar     ! Name of the variable
      REAL(sp)        , INTENT(  out), DIMENSION(:,:,:)       ::   pvar      ! read field
      INTEGER         , INTENT(in   )              , OPTIONAL ::   ktime     ! record number
      CHARACTER(len=1), INTENT(in   )              , OPTIONAL ::   cd_type   ! nature of grid-points (T, U, V, F, W)
      REAL(sp)        , INTENT(in   )              , OPTIONAL ::   psgn      ! -1.(1.) : (not) change sign across the north fold
      INTEGER         , INTENT(in   )              , OPTIONAL ::   kfill     ! value of kfillmode in lbc_lbk
      INTEGER         , INTENT(in   ), DIMENSION(3), OPTIONAL ::   kstart    ! start axis position of the reading
      INTEGER         , INTENT(in   ), DIMENSION(3), OPTIONAL ::   kcount    ! number of points in each axis
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 ) CALL iom_get_123d( kiomid, kdom, cdvar, pvsp3d = pvar,                 &
            &                                               cd_type = cd_type, psgn_sp = psgn, kfill = kfill,   &
            &                                               ktime = ktime, kstart = kstart, kcount = kcount )
      ENDIF
   END SUBROUTINE iom_g3d_sp

   SUBROUTINE iom_g3d_dp( kiomid, kdom, cdvar, pvar, ktime, cd_type, psgn, kfill, kstart, kcount )
      INTEGER         , INTENT(in   )                         ::   kiomid    ! Identifier of the file
      INTEGER         , INTENT(in   )                         ::   kdom      ! Type of domain to be read
      CHARACTER(len=*), INTENT(in   )                         ::   cdvar     ! Name of the variable
      REAL(dp)        , INTENT(  out), DIMENSION(:,:,:)       ::   pvar      ! read field
      INTEGER         , INTENT(in   )              , OPTIONAL ::   ktime     ! record number
      CHARACTER(len=1), INTENT(in   )              , OPTIONAL ::   cd_type   ! nature of grid-points (T, U, V, F, W)
      REAL(dp)        , INTENT(in   )              , OPTIONAL ::   psgn      ! -1.(1.) : (not) change sign across the north fold
      INTEGER         , INTENT(in   )              , OPTIONAL ::   kfill     ! value of kfillmode in lbc_lbk
      INTEGER         , INTENT(in   ), DIMENSION(3), OPTIONAL ::   kstart    ! start axis position of the reading
      INTEGER         , INTENT(in   ), DIMENSION(3), OPTIONAL ::   kcount    ! number of points in each axis
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 ) CALL iom_get_123d( kiomid, kdom, cdvar, pvdp3d = pvar,                 &
            &                                               cd_type = cd_type, psgn_dp = psgn, kfill = kfill,   &
            &                                               ktime = ktime, kstart = kstart, kcount = kcount )
      ENDIF
   END SUBROUTINE iom_g3d_dp

   !!----------------------------------------------------------------------

   SUBROUTINE iom_get_123d( kiomid , kdom, cdvar, pvsp1d, pvsp2d, pvsp3d, pvdp1d, pvdp2d, pvdp3d,   &
         &                  ktime , cd_type, psgn_sp, psgn_dp, kfill, kstart, kcount )
      !!-----------------------------------------------------------------------
      !!                  ***  ROUTINE  iom_get_123d  ***
      !!
      !! ** Purpose : read a 1D/2D/3D variable
      !!
      !! ** Method : read ONE record at each CALL
      !!-----------------------------------------------------------------------
      INTEGER                    , INTENT(in   )           ::   kiomid    ! Identifier of the file
      INTEGER                    , INTENT(in   )           ::   kdom      ! Type of domain to be read
      CHARACTER(len=*)           , INTENT(in   )           ::   cdvar     ! Name of the variable
      REAL(sp), DIMENSION(:)     , INTENT(  out), OPTIONAL ::   pvsp1d    ! read field (1D case), single precision
      REAL(sp), DIMENSION(:,:)   , INTENT(  out), OPTIONAL ::   pvsp2d    ! read field (2D case), single precision
      REAL(sp), DIMENSION(:,:,:) , INTENT(  out), OPTIONAL ::   pvsp3d    ! read field (3D case), single precision
      REAL(dp), DIMENSION(:)     , INTENT(  out), OPTIONAL ::   pvdp1d    ! read field (1D case), double precision
      REAL(dp), DIMENSION(:,:)   , INTENT(  out), OPTIONAL ::   pvdp2d    ! read field (2D case), double precision
      REAL(dp), DIMENSION(:,:,:) , INTENT(  out), OPTIONAL ::   pvdp3d    ! read field (3D case), double precision
      INTEGER                    , INTENT(in   ), OPTIONAL ::   ktime     ! record number
      CHARACTER(len=1)           , INTENT(in   ), OPTIONAL ::   cd_type   ! nature of grid-points (T, U, V, F, W)
      REAL(sp)                   , INTENT(in   ), OPTIONAL ::   psgn_sp   ! -1.(1.) : (not) change sign across the north fold
      REAL(dp)                   , INTENT(in   ), OPTIONAL ::   psgn_dp   ! -1.(1.) : (not) change sign across the north fold
      INTEGER                    , INTENT(in   ), OPTIONAL ::   kfill     ! value of kfillmode in lbc_lbk
      INTEGER , DIMENSION(:)     , INTENT(in   ), OPTIONAL ::   kstart    ! start position of the reading in each axis
      INTEGER , DIMENSION(:)     , INTENT(in   ), OPTIONAL ::   kcount    ! number of points to be read in each axis
      !
      LOGICAL                        ::   llok        ! true if ok!
      INTEGER                        ::   jl          ! loop on number of dimension
      INTEGER                        ::   idom        ! type of domain
      INTEGER                        ::   idvar       ! id of the variable
      INTEGER                        ::   inbdim      ! number of dimensions of the variable
      INTEGER                        ::   idmspc      ! number of spatial dimensions
      INTEGER                        ::   itime       ! record number
      INTEGER                        ::   istop       ! temporary value of nstop
      INTEGER                        ::   ix1, ix2, iy1, iy2   ! subdomain indexes
      INTEGER                        ::   ji, jj      ! loop counters
      INTEGER                        ::   irankpv     !
      INTEGER                        ::   ind1, ind2  ! substring index
      INTEGER, DIMENSION(jpmax_dims) ::   istart      ! starting point to read for each axis
      INTEGER, DIMENSION(jpmax_dims) ::   icnt        ! number of value to read along each axis
      INTEGER, DIMENSION(jpmax_dims) ::   idimsz      ! size of the dimensions of the variable
      INTEGER, DIMENSION(jpmax_dims) ::   ishape      ! size of the dimensions of the variable
      REAL(sp)                       ::   zscf_sp, zofs_sp  ! sacle_factor and add_offset, single precision
      REAL(dp)                       ::   zscf_dp, zofs_dp  ! sacle_factor and add_offset, double precision
      REAL(sp)                       ::   zsgn_sp     ! local value of psgn, single precision
      REAL(dp)                       ::   zsgn_dp     ! local value of psgn, double precision
      INTEGER                        ::   itmp        ! temporary integer
      CHARACTER(LEN=256)             ::   clinfo      ! info character
      CHARACTER(LEN=256)             ::   clname      ! file name
      CHARACTER(LEN=1)               ::   clrankpv, cldmspc      !
      CHARACTER(LEN=1)               ::   cl_type     ! local value of cd_type
      LOGICAL                        ::   ll_only3rd  ! T => if kstart, kcount present then *only* use values for 3rd spatial dimension.
      LOGICAL                        ::   llis1d, llis2d, llis3d
      LOGICAL                        ::   llsp        ! use single precision
      INTEGER                        ::   inlev       ! number of levels for 3D data
      !---------------------------------------------------------------------
      CHARACTER(LEN=lc)                               ::   context
      !
      CALL set_xios_context(kiomid, context)
      !
      llsp = PRESENT(pvsp1d) .OR. PRESENT(pvsp2d) .OR. PRESENT(pvsp3d)
      IF( llsp ) THEN 
         llis1d = PRESENT(pvsp1d)   ;   IF( llis1d )   ishape(1:1) = SHAPE(pvsp1d)
         llis2d = PRESENT(pvsp2d)   ;   IF( llis2d )   ishape(1:2) = SHAPE(pvsp2d)
         llis3d = PRESENT(pvsp3d)   ;   IF( llis3d )   ishape(1:3) = SHAPE(pvsp3d)
      ELSE
         llis1d = PRESENT(pvdp1d)   ;   IF( llis1d )   ishape(1:1) = SHAPE(pvdp1d)
         llis2d = PRESENT(pvdp2d)   ;   IF( llis2d )   ishape(1:2) = SHAPE(pvdp2d)
         llis3d = PRESENT(pvdp3d)   ;   IF( llis3d )   ishape(1:3) = SHAPE(pvdp3d)
      ENDIF
      inlev = -1
      IF( llis3d )   inlev = ishape(3)
      !
      idom = kdom
      istop = nstop
      !
      IF(context == "NONE") THEN
         clname = iom_file(kiomid)%name   !   esier to read
         clinfo = '          iom_get_123d, file: '//TRIM(clname)//', var: '//TRIM(cdvar)
         ! check kcount and kstart optionals parameters...
         IF( PRESENT(kcount) .AND. .NOT. PRESENT(kstart) ) CALL ctl_stop(TRIM(clinfo), 'kcount present needs kstart present')
         IF( PRESENT(kstart) .AND. .NOT. PRESENT(kcount) ) CALL ctl_stop(TRIM(clinfo), 'kstart present needs kcount present')
         IF( PRESENT(kstart) .AND. idom /= jpdom_unknown .AND. idom /= jpdom_auto_xy ) &
            &          CALL ctl_stop(TRIM(clinfo), 'kstart present needs idom = jpdom_unknown or idom = jpdom_auto_xy')
         IF( idom == jpdom_auto_xy .AND. .NOT. PRESENT(kstart) ) &
            &          CALL ctl_stop(TRIM(clinfo), 'idom = jpdom_auto_xy requires kstart to be present')
         !
         ! Search for the variable in the data base (eventually actualize data)
         !
         idvar = iom_varid( kiomid, cdvar )
         IF( idvar > 0 ) THEN
            !
            idimsz(:) = iom_file(kiomid)%dimsz(:, idvar)      ! to write iom_file(kiomid)%dimsz in a shorter way
            inbdim = iom_file(kiomid)%ndims(idvar)            ! number of dimensions in the file
            idmspc = inbdim                                   ! number of spatial dimensions in the file
            IF( iom_file(kiomid)%luld(idvar) )   idmspc = inbdim - 1
            IF( idmspc > 3 )   CALL ctl_stop(TRIM(clinfo), 'the file has more than 3 spatial dimensions this case is not coded...')
            !
            ! Identify the domain in case of jpdom_auto definition
            IF( idom == jpdom_auto .OR. idom == jpdom_auto_xy ) THEN
               idom = jpdom_global   ! default
               ! else: if the file name finishes with _xxxx.nc with xxxx any number
               ind1 = INDEX( clname, '_', back = .TRUE. ) + 1
               ind2 = INDEX( clname, '.', back = .TRUE. ) - 1
               IF( ind2 > ind1 ) THEN   ;   IF( VERIFY( clname(ind1:ind2), '0123456789' ) == 0 )   idom = jpdom_local   ;   ENDIF
            ENDIF
            !
            ! check the consistency between input array and data rank in the file
            !
            ! initializations
            itime = 1
            IF( PRESENT(ktime) ) itime = ktime
            !
            irankpv = 1 * COUNT( (/ llis1d /) ) + 2 * COUNT( (/ llis2d /) ) + 3 * COUNT( (/ llis3d /) )
            WRITE(clrankpv, fmt='(i1)') irankpv
            WRITE(cldmspc , fmt='(i1)') idmspc
            !
            IF(     idmspc <  irankpv ) THEN                     ! it seems we want to read more than we can...
               IF(     llis3d .AND. idmspc == 2 ) THEN           !   3D input array from 2D spatial data in the file:
                  llok = inlev == 1                              !     -> 3rd dimension must be equal to 1
               ELSEIF( llis3d .AND. idmspc == 1 ) THEN           !   3D input array from 1D spatial data in the file:
                  llok = inlev == 1 .AND. ishape(2) == 1         !     -> 2nd and 3rd dimensions must be equal to 1
               ELSEIF( llis3d .AND. idmspc == 2 ) THEN           !   2D input array from 1D spatial data in the file:
                  llok = ishape(2) == 1                          !     -> 2nd dimension must be equal to 1
               ELSE
                  llok = .FALSE.
               ENDIF
               IF( .NOT. llok )   CALL ctl_stop( TRIM(clinfo), 'The file has only '//cldmspc//' spatial dimension',   &
                  &                                            '=> cannot read a true '//clrankpv//'D array from this file...' )
            ELSEIF( idmspc == irankpv ) THEN
               IF( llis1d .AND. idom /= jpdom_unknown )   &
                  &   CALL ctl_stop( TRIM(clinfo), 'case not coded...You must use jpdom_unknown' )
            ELSEIF( idmspc >  irankpv ) THEN                     ! it seems we want to read less than we should...
                  IF( llis2d .AND. itime == 1 .AND. idimsz(3) == 1 .AND. idmspc == 3 ) THEN
                     CALL ctl_warn( TRIM(clinfo), '2D array input but 3 spatial dimensions in the file...'              ,   &
                           &         'As the size of the z dimension is 1 and as we try to read the first record, ',   &
                           &         'we accept this case, even if there is a possible mix-up between z and time dimension' )
                     idmspc = idmspc - 1
                  !!GS: possibility to read 3D ABL atmopsheric forcing and use 1st level to force BULK simulation
                  !ELSE
                  !   CALL ctl_stop( TRIM(clinfo), 'To keep iom lisibility, when reading a '//clrankpv//'D array,',   &
                  !      &                         'we do not accept data with '//cldmspc//' spatial dimensions'  ,   &
                  !      &                         'Use ncwa -a to suppress the unnecessary dimensions' )
                  ENDIF
            ENDIF
            !
            ! definition of istart and icnt
            !
            icnt  (:) = 1              ! default definition (simple way to deal with special cases listed above)
            istart(:) = 1              ! default definition (simple way to deal with special cases listed above)
            istart(idmspc+1) = itime   ! temporal dimenstion
            !
            IF( idom == jpdom_unknown ) THEN
               IF( PRESENT(kstart) .AND. idom /= jpdom_auto_xy ) THEN
                  istart(1:idmspc) = kstart(1:idmspc)
                  icnt  (1:idmspc) = kcount(1:idmspc)
               ELSE
                  icnt  (1:idmspc) = idimsz(1:idmspc)
               ENDIF
            ELSE   !   not a 1D array as pv(sd)p1d requires jpdom_unknown
               ! we do not read the overlap and the extra-halos -> from Nis0 to Nie0 and from Njs0 to Nje0
               IF( idom == jpdom_global )   istart(1:2) = (/ mig0(Nis0), mjg0(Njs0) /)
               icnt(1:2) = (/ Ni_0, Nj_0 /)
               IF( llis3d ) THEN
                  IF( idom == jpdom_auto_xy ) THEN
                     istart(3) = kstart(3)
                     icnt  (3) = kcount(3)
                  ELSE
                     icnt  (3) = inlev
                  ENDIF
               ENDIF
            ENDIF
            !
            ! check that istart and icnt can be used with this file
            !-
            DO jl = 1, jpmax_dims
               itmp = istart(jl)+icnt(jl)-1
               IF( itmp > idimsz(jl) .AND. idimsz(jl) /= 0 ) THEN
                  WRITE( ctmp1, FMT="('(istart(', i1, ') + icnt(', i1, ') - 1) = ', i5)" ) jl, jl, itmp
                  WRITE( ctmp2, FMT="(' is larger than idimsz(', i1,') = ', i5)"         ) jl, idimsz(jl)
                  CALL ctl_stop( TRIM(clinfo), 'start and count too big regarding to the size of the data, ', ctmp1, ctmp2 )
               ENDIF
            END DO
            !
            ! check that icnt matches the input array
            !-
            IF( idom == jpdom_unknown ) THEN
               ctmp1 = 'd'
            ELSE                    ! we must redefine ishape as we don't read the full array
               IF( llis2d ) THEN
                  IF( llsp ) THEN   ;   ishape(1:2) = SHAPE(pvsp2d(Nis0:Nie0,Njs0:Nje0  ))
                  ELSE              ;   ishape(1:2) = SHAPE(pvdp2d(Nis0:Nie0,Njs0:Nje0  ))
                  ENDIF
                  ctmp1 = 'd(Nis0:Nie0,Njs0:Nje0)'
               ENDIF
               IF( llis3d ) THEN
                  IF( llsp ) THEN   ;   ishape(1:3) = SHAPE(pvsp3d(Nis0:Nie0,Njs0:Nje0,:))
                  ELSE              ;   ishape(1:3) = SHAPE(pvdp3d(Nis0:Nie0,Njs0:Nje0,:))
                  ENDIF
                  ctmp1 = 'd(Nis0:Nie0,Njs0:Nje0,:)'
               ENDIF
            ENDIF
            DO jl = 1, irankpv
               WRITE( ctmp2, FMT="(', ', i1,'): ', i5,' /= icnt(', i1,'):', i5)" ) jl, ishape(jl), jl, icnt(jl)
               IF( llsp ) THEN
                  IF( ishape(jl) /= icnt(jl) )   CALL ctl_stop( TRIM(clinfo), 'size(pvsp'//clrankpv//TRIM(ctmp1)//TRIM(ctmp2) )
               ELSE
                  IF( ishape(jl) /= icnt(jl) )   CALL ctl_stop( TRIM(clinfo), 'size(pvdp'//clrankpv//TRIM(ctmp1)//TRIM(ctmp2) )
               ENDIF
            END DO

         ENDIF

         ! read the data
         !-
         IF( idvar > 0 .AND. istop == nstop ) THEN   ! no additional errors until this point...
            !
            ! find the right index of the array to be read
            IF( idom /= jpdom_unknown ) THEN   ;   ix1 = Nis0   ;   ix2 = Nie0      ;   iy1 = Njs0   ;   iy2 = Nje0
            ELSE                               ;   ix1 = 1      ;   ix2 = icnt(1)   ;   iy1 = 1      ;   iy2 = icnt(2)
            ENDIF

            CALL iom_nf90_get( kiomid, idvar, inbdim, istart, icnt, ix1, ix2, iy1, iy2,   &
               &               pvsp1d, pvsp2d, pvsp3d, pvdp1d, pvdp2d, pvdp3d )

            IF( istop == nstop ) THEN   ! no additional errors until this point...
               IF(lwp) WRITE(numout,"(10x,' read ',a,' (rec: ',i6,') in ',a,' ok')") TRIM(cdvar), itime, TRIM(iom_file(kiomid)%name)

               cl_type = 'T'
               IF( PRESENT(cd_type) )   cl_type = cd_type
               IF( llsp ) THEN
                  zsgn_sp = 1._sp
                  IF( PRESENT(psgn_sp) )   zsgn_sp = psgn_sp
                  IF(     llis2d .AND. idom /= jpdom_unknown .AND. cl_type /= 'Z' ) THEN
                     CALL lbc_lnk( 'iom', pvsp2d, cl_type, zsgn_sp, kfillmode = kfill )
                  ELSEIF( llis3d .AND. idom /= jpdom_unknown .AND. cl_type /= 'Z' ) THEN
                     CALL lbc_lnk( 'iom', pvsp3d, cl_type, zsgn_sp, kfillmode = kfill )
                  ENDIF
               ELSE
                  zsgn_dp = 1._dp
                  IF( PRESENT(psgn_dp) )   zsgn_dp = psgn_dp
                  IF(     llis2d .AND. idom /= jpdom_unknown .AND. cl_type /= 'Z' ) THEN
                     CALL lbc_lnk( 'iom', pvdp2d, cl_type, zsgn_dp, kfillmode = kfill )
                  ELSEIF( llis3d .AND. idom /= jpdom_unknown .AND. cl_type /= 'Z' ) THEN
                     CALL lbc_lnk( 'iom', pvdp3d, cl_type, zsgn_dp, kfillmode = kfill )
                  ENDIF
               ENDIF
               !--- overlap areas and extra hallows (mpp)
               !
            ELSE
               ! return if istop == nstop is false
               RETURN
            ENDIF
         ELSE
            ! return if statment idvar > 0 .AND. istop == nstop is false
            RETURN
         ENDIF
         !
      ELSE        ! read using XIOS. Only if key_xios is defined
#if defined key_xios
!would be good to be able to check which context is active and swap only if current is not restart
         idvar = iom_varid( kiomid, cdvar )
         CALL iom_swap(context)
         IF( llsp ) THEN
            zsgn_sp = 1._sp   ;   IF( PRESENT(psgn_sp) )   zsgn_sp = psgn_sp
         ELSE
            zsgn_dp = 1._dp   ;   IF( PRESENT(psgn_dp) )   zsgn_dp = psgn_dp
         ENDIF
         cl_type = 'T'
         IF( PRESENT(cd_type) )   cl_type = cd_type

         IF(     llis3d ) THEN
            IF(lwp) WRITE(numout,*) 'XIOS RST READ (3D): ',TRIM(cdvar)
            IF( llsp ) THEN   ;   CALL xios_recv_field( TRIM(cdvar), pvsp3d )
            ELSE              ;   CALL xios_recv_field( TRIM(cdvar), pvdp3d )
            ENDIF
            IF(idom /= jpdom_unknown .AND. cl_type /= 'Z' ) THEN
               IF( llsp ) THEN   ;   CALL lbc_lnk( 'iom', pvsp3d, cl_type, zsgn_sp, kfillmode = kfill)
               ELSE              ;   CALL lbc_lnk( 'iom', pvdp3d, cl_type, zsgn_dp, kfillmode = kfill)
               ENDIF
            ENDIF
         ELSEIF( llis2d ) THEN
            IF(lwp) WRITE(numout,*) 'XIOS RST READ (2D): ', TRIM(cdvar)
            IF( llsp ) THEN   ;   CALL xios_recv_field( TRIM(cdvar), pvsp2d )
            ELSE              ;   CALL xios_recv_field( TRIM(cdvar), pvdp2d )
            ENDIF
            IF(idom /= jpdom_unknown .AND. cl_type /= 'Z' ) THEN
               IF( llsp ) THEN   ;   CALL lbc_lnk('iom', pvsp2d, cl_type, zsgn_sp, kfillmode = kfill)
               ELSE              ;   CALL lbc_lnk('iom', pvdp2d, cl_type, zsgn_dp, kfillmode = kfill)
               ENDIF
            ENDIF
         ELSEIF( llis1d ) THEN
            IF(lwp) WRITE(numout,*) 'XIOS RST READ (1D): ', TRIM(cdvar)
            IF( llsp ) THEN   ;   CALL xios_recv_field( TRIM(cdvar), pvsp1d )
            ELSE              ;   CALL xios_recv_field( TRIM(cdvar), pvdp1d )
            ENDIF
         ENDIF
         CALL iom_swap(cxios_context)
#else
         istop = istop + 1
         clinfo = 'Can not use XIOS in iom_get_123d, file: '//TRIM(clname)//', var:'//TRIM(cdvar)
#endif
      ENDIF

      !--- Apply scale_factor and offset
      IF( llsp ) THEN
         zscf_sp = iom_file(kiomid)%scf(idvar)      ! scale factor
         zofs_sp = iom_file(kiomid)%ofs(idvar)      ! offset
         IF(     llis1d ) THEN
            IF( zscf_sp /= 1._sp )   pvsp1d(:    ) = pvsp1d(:    ) * zscf_sp
            IF( zofs_sp /= 0._sp )   pvsp1d(:    ) = pvsp1d(:    ) + zofs_sp
         ELSEIF( llis2d ) THEN
            IF( zscf_sp /= 1._sp )   pvsp2d(:,:  ) = pvsp2d(:,:  ) * zscf_sp
            IF( zofs_sp /= 0._sp )   pvsp2d(:,:  ) = pvsp2d(:,:  ) + zofs_sp
         ELSEIF( llis3d ) THEN
            IF( zscf_sp /= 1._sp )   pvsp3d(:,:,:) = pvsp3d(:,:,:) * zscf_sp
            IF( zofs_sp /= 0._sp )   pvsp3d(:,:,:) = pvsp3d(:,:,:) + zofs_sp
         ENDIF
      ELSE
         zscf_dp = iom_file(kiomid)%scf(idvar)      ! scale factor
         zofs_dp = iom_file(kiomid)%ofs(idvar)      ! offset
         IF(     llis1d ) THEN
            IF( zscf_dp /= 1._dp )   pvdp1d(:    ) = pvdp1d(:    ) * zscf_dp
            IF( zofs_dp /= 0._dp )   pvdp1d(:    ) = pvdp1d(:    ) + zofs_dp
         ELSEIF( llis2d ) THEN
            IF( zscf_dp /= 1._dp )   pvdp2d(:,:  ) = pvdp2d(:,:  ) * zscf_dp
            IF( zofs_dp /= 0._dp )   pvdp2d(:,:  ) = pvdp2d(:,:  ) + zofs_dp
         ELSEIF( llis3d ) THEN
            IF( zscf_dp /= 1._dp )   pvdp3d(:,:,:) = pvdp3d(:,:,:) * zscf_dp
            IF( zofs_dp /= 0._dp )   pvdp3d(:,:,:) = pvdp3d(:,:,:) + zofs_dp
         ENDIF
      ENDIF
      !
   END SUBROUTINE iom_get_123d

   SUBROUTINE iom_get_var( cdname, z2d)
      CHARACTER(LEN=*), INTENT(in ) ::   cdname
      REAL(wp), DIMENSION(jpi,jpj) ::   z2d
#if defined key_xios
      IF( xios_field_is_active( cdname, at_current_timestep_arg = .TRUE. ) ) THEN
         z2d(:,:) = 0._wp
         CALL xios_recv_field( cdname, z2d)
      ENDIF
#else
      IF( .FALSE. )   WRITE(numout,*) cdname, z2d ! useless test to avoid compilation warnings
#endif
   END SUBROUTINE iom_get_var


   FUNCTION iom_getszuld ( kiomid )
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION  iom_getszuld  ***
      !!
      !! ** Purpose : get the size of the unlimited dimension in a file
      !!              (return -1 if not found)
      !!-----------------------------------------------------------------------
      INTEGER, INTENT(in   ) ::   kiomid   ! file Identifier
      !
      INTEGER                ::   iom_getszuld
      !!-----------------------------------------------------------------------
      iom_getszuld = -1
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%iduld > 0 )   iom_getszuld = iom_file(kiomid)%lenuld
      ENDIF
   END FUNCTION iom_getszuld


   !!----------------------------------------------------------------------
   !!                   INTERFACE iom_chkatt
   !!----------------------------------------------------------------------
   SUBROUTINE iom_chkatt( kiomid, cdatt, llok, ksize, cdvar )
      INTEGER         , INTENT(in   )                 ::   kiomid    ! Identifier of the file
      CHARACTER(len=*), INTENT(in   )                 ::   cdatt     ! Name of the attribute
      LOGICAL         , INTENT(  out)                 ::   llok      ! Error code
      INTEGER         , INTENT(  out), OPTIONAL       ::   ksize     ! Size of the attribute array
      CHARACTER(len=*), INTENT(in   ), OPTIONAL       ::   cdvar     ! Name of the variable
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 )   CALL iom_nf90_chkatt( kiomid, cdatt, llok, ksize=ksize, cdvar=cdvar )
      ENDIF
      !
   END SUBROUTINE iom_chkatt

   !!----------------------------------------------------------------------
   !!                   INTERFACE iom_getatt
   !!----------------------------------------------------------------------
   SUBROUTINE iom_g0d_iatt( kiomid, cdatt, katt0d, cdvar )
      INTEGER               , INTENT(in   )           ::   kiomid    ! Identifier of the file
      CHARACTER(len=*)      , INTENT(in   )           ::   cdatt     ! Name of the attribute
      INTEGER               , INTENT(  out)           ::   katt0d    ! read field
      CHARACTER(len=*)      , INTENT(in   ), OPTIONAL ::   cdvar     ! Name of the variable
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 )   CALL iom_nf90_getatt( kiomid, cdatt,  katt0d =  katt0d, cdvar=cdvar )
      ENDIF
   END SUBROUTINE iom_g0d_iatt

   SUBROUTINE iom_g1d_iatt( kiomid, cdatt, katt1d, cdvar )
      INTEGER               , INTENT(in   )           ::   kiomid    ! Identifier of the file
      CHARACTER(len=*)      , INTENT(in   )           ::   cdatt     ! Name of the attribute
      INTEGER, DIMENSION(:) , INTENT(  out)           ::   katt1d    ! read field
      CHARACTER(len=*)      , INTENT(in   ), OPTIONAL ::   cdvar     ! Name of the variable
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 )   CALL iom_nf90_getatt( kiomid, cdatt,  katt1d =  katt1d, cdvar=cdvar )
      ENDIF
   END SUBROUTINE iom_g1d_iatt

   SUBROUTINE iom_g0d_ratt( kiomid, cdatt, patt0d, cdvar )
      INTEGER               , INTENT(in   )           ::   kiomid    ! Identifier of the file
      CHARACTER(len=*)      , INTENT(in   )           ::   cdatt     ! Name of the attribute
      REAL(wp)              , INTENT(  out)           ::   patt0d    ! read field
      CHARACTER(len=*)      , INTENT(in   ), OPTIONAL ::   cdvar     ! Name of the variable
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 )   CALL iom_nf90_getatt( kiomid, cdatt,  patt0d =  patt0d, cdvar=cdvar )
      ENDIF
   END SUBROUTINE iom_g0d_ratt

   SUBROUTINE iom_g1d_ratt( kiomid, cdatt, patt1d, cdvar )
      INTEGER               , INTENT(in   )           ::   kiomid    ! Identifier of the file
      CHARACTER(len=*)      , INTENT(in   )           ::   cdatt     ! Name of the attribute
      REAL(wp), DIMENSION(:), INTENT(  out)           ::   patt1d    ! read field
      CHARACTER(len=*)      , INTENT(in   ), OPTIONAL ::   cdvar     ! Name of the variable
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 )   CALL iom_nf90_getatt( kiomid, cdatt,  patt1d =  patt1d, cdvar=cdvar )
      ENDIF
   END SUBROUTINE iom_g1d_ratt

   SUBROUTINE iom_g0d_catt( kiomid, cdatt, cdatt0d, cdvar )
      INTEGER               , INTENT(in   )           ::   kiomid    ! Identifier of the file
      CHARACTER(len=*)      , INTENT(in   )           ::   cdatt     ! Name of the attribute
      CHARACTER(len=*)      , INTENT(  out)           ::   cdatt0d   ! read field
      CHARACTER(len=*)      , INTENT(in   ), OPTIONAL ::   cdvar     ! Name of the variable
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 )   CALL iom_nf90_getatt( kiomid, cdatt, cdatt0d = cdatt0d, cdvar=cdvar )
      ENDIF
   END SUBROUTINE iom_g0d_catt


   !!----------------------------------------------------------------------
   !!                   INTERFACE iom_putatt
   !!----------------------------------------------------------------------
   SUBROUTINE iom_p0d_iatt( kiomid, cdatt, katt0d, cdvar )
      INTEGER               , INTENT(in   )           ::   kiomid    ! Identifier of the file
      CHARACTER(len=*)      , INTENT(in   )           ::   cdatt     ! Name of the attribute
      INTEGER               , INTENT(in   )           ::   katt0d    ! written field
      CHARACTER(len=*)      , INTENT(in   ), OPTIONAL ::   cdvar     ! Name of the variable
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 )   CALL iom_nf90_putatt( kiomid, cdatt,  katt0d =  katt0d, cdvar=cdvar )
      ENDIF
   END SUBROUTINE iom_p0d_iatt

   SUBROUTINE iom_p1d_iatt( kiomid, cdatt, katt1d, cdvar )
      INTEGER               , INTENT(in   )           ::   kiomid    ! Identifier of the file
      CHARACTER(len=*)      , INTENT(in   )           ::   cdatt     ! Name of the attribute
      INTEGER, DIMENSION(:) , INTENT(in   )           ::   katt1d    ! written field
      CHARACTER(len=*)      , INTENT(in   ), OPTIONAL ::   cdvar     ! Name of the variable
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 )   CALL iom_nf90_putatt( kiomid, cdatt,  katt1d =  katt1d, cdvar=cdvar )
      ENDIF
   END SUBROUTINE iom_p1d_iatt

   SUBROUTINE iom_p0d_ratt( kiomid, cdatt, patt0d, cdvar )
      INTEGER               , INTENT(in   )           ::   kiomid    ! Identifier of the file
      CHARACTER(len=*)      , INTENT(in   )           ::   cdatt     ! Name of the attribute
      REAL(wp)              , INTENT(in   )           ::   patt0d    ! written field
      CHARACTER(len=*)      , INTENT(in   ), OPTIONAL ::   cdvar     ! Name of the variable
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 )   CALL iom_nf90_putatt( kiomid, cdatt,  patt0d =  patt0d, cdvar=cdvar )
      ENDIF
   END SUBROUTINE iom_p0d_ratt

   SUBROUTINE iom_p1d_ratt( kiomid, cdatt, patt1d, cdvar )
      INTEGER               , INTENT(in   )           ::   kiomid    ! Identifier of the file
      CHARACTER(len=*)      , INTENT(in   )           ::   cdatt     ! Name of the attribute
      REAL(wp), DIMENSION(:), INTENT(in   )           ::   patt1d    ! written field
      CHARACTER(len=*)      , INTENT(in   ), OPTIONAL ::   cdvar     ! Name of the variable
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 )   CALL iom_nf90_putatt( kiomid, cdatt,  patt1d =  patt1d, cdvar=cdvar )
      ENDIF
   END SUBROUTINE iom_p1d_ratt

   SUBROUTINE iom_p0d_catt( kiomid, cdatt, cdatt0d, cdvar )
      INTEGER               , INTENT(in   )           ::   kiomid    ! Identifier of the file
      CHARACTER(len=*)      , INTENT(in   )           ::   cdatt     ! Name of the attribute
      CHARACTER(len=*)      , INTENT(in   )           ::   cdatt0d   ! written field
      CHARACTER(len=*)      , INTENT(in   ), OPTIONAL ::   cdvar     ! Name of the variable
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 )   CALL iom_nf90_putatt( kiomid, cdatt, cdatt0d = cdatt0d, cdvar=cdvar )
      ENDIF
   END SUBROUTINE iom_p0d_catt


   !!----------------------------------------------------------------------
   !!                   INTERFACE iom_rstput
   !!----------------------------------------------------------------------
   SUBROUTINE iom_rp0d_sp( kt, kwrite, kiomid, cdvar, pvar, ktype )
      INTEGER         , INTENT(in)                         ::   kt       ! ocean time-step
      INTEGER         , INTENT(in)                         ::   kwrite   ! writing time-step
      INTEGER         , INTENT(in)                         ::   kiomid   ! Identifier of the file
      CHARACTER(len=*), INTENT(in)                         ::   cdvar    ! time axis name
      REAL(sp)        , INTENT(in)                         ::   pvar     ! written field
      INTEGER         , INTENT(in), OPTIONAL               ::   ktype    ! variable external type
      !
      LOGICAL           :: llx                ! local xios write flag
      INTEGER           :: ivid   ! variable id
      CHARACTER(LEN=lc) :: context
      !
      CALL iom_rp0123d( kt, kwrite, kiomid, cdvar, ktype, pvsp0d = pvar )
      !
   END SUBROUTINE iom_rp0d_sp

   SUBROUTINE iom_rp0d_dp( kt, kwrite, kiomid, cdvar, pvar, ktype )
      INTEGER         , INTENT(in)                         ::   kt       ! ocean time-step
      INTEGER         , INTENT(in)                         ::   kwrite   ! writing time-step
      INTEGER         , INTENT(in)                         ::   kiomid   ! Identifier of the file
      CHARACTER(len=*), INTENT(in)                         ::   cdvar    ! time axis name
      REAL(dp)        , INTENT(in)                         ::   pvar     ! written field
      INTEGER         , INTENT(in), OPTIONAL               ::   ktype    ! variable external type
      !
      LOGICAL           :: llx                ! local xios write flag
      INTEGER           :: ivid   ! variable id
      CHARACTER(LEN=lc) :: context
      !
      CALL iom_rp0123d( kt, kwrite, kiomid, cdvar, ktype, pvdp0d = pvar )
      !
   END SUBROUTINE iom_rp0d_dp


   SUBROUTINE iom_rp1d_sp( kt, kwrite, kiomid, cdvar, pvar, ktype )
      INTEGER         , INTENT(in)                         ::   kt       ! ocean time-step
      INTEGER         , INTENT(in)                         ::   kwrite   ! writing time-step
      INTEGER         , INTENT(in)                         ::   kiomid   ! Identifier of the file
      CHARACTER(len=*), INTENT(in)                         ::   cdvar    ! time axis name
      REAL(sp)        , INTENT(in), DIMENSION(          :) ::   pvar     ! written field
      INTEGER         , INTENT(in), OPTIONAL               ::   ktype    ! variable external type
      !
      LOGICAL           :: llx                ! local xios write flag
      INTEGER           :: ivid   ! variable id
      CHARACTER(LEN=lc) :: context
      !
      CALL iom_rp0123d( kt, kwrite, kiomid, cdvar, ktype, pvsp1d = pvar )
      !
   END SUBROUTINE iom_rp1d_sp

   SUBROUTINE iom_rp1d_dp( kt, kwrite, kiomid, cdvar, pvar, ktype )
      INTEGER         , INTENT(in)                         ::   kt       ! ocean time-step
      INTEGER         , INTENT(in)                         ::   kwrite   ! writing time-step
      INTEGER         , INTENT(in)                         ::   kiomid   ! Identifier of the file
      CHARACTER(len=*), INTENT(in)                         ::   cdvar    ! time axis name
      REAL(dp)        , INTENT(in), DIMENSION(          :) ::   pvar     ! written field
      INTEGER         , INTENT(in), OPTIONAL               ::   ktype    ! variable external type
      !
      LOGICAL           :: llx                ! local xios write flag
      INTEGER           :: ivid   ! variable id
      CHARACTER(LEN=lc) :: context
      !
      CALL iom_rp0123d( kt, kwrite, kiomid, cdvar, ktype, pvdp1d = pvar )
      !
   END SUBROUTINE iom_rp1d_dp

   
   SUBROUTINE iom_rp2d_sp( kt, kwrite, kiomid, cdvar, pvar, ktype )
      INTEGER         , INTENT(in)                         ::   kt       ! ocean time-step
      INTEGER         , INTENT(in)                         ::   kwrite   ! writing time-step
      INTEGER         , INTENT(in)                         ::   kiomid   ! Identifier of the file
      CHARACTER(len=*), INTENT(in)                         ::   cdvar    ! time axis name
      REAL(sp)        , INTENT(in), DIMENSION(:,    :    ) ::   pvar     ! written field
      INTEGER         , INTENT(in), OPTIONAL               ::   ktype    ! variable external type
      !
      LOGICAL            :: llx
      INTEGER            :: ivid   ! variable id
      CHARACTER(LEN=lc)  :: context
      !
      CALL iom_rp0123d( kt, kwrite, kiomid, cdvar, ktype, pvsp2d = pvar )
      !
   END SUBROUTINE iom_rp2d_sp

   SUBROUTINE iom_rp2d_dp( kt, kwrite, kiomid, cdvar, pvar, ktype )
      INTEGER         , INTENT(in)                         ::   kt       ! ocean time-step
      INTEGER         , INTENT(in)                         ::   kwrite   ! writing time-step
      INTEGER         , INTENT(in)                         ::   kiomid   ! Identifier of the file
      CHARACTER(len=*), INTENT(in)                         ::   cdvar    ! time axis name
      REAL(dp)        , INTENT(in), DIMENSION(:,    :    ) ::   pvar     ! written field
      INTEGER         , INTENT(in), OPTIONAL               ::   ktype    ! variable external type
      !
      LOGICAL           :: llx
      INTEGER           :: ivid   ! variable id
      CHARACTER(LEN=lc) :: context
      !
      CALL iom_rp0123d( kt, kwrite, kiomid, cdvar, ktype, pvdp2d = pvar )
      !
   END SUBROUTINE iom_rp2d_dp


   SUBROUTINE iom_rp3d_sp( kt, kwrite, kiomid, cdvar, pvar, ktype )
      INTEGER         , INTENT(in)                         ::   kt       ! ocean time-step
      INTEGER         , INTENT(in)                         ::   kwrite   ! writing time-step
      INTEGER         , INTENT(in)                         ::   kiomid   ! Identifier of the file
      CHARACTER(len=*), INTENT(in)                         ::   cdvar    ! time axis name
      REAL(sp)        , INTENT(in),       DIMENSION(:,:,:) ::   pvar     ! written field
      INTEGER         , INTENT(in), OPTIONAL               ::   ktype    ! variable external type
      !
      LOGICAL           :: llx                 ! local xios write flag
      INTEGER           :: ivid   ! variable id
      CHARACTER(LEN=lc) :: context
      !
      CALL iom_rp0123d( kt, kwrite, kiomid, cdvar, ktype, pvsp3d = pvar )
      !
   END SUBROUTINE iom_rp3d_sp

   SUBROUTINE iom_rp3d_dp( kt, kwrite, kiomid, cdvar, pvar, ktype )
      INTEGER         , INTENT(in)                         ::   kt       ! ocean time-step
      INTEGER         , INTENT(in)                         ::   kwrite   ! writing time-step
      INTEGER         , INTENT(in)                         ::   kiomid   ! Identifier of the file
      CHARACTER(len=*), INTENT(in)                         ::   cdvar    ! time axis name
      REAL(dp)        , INTENT(in),       DIMENSION(:,:,:) ::   pvar     ! written field
      INTEGER         , INTENT(in), OPTIONAL               ::   ktype    ! variable external type
      !
      LOGICAL           :: llx                 ! local xios write flag
      INTEGER           :: ivid   ! variable id
      CHARACTER(LEN=lc) :: context
      !
      CALL iom_rp0123d( kt, kwrite, kiomid, cdvar, ktype, pvdp3d = pvar )
      !
   END SUBROUTINE iom_rp3d_dp

   SUBROUTINE iom_rp0123d( kt, kwrite, kiomid, cdvar, ktype, pvsp0d, pvsp1d, pvsp2d, pvsp3d, pvdp0d, pvdp1d, pvdp2d, pvdp3d )
      INTEGER                    , INTENT(in)           ::   kt       ! ocean time-step
      INTEGER                    , INTENT(in)           ::   kwrite   ! writing time-step
      INTEGER                    , INTENT(in)           ::   kiomid   ! Identifier of the file
      CHARACTER(len=*)           , INTENT(in)           ::   cdvar    ! time axis name
      INTEGER                    , INTENT(in), OPTIONAL ::   ktype    ! variable external type
      REAL(sp)                   , INTENT(in), OPTIONAL ::   pvsp0d    ! read field (0D case), single precision
      REAL(sp), DIMENSION(:)     , INTENT(in), OPTIONAL ::   pvsp1d    ! read field (1D case), single precision
      REAL(sp), DIMENSION(:,:)   , INTENT(in), OPTIONAL ::   pvsp2d    ! read field (2D case), single precision
      REAL(sp), DIMENSION(:,:,:) , INTENT(in), OPTIONAL ::   pvsp3d    ! read field (3D case), single precision
      REAL(dp)                   , INTENT(in), OPTIONAL ::   pvdp0d    ! read field (0D case), double precision
      REAL(dp), DIMENSION(:)     , INTENT(in), OPTIONAL ::   pvdp1d    ! read field (1D case), double precision
      REAL(dp), DIMENSION(:,:)   , INTENT(in), OPTIONAL ::   pvdp2d    ! read field (2D case), double precision
      REAL(dp), DIMENSION(:,:,:) , INTENT(in), OPTIONAL ::   pvdp3d    ! read field (3D case), double precision
      !
      LOGICAL           :: llx                 ! local xios write flag
      INTEGER           :: ivid   ! variable id
      CHARACTER(LEN=lc) :: context
      !
      CALL set_xios_context(kiomid, context)

      llx = .NOT. (context == "NONE")

      IF( llx ) THEN
#ifdef key_xios
         IF( kt == kwrite ) THEN
            IF(lwp) write(numout,*) 'RESTART: write (XIOS) ',TRIM(cdvar)
            CALL iom_swap(context)
            IF( PRESENT(pvsp0d) )   CALL iom_put(TRIM(cdvar), pvsp0d)
            IF( PRESENT(pvsp1d) )   CALL iom_put(TRIM(cdvar), pvsp1d)
            IF( PRESENT(pvsp2d) )   CALL iom_put(TRIM(cdvar), pvsp2d)
            IF( PRESENT(pvsp3d) )   CALL iom_put(TRIM(cdvar), pvsp3d)
            IF( PRESENT(pvdp0d) )   CALL iom_put(TRIM(cdvar), pvdp0d)
            IF( PRESENT(pvdp1d) )   CALL iom_put(TRIM(cdvar), pvdp1d)
            IF( PRESENT(pvdp2d) )   CALL iom_put(TRIM(cdvar), pvdp2d)
            IF( PRESENT(pvdp3d) )   CALL iom_put(TRIM(cdvar), pvdp3d)
            CALL iom_swap(cxios_context)
         ELSE
            IF(lwp) write(numout,*) 'RESTART: define (XIOS)',TRIM(cdvar)
            CALL iom_swap(context)
            CALL iom_set_rstw_active( TRIM(cdvar), rs0 = pvsp0d, rs1 = pvsp1d, rs2 = pvsp2d, rs3 = pvsp3d   &
               &                                 , rd0 = pvdp0d, rd1 = pvdp1d, rd2 = pvdp2d, rd3 = pvdp3d )
            CALL iom_swap(cxios_context)
         ENDIF
#endif
      ELSE
         IF( kiomid > 0 ) THEN
            IF( iom_file(kiomid)%nfid > 0 ) THEN
               ivid = iom_varid( kiomid, cdvar, ldstop = .FALSE. )
               CALL iom_nf90_rstput( kt, kwrite, kiomid, cdvar, ivid, ktype, pvsp0d, pvsp1d, pvsp2d, pvsp3d,   &
                  &                                                          pvdp0d, pvdp1d, pvdp2d, pvdp3d )
            ENDIF
         ENDIF
      ENDIF
   END SUBROUTINE iom_rp0123d


  SUBROUTINE iom_delay_rst( cdaction, cdcpnt, kncid )
      !!---------------------------------------------------------------------
      !!   Routine iom_delay_rst: used read/write restart related to mpp_delay
      !!
      !!---------------------------------------------------------------------
      CHARACTER(len=*), INTENT(in   ) ::   cdaction        !
      CHARACTER(len=*), INTENT(in   ) ::   cdcpnt
      INTEGER         , INTENT(in   ) ::   kncid
      !
      INTEGER  :: ji
      INTEGER  :: indim
      LOGICAL  :: llattexist
      REAL(wp), ALLOCATABLE, DIMENSION(:) ::   zreal1d
      !!---------------------------------------------------------------------
      !
      !                                      ===================================
      IF( TRIM(cdaction) == 'READ' ) THEN   ! read restart related to mpp_delay !
         !                                   ===================================
         DO ji = 1, nbdelay
            IF ( c_delaycpnt(ji) == cdcpnt ) THEN
               CALL iom_chkatt( kncid, 'DELAY_'//c_delaylist(ji), llattexist, indim )
               IF( llattexist )  THEN
                  ALLOCATE( todelay(ji)%z1d(indim) )
                  CALL iom_getatt( kncid, 'DELAY_'//c_delaylist(ji), todelay(ji)%z1d(:) )
                  ndelayid(ji) = 0   ! set to 0 to specify that the value was read in the restart
               ENDIF
           ENDIF
         END DO
         !                                   ====================================
      ELSE                                  ! write restart related to mpp_delay !
         !                                   ====================================
         DO ji = 1, nbdelay   ! save only ocean delayed global communication variables
            IF ( c_delaycpnt(ji) == cdcpnt ) THEN
               IF( ASSOCIATED(todelay(ji)%z1d) ) THEN
                  CALL mpp_delay_rcv(ji)   ! make sure %z1d is received
                  CALL iom_putatt( kncid, 'DELAY_'//c_delaylist(ji), todelay(ji)%z1d(:) )
               ENDIF
            ENDIF
         END DO
         !
      ENDIF

   END SUBROUTINE iom_delay_rst



   !!----------------------------------------------------------------------
   !!                   INTERFACE iom_put
   !!----------------------------------------------------------------------
   SUBROUTINE iom_p0d_sp( cdname, pfield0d )
      CHARACTER(LEN=*), INTENT(in) ::   cdname
      REAL(sp)        , INTENT(in) ::   pfield0d
      !!      REAL(wp)        , DIMENSION(jpi,jpj) ::   zz     ! masson
#if defined key_xios
!!clem      zz(:,:)=pfield0d
!!clem      CALL xios_send_field(cdname, zz)
      CALL xios_send_field(cdname, (/pfield0d/))
#else
      IF( .FALSE. )   WRITE(numout,*) cdname, pfield0d   ! useless test to avoid compilation warnings
#endif
   END SUBROUTINE iom_p0d_sp

   SUBROUTINE iom_p0d_dp( cdname, pfield0d )
      CHARACTER(LEN=*), INTENT(in) ::   cdname
      REAL(dp)        , INTENT(in) ::   pfield0d
!!      REAL(wp)        , DIMENSION(jpi,jpj) ::   zz     ! masson
#if defined key_xios
!!clem      zz(:,:)=pfield0d
!!clem      CALL xios_send_field(cdname, zz)
      CALL xios_send_field(cdname, (/pfield0d/))
#else
      IF( .FALSE. )   WRITE(numout,*) cdname, pfield0d   ! useless test to avoid compilation warnings
#endif
   END SUBROUTINE iom_p0d_dp


   SUBROUTINE iom_p1d_sp( cdname, pfield1d )
      CHARACTER(LEN=*)          , INTENT(in) ::   cdname
      REAL(sp),     DIMENSION(:), INTENT(in) ::   pfield1d
#if defined key_xios
      CALL xios_send_field( cdname, RESHAPE( (/pfield1d/), (/1,1,SIZE(pfield1d)/) ) )
#else
      IF( .FALSE. )   WRITE(numout,*) cdname, pfield1d   ! useless test to avoid compilation warnings
#endif
   END SUBROUTINE iom_p1d_sp

   SUBROUTINE iom_p1d_dp( cdname, pfield1d )
      CHARACTER(LEN=*)          , INTENT(in) ::   cdname
      REAL(dp),     DIMENSION(:), INTENT(in) ::   pfield1d
#if defined key_xios
      CALL xios_send_field( cdname, RESHAPE( (/pfield1d/), (/1,1,SIZE(pfield1d)/) ) )
#else
      IF( .FALSE. )   WRITE(numout,*) cdname, pfield1d   ! useless test to avoid compilation warnings
#endif
   END SUBROUTINE iom_p1d_dp

   SUBROUTINE iom_p2d_sp( cdname, pfield2d )
      CHARACTER(LEN=*)            , INTENT(in) ::   cdname
      REAL(sp),     DIMENSION(:,:), INTENT(in) ::   pfield2d
      IF( iom_use(cdname) ) THEN
#if defined key_xios
         IF( is_tile(pfield2d) == 1 ) THEN
#if ! defined key_xios3
            CALL xios_send_field( cdname, pfield2d, ntile - 1 )
#endif
         ELSE IF( .NOT. l_istiled .OR. ntile == nijtile ) THEN
            CALL xios_send_field( cdname, pfield2d )
         ENDIF
#else
         WRITE(numout,*) pfield2d   ! iom_use(cdname) = .F. -> useless test to avoid compilation warnings
#endif
      ENDIF
   END SUBROUTINE iom_p2d_sp

   SUBROUTINE iom_p2d_dp( cdname, pfield2d )
      CHARACTER(LEN=*)            , INTENT(in) ::   cdname
      REAL(dp),     DIMENSION(:,:), INTENT(in) ::   pfield2d
      IF( iom_use(cdname) ) THEN
#if defined key_xios
         IF( is_tile(pfield2d) == 1 ) THEN
#if ! defined key_xios3
            CALL xios_send_field( cdname, pfield2d, ntile - 1 )
#endif
         ELSE IF( .NOT. l_istiled .OR. ntile == nijtile ) THEN
            CALL xios_send_field( cdname, pfield2d )
         ENDIF
#else
         WRITE(numout,*) pfield2d   ! iom_use(cdname) = .F. -> useless test to avoid compilation warnings
#endif
      ENDIF
   END SUBROUTINE iom_p2d_dp

   SUBROUTINE iom_p3d_sp( cdname, pfield3d )
      CHARACTER(LEN=*)                , INTENT(in) ::   cdname
      REAL(sp),       DIMENSION(:,:,:), INTENT(in) ::   pfield3d
      IF( iom_use(cdname) ) THEN
#if defined key_xios
         IF( is_tile(pfield3d) == 1 ) THEN
#if ! defined key_xios3
            CALL xios_send_field( cdname, pfield3d, ntile - 1 )
#endif
         ELSE IF( .NOT. l_istiled .OR. ntile == nijtile ) THEN
            CALL xios_send_field( cdname, pfield3d )
         ENDIF
#else
         WRITE(numout,*) pfield3d   ! iom_use(cdname) = .F. -> useless test to avoid compilation warnings
#endif
      ENDIF
   END SUBROUTINE iom_p3d_sp

   SUBROUTINE iom_p3d_dp( cdname, pfield3d )
      CHARACTER(LEN=*)                , INTENT(in) ::   cdname
      REAL(dp),       DIMENSION(:,:,:), INTENT(in) ::   pfield3d
      IF( iom_use(cdname) ) THEN
#if defined key_xios
         IF( is_tile(pfield3d) == 1 ) THEN
#if ! defined key_xios3
            CALL xios_send_field( cdname, pfield3d, ntile - 1 )
#endif
         ELSE IF( .NOT. l_istiled .OR. ntile == nijtile ) THEN
            CALL xios_send_field( cdname, pfield3d )
         ENDIF
#else
         WRITE(numout,*) pfield3d   ! iom_use(cdname) = .F. -> useless test to avoid compilation warnings
#endif
      ENDIF
   END SUBROUTINE iom_p3d_dp

   SUBROUTINE iom_p4d_sp( cdname, pfield4d )
      CHARACTER(LEN=*)                , INTENT(in) ::   cdname
      REAL(sp),       DIMENSION(:,:,:,:), INTENT(in) ::   pfield4d
      IF( iom_use(cdname) ) THEN
#if defined key_xios
         IF( is_tile(pfield4d) == 1 ) THEN
#if ! defined key_xios3
            CALL xios_send_field( cdname, pfield4d, ntile - 1 )
#endif
         ELSE IF( .NOT. l_istiled .OR. ntile == nijtile ) THEN
            CALL xios_send_field( cdname, pfield4d )
         ENDIF
#else
         WRITE(numout,*) pfield4d   ! iom_use(cdname) = .F. -> useless test to avoid compilation warnings
#endif
      ENDIF
   END SUBROUTINE iom_p4d_sp

   SUBROUTINE iom_p4d_dp( cdname, pfield4d )
      CHARACTER(LEN=*)                , INTENT(in) ::   cdname
      REAL(dp),       DIMENSION(:,:,:,:), INTENT(in) ::   pfield4d
      IF( iom_use(cdname) ) THEN
#if defined key_xios
         IF( is_tile(pfield4d) == 1 ) THEN
#if ! defined key_xios3
            CALL xios_send_field( cdname, pfield4d, ntile - 1 )
#endif
         ELSE IF( .NOT. l_istiled .OR. ntile == nijtile ) THEN
            CALL xios_send_field( cdname, pfield4d )
         ENDIF
#else
         WRITE(numout,*) pfield4d   ! iom_use(cdname) = .F. -> useless test to avoid compilation warnings
#endif
      ENDIF
   END SUBROUTINE iom_p4d_dp

#if defined key_xios
   !!----------------------------------------------------------------------
   !!   'key_xios'                                         XIOS interface
   !!----------------------------------------------------------------------

   SUBROUTINE iom_set_domain_attr( cdid, ni_glo, nj_glo, ibegin, jbegin, ni, nj,                                               &
      &                                    data_dim, data_ibegin, data_ni, data_jbegin, data_nj, lonvalue, latvalue, mask,     &
      &                                  ntiles, tile_ibegin, tile_jbegin, tile_ni, tile_nj,                                   &
      &                                  tile_data_ibegin, tile_data_jbegin, tile_data_ni, tile_data_nj,                       &
      &                                    nvertex, bounds_lon, bounds_lat, area )
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)                  , INTENT(in) ::   cdid
      INTEGER                 , OPTIONAL, INTENT(in) ::   ni_glo, nj_glo, ibegin, jbegin, ni, nj
      INTEGER,  DIMENSION(:)  , OPTIONAL, INTENT(in) ::   tile_ibegin, tile_jbegin, tile_ni, tile_nj
      INTEGER,  DIMENSION(:)  , OPTIONAL, INTENT(in) ::   tile_data_ibegin, tile_data_jbegin, tile_data_ni, tile_data_nj
      INTEGER                 , OPTIONAL, INTENT(in) ::   data_dim, data_ibegin, data_ni, data_jbegin, data_nj
      INTEGER                 , OPTIONAL, INTENT(in) ::   nvertex, ntiles
      REAL(dp), DIMENSION(:)  , OPTIONAL, INTENT(in) ::   lonvalue, latvalue
      REAL(dp), DIMENSION(:,:), OPTIONAL, INTENT(in) ::   bounds_lon, bounds_lat, area
      LOGICAL , DIMENSION(:)  , OPTIONAL, INTENT(in) ::   mask
      !!----------------------------------------------------------------------
      !
      IF( xios_is_valid_domain     (cdid) ) THEN
         CALL xios_set_domain_attr     ( cdid, ni_glo=ni_glo, nj_glo=nj_glo, ibegin=ibegin, jbegin=jbegin, ni=ni, nj=nj,   &
            &    data_dim=data_dim, data_ibegin=data_ibegin, data_ni=data_ni, data_jbegin=data_jbegin, data_nj=data_nj ,   &
#if ! defined key_xios3
            &    ntiles=ntiles, tile_ibegin=tile_ibegin, tile_jbegin=tile_jbegin, tile_ni=tile_ni, tile_nj=tile_nj,        &
            &    tile_data_ibegin=tile_data_ibegin, tile_data_jbegin=tile_data_jbegin,                                     &
            &    tile_data_ni=tile_data_ni, tile_data_nj=tile_data_nj,                                                     &
#endif
            &    lonvalue_1D=lonvalue, latvalue_1D=latvalue, mask_1D=mask, nvertex=nvertex, bounds_lon_1D=bounds_lon,      &
            &    bounds_lat_1D=bounds_lat, area=area, type='curvilinear')
      ENDIF
      IF( xios_is_valid_domaingroup(cdid) ) THEN
         CALL xios_set_domaingroup_attr( cdid, ni_glo=ni_glo, nj_glo=nj_glo, ibegin=ibegin, jbegin=jbegin, ni=ni, nj=nj,   &
            &    data_dim=data_dim, data_ibegin=data_ibegin, data_ni=data_ni, data_jbegin=data_jbegin, data_nj=data_nj ,   &
#if ! defined key_xios3
            &    ntiles=ntiles, tile_ibegin=tile_ibegin, tile_jbegin=tile_jbegin, tile_ni=tile_ni, tile_nj=tile_nj,        &
            &    tile_data_ibegin=tile_data_ibegin, tile_data_jbegin=tile_data_jbegin,                                     &
            &    tile_data_ni=tile_data_ni, tile_data_nj=tile_data_nj,                                                     &
#endif
            &    lonvalue_1D=lonvalue, latvalue_1D=latvalue, mask_1D=mask, nvertex=nvertex, bounds_lon_1D=bounds_lon,      &
            &    bounds_lat_1D=bounds_lat, area=area, type='curvilinear' )
      ENDIF
      !
      CALL xios_solve_inheritance()
      !
   END SUBROUTINE iom_set_domain_attr


   SUBROUTINE iom_set_zoom_domain_attr( cdid, ibegin, jbegin, ni, nj )
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*), INTENT(in) ::   cdid
      INTEGER         , INTENT(in) ::   ibegin, jbegin, ni, nj
      !
      TYPE(xios_gridgroup) :: gridgroup_hdl
      TYPE(xios_grid)      :: grid_hdl
      TYPE(xios_domain)    :: domain_hdl
      TYPE(xios_axis)      :: axis_hdl
      CHARACTER(LEN=64)    :: cldomrefid   ! domain_ref name
      CHARACTER(len=1)     :: cl1          ! last character of this name
      !!----------------------------------------------------------------------
      !
#if ! defined key_xios3
      IF( xios_is_valid_zoom_domain(cdid) ) THEN
         ! define the zoom_domain attributs
         CALL xios_set_zoom_domain_attr( cdid, ibegin=ibegin, jbegin=jbegin, ni=ni, nj=nj )
#else
      IF( xios_is_valid_extract_domain(cdid) ) THEN
         ! define the zoom_domain attributs
         CALL xios_set_extract_domain_attr( cdid, ibegin=ibegin, jbegin=jbegin, ni=ni, nj=nj )
#endif
         ! define a new 2D grid with this new domain
         CALL xios_get_handle("grid_definition", gridgroup_hdl )
         CALL xios_add_child(gridgroup_hdl, grid_hdl, TRIM(cdid)//'_2D' )   ! add a new 2D grid to grid_definition
         CALL xios_add_child(grid_hdl, domain_hdl, TRIM(cdid) )             ! add its domain
         ! define a new 3D grid with this new domain
         CALL xios_add_child(gridgroup_hdl, grid_hdl, TRIM(cdid)//'_3D' )   ! add a new 3D grid to grid_definition
         CALL xios_add_child(grid_hdl, domain_hdl, TRIM(cdid) )             ! add its domain
         ! vertical axis
         cl1 = cdid(LEN_TRIM(cdid):)                                        ! last letter of cdid
         cl1 = CHAR(ICHAR(cl1)+32)                                          ! from upper to lower case
         CALL xios_add_child(grid_hdl, axis_hdl, 'depth'//cl1)              ! add its axis
      ENDIF
      !
   END SUBROUTINE iom_set_zoom_domain_attr


   SUBROUTINE iom_set_axis_attr( cdid, paxis, bounds )
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)      , INTENT(in) ::   cdid
      REAL(wp), DIMENSION(:)  , OPTIONAL, INTENT(in) ::   paxis
      REAL(wp), DIMENSION(:,:), OPTIONAL, INTENT(in) ::   bounds
      !!----------------------------------------------------------------------
      IF( PRESENT(paxis) ) THEN
         IF( xios_is_valid_axis     (cdid) )   CALL xios_set_axis_attr     ( cdid, n_glo=SIZE(paxis), value=real(paxis, dp) )
         IF( xios_is_valid_axisgroup(cdid) )   CALL xios_set_axisgroup_attr( cdid, n_glo=SIZE(paxis), value=real(paxis, dp) )
      ENDIF
      IF( PRESENT(bounds) ) THEN
         IF( xios_is_valid_axis     (cdid) )   CALL xios_set_axis_attr     ( cdid, bounds=real(bounds, dp) )
         IF( xios_is_valid_axisgroup(cdid) )   CALL xios_set_axisgroup_attr( cdid, bounds=real(bounds, dp) )
      ELSE
         IF( xios_is_valid_axis     (cdid) )   CALL xios_set_axis_attr     ( cdid)
         IF( xios_is_valid_axisgroup(cdid) )   CALL xios_set_axisgroup_attr( cdid)
      END IF
      CALL xios_solve_inheritance()
   END SUBROUTINE iom_set_axis_attr


   SUBROUTINE iom_set_field_attr( cdid, freq_op, freq_offset )
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)             , INTENT(in) ::   cdid
      TYPE(xios_duration), OPTIONAL, INTENT(in) ::   freq_op
      TYPE(xios_duration), OPTIONAL, INTENT(in) ::   freq_offset
      !!----------------------------------------------------------------------
      IF( xios_is_valid_field     (cdid) )   CALL xios_set_field_attr     ( cdid, freq_op=freq_op, freq_offset=freq_offset )
      IF( xios_is_valid_fieldgroup(cdid) )   CALL xios_set_fieldgroup_attr( cdid, freq_op=freq_op, freq_offset=freq_offset )
      CALL xios_solve_inheritance()
   END SUBROUTINE iom_set_field_attr


   SUBROUTINE iom_set_file_attr( cdid, name, name_suffix )
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)          , INTENT(in) ::   cdid
      CHARACTER(LEN=*),OPTIONAL , INTENT(in) ::   name, name_suffix
      !!----------------------------------------------------------------------
      IF( xios_is_valid_file     (cdid) )   CALL xios_set_file_attr     ( cdid, name=name, name_suffix=name_suffix )
      IF( xios_is_valid_filegroup(cdid) )   CALL xios_set_filegroup_attr( cdid, name=name, name_suffix=name_suffix )
      CALL xios_solve_inheritance()
   END SUBROUTINE iom_set_file_attr


   SUBROUTINE iom_get_file_attr( cdid, name, name_suffix, output_freq )
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)          , INTENT(in ) ::   cdid
      CHARACTER(LEN=*),OPTIONAL , INTENT(out) ::   name, name_suffix
      TYPE(xios_duration), OPTIONAL , INTENT(out) :: output_freq
      LOGICAL                                 ::   llexist1,llexist2,llexist3
      !---------------------------------------------------------------------
      IF( PRESENT( name        ) )   name = ''          ! default values
      IF( PRESENT( name_suffix ) )   name_suffix = ''
      IF( PRESENT( output_freq ) )   output_freq = xios_duration(0,0,0,0,0,0)
      IF( xios_is_valid_file     (cdid) ) THEN
         CALL xios_solve_inheritance()
         CALL xios_is_defined_file_attr     ( cdid, name = llexist1, name_suffix = llexist2, output_freq = llexist3)
         IF(llexist1)   CALL xios_get_file_attr     ( cdid, name = name )
         IF(llexist2)   CALL xios_get_file_attr     ( cdid, name_suffix = name_suffix )
         IF(llexist3)   CALL xios_get_file_attr     ( cdid, output_freq = output_freq )
      ENDIF
      IF( xios_is_valid_filegroup(cdid) ) THEN
         CALL xios_solve_inheritance()
         CALL xios_is_defined_filegroup_attr( cdid, name = llexist1, name_suffix = llexist2, output_freq = llexist3)
         IF(llexist1)   CALL xios_get_filegroup_attr( cdid, name = name )
         IF(llexist2)   CALL xios_get_filegroup_attr( cdid, name_suffix = name_suffix )
         IF(llexist3)   CALL xios_get_filegroup_attr( cdid, output_freq = output_freq )
      ENDIF
   END SUBROUTINE iom_get_file_attr


   SUBROUTINE iom_set_grid_attr( cdid, mask )
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)                   , INTENT(in) ::   cdid
      LOGICAL, DIMENSION(:,:,:), OPTIONAL, INTENT(in) ::   mask
      !!----------------------------------------------------------------------
      IF( xios_is_valid_grid     (cdid) )   CALL xios_set_grid_attr     ( cdid, mask_3D=mask )
      IF( xios_is_valid_gridgroup(cdid) )   CALL xios_set_gridgroup_attr( cdid, mask_3D=mask )
      CALL xios_solve_inheritance()
   END SUBROUTINE iom_set_grid_attr

   SUBROUTINE iom_setkt( kt, cdname )
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      INTEGER         , INTENT(in) ::   kt
      CHARACTER(LEN=*), INTENT(in) ::   cdname
      !!----------------------------------------------------------------------
      CALL iom_swap( cdname )   ! swap to cdname context
      CALL xios_update_calendar(kt)
      IF( cdname /= TRIM(cxios_context) )   CALL iom_swap( cxios_context )   ! return back to nemo context
   END SUBROUTINE iom_setkt

   SUBROUTINE iom_context_finalize( cdname )
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*), INTENT(in) :: cdname
      CHARACTER(LEN=120)           :: clname
      CHARACTER(LEN=256)           :: cltmpn    ! tempory name to store clname (in writting mode)
      INTEGER                      :: iln
      !!----------------------------------------------------------------------
      clname = TRIM(cdname)
      IF ( .NOT. Agrif_Root() ) THEN
         iln    = INDEX(clname,'/', BACK=.TRUE.)
         cltmpn = clname(1:iln)
         clname = clname(iln+1:LEN_TRIM(clname))
         clname = TRIM(cltmpn)//TRIM(Agrif_CFixed())//'_'//TRIM(clname)
      ENDIF

      IF( xios_is_valid_context(clname) ) THEN
         CALL iom_swap( cdname )   ! swap to cdname context
         CALL xios_context_finalize() ! finalize the context
         IF( cdname /= cxios_context ) CALL iom_swap( cxios_context )   ! return back to nemo context
      ENDIF
      !
   END SUBROUTINE iom_context_finalize


   SUBROUTINE set_grid( cdgrd, plon, plat, ldxios, ldrxios )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE set_grid  ***
      !!
      !! ** Purpose :   define horizontal grids
      !!----------------------------------------------------------------------
      CHARACTER(LEN=1)            , INTENT(in) ::   cdgrd
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   plon
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   plat
      !
      REAL(wp), DIMENSION(A2D(0),jpk) ::   zmask
      INTEGER :: jn
      INTEGER, DIMENSION(nijtile) :: ini, inj, idb
      LOGICAL, INTENT(IN) :: ldxios, ldrxios
      !!----------------------------------------------------------------------
      !
      ! nn_hls halo points
      CALL iom_set_domain_attr("grid_"//cdgrd, ni_glo=Ni0glo,nj_glo=Nj0glo,ibegin=mig0(Nis0)-1,jbegin=mjg0(Njs0)-1,ni=Ni_0,nj=Nj_0)
      CALL iom_set_domain_attr("grid_"//cdgrd, data_dim=2, data_ibegin = -nn_hls, data_ni=jpi, data_jbegin = -nn_hls, data_nj=jpj)

      ! Inner domain only
      CALL iom_set_domain_attr("grid_"//cdgrd//"_inner", ni_glo = Ni0glo, nj_glo = Nj0glo,   &
         &                     ibegin = mig0(Nis0) - 1, jbegin = mjg0(Njs0) - 1, ni = Ni_0, nj = Nj_0)
      CALL iom_set_domain_attr("grid_"//cdgrd//"_inner", data_dim=2, data_ibegin = 0, data_ni=Ni_0, data_jbegin = 0, data_nj=Nj_0)

      IF( ln_tile ) THEN
         DO jn = 1, nijtile
            ini(jn) = ntei_a(jn) - ntsi_a(jn) + 1     ! Tile size in i and j
            inj(jn) = ntej_a(jn) - ntsj_a(jn) + 1
            idb(jn) = -nn_hls                         ! Tile data offset (halo size)
         END DO

         ! Data includes all halo points
         CALL iom_set_domain_attr("grid_"//cdgrd, ntiles=nijtile,                                     &
            & tile_ibegin=ntsi_a(1:nijtile) - nn_hls - 1, tile_jbegin=ntsj_a(1:nijtile) - nn_hls - 1, &
            & tile_ni=ini(:), tile_nj=inj(:),                                                         &
            & tile_data_ibegin=idb(:), tile_data_jbegin=idb(:),                                       &
            & tile_data_ni=ini(:) - 2 * idb(:), tile_data_nj=inj(:) - 2 * idb(:))
         ! Data contains no halo points
         idb(:) = 0
         CALL iom_set_domain_attr("grid_"//cdgrd//"_inner", ntiles=nijtile,                           &
            & tile_ibegin=ntsi_a(1:nijtile) - nn_hls - 1, tile_jbegin=ntsj_a(1:nijtile) - nn_hls - 1, &
            & tile_ni=ini(:), tile_nj=inj(:),                                                         &
            & tile_data_ibegin=idb(:), tile_data_jbegin=idb(:),                                       &
            & tile_data_ni=ini(:) - 2 * idb(:), tile_data_nj=inj(:) - 2 * idb(:))
      ENDIF

!don't define lon and lat for restart reading context.
      IF ( .NOT.ldrxios ) &
         CALL iom_set_domain_attr("grid_"//cdgrd, lonvalue = real(RESHAPE(plon(Nis0:Nie0, Njs0:Nje0),(/ Ni_0*Nj_0 /)),dp),   &
         &                                        latvalue = real(RESHAPE(plat(Nis0:Nie0, Njs0:Nje0),(/ Ni_0*Nj_0 /)),dp ))
      !
      IF ( ln_mskland .AND. (.NOT.ldxios) ) THEN
         ! mask land points, keep values on coast line -> specific mask for U, V and W points
         SELECT CASE ( cdgrd )
         CASE('T')   ;   zmask(:,:,:) = tmask(Nis0  :Nie0  , Njs0:Nje0,:)
         CASE('U')   ;   zmask(:,:,:) = tmask(Nis0  :Nie0  , Njs0:Nje0,:) + tmask(Nis0+1:Nie0+1, Njs0  :Nje0  ,:)
         CASE('V')   ;   zmask(:,:,:) = tmask(Nis0  :Nie0  , Njs0:Nje0,:) + tmask(Nis0  :Nie0  , Njs0+1:Nje0+1,:)
         CASE('F')   ;   zmask(:,:,:) = tmask(Nis0  :Nie0  , Njs0:Nje0,:) + tmask(Nis0  :Nie0  , Njs0+1:Nje0+1,:)   &
            &                         + tmask(Nis0+1:Nie0+1, Njs0:Nje0,:) + tmask(Nis0+1:Nie0+1, Njs0+1:Nje0+1,:)
         CASE('W')   ;   zmask(:,:,2:jpk) = tmask(Nis0:Nie0, Njs0:Nje0,1:jpkm1) + tmask(Nis0:Nie0, Njs0:Nje0,2:jpk)
                         zmask(:,:,1    ) = tmask(Nis0:Nie0, Njs0:Nje0,1)
         END SELECT
         !
         CALL iom_set_domain_attr( "grid_"//cdgrd             , mask=RESHAPE(zmask(:,:,1),(/Ni_0*Nj_0    /)) /= 0. )
         CALL iom_set_grid_attr  ( "grid_"//cdgrd//"_3D"      , mask=RESHAPE(zmask(:,:,:),(/Ni_0,Nj_0,jpk/)) /= 0. )
         CALL iom_set_domain_attr( "grid_"//cdgrd//"_inner"   , mask=RESHAPE(zmask(:,:,1),(/Ni_0*Nj_0    /)) /= 0. )
         CALL iom_set_grid_attr  ( "grid_"//cdgrd//"_3D_inner", mask=RESHAPE(zmask(:,:,:),(/Ni_0,Nj_0,jpk/)) /= 0. )
      ENDIF
      !
   END SUBROUTINE set_grid

   SUBROUTINE set_grid_bounds( cdgrd, plon_cnr, plat_cnr, plon_pnt, plat_pnt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE set_grid_bounds  ***
      !!
      !! ** Purpose :   define horizontal grid corners
      !!
      !!----------------------------------------------------------------------
      CHARACTER(LEN=1)                      , INTENT(in) :: cdgrd
      REAL(wp), DIMENSION(jpi,jpj)          , INTENT(in) :: plon_cnr, plat_cnr  ! Lat/lon coord. of a contiguous vertex of cell (i,j)
      REAL(wp), DIMENSION(jpi,jpj), OPTIONAL, INTENT(in) :: plon_pnt, plat_pnt  ! Lat/lon coord. of the point of cell (i,j)
      !
      INTEGER :: ji, jj, jn
      INTEGER :: icnr, jcnr                             ! Offset such that the vertex coordinate (i+icnr,j+jcnr)
      !                                                 ! represents the
      !                                                 bottom-left corner of
      !                                                 cell (i,j)
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) :: z_bnds      ! Lat/lon coordinates of the vertices of cell (i,j)
      REAL(wp), ALLOCATABLE, DIMENSION(:,:)     :: z_fld       ! Working array to determine where to rotate cells
      REAL(wp), ALLOCATABLE, DIMENSION(:,:)     :: z_rot       ! Lat/lon working array for rotation of cells
      !!----------------------------------------------------------------------
      !
      ALLOCATE( z_bnds(4,jpi,jpj,2), z_fld(jpi,jpj), z_rot(4,2)  )
      !
      ! Offset of coordinate representing bottom-left corner
      SELECT CASE ( TRIM(cdgrd) )
      CASE ('T', 'W')   ;   icnr = -1   ;   jcnr = -1
      CASE ('U')        ;   icnr =  0   ;   jcnr = -1
      CASE ('V')        ;   icnr = -1   ;   jcnr =  0
      CASE ('F')        ;   icnr =  0   ;   jcnr =  0
      END SELECT
      !
      z_fld(:,:) = 1._wp
      CALL lbc_lnk( 'iom', z_fld, cdgrd, -1.0_wp )    ! Working array for location of northfold
      !
      ! Cell vertices that can be defined
      DO_2D( 0, 0, 0, 0 )
         z_bnds(1,ji,jj,1) = plat_cnr(ji+icnr,  jj+jcnr  ) ! Bottom-left
         z_bnds(2,ji,jj,1) = plat_cnr(ji+icnr+1,jj+jcnr  ) ! Bottom-right
         z_bnds(3,ji,jj,1) = plat_cnr(ji+icnr+1,jj+jcnr+1) ! Top-right
         z_bnds(4,ji,jj,1) = plat_cnr(ji+icnr,  jj+jcnr+1) ! Top-left
         z_bnds(1,ji,jj,2) = plon_cnr(ji+icnr,  jj+jcnr  ) ! Bottom-left
         z_bnds(2,ji,jj,2) = plon_cnr(ji+icnr+1,jj+jcnr  ) ! Bottom-right
         z_bnds(3,ji,jj,2) = plon_cnr(ji+icnr+1,jj+jcnr+1) ! Top-right
         z_bnds(4,ji,jj,2) = plon_cnr(ji+icnr,  jj+jcnr+1) ! Top-left
      END_2D
      !
      DO_2D( 0, 0, 0, 0 )
         IF( z_fld(ji,jj) == -1. ) THEN
            z_rot(1,:) = z_bnds(3,ji,jj,:) ; z_rot(2,:) = z_bnds(4,ji,jj,:)
            z_rot(3,:) = z_bnds(1,ji,jj,:) ; z_rot(4,:) = z_bnds(2,ji,jj,:)
            z_bnds(:,ji,jj,:) = z_rot(:,:)
         ENDIF
      END_2D
      !
      CALL iom_set_domain_attr("grid_"//cdgrd, bounds_lat = real(RESHAPE(z_bnds(:,Nis0:Nie0,Njs0:Nje0,1),(/ 4,Ni_0*Nj_0 /)), dp),           &
          &                                    bounds_lon = real(RESHAPE(z_bnds(:,Nis0:Nie0,Njs0:Nje0,2),(/ 4,Ni_0*Nj_0 /)), dp), nvertex=4 )
      !
      DEALLOCATE( z_bnds, z_fld, z_rot )
      !
   END SUBROUTINE set_grid_bounds

   SUBROUTINE set_grid_znl( plat )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE set_grid_znl  ***
      !!
      !! ** Purpose :   define grids for zonal mean
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   plat
      !
      INTEGER  :: ix, iy
      REAL(wp), DIMENSION(:), ALLOCATABLE  ::   zlon
      !!----------------------------------------------------------------------
      !
      ALLOCATE( zlon(Ni_0*Nj_0) )       ;       zlon(:) = 0._wp
      !
!      CALL dom_ngb( -168.53_wp, 65.03_wp, ix, iy, 'T' ) !  i-line that passes through Bering Strait: Reference latitude (used in plots)
      CALL dom_ngb( 180.0_wp, 90.0_wp, ix, iy, 'T' ) !  i-line that passes near the North Pole : Reference latitude (used in plots)
      CALL iom_set_domain_attr("gznl", ni_glo=Ni0glo, nj_glo=Nj0glo, ibegin=mig0(Nis0)-1, jbegin=mjg0(Njs0)-1, ni=Ni_0, nj=Nj_0)
      CALL iom_set_domain_attr("gznl", data_dim=2, data_ibegin = -nn_hls, data_ni = jpi, data_jbegin = -nn_hls, data_nj = jpj)
      CALL iom_set_domain_attr("gznl", lonvalue = real(zlon, dp),   &
         &                             latvalue = real(RESHAPE(plat(Nis0:Nie0, Njs0:Nje0),(/ Ni_0*Nj_0 /)),dp))
      CALL iom_set_zoom_domain_attr("ptr", ibegin=ix-1, jbegin=0, ni=1, nj=Nj0glo)
      !
      CALL iom_update_file_name('ptr')
      !
   END SUBROUTINE set_grid_znl


   SUBROUTINE set_scalar
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE set_scalar  ***
      !!
      !! ** Purpose :   define fake grids for scalar point
      !!
      !!----------------------------------------------------------------------
      REAL(dp), DIMENSION(1)   ::   zz = 1.
      !!----------------------------------------------------------------------
      !
      CALL iom_set_domain_attr('scalarpoint', ni_glo=jpnij, nj_glo=1, ibegin=narea-1, jbegin=0, ni=1, nj=1)
      CALL iom_set_domain_attr('scalarpoint', data_dim=2, data_ibegin = 1, data_ni = 1, data_jbegin = 1, data_nj = 1)
      !
      zz = REAL( narea, wp )
      CALL iom_set_domain_attr('scalarpoint', lonvalue=zz, latvalue=zz)
      !
   END SUBROUTINE set_scalar


   SUBROUTINE set_xmlatt
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE set_xmlatt  ***
      !!
      !! ** Purpose :   automatic definitions of some of the xml attributs...
      !!
      !!----------------------------------------------------------------------
      CHARACTER(len=1),DIMENSION( 3) ::   clgrd                    ! suffix name
      CHARACTER(len=256)             ::   clsuff                   ! suffix name
      CHARACTER(len=1)               ::   cl1                      ! 1 character
      CHARACTER(len=2)               ::   cl2                      ! 2 characters
      CHARACTER(len=3)               ::   cl3                      ! 3 characters
      INTEGER                        ::   ji, jg                   ! loop counters
      INTEGER                        ::   ix, iy                   ! i-,j- index
      REAL(wp)        ,DIMENSION(11) ::   zlontao                  ! longitudes of tao    moorings
      REAL(wp)        ,DIMENSION( 7) ::   zlattao                  ! latitudes  of tao    moorings
      REAL(wp)        ,DIMENSION( 4) ::   zlonrama                 ! longitudes of rama   moorings
      REAL(wp)        ,DIMENSION(11) ::   zlatrama                 ! latitudes  of rama   moorings
      REAL(wp)        ,DIMENSION( 3) ::   zlonpira                 ! longitudes of pirata moorings
      REAL(wp)        ,DIMENSION( 9) ::   zlatpira                 ! latitudes  of pirata moorings
      TYPE(xios_duration)            ::   f_op, f_of
      !!----------------------------------------------------------------------
      !
      ! frequency of the call of iom_put (attribut: freq_op)
      f_op%timestep = 1        ;  f_of%timestep =  0  ; CALL iom_set_field_attr('field_definition', freq_op=f_op, freq_offset=f_of)
      f_op%timestep = 2        ;  f_of%timestep =  0  ; CALL iom_set_field_attr('trendT_even'     , freq_op=f_op, freq_offset=f_of)
      f_op%timestep = 2        ;  f_of%timestep = -1  ; CALL iom_set_field_attr('trendT_odd'      , freq_op=f_op, freq_offset=f_of)
      f_op%timestep = nn_fsbc  ;  f_of%timestep =  0  ; CALL iom_set_field_attr('SBC'             , freq_op=f_op, freq_offset=f_of)
      f_op%timestep = nn_fsbc  ;  f_of%timestep =  0  ; CALL iom_set_field_attr('SBC_scalar'      , freq_op=f_op, freq_offset=f_of)
      f_op%timestep = nn_fsbc  ;  f_of%timestep =  0  ; CALL iom_set_field_attr('ABL'             , freq_op=f_op, freq_offset=f_of)
      IF( ln_diadct  ) THEN
        f_op%timestep = nn_dctwri ;  f_of%timestep =  0  ; CALL iom_set_field_attr('DCT'          , freq_op=f_op, freq_offset=f_of)
      ENDIF

      ! output file names (attribut: name)
      DO ji = 1, 9
         WRITE(cl1,'(i1)') ji
         CALL iom_update_file_name('file'//cl1)
      END DO
      DO ji = 1, 99
         WRITE(cl2,'(i2.2)') ji
         CALL iom_update_file_name('file'//cl2)
      END DO
      DO ji = 1, 999
         WRITE(cl3,'(i3.3)') ji
         CALL iom_update_file_name('file'//cl3)
      END DO

      ! Zooms...
      clgrd = (/ 'T', 'U', 'W' /)
      DO jg = 1, SIZE(clgrd)                                                                   ! grid type
         cl1 = clgrd(jg)
         ! Equatorial section (attributs: jbegin, ni, name_suffix)
         CALL dom_ngb( 0.0_wp, 0.0_wp, ix, iy, cl1 )
         CALL iom_set_zoom_domain_attr('Eq'//cl1, ibegin=0, jbegin=iy-1, ni=Ni0glo, nj=1 )
         CALL iom_get_file_attr   ('Eq'//cl1, name_suffix = clsuff             )
         CALL iom_set_file_attr   ('Eq'//cl1, name_suffix = TRIM(clsuff)//'_Eq')
         CALL iom_update_file_name('Eq'//cl1)
      END DO
      ! TAO moorings (attributs: ibegin, jbegin, name_suffix)
      zlontao = (/ 137.0, 147.0, 156.0, 165.0, -180.0, -170.0, -155.0, -140.0, -125.0, -110.0, -95.0 /)
      zlattao = (/  -8.0,  -5.0,  -2.0,   0.0,    2.0,    5.0,    8.0 /)
      CALL set_mooring( zlontao, zlattao )
      ! RAMA moorings (attributs: ibegin, jbegin, name_suffix)
      zlonrama = (/  55.0,  67.0, 80.5, 90.0 /)
      zlatrama = (/ -16.0, -12.0, -8.0, -4.0, -1.5, 0.0, 1.5, 4.0, 8.0, 12.0, 15.0 /)
      CALL set_mooring( zlonrama, zlatrama )
      ! PIRATA moorings (attributs: ibegin, jbegin, name_suffix)
      zlonpira = (/ -38.0, -23.0, -10.0 /)
      zlatpira = (/ -19.0, -14.0,  -8.0, 0.0, 4.0, 8.0, 12.0, 15.0, 20.0 /)
      CALL set_mooring( zlonpira, zlatpira )
      !
   END SUBROUTINE set_xmlatt


   SUBROUTINE set_mooring( plon, plat )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE set_mooring  ***
      !!
      !! ** Purpose :   automatic definitions of moorings xml attributs...
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:), INTENT(in) ::   plon, plat   ! longitudes/latitudes oft the mooring
      !
!!$      CHARACTER(len=1),DIMENSION(4) ::   clgrd = (/ 'T', 'U', 'V', 'W' /)   ! suffix name
      CHARACTER(len=1),DIMENSION(1) ::   clgrd = (/ 'T' /)        ! suffix name
      CHARACTER(len=256)            ::   clname                   ! file name
      CHARACTER(len=256)            ::   clsuff                   ! suffix name
      CHARACTER(len=1)              ::   cl1                      ! 1 character
      CHARACTER(len=6)              ::   clon,clat                ! name of longitude, latitude
      INTEGER                       ::   ji, jj, jg               ! loop counters
      INTEGER                       ::   ix, iy                   ! i-,j- index
      REAL(wp)                      ::   zlon, zlat
      !!----------------------------------------------------------------------
      DO jg = 1, SIZE(clgrd)
         cl1 = clgrd(jg)
         DO ji = 1, SIZE(plon)
            DO jj = 1, SIZE(plat)
               zlon = plon(ji)
               zlat = plat(jj)
               ! modifications for RAMA moorings
               IF( zlon ==  67. .AND. zlat ==  15. )   zlon =  65.
               IF( zlon ==  90. .AND. zlat <=  -4. )   zlon =  95.
               IF( zlon ==  95. .AND. zlat ==  -4. )   zlat =  -5.
               ! modifications for PIRATA moorings
               IF( zlon == -38. .AND. zlat == -19. )   zlon = -34.
               IF( zlon == -38. .AND. zlat == -14. )   zlon = -32.
               IF( zlon == -38. .AND. zlat ==  -8. )   zlon = -30.
               IF( zlon == -38. .AND. zlat ==   0. )   zlon = -35.
               IF( zlon == -23. .AND. zlat ==  20. )   zlat =  21.
               IF( zlon == -10. .AND. zlat == -14. )   zlat = -10.
               IF( zlon == -10. .AND. zlat ==  -8. )   zlat =  -6.
               IF( zlon == -10. .AND. zlat ==   4. ) THEN   ;   zlon = 0.   ;   zlat = 0.   ;   ENDIF
               CALL dom_ngb( zlon, zlat, ix, iy, cl1 )
               IF( zlon >= 0. ) THEN
                  IF( zlon == REAL(NINT(zlon), wp) ) THEN   ;   WRITE(clon, '(i3,  a)') NINT( zlon), 'e'
                  ELSE                                      ;   WRITE(clon, '(f5.1,a)')       zlon , 'e'
                  ENDIF
               ELSE
                  IF( zlon == REAL(NINT(zlon), wp) ) THEN   ;   WRITE(clon, '(i3,  a)') NINT(-zlon), 'w'
                  ELSE                                      ;   WRITE(clon, '(f5.1,a)')      -zlon , 'w'
                  ENDIF
               ENDIF
               IF( zlat >= 0. ) THEN
                  IF( zlat == REAL(NINT(zlat), wp) ) THEN   ;   WRITE(clat, '(i2,  a)') NINT( zlat), 'n'
                  ELSE                                      ;   WRITE(clat, '(f4.1,a)')       zlat , 'n'
                  ENDIF
               ELSE
                  IF( zlat == REAL(NINT(zlat), wp) ) THEN   ;   WRITE(clat, '(i2,  a)') NINT(-zlat), 's'
                  ELSE                                      ;   WRITE(clat, '(f4.1,a)')      -zlat , 's'
                  ENDIF
               ENDIF
               clname = TRIM(ADJUSTL(clat))//TRIM(ADJUSTL(clon))
               CALL iom_set_zoom_domain_attr(TRIM(clname)//cl1, ibegin= ix-1, jbegin= iy-1, ni=1, nj=1)

               CALL iom_get_file_attr   (TRIM(clname)//cl1, name_suffix = clsuff                         )
               CALL iom_set_file_attr   (TRIM(clname)//cl1, name_suffix = TRIM(clsuff)//'_'//TRIM(clname))
               CALL iom_update_file_name(TRIM(clname)//cl1)
            END DO
         END DO
      END DO

   END SUBROUTINE set_mooring


   SUBROUTINE iom_update_file_name( cdid )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE iom_update_file_name  ***
      !!
      !! ** Purpose :
      !!
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)          , INTENT(in) ::   cdid
      !
      CHARACTER(LEN=256) ::   clname, cltmpn
      CHARACTER(LEN=20)  ::   clfreq
      CHARACTER(LEN=20)  ::   cldate
      INTEGER            ::   idx
      INTEGER            ::   jn, iln
      INTEGER            ::   itrlen
      INTEGER            ::   iyear, imonth, iday, isec
      REAL(dp)           ::   zsec
      LOGICAL            ::   llexist
      TYPE(xios_duration)   ::   output_freq
      !!----------------------------------------------------------------------
      !
      DO jn = 1, 2
         !
         output_freq = xios_duration(0,0,0,0,0,0)
         IF( jn == 1 )   CALL iom_get_file_attr( cdid, name        = clname, output_freq = output_freq )
         IF( jn == 2 )   CALL iom_get_file_attr( cdid, name_suffix = clname )
         !
         IF ( TRIM(clname) /= '' ) THEN
            !
            idx = INDEX(clname,'@expname@') + INDEX(clname,'@EXPNAME@')
            DO WHILE ( idx /= 0 )
               clname = clname(1:idx-1)//TRIM(cexper)//clname(idx+9:LEN_TRIM(clname))
               idx = INDEX(clname,'@expname@') + INDEX(clname,'@EXPNAME@')
            END DO
            !
            idx = INDEX(clname,'@freq@') + INDEX(clname,'@FREQ@')
            DO WHILE ( idx /= 0 )
              IF ( output_freq%timestep /= 0) THEN
                  WRITE(clfreq,'(I18,A2)')INT(output_freq%timestep),'ts'
                  itrlen = LEN_TRIM(ADJUSTL(clfreq))
              ELSE IF ( output_freq%second /= 0 ) THEN
                  WRITE(clfreq,'(I19,A1)')INT(output_freq%second),'s'
                  itrlen = LEN_TRIM(ADJUSTL(clfreq))
              ELSE IF ( output_freq%minute /= 0 ) THEN
                  WRITE(clfreq,'(I18,A2)')INT(output_freq%minute),'mi'
                  itrlen = LEN_TRIM(ADJUSTL(clfreq))
              ELSE IF ( output_freq%hour /= 0 ) THEN
                  WRITE(clfreq,'(I19,A1)')INT(output_freq%hour),'h'
                  itrlen = LEN_TRIM(ADJUSTL(clfreq))
              ELSE IF ( output_freq%day /= 0 ) THEN
                  WRITE(clfreq,'(I19,A1)')INT(output_freq%day),'d'
                  itrlen = LEN_TRIM(ADJUSTL(clfreq))
              ELSE IF ( output_freq%month /= 0 ) THEN
                  WRITE(clfreq,'(I19,A1)')INT(output_freq%month),'m'
                  itrlen = LEN_TRIM(ADJUSTL(clfreq))
              ELSE IF ( output_freq%year /= 0 ) THEN
                  WRITE(clfreq,'(I19,A1)')INT(output_freq%year),'y'
                  itrlen = LEN_TRIM(ADJUSTL(clfreq))
              ELSE
                  CALL ctl_stop('error in the name of file id '//TRIM(cdid),   &
                     & ' attribute output_freq is undefined -> cannot replace @freq@ in '//TRIM(clname) )
              ENDIF
              clname = clname(1:idx-1)//TRIM(ADJUSTL(clfreq))//clname(idx+6:LEN_TRIM(clname))
              idx = INDEX(clname,'@freq@') + INDEX(clname,'@FREQ@')
            END DO
            !
            idx = INDEX(clname,'@startdate@') + INDEX(clname,'@STARTDATE@')
            DO WHILE ( idx /= 0 )
               cldate = iom_sdate( fjulday - rn_Dt / rday )
               clname = clname(1:idx-1)//TRIM(cldate)//clname(idx+11:LEN_TRIM(clname))
               idx = INDEX(clname,'@startdate@') + INDEX(clname,'@STARTDATE@')
            END DO
            !
            idx = INDEX(clname,'@startdatefull@') + INDEX(clname,'@STARTDATEFULL@')
            DO WHILE ( idx /= 0 )
               cldate = iom_sdate( fjulday - rn_Dt / rday, ldfull = .TRUE. )
               clname = clname(1:idx-1)//TRIM(cldate)//clname(idx+15:LEN_TRIM(clname))
               idx = INDEX(clname,'@startdatefull@') + INDEX(clname,'@STARTDATEFULL@')
            END DO
            !
            idx = INDEX(clname,'@enddate@') + INDEX(clname,'@ENDDATE@')
            DO WHILE ( idx /= 0 )
               cldate = iom_sdate( fjulday + rn_Dt / rday * REAL( nitend - nit000, wp ), ld24 = .TRUE. )
               clname = clname(1:idx-1)//TRIM(cldate)//clname(idx+9:LEN_TRIM(clname))
               idx = INDEX(clname,'@enddate@') + INDEX(clname,'@ENDDATE@')
            END DO
            !
            idx = INDEX(clname,'@enddatefull@') + INDEX(clname,'@ENDDATEFULL@')
            DO WHILE ( idx /= 0 )
               cldate = iom_sdate( fjulday + rn_Dt / rday * REAL( nitend - nit000, wp ), ld24 = .TRUE., ldfull = .TRUE. )
               clname = clname(1:idx-1)//TRIM(cldate)//clname(idx+13:LEN_TRIM(clname))
               idx = INDEX(clname,'@enddatefull@') + INDEX(clname,'@ENDDATEFULL@')
            END DO
            !
            IF( (jn ==1).AND.(.NOT. Agrif_Root())) THEN
               iln    = INDEX(clname,'/', BACK=.TRUE.)
               cltmpn = clname(1:iln)
               clname = clname(iln+1:LEN_TRIM(clname))
               clname = TRIM(cltmpn)//TRIM(Agrif_CFixed())//'_'//TRIM(clname)
            ENDIF
            IF( jn == 1 )   CALL iom_set_file_attr( cdid, name        = clname )
            IF( jn == 2 )   CALL iom_set_file_attr( cdid, name_suffix = clname )
            !
         ENDIF
         !
      END DO
      !
   END SUBROUTINE iom_update_file_name


   FUNCTION iom_sdate( pjday, ld24, ldfull )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE iom_sdate  ***
      !!
      !! ** Purpose :   send back the date corresponding to the given julian day
      !!----------------------------------------------------------------------
      REAL(dp), INTENT(in   )           ::   pjday    ! julian day
      LOGICAL , INTENT(in   ), OPTIONAL ::   ld24     ! true to force 24:00 instead of 00:00
      LOGICAL , INTENT(in   ), OPTIONAL ::   ldfull   ! true to get the compleate date: yyyymmdd_hh:mm:ss
      !
      CHARACTER(LEN=20) ::   iom_sdate
      CHARACTER(LEN=50) ::   clfmt                         !  format used to write the date
      INTEGER           ::   iyear, imonth, iday, ihour, iminute, isec
      REAL(dp)          ::   zsec
      LOGICAL           ::   ll24, llfull
      !!----------------------------------------------------------------------
      !
      IF( PRESENT(ld24) ) THEN   ;   ll24 = ld24
      ELSE                       ;   ll24 = .FALSE.
      ENDIF
      !
      IF( PRESENT(ldfull) ) THEN   ;   llfull = ldfull
      ELSE                         ;   llfull = .FALSE.
      ENDIF
      !
      CALL ju2ymds( pjday, iyear, imonth, iday, zsec )
      isec = NINT(zsec)
      !
      IF ( ll24 .AND. isec == 0 ) THEN   ! 00:00 of the next day -> move to 24:00 of the current day
         CALL ju2ymds( pjday - 1.0_wp, iyear, imonth, iday, zsec )
         isec = 86400
      ENDIF
      !
      IF( iyear < 10000 ) THEN   ;   clfmt = "i4.4,2i2.2"                ! format used to write the date
      ELSE                       ;   WRITE(clfmt, "('i',i1,',2i2.2')") INT(LOG10(REAL(iyear,wp))) + 1
      ENDIF
      !
!$AGRIF_DO_NOT_TREAT
      ! needed in the conv
      IF( llfull ) THEN
         clfmt = TRIM(clfmt)//",'_',i2.2,':',i2.2,':',i2.2"
         ihour   = isec / 3600
         isec    = MOD(isec, 3600)
         iminute = isec / 60
         isec    = MOD(isec, 60)
         WRITE(iom_sdate, '('//TRIM(clfmt)//')') iyear, imonth, iday, ihour, iminute, isec    ! date of the end of run
      ELSE
         WRITE(iom_sdate, '('//TRIM(clfmt)//')') iyear, imonth, iday                          ! date of the end of run
      ENDIF
!$AGRIF_END_DO_NOT_TREAT
      !
   END FUNCTION iom_sdate

#else
   !!----------------------------------------------------------------------
   !!   NOT 'key_xios'                               a few dummy routines
   !!----------------------------------------------------------------------
   SUBROUTINE iom_setkt( kt, cdname )
      INTEGER         , INTENT(in)::   kt
      CHARACTER(LEN=*), INTENT(in) ::   cdname
      IF( .FALSE. )   WRITE(numout,*) kt, cdname   ! useless test to avoid compilation warnings
   END SUBROUTINE iom_setkt

   SUBROUTINE iom_context_finalize( cdname )
      CHARACTER(LEN=*), INTENT(in) ::   cdname
      IF( .FALSE. )   WRITE(numout,*)  cdname   ! useless test to avoid compilation warnings
   END SUBROUTINE iom_context_finalize

   SUBROUTINE iom_update_file_name( cdid )
      CHARACTER(LEN=*), INTENT(in) ::   cdid
      IF( .FALSE. )   WRITE(numout,*)  cdid   ! useless test to avoid compilation warnings
   END SUBROUTINE iom_update_file_name

#endif

   LOGICAL FUNCTION iom_use( cdname )
      CHARACTER(LEN=*), INTENT(in) ::   cdname
#if defined key_xios
      iom_use = xios_field_is_active( cdname )
#else
      iom_use = .FALSE.
#endif
   END FUNCTION iom_use

   SUBROUTINE iom_miss_val( cdname, pmiss_val )
      CHARACTER(LEN=*), INTENT(in ) ::   cdname
      REAL(wp)        , INTENT(out) ::   pmiss_val
      REAL(dp)                      ::   ztmp_pmiss_val
#if defined key_xios
      ! get missing value
      CALL xios_get_field_attr( cdname, default_value = ztmp_pmiss_val )
      pmiss_val = ztmp_pmiss_val
#else
      IF( .FALSE. )   WRITE(numout,*) cdname, pmiss_val   ! useless test to avoid compilation warnings
      IF( .FALSE. )   pmiss_val = 0._wp                   ! useless assignment to avoid compilation warnings
#endif
   END SUBROUTINE iom_miss_val

   !!======================================================================
END MODULE iom
