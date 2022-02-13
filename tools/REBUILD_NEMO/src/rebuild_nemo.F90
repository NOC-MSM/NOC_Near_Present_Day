PROGRAM rebuild_nemo
#define key_netcdf4
   !!=========================================================================
   !!                        ***  rebuild_nemo  ***
   !!=========================================================================
   !!
   !!  A routine to rebuild NEMO files from multiple processors into one file.
   !!  This routine is designed to be much quicker than the old IOIPSL rebuild
   !!  but at the cost of an increased memory usage.
   !!
   !!  NEMO rebuild has the following features:
   !!     * dynamically works out what variables require rebuilding
   !!     * does not copy subdomain halo regions
   !!     * works for 1,2,3 and 4d arrays or types for all valid NetCDF types
   !!     * utilises OMP shared memory parallelisation where applicable
   !!     * time 'slicing' for lower memory use 
   !!       (only for 4D vars with unlimited dimension)
   !!
   !!  Ed Blockley - August 2011
   !!  (based on original code by Matt Martin)
   !!  Julien Palmieri and Andrew Coward - September 2018 (add compression and chunking)
   !!
   !!-------------------------------------------------------------------------
   !!
   !!  The code reads the filestem and number of subdomains from the namelist file nam_rebuild.
   !! 
   !!  The 1st subdomain file is used to determine the dimensions and variables in all the input files. 
   !!  It is also used to find which dimensions (and hence which variables) require rebuilding 
   !!  as well as information about the global domain.
   !!
   !!  It then opens all the input files (unbuffered) and creates an array of netcdf identifiers
   !!  before looping through all the variables and updating the rebuilt output file (either by direct 
   !!  copying or looping over the number of domains and rebuilding as appropriate).
   !!  
   !!  The code looks more complicated than it is because it has lots of case statements to deal with all 
   !!  the various NetCDF data types and with various data dimensions (up to 4d).
   !!
   !!  Diagnostic output is written to numout (default 6 - stdout)
   !!  and errors are written to numerr (default 0 - stderr).
   !!
   !!  If time slicing is specified the code will use less memory but take a little longer.
   !!  It does this by breaking down the 4D input variables over their 4th dimension 
   !!  (generally time) by way of a while loop.
   !!
   !!-------------------------------------------------------------------------------
  
   USE netcdf

!$ USE omp_lib           ! Note OpenMP sentinel

   IMPLICIT NONE

   ! kind specifications
   INTEGER,PARAMETER :: i1=SELECTED_INT_KIND(2)          ! NF90_BYTE 
   INTEGER,PARAMETER :: i2=SELECTED_INT_KIND(4)          ! NF90_SHORT
   INTEGER,PARAMETER :: i4=SELECTED_INT_KIND(9)          ! NF90_INT
   INTEGER,PARAMETER :: sp=SELECTED_REAL_KIND(6,37)      ! NF90_FLOAT
   INTEGER,PARAMETER :: dp=SELECTED_REAL_KIND(12,307)    ! NF90_DOUBLE

   INTEGER,PARAMETER :: numnam = 11
   INTEGER,PARAMETER :: numout = 6
   INTEGER,PARAMETER :: numerr = 0

   LOGICAL, PARAMETER :: l_verbose = .true. 
    
   CHARACTER(LEN=nf90_max_name) :: filebase, suffix, attname, dimname, varname, time, date, zone, timestamp
   CHARACTER(LEN=nf90_max_name), ALLOCATABLE :: filenames(:), indimnames(:)
   CHARACTER(LEN=nf90_max_name), DIMENSION(2) :: dims
   CHARACTER(LEN=256) :: cnampath, cdimlst, cdim
   CHARACTER(LEN=50)  :: clibnc ! netcdf library version

   INTEGER :: ndomain, ifile, ndomain_file, nslicesize, deflate_level
   INTEGER :: ncid, outid, idim, istop
   INTEGER :: natts, attid, xtype, varid, rbdims 
   INTEGER :: jv, ndims, nvars, dimlen, dimids(4)
   INTEGER :: dimid, unlimitedDimId, di, dj, dr
   INTEGER :: nmax_unlimited, nt, ntslice 
   INTEGER :: fchunksize = 32000000  ! NetCDF global file chunk cache size
   INTEGER :: patchchunk             ! NetCDF processor-domain file chunk cache size
   INTEGER :: nthreads = 1
   INTEGER :: chunkalg = 0          ! NetCDF4 variable chunking algorithm
                                    ! Default variable chunksizes (typical ORCA025
                                    ! recommendations which can be adjusted via namelist
                                    ! or will be bounded if too large for domain.)
   INTEGER :: nc4_xchunk = 206      ! Default x (longitude) variable chunk size
   INTEGER :: nc4_ychunk = 135      ! Default y (latitude)  variable chunk size
   INTEGER :: nc4_zchunk = 1        ! Default z (depth) variable chunk size (almost always 1)
   INTEGER :: nc4_tchunk = 1        ! Default t (time)  variable chunk size (almost always 1)
   INTEGER, ALLOCATABLE  :: outdimids(:), outdimlens(:), indimlens(:), inncids(:)
   INTEGER, ALLOCATABLE  :: chunksizes(:)
   INTEGER, ALLOCATABLE  :: global_sizes(:), rebuild_dims(:)
   INTEGER, DIMENSION(2) :: halo_start, halo_end, local_sizes
   INTEGER, DIMENSION(2) :: idomain, jdomain, rdomain, start_pos
   INTEGER :: ji, jj, jk, jl, jr
   INTEGER :: nargs                 ! number of arguments
   INTEGER, EXTERNAL :: iargc
 
   REAL(sp) :: ValMin, ValMax, InMin, InMax, rmdi
   REAL(dp), ALLOCATABLE :: mdiVals(:)

   ! NF90_BYTE local data arrays
   INTEGER(i1), ALLOCATABLE, SAVE, DIMENSION(:) :: localdata_1d_i1
   INTEGER(i1), ALLOCATABLE, SAVE, DIMENSION(:,:) :: localdata_2d_i1
   INTEGER(i1), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: localdata_3d_i1
   INTEGER(i1), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) :: localdata_4d_i1

   ! NF90_SHORT local data arrays
   INTEGER(i2), ALLOCATABLE, SAVE, DIMENSION(:) :: localdata_1d_i2
   INTEGER(i2), ALLOCATABLE, SAVE, DIMENSION(:,:) :: localdata_2d_i2
   INTEGER(i2), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: localdata_3d_i2
   INTEGER(i2), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) :: localdata_4d_i2

   ! NF90_INT local data arrays
   INTEGER(i4), ALLOCATABLE, SAVE, DIMENSION(:) :: localdata_1d_i4
   INTEGER(i4), ALLOCATABLE, SAVE, DIMENSION(:,:) :: localdata_2d_i4
   INTEGER(i4), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: localdata_3d_i4
   INTEGER(i4), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) :: localdata_4d_i4

   ! NF90_FLOAT local data arrays
   REAL(sp), ALLOCATABLE, SAVE, DIMENSION(:) :: localdata_1d_sp
   REAL(sp), ALLOCATABLE, SAVE, DIMENSION(:,:) :: localdata_2d_sp
   REAL(sp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: localdata_3d_sp
   REAL(sp), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) :: localdata_4d_sp

   ! NF90_DOUBLE local data arrays
   REAL(dp), ALLOCATABLE, SAVE, DIMENSION(:) :: localdata_1d_dp
   REAL(dp), ALLOCATABLE, SAVE, DIMENSION(:,:) :: localdata_2d_dp
   REAL(dp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: localdata_3d_dp
   REAL(dp), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) :: localdata_4d_dp

   ! NF90_BYTE global data arrays
   INTEGER(i1) :: globaldata_0d_i1
   INTEGER(i1), ALLOCATABLE, DIMENSION(:) :: globaldata_1d_i1
   INTEGER(i1), ALLOCATABLE, DIMENSION(:,:) :: globaldata_2d_i1
   INTEGER(i1), ALLOCATABLE, DIMENSION(:,:,:) :: globaldata_3d_i1
   INTEGER(i1), ALLOCATABLE, DIMENSION(:,:,:,:) :: globaldata_4d_i1

   ! NF90_SHORT global data arrays
   INTEGER(i2) :: globaldata_0d_i2
   INTEGER(i2), ALLOCATABLE, DIMENSION(:) :: globaldata_1d_i2
   INTEGER(i2), ALLOCATABLE, DIMENSION(:,:) :: globaldata_2d_i2
   INTEGER(i2), ALLOCATABLE, DIMENSION(:,:,:) :: globaldata_3d_i2
   INTEGER(i2), ALLOCATABLE, DIMENSION(:,:,:,:) :: globaldata_4d_i2

   ! NF90_INT global data arrays
   INTEGER(i4) :: globaldata_0d_i4
   INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: globaldata_1d_i4
   INTEGER(i4), ALLOCATABLE, DIMENSION(:,:) :: globaldata_2d_i4
   INTEGER(i4), ALLOCATABLE, DIMENSION(:,:,:) :: globaldata_3d_i4
   INTEGER(i4), ALLOCATABLE, DIMENSION(:,:,:,:) :: globaldata_4d_i4
 
   ! NF90_FLOAT global data arrays
   REAL(sp) :: globaldata_0d_sp
   REAL(sp), ALLOCATABLE, DIMENSION(:) :: globaldata_1d_sp
   REAL(sp), ALLOCATABLE, DIMENSION(:,:) :: globaldata_2d_sp
   REAL(sp), ALLOCATABLE, DIMENSION(:,:,:) :: globaldata_3d_sp
   REAL(sp), ALLOCATABLE, DIMENSION(:,:,:,:) :: globaldata_4d_sp

   ! NF90_DOUBLE global data arrays
   REAL(dp) :: globaldata_0d_dp
   REAL(dp), ALLOCATABLE, DIMENSION(:) :: globaldata_1d_dp
   REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: globaldata_2d_dp
   REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: globaldata_3d_dp
   REAL(dp), ALLOCATABLE, DIMENSION(:,:,:,:) :: globaldata_4d_dp

   LOGICAL :: l_valid     = .false.
   LOGICAL :: l_noRebuild = .false.
   LOGICAL :: l_findDims  = .true.
   LOGICAL :: l_maskout   = .false.
   LOGICAL :: l_namexist  = .false.

   NAMELIST/nam_rebuild/ filebase, ndomain, dims, nslicesize, l_maskout, deflate_level, &
                       & nc4_xchunk, nc4_ychunk, nc4_zchunk, nc4_tchunk, fchunksize         


   !End of definitions 

!--------------------------------------------------------------------------------
!0. OMP setup

!$OMP PARALLEL DEFAULT(NONE) SHARED(nthreads)
!$OMP MASTER
!$      nthreads = omp_get_num_threads()
!$      WRITE(numout,*) 'Running OMP with ',nthreads,' thread(s).'
!$OMP END MASTER
!$OMP END PARALLEL
 
!--------------------------------------------------------------------------------
!1.0 Check netcdf version for warning
   clibnc = TRIM(nf90_inq_libvers())
   IF (ICHAR(clibnc(1:1)) <= 3) THEN
      PRINT *, '=========================================================='
      PRINT *, 'You are using old netcdf library (',TRIM(clibnc),').'
      PRINT *, 'REBUILD_NEMO support of old netcdf library will end soon'
      PRINT *, 'please consider moving to netcdf 4 or higher'
      PRINT *, '=========================================================='
   END IF

!1.1 Get the namelist path
   !Determine the number of arguments on the command line
   nargs=COMMAND_ARGUMENT_COUNT()
   !Check that the required argument is present, if it is not then set it to the default value: nam_rebuild
   IF (nargs == 0) THEN
      WRITE(numout,*)
      WRITE(numout,*) 'W A R N I N G : Namelist path not supplied as command line argument. Using default, nam_rebuild.'
      cnampath='nam_rebuild'
   ELSE IF (nargs == 1) THEN
      CALL GET_COMMAND_ARGUMENT(1, cnampath)
   ELSE 
      WRITE(numerr,*) 'E R R O R ! : Incorrect number of command line arguments. Please supply only'
      WRITE(numerr,*) '         the path to the namelist file, or no arguments to use default value'
      STOP 1
   END IF

   ! check presence of namelist
   INQUIRE(FILE=cnampath, EXIST=l_namexist)
   IF (.NOT. l_namexist) THEN
      WRITE(numout,*)
      WRITE(numout,*) 'E R R O R : Namelist '//TRIM(cnampath)//' not present.'
      STOP 42
   END IF

!1.2 Read in the namelist 

   dims(:) = ""
   nslicesize = 0
   deflate_level = 0
   OPEN( UNIT=numnam, FILE=TRIM(cnampath), FORM='FORMATTED', STATUS='OLD' )
   READ( numnam, nam_rebuild )
   CLOSE( numnam )
   IF( .NOT. ALL(dims(:) == "") ) l_findDims = .false.

!1.3 Set up the filenames and fileids

   ALLOCATE(filenames(ndomain))
   IF (l_verbose) WRITE(numout,*) 'Rebuilding the following files:'
   DO ifile = 1, ndomain
      WRITE(suffix,'(i4.4)') ifile-1
      filenames(ifile) = TRIM(filebase)//'_'//TRIM(suffix)//'.nc'
      IF (l_verbose) WRITE(numout,*) TRIM(filenames(ifile))
   END DO
   ALLOCATE(inncids(ndomain))
  
!---------------------------------------------------------------------------
!2. Read in the global dimensions from the first input file and set up the output file
 
   CALL check_nf90( nf90_open( TRIM(filenames(1)), nf90_share, ncid ) )
   CALL check_nf90( nf90_inquire( ncid, ndims, nvars, natts ) )
    
!2.0 Read in the total number of processors the file is expecting and check it's correct  
  
   CALL check_nf90( nf90_get_att( ncid, nf90_global, 'DOMAIN_number_total', ndomain_file ) )
   IF( ndomain /= ndomain_file ) THEN
      WRITE(numerr,*) 'ERROR! : number of files to rebuild in file does not agree with namelist'
      WRITE(numerr,*) 'Attribute DOMAIN_number_total is : ', ndomain_file
      WRITE(numerr,*) 'Number of files specified in namelist is: ', ndomain
      STOP 2
   ENDIF
  
!2.1 Set up the output file
#if defined key_netcdf4
   CALL check_nf90( nf90_create( TRIM(filebase)//'.nc', nf90_netcdf4, outid, chunksize=fchunksize ) )
#else
   CALL check_nf90( nf90_create( TRIM(filebase)//'.nc', nf90_64bit_offset, outid, chunksize=fchunksize ) )
#endif

!2.2 Set up dimensions in output file

!2.2.0 Find out how many dimensions are required to be rebuilt and which ones they are
   CALL check_nf90( nf90_inquire_attribute( ncid, nf90_global, 'DOMAIN_dimensions_ids', xtype, rbdims, attid ) )

   ALLOCATE(rebuild_dims(rbdims))
   CALL check_nf90( nf90_get_att( ncid, nf90_global, 'DOMAIN_dimensions_ids', rebuild_dims ) )

   ALLOCATE(global_sizes(rbdims))
   CALL check_nf90( nf90_get_att( ncid, nf90_global, 'DOMAIN_size_global', global_sizes ) )
   IF (l_verbose) WRITE(numout,*) 'Size of global arrays: ', global_sizes


!2.2.1 Copy the dimensions into the output file apart from rebuild_dims() which are dimensioned globally 
   ALLOCATE(indimlens(ndims), indimnames(ndims), outdimlens(ndims))
   CALL check_nf90( nf90_inquire( ncid, unlimitedDimId = unlimitedDimId ) )
   istop = 0
   DO idim = 1, ndims
      CALL check_nf90( nf90_inquire_dimension( ncid, idim, dimname, dimlen ) )
      CALL check_nf90( nf90_get_att( ncid, nf90_global, 'DOMAIN_size_local', local_sizes ) )
      indimlens(idim) = dimlen    
      indimnames(idim) = dimname
      IF (l_findDims) THEN
         IF( idim == rebuild_dims(1) ) THEN
            IF( dimlen == local_sizes(1) ) THEN 
               dimlen = global_sizes(1)
               dims(1) = trim(dimname)
            ELSE
               istop = 1
            ENDIF
         ENDIF
         IF( rbdims > 1 .AND. idim == rebuild_dims(2) ) THEN
            IF( dimlen == local_sizes(2) ) THEN
               dimlen = global_sizes(2)
               dims(2) = trim(dimname)
            ELSE
               istop = 1
            ENDIF
         ENDIF
      ELSE ! l_findDims = false
         IF( TRIM(dimname) == TRIM(dims(1))) THEN
            dimlen = global_sizes(1)
            rebuild_dims(1) = idim
         ENDIF
         IF( rbdims > 1 .AND. TRIM(dimname) == TRIM(dims(2))) THEN
            dimlen = global_sizes(2)
            rebuild_dims(2) = idim
         ENDIF
      ENDIF

      IF( idim == unlimitedDimId ) THEN
         CALL check_nf90( nf90_def_dim( outid, dimname, nf90_unlimited, dimid) )
         nmax_unlimited = dimlen
      ELSE
         CALL check_nf90( nf90_def_dim( outid, dimname, dimlen, dimid) )
      ENDIF
      outdimlens(idim) = dimlen
   END DO
   ! nmax_unlimited is only used for time-slicing so we set it to be at least 1 to 
   ! account for files with no record dimension or zero length record dimension(!)
   nmax_unlimited = max(nmax_unlimited,1)

   IF( istop == 1 ) THEN
      WRITE(numerr,*) 'ERROR! : DOMAIN_local_sizes attribute does not match rebuild dimension lengths in the first file'
      WRITE(numerr,*) 'Attribute DOMAIN_local_sizes is : ', local_sizes
      WRITE(numerr,*) 'Dimensions to be rebuilt are of size : ', outdimlens(rebuild_dims(1)), outdimlens(rebuild_dims(2)) 
      STOP 3
   ENDIF

   IF (l_findDims) THEN
      IF (l_verbose) WRITE(numout,*) 'Finding rebuild dimensions from the first file...'
   ELSE
      IF (l_verbose) WRITE(numout,*) 'Using rebuild dimensions given in namelist...'
   ENDIF

   IF( rbdims > 1 ) THEN
      IF (l_verbose) WRITE(numout,*) 'Rebuilding across dimensions '//TRIM(indimnames(rebuild_dims(1)))//  &
         &                      ' and '//TRIM(indimnames(rebuild_dims(2)))
   ELSE
      IF (l_verbose) WRITE(numout,*) 'Rebuilding across dimension '//TRIM(indimnames(rebuild_dims(1)))
   ENDIF
      
!2.2.2 Copy the global attributes into the output file, apart from those beginning with DOMAIN_  
!      Also need to change the file_name attribute and the TimeStamp attribute.
   DO attid = 1, natts
      CALL check_nf90( nf90_inq_attname( ncid, nf90_global, attid, attname ) )
      IF( INDEX( attname, "DOMAIN_" ) == 1 ) CYCLE
      IF( INDEX( attname, "file_name") == 1 ) CYCLE
      IF( INDEX( attname, "associate_file") == 1 ) CYCLE
      IF (l_verbose) WRITE(numout,*) 'Copying attribute '//TRIM(attname)//' into destination file...'
      CALL check_nf90( nf90_copy_att( ncid, nf90_global, attname, outid, nf90_global ) )
   END DO
   CALL check_nf90( nf90_put_att( outid, nf90_global, "file_name", TRIM(filebase)//'.nc') )
   IF (l_verbose) WRITE(numout,*) 'Writing new file_name attribute'  
   CALL DATE_AND_TIME ( date=date, time=time, zone=zone )
   timestamp = date(7:8) // "/" // date(5:6) // "/" // date(1:4) // " " // &
               time(1:2) // ":" // time(3:4) // ":" // time(5:6) // " " // &
               zone  
   CALL check_nf90( nf90_put_att( outid, nf90_global, "TimeStamp", timestamp ) )
   IF (l_verbose) WRITE(numout,*) 'Writing new TimeStamp attribute'
  
!2.2.3 Copy the variable definitions and attributes into the output file.
   ALLOCATE(mdiVals(nvars))
   mdiVals(:)=0
   DO jv = 1, nvars
      CALL check_nf90( nf90_inquire_variable( ncid, jv, varname, xtype, ndims, dimids, natts ) )
      ALLOCATE(outdimids(ndims))
      ALLOCATE(chunksizes(ndims))
      IF( ndims > 0 ) then
        DO idim = 1, ndims
           outdimids(idim) = dimids(idim)
           chunksizes(idim) = outdimlens(dimids(idim))
           cdim='|'//TRIM(indimnames(dimids(idim)))//'|'

! trick to find var in a list of suggestion (var0 and var1 : INDEX(|var0|var1|,|var|)
           cdimlst='|x|x_grid_T|x_grid_U|x_grid_V|x_grid_W|'
           if( INDEX(TRIM(cdimlst),TRIM(cdim)) > 0 ) &
    &                             chunksizes(idim) = min(outdimlens(dimids(idim)), max(nc4_xchunk,1))

           cdimlst='|y|y_grid_T|y_grid_U|y_grid_V|y_grid_W|'
           if( INDEX(TRIM(cdimlst),TRIM(cdim)) > 0 ) &
    &                             chunksizes(idim) = min(outdimlens(dimids(idim)), max(nc4_ychunk,1))

           cdimlst='|z|deptht|depthu|depthv|depthw|depth|nav_lev|'
           if( INDEX(TRIM(cdimlst),TRIM(cdim)) > 0 ) &
    &                             chunksizes(idim) = min(outdimlens(dimids(idim)), max(nc4_zchunk,1))

           cdimlst='|t|time|time_counter|'
           if( INDEX(TRIM(cdimlst),TRIM(cdim)) > 0 ) &
    &                             chunksizes(idim) = min(outdimlens(dimids(idim)), max(nc4_tchunk,1))

        END DO
#if defined key_netcdf4
        CALL check_nf90( nf90_def_var( outid, varname, xtype, outdimids, varid, &
                                       deflate_level=deflate_level ) )
        IF (l_verbose) WRITE(numout,*) 'Dims    : ',ndims, outdimids(1:ndims)
        IF (l_verbose) WRITE(numout,*) 'names   : ',(TRIM(indimnames(dimids(idim)))//' ',idim=1,ndims)
        IF (l_verbose) WRITE(numout,*) 'lens    : ',(outdimlens(dimids(idim)),idim=1,ndims)
        IF (l_verbose) WRITE(numout,*) 'Chunking: ',chunksizes
        IF (l_verbose) WRITE(numout,*) 'Deflation : ',deflate_level
        IF (l_verbose) WRITE(numout,*) 'Chunk algo: ',chunkalg
        CALL check_nf90( nf90_def_var_chunking( outid, varid, chunkalg, &
   &                                 chunksizes ) )
      ELSE
        CALL check_nf90( nf90_def_var( outid, varname, xtype, outdimids, varid ) )
#else
      CALL check_nf90( nf90_def_var( outid, varname, xtype, outdimids, varid ) )
#endif
      ENDIF
      DEALLOCATE(outdimids)
      DEALLOCATE(chunksizes)
      IF (l_verbose) WRITE(numout,*) 'Defining variable '//TRIM(varname)//'...' 
      IF( natts > 0 ) THEN
         DO attid = 1, natts
            CALL check_nf90( nf90_inq_attname( ncid, varid, attid, attname ) )
            IF ( attname == "_FillValue" ) THEN
               CALL check_nf90( nf90_get_att( ncid, varid, attname, rmdi ) )
               mdiVals(jv)=rmdi
            ENDIF
            CALL check_nf90( nf90_copy_att( ncid, varid, attname, outid, varid ) )
         END DO
      ENDIF
   END DO

!2.3 End definitions in output file and copy 1st file ncid to the inncids array

   CALL check_nf90( nf90_enddef( outid ) )
   inncids(1) = ncid
   IF (l_verbose) WRITE(numout,*) 'Finished defining output file.'
  
!---------------------------------------------------------------------------
!3. Read in data from each file for each variable 

!3.1 Open each file and store the ncid in inncids array

   IF (l_verbose) WRITE(numout,*) 'Opening input files...'

   ! Set a file chunk cache size for the processor-domain files that scales with the number of processors
   patchchunk = max(8192, fchunksize/ndomain)

   ! open files 
   DO ifile = 2, ndomain
      CALL check_nf90( nf90_open( TRIM(filenames(ifile)), nf90_share, ncid, chunksize=patchchunk ) )
      inncids(ifile) = ncid
   END DO
   IF (l_verbose) WRITE(numout,*) 'All input files open.'

   DO jv = 1, nvars

      ValMin = 1.e10
      ValMax = -1.e10
      l_valid = .false.
      istop = nf90_noerr
      nt = 1
      ntslice = nmax_unlimited
      IF( nslicesize == 0 ) nslicesize = nmax_unlimited

!3.2 Inquire variable to find out name and how many dimensions it has
!    and importantly whether it contains the dimensions in rebuild_dims()

      ncid = inncids(1)
      CALL check_nf90( nf90_inquire_variable( ncid, jv, varname, xtype, ndims, dimids, natts ) )

      l_noRebuild = .true.
      IF( ANY( dimids(1:ndims) == rebuild_dims(1) )) l_noRebuild = .false.
      IF( rbdims > 1 ) THEN
         IF( ANY( dimids(1:ndims) == rebuild_dims(2) )) l_noRebuild = .false.
      ENDIF

!3.2.0 start while loop for time slicing

      DO WHILE( nt <= nmax_unlimited )

         IF( ndims > 3 ) THEN
            ntslice = MIN( nslicesize, nmax_unlimited + 1 - nt )
         ENDIF

      IF (l_noRebuild) THEN

         IF( nslicesize == nmax_unlimited .OR. ndims <= 3 ) THEN
            IF (l_verbose) WRITE(numout,*) 'Copying data from variable '//TRIM(varname)//'...'
         ELSE
            IF (l_verbose) WRITE(numout,'(A,I3,A,I3,A)') ' Copying data from variable '  &
            &                 //TRIM(varname)//' for slices ',nt,' to ',nt+ntslice-1,' ...'
         ENDIF

!3.2.1 If rebuilding not required then just need to read in variable
!      for copying direct into output file after the OMP (files) loop.
         IF( ndims == 0 ) THEN

            SELECT CASE( xtype )
               CASE( NF90_BYTE )
                  CALL check_nf90( nf90_get_var( ncid, jv, globaldata_0d_i1 ) )
               CASE( NF90_SHORT )
                  CALL check_nf90( nf90_get_var( ncid, jv, globaldata_0d_i2 ) )
               CASE( NF90_INT )
                  CALL check_nf90( nf90_get_var( ncid, jv, globaldata_0d_i4 ) )
               CASE( NF90_FLOAT )
                  CALL check_nf90( nf90_get_var( ncid, jv, globaldata_0d_sp ) )
               CASE( NF90_DOUBLE )
                  CALL check_nf90( nf90_get_var( ncid, jv, globaldata_0d_dp ) )
               CASE DEFAULT
                  WRITE(numerr,*) 'Unknown nf90 type: ', xtype
                  STOP 4
            END SELECT

         ELSEIF( ndims == 1 ) THEN

            SELECT CASE( xtype )
               CASE( NF90_BYTE )
                  ALLOCATE(globaldata_1d_i1(indimlens(dimids(1))))
                  CALL check_nf90( nf90_get_var( ncid, jv, globaldata_1d_i1 ) )
               CASE( NF90_SHORT )
                  ALLOCATE(globaldata_1d_i2(indimlens(dimids(1))))
                  CALL check_nf90( nf90_get_var( ncid, jv, globaldata_1d_i2 ) )
               CASE( NF90_INT )
                  ALLOCATE(globaldata_1d_i4(indimlens(dimids(1))))
                  CALL check_nf90( nf90_get_var( ncid, jv, globaldata_1d_i4 ) )
               CASE( NF90_FLOAT )
                  ALLOCATE(globaldata_1d_sp(indimlens(dimids(1))))
                  CALL check_nf90( nf90_get_var( ncid, jv, globaldata_1d_sp ) )
               CASE( NF90_DOUBLE )
                  ALLOCATE(globaldata_1d_dp(indimlens(dimids(1))))
                  CALL check_nf90( nf90_get_var( ncid, jv, globaldata_1d_dp ) )
               CASE DEFAULT
                  WRITE(numerr,*) 'Unknown nf90 type: ', xtype
                  STOP 4
            END SELECT

         ELSEIF( ndims == 2 ) THEN

            SELECT CASE( xtype )
               CASE( NF90_BYTE )
                  ALLOCATE(globaldata_2d_i1(indimlens(dimids(1)),indimlens(dimids(2))))
                  CALL check_nf90( nf90_get_var( ncid, jv, globaldata_2d_i1 ) )
               CASE( NF90_SHORT )
                  ALLOCATE(globaldata_2d_i2(indimlens(dimids(1)),indimlens(dimids(2))))
                  CALL check_nf90( nf90_get_var( ncid, jv, globaldata_2d_i2 ) )
               CASE( NF90_INT )
                  ALLOCATE(globaldata_2d_i4(indimlens(dimids(1)),indimlens(dimids(2))))
                  CALL check_nf90( nf90_get_var( ncid, jv, globaldata_2d_i4 ) )
               CASE( NF90_FLOAT )
                  ALLOCATE(globaldata_2d_sp(indimlens(dimids(1)),indimlens(dimids(2))))
                  CALL check_nf90( nf90_get_var( ncid, jv, globaldata_2d_sp ) )
               CASE( NF90_DOUBLE )
                  ALLOCATE(globaldata_2d_dp(indimlens(dimids(1)),indimlens(dimids(2))))
                  CALL check_nf90( nf90_get_var( ncid, jv, globaldata_2d_dp ) )
               CASE DEFAULT
                  WRITE(numerr,*) 'Unknown nf90 type: ', xtype
                  STOP 4
            END SELECT

         ELSEIF( ndims == 3 ) THEN

            SELECT CASE( xtype )
               CASE( NF90_BYTE )
                  ALLOCATE(globaldata_3d_i1(indimlens(dimids(1)),indimlens(dimids(2)),       &
                     &                      indimlens(dimids(3))))
                  CALL check_nf90( nf90_get_var( ncid, jv, globaldata_3d_i1 ) )
               CASE( NF90_SHORT )
                  ALLOCATE(globaldata_3d_i2(indimlens(dimids(1)),indimlens(dimids(2)),       &
                     &                      indimlens(dimids(3))))
                  CALL check_nf90( nf90_get_var( ncid, jv, globaldata_3d_i2 ) )
               CASE( NF90_INT )
                  ALLOCATE(globaldata_3d_i4(indimlens(dimids(1)),indimlens(dimids(2)),       &
                     &                      indimlens(dimids(3))))
                  CALL check_nf90( nf90_get_var( ncid, jv, globaldata_3d_i4 ) )
               CASE( NF90_FLOAT )
                  ALLOCATE(globaldata_3d_sp(indimlens(dimids(1)),indimlens(dimids(2)),       &
                     &                      indimlens(dimids(3))))
                  CALL check_nf90( nf90_get_var( ncid, jv, globaldata_3d_sp ) )
               CASE( NF90_DOUBLE )
                  ALLOCATE(globaldata_3d_dp(indimlens(dimids(1)),indimlens(dimids(2)),       &
                     &                      indimlens(dimids(3))))
                  CALL check_nf90( nf90_get_var( ncid, jv, globaldata_3d_dp ) )
               CASE DEFAULT
                  WRITE(numerr,*) 'Unknown nf90 type: ', xtype
                  STOP 4
            END SELECT

         ELSEIF( ndims == 4 ) THEN

            SELECT CASE( xtype )
               CASE( NF90_BYTE )
                  ALLOCATE(globaldata_4d_i1(indimlens(dimids(1)),indimlens(dimids(2)),       &
                     &                      indimlens(dimids(3)),ntslice))
                  CALL check_nf90( nf90_get_var( ncid, jv, globaldata_4d_i1, start=(/1,1,1,nt/) ) )
               CASE( NF90_SHORT )
                  ALLOCATE(globaldata_4d_i2(indimlens(dimids(1)),indimlens(dimids(2)),       &
                     &                      indimlens(dimids(3)),ntslice))
                  CALL check_nf90( nf90_get_var( ncid, jv, globaldata_4d_i2, start=(/1,1,1,nt/) ) )
               CASE( NF90_INT )
                  ALLOCATE(globaldata_4d_i4(indimlens(dimids(1)),indimlens(dimids(2)),       &
                     &                      indimlens(dimids(3)),ntslice))
                  CALL check_nf90( nf90_get_var( ncid, jv, globaldata_4d_i4, start=(/1,1,1,nt/) ) )
               CASE( NF90_FLOAT )
                  ALLOCATE(globaldata_4d_sp(indimlens(dimids(1)),indimlens(dimids(2)),       &
                     &                      indimlens(dimids(3)),ntslice))
                  CALL check_nf90( nf90_get_var( ncid, jv, globaldata_4d_sp, start=(/1,1,1,nt/) ) )
               CASE( NF90_DOUBLE )
                  ALLOCATE(globaldata_4d_dp(indimlens(dimids(1)),indimlens(dimids(2)),       &
                     &                      indimlens(dimids(3)),ntslice))
                  CALL check_nf90( nf90_get_var( ncid, jv, globaldata_4d_dp, start=(/1,1,1,nt/) ) )
               CASE DEFAULT
                  WRITE(numerr,*) 'Unknown nf90 type: ', xtype
                  STOP 4
            END SELECT

         ENDIF

      ELSE  ! l_noRebuild = .false.

!3.2.2 For variables that require rebuilding we need to read in from all ndomain files
!      Here we allocate global variables ahead of looping over files
         IF( nslicesize == nmax_unlimited .OR. ndims <= 3 ) THEN
            IF (l_verbose) WRITE(numout,*) 'Rebuilding data from variable '//TRIM(varname)//'...'
         ELSE
            IF (l_verbose) WRITE(numout,'(A,I3,A,I3,A)') ' Rebuilding data from variable '  &
            &                 //TRIM(varname)//' for slices ',nt,' to ',nt+ntslice-1,' ...'
         ENDIF
         IF( ndims == 1 ) THEN

            SELECT CASE( xtype )
               CASE( NF90_BYTE )
                  ALLOCATE(globaldata_1d_i1(outdimlens(dimids(1))))
                  IF (l_maskout) globaldata_1d_i1(:)=mdiVals(jv)
               CASE( NF90_SHORT )
                  ALLOCATE(globaldata_1d_i2(outdimlens(dimids(1))))
                  IF (l_maskout) globaldata_1d_i2(:)=mdiVals(jv)
               CASE( NF90_INT )
                  ALLOCATE(globaldata_1d_i4(outdimlens(dimids(1))))
                  IF (l_maskout) globaldata_1d_i4(:)=mdiVals(jv)
               CASE( NF90_FLOAT )
                  ALLOCATE(globaldata_1d_sp(outdimlens(dimids(1))))
                  IF (l_maskout) globaldata_1d_sp(:)=mdiVals(jv)
               CASE( NF90_DOUBLE )
                  ALLOCATE(globaldata_1d_dp(outdimlens(dimids(1))))
                  IF (l_maskout) globaldata_1d_dp(:)=mdiVals(jv)
               CASE DEFAULT
                  WRITE(numerr,*) 'Unknown nf90 type: ', xtype
                  STOP 4
            END SELECT

         ELSEIF( ndims == 2 ) THEN

            SELECT CASE( xtype )
               CASE( NF90_BYTE )
                  ALLOCATE(globaldata_2d_i1(outdimlens(dimids(1)),outdimlens(dimids(2))))
                  IF (l_maskout) globaldata_2d_i1(:,:)=mdiVals(jv)
               CASE( NF90_SHORT )
                  ALLOCATE(globaldata_2d_i2(outdimlens(dimids(1)),outdimlens(dimids(2))))
                  IF (l_maskout) globaldata_2d_i2(:,:)=mdiVals(jv)
               CASE( NF90_INT )
                  ALLOCATE(globaldata_2d_i4(outdimlens(dimids(1)),outdimlens(dimids(2))))
                  IF (l_maskout) globaldata_2d_i4(:,:)=mdiVals(jv)
               CASE( NF90_FLOAT )
                  ALLOCATE(globaldata_2d_sp(outdimlens(dimids(1)),outdimlens(dimids(2))))
                  IF (l_maskout) globaldata_2d_sp(:,:)=mdiVals(jv)
               CASE( NF90_DOUBLE )
                  ALLOCATE(globaldata_2d_dp(outdimlens(dimids(1)),outdimlens(dimids(2))))
                  IF (l_maskout) globaldata_2d_dp(:,:)=mdiVals(jv)
               CASE DEFAULT
                  WRITE(numerr,*) 'Unknown nf90 type: ', xtype
                  STOP 4
            END SELECT

         ELSEIF( ndims == 3 ) THEN

            SELECT CASE( xtype )
               CASE( NF90_BYTE )
                  ALLOCATE(globaldata_3d_i1(outdimlens(dimids(1)),outdimlens(dimids(2)),     &
                     &                      outdimlens(dimids(3))))
                  IF (l_maskout) globaldata_3d_i1(:,:,:)=mdiVals(jv)
               CASE( NF90_SHORT )
                  ALLOCATE(globaldata_3d_i2(outdimlens(dimids(1)),outdimlens(dimids(2)),     &
                     &                      outdimlens(dimids(3))))
                  IF (l_maskout) globaldata_3d_i2(:,:,:)=mdiVals(jv)
               CASE( NF90_INT )
                  ALLOCATE(globaldata_3d_i4(outdimlens(dimids(1)),outdimlens(dimids(2)),     &
                     &                      outdimlens(dimids(3))))
                  IF (l_maskout) globaldata_3d_i4(:,:,:)=mdiVals(jv)
               CASE( NF90_FLOAT )
                  ALLOCATE(globaldata_3d_sp(outdimlens(dimids(1)),outdimlens(dimids(2)),     &
                     &                      outdimlens(dimids(3))))
                  IF (l_maskout) globaldata_3d_sp(:,:,:)=mdiVals(jv)
               CASE( NF90_DOUBLE )
                  ALLOCATE(globaldata_3d_dp(outdimlens(dimids(1)),outdimlens(dimids(2)),     &
                     &                      outdimlens(dimids(3))))
                  IF (l_maskout) globaldata_3d_dp(:,:,:)=mdiVals(jv)
               CASE DEFAULT
                  WRITE(numerr,*) 'Unknown nf90 type: ', xtype
                  STOP 4
            END SELECT

         ELSEIF( ndims == 4 ) THEN

            SELECT CASE( xtype )
               CASE( NF90_BYTE )
                  ALLOCATE(globaldata_4d_i1(outdimlens(dimids(1)),outdimlens(dimids(2)),     &
                     &                      outdimlens(dimids(3)),ntslice))
                  IF (l_maskout) globaldata_4d_i1(:,:,:,:)=mdiVals(jv)
               CASE( NF90_SHORT )
                  ALLOCATE(globaldata_4d_i2(outdimlens(dimids(1)),outdimlens(dimids(2)),     &
                     &                      outdimlens(dimids(3)),ntslice))
                  IF (l_maskout) globaldata_4d_i2(:,:,:,:)=mdiVals(jv)
               CASE( NF90_INT )
                  ALLOCATE(globaldata_4d_i4(outdimlens(dimids(1)),outdimlens(dimids(2)),     &
                     &                      outdimlens(dimids(3)),ntslice))
                  IF (l_maskout) globaldata_4d_i4(:,:,:,:)=mdiVals(jv)
               CASE( NF90_FLOAT )
                  ALLOCATE(globaldata_4d_sp(outdimlens(dimids(1)),outdimlens(dimids(2)),     &
                     &                      outdimlens(dimids(3)),ntslice))
                  IF (l_maskout) globaldata_4d_sp(:,:,:,:)=mdiVals(jv)
               CASE( NF90_DOUBLE )
                  ALLOCATE(globaldata_4d_dp(outdimlens(dimids(1)),outdimlens(dimids(2)),     &
                     &                      outdimlens(dimids(3)),ntslice))
                  IF (l_maskout) globaldata_4d_dp(:,:,:,:)=mdiVals(jv)
               CASE DEFAULT
                  WRITE(numerr,*) 'Unknown nf90 type: ', xtype
                  STOP 4
            END SELECT
         ELSE
            WRITE(numerr,*) 'ERROR! : A netcdf variable has more than 4 dimensions which is not taken into account'
            STOP 4
         ENDIF

!$OMP  PARALLEL DO DEFAULT(NONE)                                                          &
!$OMP& PRIVATE(ifile,ncid,xtype,start_pos,local_sizes,InMin,InMax,natts,                  &
!$OMP&         ndims,attid,attname,dimids,idim,dimname,dimlen,unlimitedDimId,             &
!$OMP&         halo_start,halo_end,idomain,jdomain,rdomain,di,dj,dr,                      &
!$OMP&         localdata_1d_i2,localdata_1d_i4,localdata_1d_sp,localdata_1d_dp,           &
!$OMP&         localdata_2d_i2,localdata_2d_i4,localdata_2d_sp,localdata_2d_dp,           &
!$OMP&         localdata_3d_i2,localdata_3d_i4,localdata_3d_sp,localdata_3d_dp,           &
!$OMP&         localdata_4d_i2,localdata_4d_i4,localdata_4d_sp,localdata_4d_dp,           &
!$OMP&         localdata_1d_i1,localdata_2d_i1,localdata_3d_i1,localdata_4d_i1)           &
!$OMP& SHARED(jv,nvars,varname,filenames,ValMin,ValMax,indimlens,outdimlens,rbdims,       &
!$OMP&        ndomain,outid,fchunksize,istop,l_valid,nthreads,inncids,rebuild_dims,       &
!$OMP&        globaldata_1d_i2,globaldata_1d_i4,globaldata_1d_sp,globaldata_1d_dp,        &
!$OMP&        globaldata_2d_i2,globaldata_2d_i4,globaldata_2d_sp,globaldata_2d_dp,        &
!$OMP&        globaldata_3d_i2,globaldata_3d_i4,globaldata_3d_sp,globaldata_3d_dp,        &
!$OMP&        globaldata_4d_i2,globaldata_4d_i4,globaldata_4d_sp,globaldata_4d_dp,        &
!$OMP&        globaldata_1d_i1,globaldata_2d_i1,globaldata_3d_i1,globaldata_4d_i1,        &
!$OMP&        ntslice,nt,nmax_unlimited,indimnames,dims,patchchunk)

         DO ifile = 1, ndomain

            ncid = inncids(ifile)
!$OMP CRITICAL
            CALL check_nf90( nf90_get_att( ncid, nf90_global, 'DOMAIN_size_local', local_sizes ), istop )
            CALL check_nf90( nf90_get_att( ncid, nf90_global, 'DOMAIN_position_first', start_pos ), istop )
            CALL check_nf90( nf90_get_att( ncid, nf90_global, 'DOMAIN_halo_size_start', halo_start ), istop )
            CALL check_nf90( nf90_get_att( ncid, nf90_global, 'DOMAIN_halo_size_end', halo_end ), istop )
            CALL check_nf90( nf90_inquire_variable( ncid, jv, varname, xtype, ndims, dimids, natts ), istop )
            CALL check_nf90( nf90_inquire( ncid, unlimitedDimId = unlimitedDimId ), istop )
!$OMP END CRITICAL

            ! set defaults for rebuilding so that i is 1st, j 2nd
            di=1
            dj=2

            IF( rbdims == 1 ) THEN
               ! override defaults above and set other variables
               start_pos(2) = 1
               local_sizes(2) = outdimlens(3-dimids(2))
               halo_end(2) = 0
               halo_start(2) = 0
               di=rebuild_dims(1)
               dj=3-di
            ENDIF

!3.3.1 Generate local domain interior sizes from local_sizes and halo sizes
!      idomain defines the 1st and last interior points in the i direction and
!      jdomain defines the 1st and last interior points in the j direction

            idomain(1) = 1 + halo_start(di)
            idomain(2) = local_sizes(di) - halo_end(di)
            jdomain(1) = 1 + halo_start(dj)
            jdomain(2) = local_sizes(dj) - halo_end(dj)

!3.3.2 For rbdims or more dimensions put the data array from this input file into the correct
!      part of the output data array. Assume the first dimensions are those to be rebuilt.

            IF( ndims == 1 ) THEN

               IF( rebuild_dims(1) == 1 ) THEN
                  dr = di
                  rdomain = idomain
               ELSE
                  dr = dj
                  rdomain = jdomain
               ENDIF

              SELECT CASE( xtype )
                  CASE( NF90_BYTE )
                     ALLOCATE(localdata_1d_i1(local_sizes(dr)))
                     CALL check_nf90( nf90_get_var( ncid, jv, localdata_1d_i1 ), istop )
                     DO jr = rdomain(1), rdomain(2)
                        globaldata_1d_i1(start_pos(dr) + jr - 1) = localdata_1d_i1(jr)
                     END DO
                     DEALLOCATE(localdata_1d_i1)
                  CASE( NF90_SHORT )
                     ALLOCATE(localdata_1d_i2(local_sizes(dr)))
                     CALL check_nf90( nf90_get_var( ncid, jv, localdata_1d_i2 ), istop )
                     DO jr = rdomain(1), rdomain(2)
                        globaldata_1d_i2(start_pos(dr) + jr - 1) = localdata_1d_i2(jr)
                     END DO
                     DEALLOCATE(localdata_1d_i2)
                  CASE( NF90_INT )
                     ALLOCATE(localdata_1d_i4(local_sizes(dr)))
                     CALL check_nf90( nf90_get_var( ncid, jv, localdata_1d_i4 ), istop )
                     DO jr = rdomain(1), rdomain(2)
                        globaldata_1d_i4(start_pos(dr) + jr - 1) = localdata_1d_i4(jr)
                     END DO
                     DEALLOCATE(localdata_1d_i4)
                  CASE( NF90_FLOAT )
                     ALLOCATE(localdata_1d_sp(local_sizes(dr)))
                     CALL check_nf90( nf90_get_var( ncid, jv, localdata_1d_sp ), istop )
                     DO jr = rdomain(1), rdomain(2)
                        globaldata_1d_sp(start_pos(dr) + jr - 1) = localdata_1d_sp(jr)
                     END DO
                     DEALLOCATE(localdata_1d_sp)
                  CASE( NF90_DOUBLE )
                     ALLOCATE(localdata_1d_dp(local_sizes(dr)))
                     CALL check_nf90( nf90_get_var( ncid, jv, localdata_1d_dp ), istop )
                     DO jr = rdomain(1), rdomain(2)
                        globaldata_1d_dp(start_pos(dr) + jr - 1) = localdata_1d_dp(jr)
                     END DO
                     DEALLOCATE(localdata_1d_dp)
                  CASE DEFAULT
                     WRITE(numerr,*) 'Unknown nf90 type: ', xtype
                     istop = istop + 1
               END SELECT

            ELSEIF( ndims == 2 ) THEN
        
               SELECT CASE( xtype )
                  CASE( NF90_BYTE )
                     ALLOCATE(localdata_2d_i1(local_sizes(di),local_sizes(dj)))
                     CALL check_nf90( nf90_get_var( ncid, jv, localdata_2d_i1 ), istop )
                        DO jj = jdomain(1), jdomain(2)
                           DO ji = idomain(1), idomain(2)
                              globaldata_2d_i1(start_pos(di) + ji - 1, start_pos(dj) + jj - 1) = localdata_2d_i1(ji,jj)
                        END DO
                     END DO
                     DEALLOCATE(localdata_2d_i1)
                  CASE( NF90_SHORT )
                     ALLOCATE(localdata_2d_i2(local_sizes(di),local_sizes(dj)))
                     CALL check_nf90( nf90_get_var( ncid, jv, localdata_2d_i2 ), istop )
                     DO jj = jdomain(1), jdomain(2)
                        DO ji = idomain(1), idomain(2)
                           globaldata_2d_i2(start_pos(di) + ji - 1, start_pos(dj) + jj - 1) = localdata_2d_i2(ji,jj)
                        END DO
                     END DO
                     DEALLOCATE(localdata_2d_i2)
                  CASE( NF90_INT )
                     ALLOCATE(localdata_2d_i4(local_sizes(di),local_sizes(dj)))
                     CALL check_nf90( nf90_get_var( ncid, jv, localdata_2d_i4 ), istop )
                     DO jj = jdomain(1), jdomain(2)
                        DO ji = idomain(1), idomain(2)
                           globaldata_2d_i4(start_pos(di) + ji - 1, start_pos(dj) + jj - 1) = localdata_2d_i4(ji,jj)
                        END DO
                     END DO
                     DEALLOCATE(localdata_2d_i4)
                  CASE( NF90_FLOAT )
                     ALLOCATE(localdata_2d_sp(local_sizes(di),local_sizes(dj)))
                     CALL check_nf90( nf90_get_var( ncid, jv, localdata_2d_sp ), istop )
                     DO jj = jdomain(1), jdomain(2)
                        DO ji = idomain(1), idomain(2)
                           globaldata_2d_sp(start_pos(di) + ji - 1, start_pos(dj) + jj - 1) = localdata_2d_sp(ji,jj)
                        END DO
                     END DO
                     DEALLOCATE(localdata_2d_sp)
                  CASE( NF90_DOUBLE )
                     ALLOCATE(localdata_2d_dp(local_sizes(di),local_sizes(dj)))
                     CALL check_nf90( nf90_get_var( ncid, jv, localdata_2d_dp ), istop )
                     DO jj = jdomain(1), jdomain(2)
                        DO ji = idomain(1), idomain(2)
                           globaldata_2d_dp(start_pos(di) + ji - 1, start_pos(dj) + jj - 1) = localdata_2d_dp(ji,jj)
                        END DO
                     END DO
                     DEALLOCATE(localdata_2d_dp) 
                  CASE DEFAULT
                     WRITE(numerr,*) 'Unknown nf90 type: ', xtype
                     istop = istop + 1
               END SELECT

            ELSEIF( ndims == 3 ) THEN

               SELECT CASE( xtype )
                  CASE( NF90_BYTE )
                     ALLOCATE(localdata_3d_i1(local_sizes(di),local_sizes(dj),indimlens(dimids(3))))
                     CALL check_nf90( nf90_get_var( ncid, jv, localdata_3d_i1 ), istop )
!$OMP  PARALLEL DO DEFAULT(NONE) PRIVATE(ji,jj,jk)   &
!$OMP& SHARED(idomain,jdomain,indimlens,dimids,start_pos,globaldata_3d_i1,localdata_3d_i1,di,dj) 
                     DO jk = 1, indimlens(dimids(3))
                        DO jj = jdomain(1), jdomain(2)
                           DO ji = idomain(1), idomain(2)
                              globaldata_3d_i1(start_pos(di) + ji - 1, start_pos(dj) + jj - 1, jk) = localdata_3d_i1(ji,jj,jk)
                           END DO
                        END DO
                     END DO
!$OMP END PARALLEL DO
                     DEALLOCATE(localdata_3d_i1)
                  CASE( NF90_SHORT )
                     ALLOCATE(localdata_3d_i2(local_sizes(di),local_sizes(dj),indimlens(dimids(3))))
                     CALL check_nf90( nf90_get_var( ncid, jv, localdata_3d_i2 ), istop )
!$OMP  PARALLEL DO DEFAULT(NONE) PRIVATE(ji,jj,jk)   &
!$OMP& SHARED(idomain,jdomain,indimlens,dimids,start_pos,globaldata_3d_i2,localdata_3d_i2,di,dj) 
                     DO jk = 1, indimlens(dimids(3))
                        DO jj = jdomain(1), jdomain(2)
                           DO ji = idomain(1), idomain(2)
                              globaldata_3d_i2(start_pos(di) + ji - 1, start_pos(dj) + jj - 1, jk) = localdata_3d_i2(ji,jj,jk)
                           END DO
                        END DO
                     END DO
!$OMP END PARALLEL DO
                     DEALLOCATE(localdata_3d_i2)
                  CASE( NF90_INT )
                     ALLOCATE(localdata_3d_i4(local_sizes(di),local_sizes(dj),indimlens(dimids(3))))
                     CALL check_nf90( nf90_get_var( ncid, jv, localdata_3d_i4 ), istop )
!$OMP  PARALLEL DO DEFAULT(NONE) PRIVATE(ji,jj,jk)   &
!$OMP& SHARED(idomain,jdomain,indimlens,dimids,start_pos,globaldata_3d_i4,localdata_3d_i4,di,dj) 
                     DO jk = 1, indimlens(dimids(3))
                        DO jj = jdomain(1), jdomain(2)
                           DO ji = idomain(1), idomain(2)
                              globaldata_3d_i4(start_pos(di) + ji - 1, start_pos(dj) + jj - 1, jk) = localdata_3d_i4(ji,jj,jk)
                           END DO
                        END DO
                     END DO
!$OMP END PARALLEL DO
                     DEALLOCATE(localdata_3d_i4)
                  CASE( NF90_FLOAT )
                     ! TG: This if statement is added to check if the 1st dimension is the corners (for lon_bounds) variables
                     ! TG: Had to add the unsatisfactory check for 'lon' as it failed for diaptr files
                     ! TG: Would like to find a better assumption for this.
                     IF ( trim(indimnames(dimids(1))) /= dims(1) .AND. indimnames(dimids(1)) .NE. 'lon' ) THEN 
                        ALLOCATE(localdata_3d_sp(indimlens(dimids(1)),local_sizes(di),local_sizes(dj)))
                        WRITE(*,*) 'test', ifile, jv, indimlens(dimids(1)),local_sizes(di),local_sizes(dj)
                        CALL check_nf90( nf90_get_var( ncid, jv, localdata_3d_sp ), istop )
                        WRITE(*,*) 'test2'
                        DO jj = jdomain(1), jdomain(2)
                           DO ji = idomain(1), idomain(2)
                              DO jk = 1, indimlens(dimids(1))
                                 globaldata_3d_sp(jk, start_pos(di) + ji - 1, start_pos(dj) + jj - 1) = localdata_3d_sp(jk,ji,jj)
                              END DO
                           END DO
                        END DO
                     ELSE
                        ALLOCATE(localdata_3d_sp(local_sizes(di),local_sizes(dj),indimlens(dimids(3))))
                        CALL check_nf90( nf90_get_var( ncid, jv, localdata_3d_sp ), istop ) 
!$OMP  PARALLEL DO DEFAULT(NONE) PRIVATE(ji,jj,jk)   &
!$OMP& SHARED(idomain,jdomain,indimlens,dimids,start_pos,globaldata_3d_sp,localdata_3d_sp,di,dj) 
                        DO jk = 1, indimlens(dimids(3))
                           DO jj = jdomain(1), jdomain(2)
                              DO ji = idomain(1), idomain(2)
                                 globaldata_3d_sp(start_pos(di) + ji - 1, start_pos(dj) + jj - 1, jk) = localdata_3d_sp(ji,jj,jk)
                              END DO
                           END DO
                        END DO
!$OMP END PARALLEL DO
                     ENDIF
                     DEALLOCATE(localdata_3d_sp)
                  CASE( NF90_DOUBLE )
                     IF ( trim(indimnames(dimids(1))) /= dims(1) ) THEN
                        ALLOCATE(localdata_3d_dp(indimlens(dimids(1)),local_sizes(di),local_sizes(dj)))
                        CALL check_nf90( nf90_get_var( ncid, jv, localdata_3d_dp ), istop ) 
                        DO jj = jdomain(1), jdomain(2)
                           DO ji = idomain(1), idomain(2)
                              DO jk = 1, indimlens(dimids(1))
                                 globaldata_3d_dp(jk, start_pos(di) + ji - 1, start_pos(dj) + jj - 1) = localdata_3d_dp(jk,ji,jj)
                              END DO
                           END DO
                        END DO
                     ELSE
                        ALLOCATE(localdata_3d_dp(local_sizes(di),local_sizes(dj),indimlens(dimids(3))))
                        CALL check_nf90( nf90_get_var( ncid, jv, localdata_3d_dp ), istop ) 
!$OMP  PARALLEL DO DEFAULT(NONE) PRIVATE(ji,jj,jk)   &
!$OMP& SHARED(idomain,jdomain,indimlens,dimids,start_pos,globaldata_3d_dp,localdata_3d_dp,di,dj)
                        DO jk = 1, indimlens(dimids(3))
                           DO jj = jdomain(1), jdomain(2)
                              DO ji = idomain(1), idomain(2)
                                 globaldata_3d_dp(start_pos(di) + ji - 1, start_pos(dj) + jj - 1, jk) = localdata_3d_dp(ji,jj,jk)
                              END DO
                           END DO
                        END DO
!$OMP END PARALLEL DO
                     ENDIF
                     DEALLOCATE(localdata_3d_dp)
                  CASE DEFAULT
                     WRITE(numerr,*) 'Unknown nf90 type: ', xtype
                     istop = istop + 1
               END SELECT
      
            ELSEIF (ndims == 4) THEN
        
               SELECT CASE( xtype )
                  CASE( NF90_BYTE )
                     ALLOCATE(localdata_4d_i1(local_sizes(di),local_sizes(dj),               &
                         &                     indimlens(dimids(3)),ntslice))
                     CALL check_nf90( nf90_get_var( ncid, jv, localdata_4d_i1, start=(/1,1,1,nt/) ), istop )
!$OMP  PARALLEL DEFAULT(NONE) PRIVATE(ji,jj,jk,jl)   &
!$OMP& SHARED(idomain,jdomain,indimlens,dimids,start_pos,globaldata_4d_i1,localdata_4d_i1,di,dj,nt,ntslice)
                     DO jl = 1, ntslice
!$OMP DO 
                        DO jk = 1, indimlens(dimids(3))
                           DO jj = jdomain(1), jdomain(2)
                              DO ji = idomain(1), idomain(2)
                                 globaldata_4d_i1(start_pos(di) + ji - 1, start_pos(dj) + jj - 1, jk, jl) = localdata_4d_i1(ji,jj,jk,jl)
                              END DO
                           END DO
                        END DO
!$OMP END DO nowait
                     END DO
!$OMP END PARALLEL
                     DEALLOCATE(localdata_4d_i1)
                  CASE( NF90_SHORT )
                     ALLOCATE(localdata_4d_i2(local_sizes(di),local_sizes(dj),               &
                        &                     indimlens(dimids(3)),ntslice))
                     CALL check_nf90( nf90_get_var( ncid, jv, localdata_4d_i2, start=(/1,1,1,nt/) ), istop )
!$OMP  PARALLEL DEFAULT(NONE) PRIVATE(ji,jj,jk,jl)   &
!$OMP& SHARED(idomain,jdomain,indimlens,dimids,start_pos,globaldata_4d_i2,localdata_4d_i2,di,dj,nt,ntslice)
                     DO jl = 1, ntslice
!$OMP DO 
                        DO jk = 1, indimlens(dimids(3))
                           DO jj = jdomain(1), jdomain(2) 
                              DO ji = idomain(1), idomain(2)
                                 globaldata_4d_i2(start_pos(di) + ji - 1, start_pos(dj) + jj - 1, jk, jl) = localdata_4d_i2(ji,jj,jk,jl)
                              END DO
                           END DO
                        END DO
!$OMP END DO nowait
                     END DO
!$OMP END PARALLEL
                     DEALLOCATE(localdata_4d_i2)
                  CASE( NF90_INT )
                     ALLOCATE(localdata_4d_i4(local_sizes(di),local_sizes(dj),               &
                        &                     indimlens(dimids(3)),ntslice))
                     CALL check_nf90( nf90_get_var( ncid, jv, localdata_4d_i4, start=(/1,1,1,nt/) ), istop )
!$OMP  PARALLEL DEFAULT(NONE) PRIVATE(ji,jj,jk,jl)   &
!$OMP& SHARED(idomain,jdomain,indimlens,dimids,start_pos,globaldata_4d_i4,localdata_4d_i4,di,dj,nt,ntslice)
                     DO jl = 1, ntslice
!$OMP DO
                        DO jk = 1, indimlens(dimids(3))
                           DO jj = jdomain(1), jdomain(2) 
                              DO ji = idomain(1), idomain(2)
                                 globaldata_4d_i4(start_pos(di) + ji - 1, start_pos(dj) + jj - 1, jk, jl) = localdata_4d_i4(ji,jj,jk,jl)
                              END DO
                           END DO
                        END DO
!$OMP END DO nowait
                     END DO
!$OMP END PARALLEL
                     DEALLOCATE(localdata_4d_i4)
                  CASE( NF90_FLOAT )
                     ALLOCATE(localdata_4d_sp(local_sizes(di),local_sizes(dj),               &
                        &                     indimlens(dimids(3)),ntslice))
                     CALL check_nf90( nf90_get_var( ncid, jv, localdata_4d_sp, start=(/1,1,1,nt/) ), istop )
!$OMP  PARALLEL DEFAULT(NONE) PRIVATE(ji,jj,jk,jl)   &
!$OMP& SHARED(idomain,jdomain,indimlens,dimids,start_pos,globaldata_4d_sp,localdata_4d_sp,di,dj,nt,ntslice) 
                     DO jl = 1, ntslice
!$OMP DO
                        DO jk = 1, indimlens(dimids(3))
                           DO jj = jdomain(1), jdomain(2) 
                              DO ji = idomain(1), idomain(2)
                                 globaldata_4d_sp(start_pos(di) + ji - 1, start_pos(dj) + jj - 1, jk, jl) = localdata_4d_sp(ji,jj,jk,jl)
                              END DO
                           END DO
                        END DO
!$OMP END DO nowait
                     END DO
!$OMP END PARALLEL
                     DEALLOCATE(localdata_4d_sp)
                  CASE( NF90_DOUBLE )
                     ALLOCATE(localdata_4d_dp(local_sizes(di),local_sizes(dj),               &
                        &                     indimlens(dimids(3)),ntslice))
                     CALL check_nf90( nf90_get_var( ncid, jv, localdata_4d_dp, start=(/1,1,1,nt/) ), istop )
!$OMP  PARALLEL DEFAULT(NONE) PRIVATE(ji,jj,jk,jl)   &
!$OMP& SHARED(idomain,jdomain,indimlens,dimids,start_pos,globaldata_4d_dp,localdata_4d_dp,di,dj,nt,ntslice) 
                     DO jl = 1, ntslice
!$OMP DO
                        DO jk = 1, indimlens(dimids(3))
                           DO jj = jdomain(1), jdomain(2) 
                              DO ji = idomain(1), idomain(2)
                                 globaldata_4d_dp(start_pos(di) + ji - 1, start_pos(dj) + jj - 1, jk, jl) = localdata_4d_dp(ji,jj,jk,jl)
                              END DO
                           END DO
                        END DO
!$OMP END DO nowait
                     END DO
!$OMP END PARALLEL
                     DEALLOCATE(localdata_4d_dp)
                  CASE DEFAULT
                     WRITE(numerr,*) 'Unknown nf90 type: ', xtype
                     istop = istop + 1
               END SELECT

            ENDIF ! l_noRebuild false

!3.4 Work out if the valid_min and valid_max attributes exist for this variable.
!    If they do then calculate the extrema over all input files.

            DO attid = 1, natts
               CALL check_nf90( nf90_inq_attname( ncid, jv, attid, attname ), istop )
               IF( INDEX( attname, "valid_min" ) == 1 ) THEN
                  CALL check_nf90( nf90_get_att( ncid, jv, attname, InMin), istop )
                  l_valid = .true.
               ENDIF
               IF( INDEX( attname, "valid_max" ) == 1 ) THEN
                  CALL check_nf90( nf90_get_att( ncid, jv, attname, InMax ), istop )
                  l_valid = .true.
               ENDIF
            END DO

            IF (l_valid) THEN
!$OMP CRITICAL
               IF( InMin < ValMin ) ValMin = InMin
               IF( InMax > ValMax ) ValMax = InMax
!$OMP END CRITICAL
            ENDIF

!3.5 Abort if failure and only 1 thread 

            IF( nthreads == 1 .AND. istop /= nf90_noerr )  THEN
               WRITE(numerr,*) '*** NEMO rebuild failed! ***'
               STOP 5
            ENDIF
        
         END DO  ! loop over files
!$OMP END PARALLEL DO

!3.6 Abort if any of the OMP threads failed
         IF( istop /= nf90_noerr )  THEN
            WRITE(numerr,*) '*** NEMO rebuild failed! ***'
            STOP 5
         ENDIF

      ENDIF ! ndims > 2

!---------------------------------------------------------------------------
!4. Write data to output file

      IF (l_verbose) WRITE(numout,*) 'Writing variable '//TRIM(varname)//'...'

!4.1 If the valid min and max attributes exist then update them in the file

      IF( l_valid ) THEN
         CALL check_nf90( nf90_put_att( outid, jv, "valid_min", ValMin ) )
         CALL check_nf90( nf90_put_att( outid, jv, "valid_max", ValMax ) )
      ENDIF

!4.2 Write the data to the output file depending on how many dimensions it has

      IF( ndims == 0 ) THEN

         SELECT CASE( xtype )
            CASE( NF90_BYTE )
               CALL check_nf90( nf90_put_var( outid, jv, globaldata_0d_i1 ) )
            CASE( NF90_SHORT )
               CALL check_nf90( nf90_put_var( outid, jv, globaldata_0d_i2 ) )
            CASE( NF90_INT )
               CALL check_nf90( nf90_put_var( outid, jv, globaldata_0d_i4 ) )
            CASE( NF90_FLOAT )
               CALL check_nf90( nf90_put_var( outid, jv, globaldata_0d_sp ) )
            CASE( NF90_DOUBLE )
               CALL check_nf90( nf90_put_var( outid, jv, globaldata_0d_dp ) )
            CASE DEFAULT   
               WRITE(numerr,*) '0d Unknown nf90 type: ', xtype
               STOP 4
         END SELECT

      ELSEIF( ndims == 1 ) THEN

         SELECT CASE( xtype )
            CASE( NF90_BYTE )
               CALL check_nf90( nf90_put_var( outid, jv, globaldata_1d_i1 ) )
               DEALLOCATE(globaldata_1d_i1)
            CASE( NF90_SHORT )
               CALL check_nf90( nf90_put_var( outid, jv, globaldata_1d_i2 ) )
               DEALLOCATE(globaldata_1d_i2)
            CASE( NF90_INT )
               CALL check_nf90( nf90_put_var( outid, jv, globaldata_1d_i4 ) )
               DEALLOCATE(globaldata_1d_i4)
            CASE( NF90_FLOAT )
               CALL check_nf90( nf90_put_var( outid, jv, globaldata_1d_sp ) )
               DEALLOCATE(globaldata_1d_sp)
            CASE( NF90_DOUBLE )
               CALL check_nf90( nf90_put_var( outid, jv, globaldata_1d_dp ) )
               DEALLOCATE(globaldata_1d_dp)
            CASE DEFAULT   
               WRITE(numerr,*) '1d Unknown nf90 type: ', xtype
               STOP 4
         END SELECT

      ELSEIF( ndims == 2 ) THEN  
     
         SELECT CASE( xtype )   
            CASE( NF90_BYTE )                   
               CALL check_nf90( nf90_put_var( outid, jv, globaldata_2d_i1 ) )
               DEALLOCATE(globaldata_2d_i1)
            CASE( NF90_SHORT )                   
               CALL check_nf90( nf90_put_var( outid, jv, globaldata_2d_i2 ) )
               DEALLOCATE(globaldata_2d_i2)
            CASE( NF90_INT )                              
               CALL check_nf90( nf90_put_var( outid, jv, globaldata_2d_i4 ) )
               DEALLOCATE(globaldata_2d_i4)
            CASE( NF90_FLOAT )                              
               CALL check_nf90( nf90_put_var( outid, jv, globaldata_2d_sp ) )
               DEALLOCATE(globaldata_2d_sp)
            CASE( NF90_DOUBLE )                                         
               CALL check_nf90( nf90_put_var( outid, jv, globaldata_2d_dp ) )
               DEALLOCATE(globaldata_2d_dp)
            CASE DEFAULT   
               WRITE(numerr,*) '2d Unknown nf90 type: ', xtype
               STOP 4
         END SELECT     
                      
      ELSEIF( ndims == 3 ) THEN
      
         SELECT CASE( xtype ) 
            CASE( NF90_BYTE )                   
               CALL check_nf90( nf90_put_var( outid, jv, globaldata_3d_i1 ) )
               DEALLOCATE(globaldata_3d_i1)
            CASE( NF90_SHORT )                   
               CALL check_nf90( nf90_put_var( outid, jv, globaldata_3d_i2 ) )
               DEALLOCATE(globaldata_3d_i2)
            CASE( NF90_INT )                              
               CALL check_nf90( nf90_put_var( outid, jv, globaldata_3d_i4 ) )
               DEALLOCATE(globaldata_3d_i4)
            CASE( NF90_FLOAT )                              
               CALL check_nf90( nf90_put_var( outid, jv, globaldata_3d_sp ) )
               DEALLOCATE(globaldata_3d_sp)
            CASE( NF90_DOUBLE )                                         
               CALL check_nf90( nf90_put_var( outid, jv, globaldata_3d_dp ) )
               DEALLOCATE(globaldata_3d_dp)
            CASE DEFAULT   
               WRITE(numerr,*) '3d Unknown nf90 type: ', xtype
               STOP 4
         END SELECT     
    
      ELSEIF( ndims == 4 ) THEN

         SELECT CASE( xtype )   
            CASE( NF90_BYTE )                   
               CALL check_nf90( nf90_put_var( outid, jv, globaldata_4d_i1, start=(/1,1,1,nt/) ) )
               DEALLOCATE(globaldata_4d_i1)
            CASE( NF90_SHORT )                   
               CALL check_nf90( nf90_put_var( outid, jv, globaldata_4d_i2, start=(/1,1,1,nt/) ) )
               DEALLOCATE(globaldata_4d_i2)
            CASE( NF90_INT )                              
               CALL check_nf90( nf90_put_var( outid, jv, globaldata_4d_i4, start=(/1,1,1,nt/) ) )
               DEALLOCATE(globaldata_4d_i4)
            CASE( NF90_FLOAT )                              
               CALL check_nf90( nf90_put_var( outid, jv, globaldata_4d_sp, start=(/1,1,1,nt/) ) )
               DEALLOCATE(globaldata_4d_sp)
            CASE( NF90_DOUBLE )                                         
               CALL check_nf90( nf90_put_var( outid, jv, globaldata_4d_dp, start=(/1,1,1,nt/) ) )
               DEALLOCATE(globaldata_4d_dp)
            CASE DEFAULT   
               WRITE(numerr,*) '4d Unknown nf90 type: ', xtype
               STOP 4
         END SELECT     
         ! why only for big data set, test the cost.
         CALL check_nf90( nf90_sync( outid ) )    ! flush buffers to disk after writing big 4D datasets
    
      ENDIF

         nt = nt + ntslice

      END DO ! WHILE loop
    
   END DO  ! loop over variables

!---------------------------------------------------------------------------
!5. Close files

!5.1 Close input files

   IF (l_verbose) WRITE(numout,*) 'Closing input files...'
   DO ifile = 1, ndomain
      ncid = inncids(ifile)
      CALL check_nf90( nf90_close( ncid ) )
   END DO

!5.2 Close output file

   IF (l_verbose) WRITE(numout,*) 'Closing output file...'
   CALL check_nf90( nf90_close( outid ) )

   IF (l_verbose) WRITE(numout,*) 'NEMO rebuild completed successfully'
   IF (l_verbose) WRITE(numout,*)

CONTAINS


   SUBROUTINE check_nf90(status, errorFlag)
   !---------------------------------------------------------------------
   !  Checks return code from nf90 library calls and warns if needed
   !  If errorFlag is present then it just increments this flag (OMP use)
   !
   !---------------------------------------------------------------------
      INTEGER, INTENT(IN   ) :: status
      INTEGER, INTENT(INOUT), OPTIONAL :: errorFlag
   !---------------------------------------------------------------------

      IF( status /= nf90_noerr ) THEN
         WRITE(numerr,*) 'ERROR! : '//TRIM(nf90_strerror(status))
         IF( PRESENT( errorFlag ) ) THEN
            errorFlag = errorFlag + status
         ELSE
            WRITE(numerr,*) "*** NEMO rebuild failed ***"
            WRITE(numerr,*)
            STOP 5
         ENDIF
      ENDIF

   END SUBROUTINE check_nf90


END PROGRAM rebuild_nemo
