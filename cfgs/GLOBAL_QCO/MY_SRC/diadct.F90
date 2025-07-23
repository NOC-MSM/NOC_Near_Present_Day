MODULE diadct
   !!======================================================================
   !!                       ***  MODULE  diadct  ***
   !! Ocean diagnostics: Compute the transport trough a sec.
   !!======================================================================
   !! History :  OPA  ! 02/1999 (Y Drillet)  original code
   !!                 ! 10/2001 (Y Drillet, R Bourdalle Badie)
   !!   NEMO     1.0  ! 10/2005 (M Laborie) F90
   !!            3.0  ! 04/2007 (G Garric) Ice sections
   !!             -   ! 04/2007 (C Bricaud) test on sec%nb_point, initialisation of ztransp1,ztransp2,...
   !!            3.4  ! 09/2011 (C Bricaud)
   !!----------------------------------------------------------------------
#if ! defined key_agrif
   !!                        ==>>  CAUTION: does not work with agrif
   !!----------------------------------------------------------------------
   !!   dia_dct      :  Compute the transport through a sec.
   !!   dia_dct_init :  Read namelist.
   !!   readsec      :  Read sections description and pathway
   !!   removepoints :  Remove points which are common to 2 procs
   !!   transport    :  Compute transport for each sections
   !!   dia_dct_wri  :  Write tranports results in ascii files
   !!   interp       :  Compute temperature/salinity/density at U-point or V-point
   !!   
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE in_out_manager  ! I/O manager
   USE iom
   USE daymod          ! calendar
   USE dianam          ! build name of file
   USE lib_mpp         ! distributed memory computing library
#if defined key_si3
   USE ice
#endif
   USE domvvl
   USE timing          ! preformance summary

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dia_dct      ! routine called by step.F90
   PUBLIC   dia_dct_init ! routine called by nemogcm.F90

   !                         !!** namelist variables **
   INTEGER         ::   nn_dct        !  Frequency of computation
   INTEGER         ::   nn_secdebug   !  Number of the section to debug
   
   INTEGER, PARAMETER :: nb_class_max  = 34
   INTEGER, PARAMETER :: nb_sec_max    = 150
   INTEGER, PARAMETER :: nb_point_max  = 2000
   INTEGER, PARAMETER :: nb_type_class = 10
   INTEGER, PARAMETER :: nb_3d_vars    = 3 
   INTEGER, PARAMETER :: nb_2d_vars    = 2 

   TYPE POINT_SECTION
      INTEGER :: I,J
   END TYPE POINT_SECTION

   TYPE COORD_SECTION
      REAL(wp) :: lon,lat
   END TYPE COORD_SECTION

   TYPE SECTION
      CHARACTER(len=60)                            :: name              ! name of the sec
      LOGICAL                                      :: llstrpond         ! true if you want the computation of salt and heat transports
      LOGICAL                                      :: ll_ice_section    ! ice surface and ice volume computation
      LOGICAL                                      :: ll_date_line      ! = T if the section crosses the date-line
      TYPE(COORD_SECTION), DIMENSION(2)            :: coordSec          ! longitude and latitude of the extremities of the sec
      INTEGER                                      :: nb_class          ! number of boundaries for density classes
      INTEGER, DIMENSION(nb_point_max)             :: direction         ! vector direction of the point in the section
      CHARACTER(len=40),DIMENSION(nb_class_max)    :: classname         ! characteristics of the class
      REAL(wp), DIMENSION(nb_class_max)            :: zsigi             ! in-situ   density classes    (99 if you don't want)
      REAL(wp), DIMENSION(nb_class_max)            :: zsigp             ! potential density classes    (99 if you don't want)
      REAL(wp), DIMENSION(nb_class_max)            :: zsal              ! salinity classes   (99 if you don't want)
      REAL(wp), DIMENSION(nb_class_max)            :: ztem              ! temperature classes(99 if you don't want)
      REAL(wp), DIMENSION(nb_class_max)            :: zlay              ! level classes      (99 if you don't want)
      REAL(wp), DIMENSION(nb_type_class,nb_class_max)  :: transport     ! transport output
      REAL(wp)                                         :: slopeSection  ! slope of the section
      INTEGER                                          :: nb_point      ! number of points in the section
      TYPE(POINT_SECTION),DIMENSION(nb_point_max)      :: listPoint     ! list of points in the sections
   END TYPE SECTION

   TYPE(SECTION),DIMENSION(nb_sec_max) :: secs ! Array of sections
 
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) ::  transports_3d 
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   ::  transports_2d  


   !! * Substitutions
#  include "single_precision_substitute.h90"  
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: diadct.F90 13286 2020-07-09 15:48:29Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS
 
   INTEGER FUNCTION diadct_alloc() 
      !!---------------------------------------------------------------------- 
      !!                   ***  FUNCTION diadct_alloc  *** 
      !!---------------------------------------------------------------------- 

      ALLOCATE( transports_3d(nb_3d_vars,nb_sec_max,nb_point_max,jpk), &
         &      transports_2d(nb_2d_vars,nb_sec_max,nb_point_max)    , STAT=diadct_alloc ) 

      CALL mpp_sum( 'diadct', diadct_alloc ) 
      IF( diadct_alloc /= 0 )   CALL ctl_stop( 'STOP', 'diadct_alloc: failed to allocate arrays' ) 

   END FUNCTION diadct_alloc

   SUBROUTINE dia_dct_init
      !!---------------------------------------------------------------------
      !!               ***  ROUTINE diadct  ***  
      !!
      !!  ** Purpose: Read the namelist parameters
      !!              Open output files
      !!
      !!---------------------------------------------------------------------
      INTEGER  ::   ios                 ! Local integer output status for namelist read
      !!
      NAMELIST/nam_diadct/ln_diadct, nn_dct, nn_dctwri, nn_secdebug
      !!---------------------------------------------------------------------

     READ  ( numnam_ref, nam_diadct, IOSTAT = ios, ERR = 901)
901  IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_diadct in reference namelist' )

     READ  ( numnam_cfg, nam_diadct, IOSTAT = ios, ERR = 902 )
902  IF( ios >  0 ) CALL ctl_nam ( ios , 'nam_diadct in configuration namelist' )
     IF(lwm) WRITE ( numond, nam_diadct )

     IF( lwp ) THEN
        WRITE(numout,*) " "
        WRITE(numout,*) "diadct_init: compute transports through sections "
        WRITE(numout,*) "~~~~~~~~~~~~~~~~~~~~~"
        WRITE(numout,*) "       Calculate transport thru sections: ln_diadct = ", ln_diadct
        WRITE(numout,*) "       Frequency of computation:          nn_dct    = ", nn_dct
        WRITE(numout,*) "       Frequency of write:                nn_dctwri = ", nn_dctwri

        IF      ( nn_secdebug .GE. 1 .AND. nn_secdebug .LE. nb_sec_max )THEN
                                            WRITE(numout,*)"       Debug section number: ", nn_secdebug 
        ELSE IF ( nn_secdebug ==  0 )THEN ; WRITE(numout,*)"       No section to debug"
        ELSE IF ( nn_secdebug == -1 )THEN ; WRITE(numout,*)"       Debug all sections"
        ELSE                              ; WRITE(numout,*)"       Wrong value for nn_secdebug : ",nn_secdebug
        ENDIF
     ENDIF

     IF( ln_diadct ) THEN
        ! control
        IF(nn_dct .GE. nn_dctwri .AND. MOD(nn_dct,nn_dctwri) .NE. 0)  &
           &  CALL ctl_stop( 'diadct: nn_dct should be smaller and a multiple of nn_dctwri' )

        ! allocate dia_dct arrays
        IF( diadct_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'diadct_alloc: failed to allocate arrays' )

        !Read section_ijglobal.diadct
        CALL readsec

        !open output file
        IF( lwm ) THEN
           CALL ctl_opn( numdct_vol,  'volume_transport', 'NEW', 'FORMATTED', 'SEQUENTIAL', -1, numout,  .FALSE. )
           CALL ctl_opn( numdct_heat, 'heat_transport'  , 'NEW', 'FORMATTED', 'SEQUENTIAL', -1, numout,  .FALSE. )
           CALL ctl_opn( numdct_salt, 'salt_transport'  , 'NEW', 'FORMATTED', 'SEQUENTIAL', -1, numout,  .FALSE. )
        ENDIF
        ! Initialise arrays to zero 
        transports_3d(:,:,:,:)=0.0 
        transports_2d(:,:,:)  =0.0 
        !
     ENDIF
     !
  END SUBROUTINE dia_dct_init
 
 
  SUBROUTINE dia_dct( kt, Kmm )
     !!---------------------------------------------------------------------
     !!               ***  ROUTINE diadct  ***  
     !!
     !!  Purpose :: Compute section transports and write it in numdct files 
     !!   
     !!  Method  :: All arrays initialised to zero in dct_init 
     !!             Each nn_dct time step call subroutine 'transports' for 
     !!               each section to sum the transports over each grid cell. 
     !!             Each nn_dctwri time step: 
     !!               Divide the arrays by the number of summations to gain 
     !!               an average value 
     !!               Call dia_dct_sum to sum relevant grid boxes to obtain 
     !!               totals for each class (density, depth, temp or sal) 
     !!               Call dia_dct_wri to write the transports into file 
     !!               Reinitialise all relevant arrays to zero 
     !!---------------------------------------------------------------------
     INTEGER, INTENT(in) ::   kt    ! ocean time step
     INTEGER, INTENT(in) ::   Kmm   ! time level index
     !
     INTEGER ::   jsec              ! loop on sections
     INTEGER ::   itotal            ! nb_sec_max*nb_type_class*nb_class_max
     LOGICAL ::   lldebug =.FALSE.  ! debug a section  
     INTEGER              , DIMENSION(1)    ::   ish     ! work array for mpp_sum
     INTEGER              , DIMENSION(3)    ::   ish2    !   "
     REAL(wp), ALLOCATABLE, DIMENSION(:)    ::   zwork   !   "  
     REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)::   zsum    !   "
     !!---------------------------------------------------------------------    
     !
     IF( ln_timing )   CALL timing_start('dia_dct')

     IF( lk_mpp )THEN
        itotal = nb_sec_max*nb_type_class*nb_class_max
        ALLOCATE( zwork(itotal) , zsum(nb_sec_max,nb_type_class,nb_class_max) )
     ENDIF    
 
     ! Initialise arrays
     zwork(:) = 0.0 
     zsum(:,:,:) = 0.0

     IF( lwp .AND. kt==nit000+nn_dct-1 ) THEN
         WRITE(numout,*) " "
         WRITE(numout,*) "diadct: compute transport"
         WRITE(numout,*) "~~~~~~~~~~~~~~~~~~~~~~~~~"
         WRITE(numout,*) "nb_sec = ",nb_sec
     ENDIF

 
     ! Compute transport and write only at nn_dctwri
     IF( MOD(kt,nn_dct)==0 ) THEN 

        DO jsec=1,nb_sec

           !debug this section computing ?
           lldebug=.FALSE.
           IF( (jsec==nn_secdebug .OR. nn_secdebug==-1) .AND.  kt==nit000+nn_dct-1 ) lldebug=.TRUE. 

           !Compute transport through section  
           CALL transport(Kmm,secs(jsec),lldebug,jsec) 

        ENDDO
             
        IF( MOD(kt,nn_dctwri)==0 )THEN

           IF( kt==nit000+nn_dctwri-1 )WRITE(numout,*)"      diadct: average transports and write at kt = ",kt         
  
           !! divide arrays by nn_dctwri/nn_dct to obtain average 
           transports_3d(:,:,:,:)=transports_3d(:,:,:,:)/(nn_dctwri/nn_dct) 
           transports_2d(:,:,:)  =transports_2d(:,:,:)  /(nn_dctwri/nn_dct) 
 
           ! Sum over each class 
           DO jsec=1,nb_sec 
              CALL dia_dct_sum(Kmm,secs(jsec),jsec) 
           ENDDO 

           !Sum on all procs 
!           IF( lk_mpp )THEN 
! removed by AV do not know the significance of this test, following test taken from diaprt 
#if ! defined key_mpi_off
           IF( .NOT. l_istiled .OR. ntile == nijtile ) THEN
              ish(1)  =  nb_sec_max*nb_type_class*nb_class_max 
              ish2    = (/nb_sec_max,nb_type_class,nb_class_max/)
              DO jsec=1,nb_sec ; zsum(jsec,:,:) = secs(jsec)%transport(:,:) ; ENDDO
              zwork(:)= RESHAPE(zsum(:,:,:), ish )
              CALL mpp_sum('diadct', zwork, ish(1))
              zsum(:,:,:)= RESHAPE(zwork,ish2)
              DO jsec=1,nb_sec ; secs(jsec)%transport(:,:) = zsum(jsec,:,:) ; ENDDO
           ENDIF
#endif

           !Write the transport
           DO jsec=1,nb_sec

              IF( lwm )CALL dia_dct_wri(kt,jsec,secs(jsec))
            
              !nullify transports values after writing
              transports_3d(:,jsec,:,:)=0.
              transports_2d(:,jsec,:  )=0.
              secs(jsec)%transport(:,:)=0.  

           ENDDO

        ENDIF 

     ENDIF

     IF( lk_mpp )THEN
        itotal = nb_sec_max*nb_type_class*nb_class_max
        DEALLOCATE( zwork , zsum  )
     ENDIF    

     IF( ln_timing )   CALL timing_stop('dia_dct')
     !
  END SUBROUTINE dia_dct


  SUBROUTINE readsec 
     !!---------------------------------------------------------------------
     !!               ***  ROUTINE readsec  ***
     !!
     !!  ** Purpose:
     !!            Read a binary file(section_ijglobal.diadct) 
     !!            generated by the tools "NEMOGCM/TOOLS/SECTIONS_DIADCT"
     !!
     !!
     !!---------------------------------------------------------------------
     INTEGER :: iptglo , iptloc                               ! Global and local number of points for a section
     INTEGER :: isec, iiglo, ijglo, iiloc, ijloc,iost,i1 ,i2  ! temporary  integer
     INTEGER :: jsec, jpt                                     ! dummy loop indices
     INTEGER, DIMENSION(2) :: icoord 
     LOGICAL               :: llbon, lldebug   ! local logical
     CHARACTER(len=160)    :: clname           ! filename
     CHARACTER(len=200)    :: cltmp
     CHARACTER(len=200)    :: clformat                          !automatic format
     TYPE(POINT_SECTION),DIMENSION(nb_point_max)  ::coordtemp   !contains listpoints coordinates read in the file
     INTEGER, DIMENSION(nb_point_max) :: directemp              !contains listpoints directions read in the files
     !!-------------------------------------------------------------------------------------

     !open input file
     !---------------
     CALL ctl_opn( numdct_in, 'section_ijglobal.diadct', 'OLD', 'UNFORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
 
     !---------------
     !Read input file
     !---------------
     
     DO jsec=1,nb_sec_max      !loop on the nb_sec sections

        IF (  jsec==nn_secdebug .OR. nn_secdebug==-1  ) &
           & WRITE(numout,*)'debuging for section number: ',jsec 

        !initialization
        !---------------
        secs(jsec)%name=''
        secs(jsec)%llstrpond    = .FALSE.  ; secs(jsec)%ll_ice_section = .FALSE.
        secs(jsec)%ll_date_line = .FALSE.  ; secs(jsec)%nb_class       = 0
        secs(jsec)%zsigi        = 99._wp   ; secs(jsec)%zsigp          = 99._wp
        secs(jsec)%zsal         = 99._wp   ; secs(jsec)%ztem           = 99._wp
        secs(jsec)%zlay         = 99._wp         
        secs(jsec)%transport    =  0._wp   ; secs(jsec)%nb_point       = 0

        !read section's number / name / computing choices / classes / slopeSection / points number
        !-----------------------------------------------------------------------------------------
        READ(numdct_in,iostat=iost)isec
        IF (iost .NE. 0 )EXIT !end of file 
        WRITE(cltmp,'(a,i4.4,a,i4.4)')'diadct: read sections : Problem of section number: isec= ',isec,' and jsec= ',jsec
        IF( jsec .NE. isec )  CALL ctl_stop( cltmp )

        IF( jsec==nn_secdebug .OR. nn_secdebug==-1 )WRITE(numout,*)"isec ",isec 

        READ(numdct_in)secs(jsec)%name
        READ(numdct_in)secs(jsec)%llstrpond
        READ(numdct_in)secs(jsec)%ll_ice_section
        READ(numdct_in)secs(jsec)%ll_date_line
        READ(numdct_in)secs(jsec)%coordSec
        READ(numdct_in)secs(jsec)%nb_class
        READ(numdct_in)secs(jsec)%zsigi
        READ(numdct_in)secs(jsec)%zsigp
        READ(numdct_in)secs(jsec)%zsal
        READ(numdct_in)secs(jsec)%ztem
        READ(numdct_in)secs(jsec)%zlay
        READ(numdct_in)secs(jsec)%slopeSection
        READ(numdct_in)iptglo

        !debug
        !-----

        IF( jsec==nn_secdebug .OR. nn_secdebug==-1 )THEN
          
            WRITE(clformat,'(a,i2,a)') '(A40,', nb_class_max,'(f8.3,1X))' 

            WRITE(numout,*)       "   Section name :                       ",TRIM(secs(jsec)%name)
            WRITE(numout,*)       "      Compute heat and salt transport ? ",secs(jsec)%llstrpond
            WRITE(numout,*)       "      Compute ice transport ?           ",secs(jsec)%ll_ice_section
            WRITE(numout,*)       "      Section crosses date-line ?       ",secs(jsec)%ll_date_line
            WRITE(numout,*)       "      Slope section :                   ",secs(jsec)%slopeSection
            WRITE(numout,*)       "      Number of points in the section:  ",iptglo
            WRITE(numout,*)       "      Number of classes                 ",secs(jsec)%nb_class
            WRITE(numout,clformat)"      Insitu density classes :          ",secs(jsec)%zsigi
            WRITE(numout,clformat)"      Potential density classes :       ",secs(jsec)%zsigp
            WRITE(numout,clformat)"      Salinity classes :                ",secs(jsec)%zsal
            WRITE(numout,clformat)"      Temperature classes :             ",secs(jsec)%ztem
            WRITE(numout,clformat)"      Depth classes :                   ",secs(jsec)%zlay
        ENDIF               

        IF( iptglo /= 0 )THEN
             
           !read points'coordinates and directions 
           !--------------------------------------
           coordtemp(:) = POINT_SECTION(0,0) !list of points read
           directemp(:) = 0                  !value of directions of each points
           DO jpt=1,iptglo
              READ(numdct_in) i1, i2
              coordtemp(jpt)%I = i1 
              coordtemp(jpt)%J = i2
           ENDDO
           READ(numdct_in) directemp(1:iptglo)
    
           !debug
           !-----
           IF( jsec==nn_secdebug .OR. nn_secdebug==-1 )THEN
              WRITE(numout,*)"      List of points in global domain:"
              DO jpt=1,iptglo
                 WRITE(numout,*)'        # I J ',jpt,coordtemp(jpt),directemp(jpt)
              ENDDO                  
           ENDIF
 
           !Now each proc selects only points that are in its domain:
           !--------------------------------------------------------
           iptloc = 0                    ! initialize number of points selected
           DO jpt = 1, iptglo            ! loop on listpoint read in the file
              !      
              iiglo=coordtemp(jpt)%I          ! global coordinates of the point
              ijglo=coordtemp(jpt)%J          !  " 

              IF( iiglo==jpiglo .AND. nimpp==1 )   iiglo = 2         !!gm BUG: Hard coded periodicity !

              iiloc=iiglo-nimpp+1   ! local coordinates of the point
              ijloc=ijglo-njmpp+1   !  "

              !verify if the point is on the local domain:(1,Nie0)*(1,Nje0)
              IF( iiloc >= 1 .AND. iiloc <= Nie0 .AND. &
                  ijloc >= 1 .AND. ijloc <= Nje0       )THEN
                 iptloc = iptloc + 1                                                 ! count local points
                 secs(jsec)%listPoint(iptloc) = POINT_SECTION(mi0(iiglo),mj0(ijglo)) ! store local coordinates
                 secs(jsec)%direction(iptloc) = directemp(jpt)                       ! store local direction
              ENDIF
              !
           END DO
     
           secs(jsec)%nb_point=iptloc !store number of section's points

           !debug
           !-----
           IF( jsec==nn_secdebug .OR. nn_secdebug==-1 )THEN
              WRITE(numout,*)"      List of points selected by the proc:"
              DO jpt = 1,iptloc
                 iiglo = secs(jsec)%listPoint(jpt)%I + nimpp - 1
                 ijglo = secs(jsec)%listPoint(jpt)%J + njmpp - 1
                 WRITE(numout,*)'         # I J : ',iiglo,ijglo
              ENDDO
           ENDIF

              IF(jsec==nn_secdebug .AND. secs(jsec)%nb_point .NE. 0)THEN
              DO jpt = 1,iptloc
                 iiglo = secs(jsec)%listPoint(jpt)%I + nimpp - 1
                 ijglo = secs(jsec)%listPoint(jpt)%J + njmpp - 1
              ENDDO
              ENDIF

           !remove redundant points between processors
           !------------------------------------------
           lldebug = .FALSE. ; IF ( jsec==nn_secdebug .OR. nn_secdebug==-1 ) lldebug = .TRUE.
           IF( iptloc .NE. 0 )THEN
              CALL removepoints(secs(jsec),'I','top_list',lldebug)
              CALL removepoints(secs(jsec),'I','bot_list',lldebug)
              CALL removepoints(secs(jsec),'J','top_list',lldebug)
              CALL removepoints(secs(jsec),'J','bot_list',lldebug)
           ENDIF
           IF(jsec==nn_secdebug .AND. secs(jsec)%nb_point .NE. 0)THEN
              DO jpt = 1,secs(jsec)%nb_point
                 iiglo = secs(jsec)%listPoint(jpt)%I + nimpp - 1
                 ijglo = secs(jsec)%listPoint(jpt)%J + njmpp - 1
              ENDDO
           ENDIF

           !debug
           !-----
           IF( jsec==nn_secdebug .OR. nn_secdebug==-1 )THEN
              WRITE(numout,*)"      List of points after removepoints:"
              iptloc = secs(jsec)%nb_point
              DO jpt = 1,iptloc
                 iiglo = secs(jsec)%listPoint(jpt)%I + nimpp - 1
                 ijglo = secs(jsec)%listPoint(jpt)%J + njmpp - 1
                 WRITE(numout,*)'         # I J : ',iiglo,ijglo
                 CALL FLUSH(numout)
              ENDDO
           ENDIF

        ELSE  ! iptglo = 0
           IF( jsec==nn_secdebug .OR. nn_secdebug==-1 )&
              WRITE(numout,*)'   No points for this section.'
        ENDIF

     ENDDO !end of the loop on jsec
 
     nb_sec = jsec-1   !number of section read in the file
     !
  END SUBROUTINE readsec


  SUBROUTINE removepoints(sec,cdind,cdextr,ld_debug)
     !!---------------------------------------------------------------------------
     !!             *** function removepoints
     !!
     !!   ** Purpose :: Remove points which are common to 2 procs
     !!
     !----------------------------------------------------------------------------
     !! * arguments
     TYPE(SECTION),INTENT(INOUT) :: sec
     CHARACTER(len=1),INTENT(IN) :: cdind   ! = 'I'/'J'
     CHARACTER(len=8),INTENT(IN) :: cdextr  ! = 'top_list'/'bot_list'
     LOGICAL,INTENT(IN)          :: ld_debug                     

     !! * Local variables
     INTEGER :: iextr         ,& !extremity of listpoint that we verify
                iind          ,& !coord     of listpoint that we verify
                itest         ,& !indice value of the side of the domain 
                                 !where points could be redundant
                isgn          ,& ! isgn= 1 : scan listpoint from start to end
                                 ! isgn=-1 : scan listpoint from end to start 
                istart,iend      !first and last points selected in listpoint
     INTEGER :: jpoint           !loop on list points
     INTEGER, DIMENSION(nb_point_max)   :: idirec !contains temporary sec%direction
     INTEGER, DIMENSION(2,nb_point_max) :: icoord !contains temporary sec%listpoint
     !----------------------------------------------------------------------------
     !
     IF( ld_debug )WRITE(numout,*)'      -------------------------'
     IF( ld_debug )WRITE(numout,*)'      removepoints in listpoint'

     !iextr=extremity of list_point that we verify
     IF      ( cdextr=='bot_list' )THEN ; iextr=1            ; isgn=1
     ELSE IF ( cdextr=='top_list' )THEN ; iextr=sec%nb_point ; isgn=-1
     ELSE    ; CALL ctl_stop("removepoints :Wrong value for cdextr")
     ENDIF
 
     !which coordinate shall we verify ?
     IF      ( cdind=='I' )THEN   ; itest=Nie0 ; iind=1
     ELSE IF ( cdind=='J' )THEN   ; itest=Nje0 ; iind=2
     ELSE    ; CALL ctl_stop("removepoints :Wrong value for cdind") 
     ENDIF

     IF( ld_debug )THEN
        WRITE(numout,*)'      case: coord/list extr/domain side'
        WRITE(numout,*)'      ', cdind,' ',cdextr,' ',itest
        WRITE(numout,*)'      Actual number of points: ',sec%nb_point
     ENDIF

     icoord(1,1:nb_point_max) = sec%listPoint%I
     icoord(2,1:nb_point_max) = sec%listPoint%J
     idirec                   = sec%direction
     sec%listPoint            = POINT_SECTION(0,0)
     sec%direction            = 0

     jpoint=iextr+isgn
     DO WHILE( jpoint .GE. 1 .AND. jpoint .LE. sec%nb_point )
         IF( icoord( iind,jpoint-isgn ) == itest .AND. icoord( iind,jpoint ) == itest )THEN ; jpoint=jpoint+isgn
         ELSE                                                                               ; EXIT
         ENDIF
     ENDDO 

     IF( cdextr=='bot_list')THEN ; istart=jpoint-1 ; iend=sec%nb_point
     ELSE                        ; istart=1        ; iend=jpoint+1
     ENDIF

     sec%listPoint(1:1+iend-istart)%I = icoord(1,istart:iend)
     sec%listPoint(1:1+iend-istart)%J = icoord(2,istart:iend)
     sec%direction(1:1+iend-istart)   = idirec(istart:iend)
     sec%nb_point                     = iend-istart+1
     
     IF( ld_debug )THEN
        WRITE(numout,*)'      Number of points after removepoints :',sec%nb_point
        WRITE(numout,*)'      sec%direction after removepoints :',sec%direction(1:sec%nb_point)
     ENDIF
      !
   END SUBROUTINE removepoints


   SUBROUTINE transport(Kmm,sec,ld_debug,jsec)
     !!-------------------------------------------------------------------------------------------
     !!                     ***  ROUTINE transport  ***
     !!
     !!  Purpose ::  Compute the transport for each point in a section 
     !! 
     !!  Method  ::  Loop over each segment, and each vertical level and add the transport 
     !!              Be aware :           
     !!              One section is a sum of segments 
     !!              One segment is defined by 2 consecutive points in sec%listPoint 
     !!              All points of sec%listPoint are positioned on the F-point of the cell 
     !! 
     !!              There are two loops:                  
     !!              loop on the segment between 2 nodes 
     !!              loop on the level jk !!
     !! 
     !!  Output  ::  Arrays containing the volume,density,heat,salt transports for each i
     !!              point in a section, summed over each nn_dct. 
     !!
     !!-------------------------------------------------------------------------------------------
     INTEGER      ,INTENT(IN)    :: Kmm         ! time level index
     TYPE(SECTION),INTENT(INOUT) :: sec
     LOGICAL      ,INTENT(IN)    :: ld_debug
     INTEGER      ,INTENT(IN)    :: jsec        ! numeric identifier of section
     !
     INTEGER ::   jk, jseg, jclass,jl, isgnu, isgnv    ! loop on level/segment/classes/ice categories
     REAL(wp)::   zumid, zvmid, zumid_ice, zvmid_ice   ! U/V ocean & ice velocity on a cell segment 
     REAL(wp)::   zTnorm                               ! transport of velocity through one cell's sides 
     REAL(wp)::   ztn, zsn, zrhoi, zrhop, zsshn, zdep  ! temperature/salinity/potential density/ssh/depth at u/v point
     TYPE(POINT_SECTION) ::   k
      !!--------------------------------------------------------
      !
      IF( ld_debug )WRITE(numout,*)'      Compute transport'

      !---------------------------!
      !  COMPUTE TRANSPORT        !
      !---------------------------!
      IF(sec%nb_point .NE. 0)THEN   

         !----------------------------------------------------------------------------------------------------
         !Compute sign for velocities:
         !
         !convention:
         !   non horizontal section: direction + is toward left hand of section
         !       horizontal section: direction + is toward north of section
         !
         !
         !       slopeSection < 0     slopeSection > 0       slopeSection=inf            slopeSection=0
         !       ----------------      -----------------     ---------------             --------------
         !
         !   isgnv=1         direction +      
         !  ______         _____             ______                                                   
         !        |           //|            |                  |                         direction +   
         !        | isgnu=1  // |            |isgnu=1           |isgnu=1                     /|\
         !        |_______  //         ______|    \\            | ---\                        |
         !               |             | isgnv=-1  \\ |         | ---/ direction +       ____________
         !               |             |          __\\|         |                    
         !               |             |     direction +        |                      isgnv=1                                 
         !                                                      
         !----------------------------------------------------------------------------------------------------
         isgnu = 1
         IF( sec%slopeSection .GT. 0 ) THEN  ; isgnv = -1 
         ELSE                                ; isgnv =  1
         ENDIF
         IF( sec%slopeSection .GE. 9999. )     isgnv =  1

         IF( ld_debug )write(numout,*)"sec%slopeSection isgnu isgnv ",sec%slopeSection,isgnu,isgnv

         !--------------------------------------!
         ! LOOP ON THE SEGMENT BETWEEN 2 NODES  !
         !--------------------------------------!
         DO jseg=1,MAX(sec%nb_point-1,0)
              
            !-------------------------------------------------------------------------------------------
            ! Select the appropriate coordinate for computing the velocity of the segment
            !
            !                      CASE(0)                                    Case (2)
            !                      -------                                    --------
            !  listPoint(jseg)                 listPoint(jseg+1)       listPoint(jseg)  F(i,j)      
            !      F(i,j)----------V(i+1,j)-------F(i+1,j)                               |
            !                                                                            |
            !                                                                            |
            !                                                                            |
            !                      Case (3)                                            U(i,j)
            !                      --------                                              |
            !                                                                            |
            !  listPoint(jseg+1) F(i,j+1)                                                |
            !                        |                                                   |
            !                        |                                                   |
            !                        |                                 listPoint(jseg+1) F(i,j-1)
            !                        |                                            
            !                        |                                            
            !                     U(i,j+1)                                            
            !                        |                                       Case(1)     
            !                        |                                       ------      
            !                        |                                            
            !                        |                 listPoint(jseg+1)             listPoint(jseg)                           
            !                        |                 F(i-1,j)-----------V(i,j) -------f(jseg)                           
            ! listPoint(jseg)     F(i,j)
            ! 
            !-------------------------------------------------------------------------------------------

            SELECT CASE( sec%direction(jseg) )
            CASE(0)   ;    k = sec%listPoint(jseg)
            CASE(1)   ;    k = POINT_SECTION(sec%listPoint(jseg)%I+1,sec%listPoint(jseg)%J)
            CASE(2)   ;    k = sec%listPoint(jseg)
            CASE(3)   ;    k = POINT_SECTION(sec%listPoint(jseg)%I,sec%listPoint(jseg)%J+1)
            END SELECT

            !---------------------------| 
            !     LOOP ON THE LEVEL     | 
            !---------------------------| 
            DO jk = 1, mbkt(k%I,k%J)            !Sum of the transport on the vertical
            !           ! compute temperature, salinity, insitu & potential density, ssh and depth at U/V point 
            SELECT CASE( sec%direction(jseg) )
               CASE(0,1) 
                  ztn   = interp(Kmm,k%I,k%J,jk,'V',ts(:,:,:,jp_tem,Kmm) ) 
                  zsn   = interp(Kmm,k%I,k%J,jk,'V',ts(:,:,:,jp_sal,Kmm) ) 
                  zrhop = interp(Kmm,k%I,k%J,jk,'V',rhop) 
                  zrhoi =interp(Kmm,k%I,k%J,jk,'V',CASTDP(rhd*rho0+rho0))
                  zsshn =  0.5*( ssh(k%I,k%J,Kmm) + ssh(k%I,k%J+1,Kmm)    ) * vmask(k%I,k%J,1) 
               CASE(2,3) 
                  ztn   = interp(Kmm,k%I,k%J,jk,'U',ts(:,:,:,jp_tem,Kmm) ) 
                  zsn   = interp(Kmm,k%I,k%J,jk,'U',ts(:,:,:,jp_sal,Kmm) ) 
                  zrhop = interp(Kmm,k%I,k%J,jk,'U',rhop) 
                  zrhoi =interp(Kmm,k%I,k%J,jk,'U',CASTDP(rhd*rho0+rho0))
                  zsshn =  0.5*( ssh(k%I,k%J,Kmm) + ssh(k%I+1,k%J,Kmm)    ) * umask(k%I,k%J,1)  
               END SELECT 
               !
               zdep= gdept(k%I,k%J,jk,Kmm) 
  
               SELECT CASE( sec%direction(jseg) )                !compute velocity with the correct direction 
               CASE(0,1)   
                  zumid=0._wp
                  zvmid=isgnv*vv(k%I,k%J,jk,Kmm)*vmask(k%I,k%J,jk) 
               CASE(2,3) 
                  zumid=isgnu*uu(k%I,k%J,jk,Kmm)*umask(k%I,k%J,jk) 
                  zvmid=0._wp
               END SELECT 
 
               !zTnorm=transport through one cell; 
               !velocity* cell's length * cell's thickness 
               zTnorm = zumid*e2u(k%I,k%J) * e3u(k%I,k%J,jk,Kmm)     & 
                  &   + zvmid*e1v(k%I,k%J) * e3v(k%I,k%J,jk,Kmm) 

!!gm  THIS is WRONG  no transport due to ssh in linear free surface case !!!!!
               IF( ln_linssh ) THEN              !add transport due to free surface 
                  IF( jk==1 ) THEN 
                     zTnorm = zTnorm + zumid* e2u(k%I,k%J) * zsshn * umask(k%I,k%J,jk)   & 
                        &            + zvmid* e1v(k%I,k%J) * zsshn * vmask(k%I,k%J,jk) 
                  ENDIF 
               ENDIF
!!gm end
              !COMPUTE TRANSPORT  
 
              transports_3d(1,jsec,jseg,jk) = transports_3d(1,jsec,jseg,jk) + zTnorm 
  
              IF( sec%llstrpond ) THEN 
                 transports_3d(2,jsec,jseg,jk) = transports_3d(2,jsec,jseg,jk)  + zTnorm * ztn * zrhop * rcp
                 transports_3d(3,jsec,jseg,jk) = transports_3d(3,jsec,jseg,jk)  + zTnorm * zsn * zrhop * 0.001
              ENDIF
   
           END DO !end of loop on the level

#if defined key_si3

           !ICE CASE    
           !------------
           IF( sec%ll_ice_section )THEN
              SELECT CASE (sec%direction(jseg))
              CASE(0)
                 zumid_ice = 0
                 zvmid_ice =  isgnv*0.5*(v_ice(k%I,k%J+1)+v_ice(k%I+1,k%J+1))
              CASE(1)
                 zumid_ice = 0
                 zvmid_ice =  isgnv*0.5*(v_ice(k%I,k%J+1)+v_ice(k%I+1,k%J+1))
              CASE(2)
                 zvmid_ice = 0
                 zumid_ice =  isgnu*0.5*(u_ice(k%I+1,k%J)+u_ice(k%I+1,k%J+1))
              CASE(3)
                 zvmid_ice = 0
                 zumid_ice =  isgnu*0.5*(u_ice(k%I+1,k%J)+u_ice(k%I+1,k%J+1))
              END SELECT
   
              zTnorm=zumid_ice*e2u(k%I,k%J)+zvmid_ice*e1v(k%I,k%J)

#if defined key_si3
              DO jl=1,jpl
                 transports_2d(1,jsec,jseg) = transports_2d(1,jsec,jseg) + (zTnorm)*       &
                                    a_i(sec%listPoint(jseg)%I,sec%listPoint(jseg)%J,jl) *  &
                                  ( h_i(sec%listPoint(jseg)%I,sec%listPoint(jseg)%J,jl) +  &
                                    h_s(sec%listPoint(jseg)%I,sec%listPoint(jseg)%J,jl) )
                                   
                 transports_2d(2,jsec,jseg) = transports_2d(2,jsec,jseg) + (zTnorm)*   &
                                    a_i(sec%listPoint(jseg)%I,sec%listPoint(jseg)%J,jl)
              END DO
#endif
   
           ENDIF !end of ice case
#endif
 
        END DO !end of loop on the segment

     ENDIF !end of sec%nb_point =0 case
     !
  END SUBROUTINE transport


  SUBROUTINE dia_dct_sum(Kmm,sec,jsec) 
     !!------------------------------------------------------------- 
     !! Purpose: Average the transport over nn_dctwri time steps  
     !! and sum over the density/salinity/temperature/depth classes 
     !! 
     !! Method:   Sum over relevant grid cells to obtain values  
     !!           for each class
     !!              There are several loops:                  
     !!              loop on the segment between 2 nodes 
     !!              loop on the level jk 
     !!              loop on the density/temperature/salinity/level classes 
     !!              test on the density/temperature/salinity/level 
     !! 
     !!  Note:    Transport through a given section is equal to the sum of transports 
     !!           computed on each proc. 
     !!           On each proc,transport is equal to the sum of transport computed through 
     !!           segments linking each point of sec%listPoint  with the next one.    
     !! 
     !!------------------------------------------------------------- 
     INTEGER      ,INTENT(IN)    :: Kmm         ! time level index
     TYPE(SECTION),INTENT(INOUT) :: sec 
     INTEGER      ,INTENT(IN)    :: jsec        ! numeric identifier of section 
 
     TYPE(POINT_SECTION) :: k 
     INTEGER  :: jk,jseg,jclass                        ! dummy variables for looping on level/segment/classes  
     REAL(wp) :: ztn, zsn, zrhoi, zrhop, zsshn, zdep ! temperature/salinity/ssh/potential density /depth at u/v point 
     !!------------------------------------------------------------- 
 
     !! Sum the relevant segments to obtain values for each class 
     IF(sec%nb_point .NE. 0)THEN    
 
        !--------------------------------------! 
        ! LOOP ON THE SEGMENT BETWEEN 2 NODES  ! 
        !--------------------------------------! 
        DO jseg=1,MAX(sec%nb_point-1,0) 
            
           !------------------------------------------------------------------------------------------- 
           ! Select the appropriate coordinate for computing the velocity of the segment 
           ! 
           !                      CASE(0)                                    Case (2) 
           !                      -------                                    -------- 
           !  listPoint(jseg)                 listPoint(jseg+1)       listPoint(jseg)  F(i,j)       
           !      F(i,j)----------V(i+1,j)-------F(i+1,j)                               | 
           !                                                                            | 
           !                                                                            | 
           !                                                                            | 
           !                      Case (3)                                            U(i,j) 
           !                      --------                                              | 
           !                                                                            | 
           !  listPoint(jseg+1) F(i,j+1)                                                | 
           !                        |                                                   | 
           !                        |                                                   | 
           !                        |                                 listPoint(jseg+1) F(i,j-1) 
           !                        |                                             
           !                        |                                             
           !                     U(i,j+1)                                             
           !                        |                                       Case(1)      
           !                        |                                       ------       
           !                        |                                             
           !                        |                 listPoint(jseg+1)             listPoint(jseg)                            
           !                        |                 F(i-1,j)-----------V(i,j) -------f(jseg)                            
           ! listPoint(jseg)     F(i,j) 
           !  
           !------------------------------------------------------------------------------------------- 
 
           SELECT CASE( sec%direction(jseg) ) 
           CASE(0)  ;   k = sec%listPoint(jseg) 
           CASE(1)  ;   k = POINT_SECTION(sec%listPoint(jseg)%I+1,sec%listPoint(jseg)%J) 
           CASE(2)  ;   k = sec%listPoint(jseg) 
           CASE(3)  ;   k = POINT_SECTION(sec%listPoint(jseg)%I,sec%listPoint(jseg)%J+1) 
           END SELECT 
 
           !---------------------------| 
           !     LOOP ON THE LEVEL     | 
           !---------------------------| 
           !Sum of the transport on the vertical  
           DO jk=1,mbkt(k%I,k%J) 
 
              ! compute temperature, salinity, insitu & potential density, ssh and depth at U/V point 
              SELECT CASE( sec%direction(jseg) ) 
              CASE(0,1) 
                 ztn   = interp(Kmm,k%I,k%J,jk,'V',ts(:,:,:,jp_tem,Kmm) ) 
                 zsn   = interp(Kmm,k%I,k%J,jk,'V',ts(:,:,:,jp_sal,Kmm) ) 
                 zrhop = interp(Kmm,k%I,k%J,jk,'V',rhop) 
                 zrhoi =interp(Kmm,k%I,k%J,jk,'V',CASTDP(rhd*rho0+rho0))

              CASE(2,3) 
                 ztn   = interp(Kmm,k%I,k%J,jk,'U',ts(:,:,:,jp_tem,Kmm) ) 
                 zsn   = interp(Kmm,k%I,k%J,jk,'U',ts(:,:,:,jp_sal,Kmm) ) 
                 zrhop = interp(Kmm,k%I,k%J,jk,'U',rhop) 
                 zrhoi =interp(Kmm,k%I,k%J,jk,'U',CASTDP(rhd*rho0+rho0))
                 zsshn =  0.5*( ssh(k%I,k%J,Kmm)    + ssh(k%I+1,k%J,Kmm)    ) * umask(k%I,k%J,1)  
              END SELECT 
 
              zdep= gdept(k%I,k%J,jk,Kmm) 
  
              !------------------------------- 
              !  LOOP ON THE DENSITY CLASSES | 
              !------------------------------- 
              !The computation is made for each density/temperature/salinity/depth class 
              DO jclass=1,MAX(1,sec%nb_class-1) 
 
                 !----------------------------------------------! 
                 !TEST ON THE DENSITY/SALINITY/TEMPERATURE/LEVEL!  
                 !----------------------------------------------! 

                 IF ( (                                                    & 
                    ((( zrhop .GE. (sec%zsigp(jclass)+1000.  )) .AND.      & 
                    (   zrhop .LE. (sec%zsigp(jclass+1)+1000. ))) .OR.     & 
                    ( sec%zsigp(jclass) .EQ. 99.)) .AND.                   & 
 
                    ((( zrhoi .GE. (sec%zsigi(jclass) + 1000.  )) .AND.    & 
                    (   zrhoi .LE. (sec%zsigi(jclass+1)+1000. ))) .OR.     & 
                    ( sec%zsigi(jclass) .EQ. 99.)) .AND.                   & 
 
                    ((( zsn .GT. sec%zsal(jclass)) .AND.                   & 
                    (   zsn .LE. sec%zsal(jclass+1))) .OR.                 & 
                    ( sec%zsal(jclass) .EQ. 99.)) .AND.                    & 
 
                    ((( ztn .GE. sec%ztem(jclass)) .AND.                   & 
                    (   ztn .LE. sec%ztem(jclass+1))) .OR.                 & 
                    ( sec%ztem(jclass) .EQ.99.)) .AND.                     & 
 
                    ((( zdep .GE. sec%zlay(jclass)) .AND.                & 
                    (   zdep .LE. sec%zlay(jclass+1))) .OR.              & 
                    ( sec%zlay(jclass) .EQ. 99. ))                         & 
                                                                   ))   THEN 
 
                    !SUM THE TRANSPORTS FOR EACH CLASSES FOR THE POSITIVE AND NEGATIVE DIRECTIONS 
                    !---------------------------------------------------------------------------- 
                    IF (transports_3d(1,jsec,jseg,jk) .GE. 0.0) THEN  
                       sec%transport(1,jclass) = sec%transport(1,jclass)+transports_3d(1,jsec,jseg,jk)*1.E-6 
                    ELSE 
                       sec%transport(2,jclass) = sec%transport(2,jclass)+transports_3d(1,jsec,jseg,jk)*1.E-6 
                    ENDIF 
                    IF( sec%llstrpond )THEN 
 
                       IF ( transports_3d(2,jsec,jseg,jk) .GE. 0.0 ) THEN 
                          sec%transport(3,jclass) = sec%transport(3,jclass)+transports_3d(2,jsec,jseg,jk) 
                       ELSE 
                          sec%transport(4,jclass) = sec%transport(4,jclass)+transports_3d(2,jsec,jseg,jk) 
                       ENDIF 
 
                       IF ( transports_3d(3,jsec,jseg,jk) .GE. 0.0 ) THEN 
                          sec%transport(5,jclass) = sec%transport(5,jclass)+transports_3d(3,jsec,jseg,jk) 
                       ELSE 
                          sec%transport(6,jclass) = sec%transport(6,jclass)+transports_3d(3,jsec,jseg,jk) 
                       ENDIF 
 
                    ELSE 
                       sec%transport( 3,jclass) = 0._wp 
                       sec%transport( 4,jclass) = 0._wp 
                       sec%transport( 5,jclass) = 0._wp 
                       sec%transport( 6,jclass) = 0._wp 
                    ENDIF 
 
                 ENDIF ! end of test if point is in class 
    
              END DO ! end of loop on the classes 
 
           END DO ! loop over jk 
 
#if defined key_si3 
 
           !ICE CASE     
           IF( sec%ll_ice_section )THEN 
 
              IF ( transports_2d(1,jsec,jseg) .GE. 0.0 ) THEN 
                 sec%transport( 7,1) = sec%transport( 7,1)+transports_2d(1,jsec,jseg)*1.E-6 
              ELSE 
                 sec%transport( 8,1) = sec%transport( 8,1)+transports_2d(1,jsec,jseg)*1.E-6 
              ENDIF 
 
              IF ( transports_2d(3,jsec,jseg) .GE. 0.0 ) THEN 
                 sec%transport( 9,1) = sec%transport( 9,1)+transports_2d(2,jsec,jseg)*1.E-6 
              ELSE 
                 sec%transport(10,1) = sec%transport(10,1)+transports_2d(2,jsec,jseg)*1.E-6 
              ENDIF 
 
           ENDIF !end of ice case 
#endif 
  
        END DO !end of loop on the segment 
 
     ELSE  !if sec%nb_point =0 
        sec%transport(1:2,:)=0. 
        IF (sec%llstrpond) sec%transport(3:6,:)=0. 
        IF (sec%ll_ice_section) sec%transport(7:10,:)=0. 
     ENDIF !end of sec%nb_point =0 case 
 
  END SUBROUTINE dia_dct_sum 


  SUBROUTINE dia_dct_wri(kt,ksec,sec)
     !!-------------------------------------------------------------
     !! Write transport output in numdct 
     !! 
     !! Purpose: Write  transports in ascii files
     !! 
     !! Method:
     !!        1. Write volume transports in "volume_transport"
     !!           Unit: Sv : area * Velocity / 1.e6 
     !! 
     !!        2. Write heat transports in "heat_transport"
     !!           Unit: Peta W : area * Velocity * T * rhop * Cp * 1.e-15
     !! 
     !!        3. Write salt transports in "salt_transport"
     !!           Unit: 10^9 Kg/m^2/s : area * Velocity * S * rhop * 1.e-9 
     !!
     !!------------------------------------------------------------- 
     !!arguments
     INTEGER, INTENT(IN)          :: kt          ! time-step
     TYPE(SECTION), INTENT(INOUT) :: sec         ! section to write   
     INTEGER ,INTENT(IN)          :: ksec        ! section number

     !!local declarations
     INTEGER               :: jclass             ! Dummy loop
     CHARACTER(len=2)      :: classe             ! Classname 
     REAL(wp)              :: zbnd1,zbnd2        ! Class bounds
     REAL(wp)              :: zslope             ! section's slope coeff
     !
     REAL(wp), DIMENSION(nb_type_class)::   zsumclasses   ! 1D workspace 
     !!------------------------------------------------------------- 

     zsumclasses(:)=0._wp
     zslope = sec%slopeSection       

 
     DO jclass=1,MAX(1,sec%nb_class-1)

        classe   = 'N       '
        zbnd1   = 0._wp
        zbnd2   = 0._wp
        zsumclasses(1:nb_type_class)=zsumclasses(1:nb_type_class)+sec%transport(1:nb_type_class,jclass)

   
        !insitu density classes transports
        IF( ( sec%zsigi(jclass)   .NE. 99._wp ) .AND. &
            ( sec%zsigi(jclass+1) .NE. 99._wp )       )THEN
           classe = 'DI       '
           zbnd1 = sec%zsigi(jclass)
           zbnd2 = sec%zsigi(jclass+1)
        ENDIF
        !potential density classes transports
        IF( ( sec%zsigp(jclass)   .NE. 99._wp ) .AND. &
            ( sec%zsigp(jclass+1) .NE. 99._wp )       )THEN
           classe = 'DP      '
           zbnd1 = sec%zsigp(jclass)
           zbnd2 = sec%zsigp(jclass+1)
        ENDIF
        !depth classes transports
        IF( ( sec%zlay(jclass)    .NE. 99._wp ) .AND. &
            ( sec%zlay(jclass+1)  .NE. 99._wp )       )THEN 
           classe = 'Z       '
           zbnd1 = sec%zlay(jclass)
           zbnd2 = sec%zlay(jclass+1)
        ENDIF
        !salinity classes transports
        IF( ( sec%zsal(jclass) .NE. 99._wp    ) .AND. &
            ( sec%zsal(jclass+1) .NE. 99._wp  )       )THEN
           classe = 'S       '
           zbnd1 = sec%zsal(jclass)
           zbnd2 = sec%zsal(jclass+1)   
        ENDIF
        !temperature classes transports
        IF( ( sec%ztem(jclass) .NE. 99._wp     ) .AND. &
            ( sec%ztem(jclass+1) .NE. 99._wp     )       ) THEN
           classe = 'T       '
           zbnd1 = sec%ztem(jclass)
           zbnd2 = sec%ztem(jclass+1)
        ENDIF
                  
        !write volume transport per class
        WRITE(numdct_vol,118) ndastp,kt,ksec,sec%name,zslope, &
                              jclass,classe,zbnd1,zbnd2,&
                              sec%transport(1,jclass),sec%transport(2,jclass), &
                              sec%transport(1,jclass)+sec%transport(2,jclass)

        IF( sec%llstrpond )THEN

           !write heat transport per class:
           WRITE(numdct_heat,119) ndastp,kt,ksec,sec%name,zslope,  &
                              jclass,classe,zbnd1,zbnd2,&
                              sec%transport(3,jclass)*1.e-15,sec%transport(4,jclass)*1.e-15, &
                              ( sec%transport(3,jclass)+sec%transport(4,jclass) )*1.e-15
           !write salt transport per class
           WRITE(numdct_salt,119) ndastp,kt,ksec,sec%name,zslope,  &
                              jclass,classe,zbnd1,zbnd2,&
                              sec%transport(5,jclass)*1.e-9,sec%transport(6,jclass)*1.e-9,&
                              (sec%transport(5,jclass)+sec%transport(6,jclass))*1.e-9
        ENDIF
     ENDDO

     zbnd1 = 0._wp
     zbnd2 = 0._wp
     jclass=0

     !write total volume transport
     WRITE(numdct_vol,118) ndastp,kt,ksec,sec%name,zslope, &
                           jclass,"total",zbnd1,zbnd2,&
                           zsumclasses(1),zsumclasses(2),zsumclasses(1)+zsumclasses(2)

     IF( sec%llstrpond )THEN

        !write total heat transport
        WRITE(numdct_heat,119) ndastp,kt,ksec,sec%name,zslope, &
                           jclass,"total",zbnd1,zbnd2,&
                           zsumclasses(3)*1.e-15,zsumclasses(4)*1.e-15,&
                           (zsumclasses(3)+zsumclasses(4) )*1.e-15
        !write total salt transport
        WRITE(numdct_salt,119) ndastp,kt,ksec,sec%name,zslope, &
                           jclass,"total",zbnd1,zbnd2,&
                           zsumclasses(5)*1.e-9,zsumclasses(6)*1.e-9,&
                           (zsumclasses(5)+zsumclasses(6))*1.e-9
     ENDIF

      
     IF ( sec%ll_ice_section) THEN
        !write total ice volume transport
        WRITE(numdct_vol,118) ndastp,kt,ksec,sec%name,zslope,&
                              jclass,"ice_vol",zbnd1,zbnd2,&
                              sec%transport(7,1),sec%transport(8,1),&
                              sec%transport(7,1)+sec%transport(8,1)
        !write total ice surface transport
        WRITE(numdct_vol,118) ndastp,kt,ksec,sec%name,zslope,&
                              jclass,"ice_surf",zbnd1,zbnd2,&
                              sec%transport(9,1),sec%transport(10,1), &
                              sec%transport(9,1)+sec%transport(10,1) 
      ENDIF
                                              
118   FORMAT(I8,1X,I8,1X,I4,1X,A30,1X,f9.2,1X,I4,3X,A8,1X,2F12.4,5X,3F12.4)
119   FORMAT(I8,1X,I8,1X,I4,1X,A30,1X,f9.2,1X,I4,3X,A8,1X,2F12.4,5X,3E15.6)
   END SUBROUTINE dia_dct_wri


   FUNCTION interp(Kmm, ki, kj, kk, cd_point, ptab)
  !!----------------------------------------------------------------------
  !!
  !!   Purpose: compute temperature/salinity/density at U-point or V-point
  !!   --------
  !!
  !!   Method:
  !!   ------
  !!
  !!   ====> full step and partial step
  !! 
  !!
  !!    |    I          |    I+1           |    Z=temperature/salinity/density at U-poinT
  !!    |               |                  |
  !!  ----------------------------------------  1. Veritcal interpolation: compute zbis
  !!    |               |                  |       interpolation between ptab(I,J,K) and ptab(I,J,K+1)
  !!    |               |                  |       zbis = 
  !!    |               |                  |      [ e3w_n(I+1,J,K,NOW)*ptab(I,J,K) + ( e3w_n(I,J,K,NOW) - e3w_n(I+1,J,K,NOW) ) * ptab(I,J,K-1) ]
  !!    |               |                  |     /[ e3w_n(I+1,J,K,NOW)             +   e3w_n(I,J,K,NOW) - e3w_n(I+1,J,K,NOW) ] 
  !!    |               |                  | 
  !!    |               |                  |    2. Horizontal interpolation: compute value at U/V point
  !!K-1 | ptab(I,J,K-1) |                  |       interpolation between zbis and ptab(I+1,J,K)  
  !!    |     .         |                  |
  !!    |     .         |                  |       interp = ( 0.5*zet2*zbis + 0.5*zet1*ptab(I+1,J,K) )/(0.5*zet2+0.5*zet1) 
  !!    |     .         |                  |
  !!  ------------------------------------------
  !!    |     .         |                  |
  !!    |     .         |                  |
  !!    |     .         |                  |
  !!K   |    zbis.......U...ptab(I+1,J,K)  |
  !!    |     .         |                  |
  !!    | ptab(I,J,K)   |                  |
  !!    |               |------------------|
  !!    |               | partials         |
  !!    |               |  steps           |
  !!  -------------------------------------------
  !!    <----zet1------><----zet2--------->
  !!
  !!
  !!   ====>  s-coordinate
  !!     
  !!    |                |                  |   1. Compute distance between T1 and U points: SQRT( zdep1^2 + (0.5 * zet1 )^2
  !!    |                |                  |      Compute distance between T2 and U points: SQRT( zdep2^2 + (0.5 * zet2 )^2
  !!    |                | ptab(I+1,J,K)    | 
  !!    |                |      T2          |   2. Interpolation between  T1 and T2 values at U point 
  !!    |                |      ^           |    
  !!    |                |      | zdep2     |    
  !!    |                |      |           |    
  !!    |       ^        U      v           |
  !!    |       |        |                  |
  !!    |       | zdep1  |                  |    
  !!    |       v        |                  |
  !!    |      T1        |                  |
  !!    | ptab(I,J,K)    |                  | 
  !!    |                |                  | 
  !!    |                |                  | 
  !!
  !!    <----zet1--------><----zet2--------->
  !!
  !!----------------------------------------------------------------------
  !*arguments
  INTEGER, INTENT(IN)                          :: Kmm          ! time level index
  INTEGER, INTENT(IN)                          :: ki, kj, kk   ! coordinate of point
  CHARACTER(len=1), INTENT(IN)                 :: cd_point     ! type of point (U, V)
  REAL(dp), DIMENSION(jpi,jpj,jpk), INTENT(IN) :: ptab         ! variable to compute at (ki, kj, kk )
  REAL(wp)                                     :: interp       ! interpolated variable 

  !*local declations
  INTEGER :: ii1, ij1, ii2, ij2                                ! local integer
  REAL(wp):: ze3t, ze3, zwgt1, zwgt2, zbis, zdepu            ! local real
  REAL(wp):: zet1, zet2                                        ! weight for interpolation 
  REAL(wp):: zdep1,zdep2                                       ! differences of depth
  REAL(wp):: zmsk                                              ! mask value
  !!----------------------------------------------------------------------

  IF( cd_point=='U' )THEN 
     ii1 = ki    ; ij1 = kj 
     ii2 = ki+1  ; ij2 = kj 

     zet1=e1t(ii1,ij1)
     zet2=e1t(ii2,ij2)
     zmsk=umask(ii1,ij1,kk)
  

  ELSE ! cd_point=='V' 
     ii1 = ki    ; ij1 = kj 
     ii2 = ki    ; ij2 = kj+1  

     zet1=e2t(ii1,ij1)
     zet2=e2t(ii2,ij2)
     zmsk=vmask(ii1,ij1,kk)

  ENDIF

  IF( ln_sco )THEN   ! s-coordinate case

     zdepu = ( gdept(ii1,ij1,kk,Kmm) +  gdept(ii2,ij2,kk,Kmm) ) * 0.5_wp 
     zdep1 = gdept(ii1,ij1,kk,Kmm) - zdepu
     zdep2 = gdept(ii2,ij2,kk,Kmm) - zdepu

     ! weights
     zwgt1 = SQRT( ( 0.5 * zet1 ) * ( 0.5 * zet1 ) + ( zdep1 * zdep1 ) )
     zwgt2 = SQRT( ( 0.5 * zet2 ) * ( 0.5 * zet2 ) + ( zdep2 * zdep2 ) )
  
     ! result
     interp = zmsk * ( zwgt2 *  ptab(ii1,ij1,kk) + zwgt1 *  ptab(ii1,ij1,kk) ) / ( zwgt2 + zwgt1 )   


  ELSE       ! full step or partial step case 

     ze3t  = e3t(ii2,ij2,kk,Kmm) - e3t(ii1,ij1,kk,Kmm) 
     zwgt1 = ( e3w(ii2,ij2,kk,Kmm) - e3w(ii1,ij1,kk,Kmm) )   &
        &    / e3w(ii2,ij2,kk,Kmm)
     zwgt2 = ( e3w(ii1,ij1,kk,Kmm) - e3w(ii2,ij2,kk,Kmm) )   &
        &    / e3w(ii1,ij1,kk,Kmm)

     IF(kk .NE. 1)THEN

        IF( ze3t >= 0. )THEN 
           ! zbis
           zbis = ptab(ii2,ij2,kk) + zwgt1 * ( ptab(ii2,ij2,kk-1) - ptab(ii2,ij2,kk) ) 
           ! result
            interp = zmsk * ( zet2 * ptab(ii1,ij1,kk) + zet1 * zbis )/( zet1 + zet2 )
        ELSE
           ! zbis
           zbis = ptab(ii1,ij1,kk) + zwgt2 * ( ptab(ii1,ij1,kk-1) - ptab(ii1,ij2,kk) )
           ! result
           interp = zmsk * ( zet2 * zbis + zet1 * ptab(ii2,ij2,kk) )/( zet1 + zet2 )
        ENDIF    

     ELSE
        interp = zmsk * (  zet2 * ptab(ii1,ij1,kk) + zet1 * ptab(ii2,ij2,kk) )/( zet1 + zet2 )
     ENDIF

  ENDIF
      !
   END FUNCTION interp

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                                             
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC ::   ln_diadct = .FALSE.
CONTAINS
   SUBROUTINE dia_dct_init
      IMPLICIT NONE
   END SUBROUTINE dia_dct_init

   SUBROUTINE dia_dct( kt, Kmm )         ! Dummy routine
      IMPLICIT NONE
      INTEGER, INTENT( in ) :: kt   ! ocean time-step index
      INTEGER, INTENT( in ) :: Kmm  ! ocean time level index
      WRITE(*,*) 'dia_dct: You should not have seen this print! error?', kt
   END SUBROUTINE dia_dct
   !
#endif

   !!======================================================================
END MODULE diadct
