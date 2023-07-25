MODULE domain
   !!==============================================================================
   !!                       ***  MODULE domain   ***
   !! Ocean initialization : domain initialization
   !!==============================================================================
   !! History :  OPA  !  1990-10  (C. Levy - G. Madec)  Original code
   !!                 !  1992-01  (M. Imbard) insert time step initialization
   !!                 !  1996-06  (G. Madec) generalized vertical coordinate
   !!                 !  1997-02  (G. Madec) creation of domwri.F
   !!                 !  2001-05  (E.Durand - G. Madec) insert closed sea
   !!   NEMO     1.0  !  2002-08  (G. Madec)  F90: Free form and module
   !!            2.0  !  2005-11  (V. Garnier) Surface pressure gradient organization
   !!            3.3  !  2010-11  (G. Madec)  initialisation in C1D configuration 
   !!            3.6  !  2013     ( J. Simeon, C. Calone, G. Madec, C. Ethe ) Online coarsening of outputs
   !!            3.7  !  2015-11  (G. Madec, A. Coward)  time varying zgr by default
   !!            4.0  !  2016-10  (G. Madec, S. Flavoni)  domain configuration / user defined interface
   !!            4.1  !  2020-02  (G. Madec, S. Techene)  introduce ssh to h0 ratio
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dom_init      : initialize the space and time domain
   !!   dom_nam       : read and contral domain namelists
   !!   dom_ctl       : control print for the ocean domain
   !!   domain_cfg    : read the global domain size in domain configuration file
   !!   cfg_write     : create the domain configuration file
   !!----------------------------------------------------------------------
   USE oce            ! ocean variables
   USE dom_oce        ! domain: ocean
   USE domtile        ! tiling utilities
#if defined key_qco
   USE domqco         ! quasi-eulerian coord.
#elif defined key_linssh
   !                  ! fix in time coord.
#else
   USE domvvl         ! variable volume coord.
#endif
#if defined key_agrif
   USE agrif_oce_interp, ONLY : Agrif_istate_ssh ! ssh interpolated from parent
#if defined key_si3
   USE agrif_ice_interp, ONLY : agrif_istate_icevol ! ssh increment from ice
#endif
#endif
   USE sbc_oce        ! surface boundary condition: ocean
   USE trc_oce        ! shared ocean & passive tracers variab
   USE phycst         ! physical constants
   USE domhgr         ! domain: set the horizontal mesh
   USE domzgr         ! domain: set the vertical mesh
   USE dommsk         ! domain: set the mask system
   USE domwri         ! domain: write the meshmask file
   USE wet_dry , ONLY : ll_wd     ! wet & drying flag
   USE closea  , ONLY : dom_clo   ! closed seas routine
   USE c1d
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O library
   USE lbclnk         ! ocean lateral boundary condition (or mpp link)
   USE lib_mpp        ! distributed memory computing library
   USE restart        ! only for lrst_oce and rst_read_ssh

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dom_init     ! called by nemogcm.F90
   PUBLIC   domain_cfg   ! called by nemogcm.F90

   !! * Substitutions
#  include "single_precision_substitute.h90"
#  include "do_loop_substitute.h90"
   !!-------------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: domain.F90 15270 2021-09-17 14:27:55Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!-------------------------------------------------------------------------
CONTAINS

   SUBROUTINE dom_init( Kbb, Kmm, Kaa )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dom_init  ***
      !!
      !! ** Purpose :   Domain initialization. Call the routines that are
      !!              required to create the arrays which define the space
      !!              and time domain of the ocean model.
      !!
      !! ** Method  : - dom_msk: compute the masks from the bathymetry file
      !!              - dom_hgr: compute or read the horizontal grid-point position
      !!                         and scale factors, and the coriolis factor
      !!              - dom_zgr: define the vertical coordinate and the bathymetry
      !!              - dom_wri: create the meshmask file (ln_meshmask=T)
      !!              - 1D configuration, move Coriolis, u and v at T-point
      !!----------------------------------------------------------------------
      INTEGER          , INTENT(in) :: Kbb, Kmm, Kaa          ! ocean time level indices
      !
      INTEGER ::   ji, jj, jk, jt   ! dummy loop indices
      INTEGER ::   iconf = 0    ! local integers
      REAL(wp)::   zrdt
      CHARACTER (len=64) ::   cform = "(A12, 3(A13, I7))"
      INTEGER , DIMENSION(jpi,jpj) ::   ik_top , ik_bot       ! top and bottom ocean level
      REAL(wp), DIMENSION(jpi,jpj) ::   z1_hu_0, z1_hv_0
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN         ! Ocean domain Parameters (control print)
         WRITE(numout,*)
         WRITE(numout,*) 'dom_init : domain initialization'
         WRITE(numout,*) '~~~~~~~~'
         !
         WRITE(numout,*)     '   Domain info'
         WRITE(numout,*)     '      dimension of model:'
         WRITE(numout,*)     '             Local domain      Global domain       Data domain '
         WRITE(numout,cform) '        ','   jpi     : ', jpi, '   jpiglo  : ', jpiglo
         WRITE(numout,cform) '        ','   jpj     : ', jpj, '   jpjglo  : ', jpjglo
         WRITE(numout,cform) '        ','   jpk     : ', jpk, '   jpkglo  : ', jpkglo
         WRITE(numout,cform) '       ' ,'   jpij    : ', jpij
         WRITE(numout,*)     '      mpp local domain info (mpp):'
         WRITE(numout,*)     '              jpni    : ', jpni, '   nn_hls  : ', nn_hls
         WRITE(numout,*)     '              jpnj    : ', jpnj, '   nn_hls  : ', nn_hls
         WRITE(numout,*)     '              jpnij   : ', jpnij
         WRITE(numout,*)     '      lateral boundary of the Global domain:'
         WRITE(numout,*)     '              cyclic east-west             :', l_Iperio
         WRITE(numout,*)     '              cyclic north-south           :', l_Jperio
         WRITE(numout,*)     '              North Pole folding           :', l_NFold
         WRITE(numout,*)     '                 type of North pole Folding:', c_NFtype
         WRITE(numout,*)     '      Ocean model configuration used:'
         WRITE(numout,*)     '              cn_cfg = ', TRIM( cn_cfg ), '   nn_cfg = ', nn_cfg
      ENDIF

      !
      !           !==  Reference coordinate system  ==!
      !
      CALL dom_nam                      ! read namelist ( namrun, namdom )
      CALL dom_tile_init                ! Tile domain

      IF( ln_c1d )   CALL     c1d_init                    ! 1D column configuration
      !
      CALL dom_hgr                      ! Horizontal mesh

      IF( ln_closea ) CALL dom_clo      ! Read in masks to define closed seas and lakes

      CALL dom_zgr( ik_top, ik_bot )    ! Vertical mesh and bathymetry (return top and bottom ocean t-level indices)

      CALL dom_msk( ik_top, ik_bot )    ! Masks
      !
      ht_0(:,:) = 0._wp  ! Reference ocean thickness
      hu_0(:,:) = 0._wp
      hv_0(:,:) = 0._wp
      hf_0(:,:) = 0._wp
      DO jk = 1, jpkm1
         ht_0(:,:) = ht_0(:,:) + e3t_0(:,:,jk) * tmask(:,:,jk)
         hu_0(:,:) = hu_0(:,:) + e3u_0(:,:,jk) * umask(:,:,jk)
         hv_0(:,:) = hv_0(:,:) + e3v_0(:,:,jk) * vmask(:,:,jk)
      END DO
      !
      DO_3D( 0, 0, 0, 0, 1, jpkm1 )
         hf_0(ji,jj) = hf_0(ji,jj) + e3f_0(ji,jj,jk)*vmask(ji,jj,jk)*vmask(ji+1,jj,jk)
      END_3D
      CALL lbc_lnk('domain', hf_0, 'F', 1._wp)
      !
      IF( lk_SWE ) THEN      ! SWE case redefine hf_0
         hf_0(:,:) = hf_0(:,:) + e3f_0(:,:,1) * ssfmask(:,:)
      ENDIF
      !
      r1_ht_0(:,:) = ssmask (:,:) / ( ht_0(:,:) + 1._wp -  ssmask (:,:) )
      r1_hu_0(:,:) = ssumask(:,:) / ( hu_0(:,:) + 1._wp -  ssumask(:,:) )
      r1_hv_0(:,:) = ssvmask(:,:) / ( hv_0(:,:) + 1._wp -  ssvmask(:,:) )
      r1_hf_0(:,:) = ssfmask(:,:) / ( hf_0(:,:) + 1._wp -  ssfmask(:,:) )
      !
      IF( ll_wd ) THEN       ! wet and drying (check ht_0 >= 0)
         DO_2D( 1, 1, 1, 1 )
            IF( ht_0(ji,jj) < 0._wp .AND. ssmask(ji,jj) == 1._wp ) THEN
               CALL ctl_stop( 'dom_init : ht_0 must be positive at potentially wet points' )
            ENDIF
         END_2D
      ENDIF
      !
      !           !==  initialisation of time varying coordinate  ==!
      !
      !                                 != ssh initialization
      !
      IF( l_SAS ) THEN        !* No ocean dynamics calculation : set to 0
         ssh(:,:,:) = 0._wp
#if defined key_agrif
      ELSEIF( .NOT.Agrif_root() .AND.    &
         &     ln_init_chfrpar ) THEN        !* Interpolate initial ssh from parent
         CALL Agrif_istate_ssh( Kbb, Kmm, Kaa )
#if defined key_si3
         ! Possibly add ssh increment from parent grid
         ! only if there is no ice model in the child grid
         CALL Agrif_istate_icevol( Kbb, Kmm, Kaa ) 
#endif
#endif
      ELSE                                   !* Read in restart file or set by user
         CALL rst_read_ssh( Kbb, Kmm, Kaa )
      ENDIF
      !     
#if defined key_qco
      !                                 != Quasi-Euerian coordinate case
      !
      IF( .NOT.l_offline )   CALL dom_qco_init( Kbb, Kmm, Kaa )
#elif defined key_linssh
      !                                 != Fix in time : key_linssh case, set through domzgr_substitute.h90
#else
      !
      IF( ln_linssh ) THEN              != Fix in time : set to the reference one for all
         !
         DO jt = 1, jpt                         ! depth of t- and w-grid-points
            gdept(:,:,:,jt) = gdept_0(:,:,:)
            gdepw(:,:,:,jt) = gdepw_0(:,:,:)
         END DO
            gde3w(:,:,:)    = gde3w_0(:,:,:)    ! = gdept as the sum of e3t
         !
         DO jt = 1, jpt                         ! vertical scale factors
            e3t (:,:,:,jt) =  e3t_0(:,:,:)
            e3u (:,:,:,jt) =  e3u_0(:,:,:)
            e3v (:,:,:,jt) =  e3v_0(:,:,:)
            e3w (:,:,:,jt) =  e3w_0(:,:,:)
            e3uw(:,:,:,jt) = e3uw_0(:,:,:)
            e3vw(:,:,:,jt) = e3vw_0(:,:,:)
         END DO
            e3f (:,:,:)    =  e3f_0(:,:,:)
         !
         DO jt = 1, jpt                         ! water column thickness and its inverse
               hu(:,:,jt) =    hu_0(:,:)
               hv(:,:,jt) =    hv_0(:,:)
            r1_hu(:,:,jt) = r1_hu_0(:,:)
            r1_hv(:,:,jt) = r1_hv_0(:,:)
         END DO
               ht   (:,:) =    ht_0(:,:)
         !
      ELSE                              != Time varying : initialize before/now/after variables
         !
         IF( .NOT.l_offline )   CALL dom_vvl_init( Kbb, Kmm, Kaa )
         !
      ENDIF
#endif

      !

#if defined key_agrif
      IF( .NOT. Agrif_Root() ) CALL Agrif_Init_Domain( Kbb, Kmm, Kaa )
#endif
      IF( ln_meshmask    )   CALL dom_wri       ! Create a domain file
      IF( .NOT.ln_rstart )   CALL dom_ctl       ! Domain control
      !
      IF( ln_write_cfg   )   CALL cfg_write     ! create the configuration file
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dom_init :   ==>>>   END of domain initialization'
         WRITE(numout,*) '~~~~~~~~'
         WRITE(numout,*)
      ENDIF
      !
   END SUBROUTINE dom_init


   SUBROUTINE dom_nam
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dom_nam  ***
      !!
      !! ** Purpose :   read domaine namelists and print the variables.
      !!
      !! ** input   : - namrun namelist
      !!              - namdom namelist
      !!              - namtile namelist
      !!              - namnc4 namelist   ! "key_netcdf4" only
      !!----------------------------------------------------------------------
      USE ioipsl
      !!
      INTEGER ::   ios   ! Local integer
      REAL(wp)::   zrdt
      !!----------------------------------------------------------------------
      !
      NAMELIST/namrun/ cn_ocerst_indir, cn_ocerst_outdir, nn_stocklist, ln_rst_list,                 &
         &             nn_no   , cn_exp   , cn_ocerst_in, cn_ocerst_out, ln_rstart , ln_reset_ts ,   &
         &             ln_rstdate, nn_rstctl,                                                        &
         &             nn_it000, nn_itend , nn_date0    , nn_time0     , nn_leapy  , nn_istate ,     &
         &             nn_stock, nn_write , ln_mskland  , ln_clobber   , nn_chunksz, ln_1st_euler  , &
         &             ln_cfmeta, ln_xios_read, nn_wxios, ln_rst_eos
      NAMELIST/namdom/ ln_linssh, rn_Dt, rn_atfp, ln_crs, ln_c1d, ln_meshmask
      NAMELIST/namtile/ ln_tile, nn_ltile_i, nn_ltile_j
#if defined key_netcdf4
      NAMELIST/namnc4/ nn_nchunks_i, nn_nchunks_j, nn_nchunks_k, ln_nc4zip
#endif
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dom_nam : domain initialization through namelist read'
         WRITE(numout,*) '~~~~~~~ '
      ENDIF
      !
      !                       !=======================!
      !                       !==  namelist namdom  ==!
      !                       !=======================!
      !
      READ  ( numnam_ref, namdom, IOSTAT = ios, ERR = 903)
903   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namdom in reference namelist' )
      READ  ( numnam_cfg, namdom, IOSTAT = ios, ERR = 904 )
904   IF( ios >  0 )   CALL ctl_nam ( ios , 'namdom in configuration namelist' )
      IF(lwm) WRITE( numond, namdom )
      !
#if defined key_linssh
      ln_linssh = lk_linssh      ! overwrite ln_linssh with the logical associated with key_linssh
#endif
      !
#if defined key_agrif
      IF( .NOT. Agrif_Root() ) THEN    ! AGRIF child, subdivide the Parent timestep
         rn_Dt = Agrif_Parent (rn_Dt ) / Agrif_Rhot()
      ENDIF
#endif
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) '   Namelist : namdom   ---   space & time domain'
         WRITE(numout,*) '      linear free surface (=T)                ln_linssh   = ', ln_linssh
         WRITE(numout,*) '      create mesh/mask file                   ln_meshmask = ', ln_meshmask
         WRITE(numout,*) '      ocean time step                         rn_Dt       = ', rn_Dt
         WRITE(numout,*) '      asselin time filter parameter           rn_atfp     = ', rn_atfp
         WRITE(numout,*) '      online coarsening of dynamical fields   ln_crs      = ', ln_crs
         WRITE(numout,*) '      single column domain (1x1pt)            ln_c1d      = ', ln_c1d
      ENDIF
      !
      ! set current model timestep rDt = 2*rn_Dt if MLF or rDt = rn_Dt if RK3
      rDt   = 2._wp * rn_Dt
      r1_Dt = 1._wp / rDt
      !
      IF( l_SAS .AND. .NOT.ln_linssh ) THEN
         CALL ctl_warn( 'SAS requires linear ssh : force ln_linssh = T' )
         ln_linssh = .TRUE.
      ENDIF
      !
#if defined key_qco
      IF( ln_linssh )   CALL ctl_stop( 'STOP','domain: key_qco and ln_linssh=T or key_linssh are incompatible' )
#endif
      !
      !                       !=======================!
      !                       !==  namelist namrun  ==!
      !                       !=======================!
      !
      READ  ( numnam_ref, namrun, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namrun in reference namelist' )
      READ  ( numnam_cfg, namrun, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namrun in configuration namelist' )
      IF(lwm) WRITE ( numond, namrun )

#if defined key_agrif
      IF( .NOT. Agrif_Root() ) THEN
            nn_it000 = (Agrif_Parent(nn_it000)-1)*Agrif_IRhot() + 1
            nn_itend =  Agrif_Parent(nn_itend)   *Agrif_IRhot()
            nn_date0 =  Agrif_Parent(nn_date0)
            nn_time0 =  Agrif_Parent(nn_time0)
            nn_leapy =  Agrif_Parent(nn_leapy)
      ENDIF
#endif
      !
      IF(lwp) THEN                  ! control print
         WRITE(numout,*) '   Namelist : namrun   ---   run parameters'
         WRITE(numout,*) '      Assimilation cycle              nn_no           = ', nn_no
         WRITE(numout,*) '      experiment name for output      cn_exp          = ', TRIM( cn_exp           )
         WRITE(numout,*) '      file prefix restart input       cn_ocerst_in    = ', TRIM( cn_ocerst_in     )
         WRITE(numout,*) '      restart input directory         cn_ocerst_indir = ', TRIM( cn_ocerst_indir  )
         WRITE(numout,*) '      file prefix restart output      cn_ocerst_out   = ', TRIM( cn_ocerst_out    )
         WRITE(numout,*) '      restart output directory        cn_ocerst_outdir= ', TRIM( cn_ocerst_outdir )
         WRITE(numout,*) '      restart logical                 ln_rstart       = ', ln_rstart
         WRITE(numout,*) '      reset TS from inital TS file    ln_reset_ts     = ', ln_reset_ts
         WRITE(numout,*) '      start with forward time step    ln_1st_euler    = ', ln_1st_euler
         WRITE(numout,*) '      control of time step            nn_rstctl       = ', nn_rstctl
         WRITE(numout,*) '      number of the first time step   nn_it000        = ', nn_it000
         WRITE(numout,*) '      number of the last time step    nn_itend        = ', nn_itend
         WRITE(numout,*) '      initial calendar date aammjj    nn_date0        = ', nn_date0
         WRITE(numout,*) '      initial time of day in hhmm     nn_time0        = ', nn_time0
         WRITE(numout,*) '      leap year calendar (0/1)        nn_leapy        = ', nn_leapy
         WRITE(numout,*) '      initial state output            nn_istate       = ', nn_istate
         IF( ln_rst_list ) THEN
            WRITE(numout,*) '      list of restart dump times      nn_stocklist    =', nn_stocklist
         ELSE
            WRITE(numout,*) '      frequency of restart file       nn_stock        = ', nn_stock
         ENDIF
#if ! defined key_xios
         WRITE(numout,*) '      frequency of output file        nn_write        = ', nn_write
#endif
         WRITE(numout,*) '      mask land points                ln_mskland      = ', ln_mskland
         WRITE(numout,*) '      date-stamp restart files        ln_rstdate      = ', ln_rstdate
         WRITE(numout,*) '      additional CF standard metadata ln_cfmeta       = ', ln_cfmeta
         WRITE(numout,*) '      overwrite an existing file      ln_clobber      = ', ln_clobber
         WRITE(numout,*) '      NetCDF chunksize (bytes)        nn_chunksz      = ', nn_chunksz
         WRITE(numout,*) '      check restart equation of state ln_rst_eos      = ', ln_rst_eos
         IF( TRIM(Agrif_CFixed()) == '0' ) THEN
            WRITE(numout,*) '      READ restart for a single file using XIOS ln_xios_read =', ln_xios_read
            WRITE(numout,*) '      Write restart using XIOS        nn_wxios   = ', nn_wxios
         ELSE
            WRITE(numout,*) "      AGRIF: nn_wxios will be ingored. See setting for parent"
            WRITE(numout,*) "      AGRIF: ln_xios_read will be ingored. See setting for parent"
         ENDIF
      ENDIF

      cexper = cn_exp         ! conversion DOCTOR names into model names (this should disappear soon)
      nrstdt = nn_rstctl
      nit000 = nn_it000
      nitend = nn_itend
      ndate0 = nn_date0
      nleapy = nn_leapy
      ninist = nn_istate
      !
      !                                        !==  Set parameters for restart reading using xIOS  ==!
      !
      IF( TRIM(Agrif_CFixed()) == '0' ) THEN
         lrxios = ln_xios_read .AND. ln_rstart
         IF( nn_wxios > 0 )   lwxios = .TRUE.           !* set output file type for XIOS based on NEMO namelist
         nxioso = nn_wxios
      ENDIF
      !                                        !==  Check consistency between ln_rstart and ln_1st_euler  ==!   (i.e. set l_1st_euler)
      l_1st_euler = ln_1st_euler
      !
      IF( ln_rstart ) THEN                              !*  Restart case
         !
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   open the restart file'
         CALL rst_read_open                                              !- Open the restart file
         !
         IF( iom_varid( numror, 'rdt', ldstop = .FALSE. ) > 0 ) THEN     !- Check time-step consistency and force Euler restart if changed
            CALL iom_get( numror, 'rdt', zrdt )
            IF( zrdt /= rn_Dt ) THEN
               IF(lwp) WRITE( numout,*)
               IF(lwp) WRITE( numout,*) '   rn_Dt = ', rn_Dt,' not equal to the READ one rdt = ', zrdt
               IF(lwp) WRITE( numout,*)
               IF(lwp) WRITE( numout,*) '      ==>>>   forced euler first time-step'
               l_1st_euler =  .TRUE.
            ENDIF
         ENDIF
         !
         IF( .NOT.l_SAS .AND. iom_varid( numror, 'sshb', ldstop = .FALSE. ) <= 0 ) THEN   !- Check absence of one of the Kbb field (here sshb)
            !                                                                             !  (any Kbb field is missing ==> all Kbb fields are missing)
            IF( .NOT.l_1st_euler ) THEN
               CALL ctl_warn('dom_nam : ssh at Kbb not found in restart files ',   &
                  &                        'l_1st_euler forced to .true. and ' ,   &
                  &                        'ssh(Kbb) = ssh(Kmm) '                  )
               l_1st_euler = .TRUE.
            ENDIF
         ENDIF
      ELSEIF( .NOT.l_1st_euler ) THEN                   !*  Initialization case
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*)'   ==>>>   Start from rest (ln_rstart=F)'
         IF(lwp) WRITE(numout,*)'           an Euler initial time step is used : l_1st_euler is forced to .true. '
         l_1st_euler = .TRUE.
      ENDIF
      !
      !                                        !==  control of output frequency  ==!
      !
      IF( .NOT. ln_rst_list ) THEN   ! we use nn_stock
         IF( nn_stock == -1 )   CALL ctl_warn( 'nn_stock = -1 --> no restart will be done' )
         IF( nn_stock == 0 .OR. nn_stock > nitend ) THEN
            WRITE(ctmp1,*) 'nn_stock = ', nn_stock, ' it is forced to ', nitend
            CALL ctl_warn( ctmp1 )
            nn_stock = nitend
         ENDIF
      ENDIF
#if ! defined key_xios
      IF( nn_write == -1 )   CALL ctl_warn( 'nn_write = -1 --> no output files will be done' )
      IF ( nn_write == 0 ) THEN
         WRITE(ctmp1,*) 'nn_write = ', nn_write, ' it is forced to ', nitend
         CALL ctl_warn( ctmp1 )
         nn_write = nitend
      ENDIF
#endif

      IF( Agrif_Root() ) THEN
         IF(lwp) WRITE(numout,*)
         SELECT CASE ( nleapy )                !==  Choose calendar for IOIPSL  ==!
         CASE (  1 )
            CALL ioconf_calendar('gregorian')
            IF(lwp) WRITE(numout,*) '   ==>>>   The IOIPSL calendar is "gregorian", i.e. leap year'
         CASE (  0 )
            CALL ioconf_calendar('noleap')
            IF(lwp) WRITE(numout,*) '   ==>>>   The IOIPSL calendar is "noleap", i.e. no leap year'
         CASE ( 30 )
            CALL ioconf_calendar('360d')
            IF(lwp) WRITE(numout,*) '   ==>>>   The IOIPSL calendar is "360d", i.e. 360 days in a year'
         END SELECT
      ENDIF
      !
      !                       !========================!
      !                       !==  namelist namtile  ==!
      !                       !========================!
      !
      READ  ( numnam_ref, namtile, IOSTAT = ios, ERR = 905 )
905   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namtile in reference namelist' )
      READ  ( numnam_cfg, namtile, IOSTAT = ios, ERR = 906 )
906   IF( ios >  0 )   CALL ctl_nam ( ios , 'namtile in configuration namelist' )
      IF(lwm) WRITE( numond, namtile )

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*)    '   Namelist : namtile   ---   Domain tiling decomposition'
         WRITE(numout,*)    '      Tiling (T) or not (F)                ln_tile    = ', ln_tile
         WRITE(numout,*)    '      Length of tile in i                  nn_ltile_i = ', nn_ltile_i
         WRITE(numout,*)    '      Length of tile in j                  nn_ltile_j = ', nn_ltile_j
         WRITE(numout,*)
         IF( ln_tile ) THEN
            WRITE(numout,*) '      The domain will be decomposed into tiles of size', nn_ltile_i, 'x', nn_ltile_j
         ELSE
            WRITE(numout,*) '      Domain tiling will NOT be used'
         ENDIF
      ENDIF
      !
#if defined key_netcdf4
      !                       !=======================!
      !                       !==  namelist namnc4  ==!   NetCDF 4 case   ("key_netcdf4" defined)
      !                       !=======================!
      !
      READ  ( numnam_ref, namnc4, IOSTAT = ios, ERR = 907)
907   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namnc4 in reference namelist' )
      READ  ( numnam_cfg, namnc4, IOSTAT = ios, ERR = 908 )
908   IF( ios >  0 )   CALL ctl_nam ( ios , 'namnc4 in configuration namelist' )
      IF(lwm) WRITE( numond, namnc4 )

      IF(lwp) THEN                        ! control print
         WRITE(numout,*)
         WRITE(numout,*) '   Namelist namnc4 - Netcdf4 chunking parameters ("key_netcdf4" defined)'
         WRITE(numout,*) '      number of chunks in i-dimension             nn_nchunks_i = ', nn_nchunks_i
         WRITE(numout,*) '      number of chunks in j-dimension             nn_nchunks_j = ', nn_nchunks_j
         WRITE(numout,*) '      number of chunks in k-dimension             nn_nchunks_k = ', nn_nchunks_k
         WRITE(numout,*) '      apply netcdf4/hdf5 chunking & compression   ln_nc4zip    = ', ln_nc4zip
      ENDIF

      ! Put the netcdf4 settings into a simple structure (snc4set, defined in in_out_manager module)
      ! Note the chunk size in the unlimited (time) dimension will be fixed at 1
      snc4set%ni   = nn_nchunks_i
      snc4set%nj   = nn_nchunks_j
      snc4set%nk   = nn_nchunks_k
      snc4set%luse = ln_nc4zip
#else
      snc4set%luse = .FALSE.        ! No NetCDF 4 case
#endif
      !
   END SUBROUTINE dom_nam


   SUBROUTINE dom_ctl
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dom_ctl  ***
      !!
      !! ** Purpose :   Domain control.
      !!
      !! ** Method  :   compute and print extrema of masked scale factors
      !!----------------------------------------------------------------------
      LOGICAL, DIMENSION(jpi,jpj) ::   llmsk
      INTEGER, DIMENSION(2)       ::   imil, imip, imi1, imi2, imal, imap, ima1, ima2
      REAL(wp)                    ::   zglmin, zglmax, zgpmin, zgpmax, ze1min, ze1max, ze2min, ze2max
      !!----------------------------------------------------------------------
      !
      llmsk = tmask_i(:,:) == 1._wp
      !
      CALL mpp_minloc( 'domain', CASTDP(glamt(:,:)), llmsk, zglmin, imil )
      CALL mpp_minloc( 'domain', CASTDP(gphit(:,:)), llmsk, zgpmin, imip )
      CALL mpp_minloc( 'domain',   CASTDP(e1t(:,:)), llmsk, ze1min, imi1 )
      CALL mpp_minloc( 'domain',   CASTDP(e2t(:,:)), llmsk, ze2min, imi2 )
      CALL mpp_maxloc( 'domain', CASTDP(glamt(:,:)), llmsk, zglmax, imal )
      CALL mpp_maxloc( 'domain', CASTDP(gphit(:,:)), llmsk, zgpmax, imap )
      CALL mpp_maxloc( 'domain',   e1t(:,:), llmsk, ze1max, ima1 )
      CALL mpp_maxloc( 'domain',   e2t(:,:), llmsk, ze2max, ima2 )
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dom_ctl : extrema of the masked scale factors'
         WRITE(numout,*) '~~~~~~~'
         WRITE(numout,"(14x,'glamt mini: ',1f10.2,' at i = ',i5,' j= ',i5)") zglmin, imil(1), imil(2)
         WRITE(numout,"(14x,'glamt maxi: ',1f10.2,' at i = ',i5,' j= ',i5)") zglmax, imal(1), imal(2)
         WRITE(numout,"(14x,'gphit mini: ',1f10.2,' at i = ',i5,' j= ',i5)") zgpmin, imip(1), imip(2)
         WRITE(numout,"(14x,'gphit maxi: ',1f10.2,' at i = ',i5,' j= ',i5)") zgpmax, imap(1), imap(2)
         WRITE(numout,"(14x,'  e1t mini: ',1f10.2,' at i = ',i5,' j= ',i5)") ze1min, imi1(1), imi1(2)
         WRITE(numout,"(14x,'  e1t maxi: ',1f10.2,' at i = ',i5,' j= ',i5)") ze1max, ima1(1), ima1(2)
         WRITE(numout,"(14x,'  e2t mini: ',1f10.2,' at i = ',i5,' j= ',i5)") ze2min, imi2(1), imi2(2)
         WRITE(numout,"(14x,'  e2t maxi: ',1f10.2,' at i = ',i5,' j= ',i5)") ze2max, ima2(1), ima2(2)
      ENDIF
      !
   END SUBROUTINE dom_ctl


   SUBROUTINE domain_cfg( cd_cfg, kk_cfg, kpi, kpj, kpk, ldIperio, ldJperio, ldNFold, cdNFtype )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE domain_cfg  ***
      !!
      !! ** Purpose :   read the domain size in domain configuration file
      !!
      !! ** Method  :   read the cn_domcfg NetCDF file
      !!----------------------------------------------------------------------
      CHARACTER(len=*), INTENT(out) ::   cd_cfg               ! configuration name
      INTEGER         , INTENT(out) ::   kk_cfg               ! configuration resolution
      INTEGER         , INTENT(out) ::   kpi, kpj, kpk        ! global domain sizes
      LOGICAL         , INTENT(out) ::   ldIperio, ldJperio   ! i- and j- periodicity
      LOGICAL         , INTENT(out) ::   ldNFold              ! North pole folding
      CHARACTER(len=1), INTENT(out) ::   cdNFtype             ! Folding type: T or F
      !
      CHARACTER(len=7) ::   catt                  ! 'T', 'F', '-' or 'UNKNOWN'
      INTEGER ::   inum, iperio, iatt             ! local integer
      REAL(wp) ::   zorca_res                     ! local scalars
      REAL(wp) ::   zperio                        !   -      -
      INTEGER, DIMENSION(4) ::   idvar, idimsz    ! size   of dimensions
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*) '           '
         WRITE(numout,*) 'domain_cfg : domain size read in ', TRIM( cn_domcfg ), ' file'
         WRITE(numout,*) '~~~~~~~~~~ '
      ENDIF
      !
      CALL iom_open( cn_domcfg, inum )
      !
      CALL iom_getatt( inum,  'CfgName', cd_cfg )   ! returns 'UNKNOWN' if not found
      CALL iom_getatt( inum, 'CfgIndex', kk_cfg )   ! returns      -999 if not found
      !
      ! ------- keep compatibility with OLD VERSION... start -------
      IF( cd_cfg == 'UNKNOWN' .AND. kk_cfg == -999 ) THEN
         IF(  iom_varid( inum, 'ORCA'       , ldstop = .FALSE. ) > 0  .AND.  &
            & iom_varid( inum, 'ORCA_index' , ldstop = .FALSE. ) > 0    ) THEN
            !
            cd_cfg = 'ORCA'
            CALL iom_get( inum, 'ORCA_index', zorca_res )   ;   kk_cfg = NINT( zorca_res )
            !
         ELSE
            CALL iom_getatt( inum, 'cn_cfg', cd_cfg )  ! returns 'UNKNOWN' if not found
            CALL iom_getatt( inum, 'nn_cfg', kk_cfg )  ! returns      -999 if not found
         ENDIF
      ENDIF
      ! ------- keep compatibility with OLD VERSION... end -------
      !
      idvar = iom_varid( inum, 'e3t_0', kdimsz = idimsz )   ! use e3t_0, that must exist, to get jp(ijk)glo
      kpi = idimsz(1)
      kpj = idimsz(2)
      kpk = idimsz(3)
      !
      CALL iom_getatt( inum, 'Iperio', iatt )   ;   ldIperio = iatt == 1   ! returns      -999 if not found -> default = .false.
      CALL iom_getatt( inum, 'Jperio', iatt )   ;   ldJperio = iatt == 1   ! returns      -999 if not found -> default = .false.
      CALL iom_getatt( inum,  'NFold', iatt )   ;   ldNFold  = iatt == 1   ! returns      -999 if not found -> default = .false.
      CALL iom_getatt( inum, 'NFtype', catt )                              ! returns 'UNKNOWN' if not found
      IF( LEN_TRIM(catt) == 1 ) THEN   ;   cdNFtype = TRIM(catt)
      ELSE                             ;   cdNFtype = '-'
      ENDIF
      !
      ! ------- keep compatibility with OLD VERSION... start -------
      IF( iatt == -999 .AND. catt == 'UNKNOWN' .AND. iom_varid( inum, 'jperio', ldstop = .FALSE. ) > 0 ) THEN                    
         CALL iom_get( inum, 'jperio', zperio )   ;   iperio = NINT( zperio )
         ldIperio = iperio == 1  .OR. iperio == 4 .OR. iperio == 6 .OR. iperio == 7   ! i-periodicity
         ldJperio = iperio == 2  .OR. iperio == 7                                     ! j-periodicity
         ldNFold  = iperio >= 3 .AND. iperio <= 6                                     ! North pole folding
         IF(     iperio == 3 .OR. iperio == 4 ) THEN   ;   cdNFtype = 'T'             !    folding at T point
         ELSEIF( iperio == 5 .OR. iperio == 6 ) THEN   ;   cdNFtype = 'F'             !    folding at F point
         ELSE                                          ;   cdNFtype = '-'             !    default value
         ENDIF
      ENDIF
      ! ------- keep compatibility with OLD VERSION... end -------
      !
      CALL iom_close( inum )
      !
      IF(lwp) THEN
         WRITE(numout,*) '   .'
         WRITE(numout,*) '   ==>>>   ', TRIM(cn_cfg), ' configuration '
         WRITE(numout,*) '   .'
         WRITE(numout,*) '      nn_cfg = ', kk_cfg
         WRITE(numout,*) '      Ni0glo = ', kpi
         WRITE(numout,*) '      Nj0glo = ', kpj
         WRITE(numout,*) '      jpkglo = ', kpk
      ENDIF
      !
   END SUBROUTINE domain_cfg


   SUBROUTINE cfg_write
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE cfg_write  ***
      !!
      !! ** Purpose :   Create the "cn_domcfg_out" file, a NetCDF file which
      !!              contains all the ocean domain informations required to
      !!              define an ocean configuration.
      !!
      !! ** Method  :   Write in a file all the arrays required to set up an
      !!              ocean configuration.
      !!
      !! ** output file :   domcfg_out.nc : domain size, characteristics, horizontal
      !!                       mesh, Coriolis parameter, and vertical scale factors
      !!                    NB: also contain ORCA family information
      !!----------------------------------------------------------------------
      INTEGER           ::   ji, jj, jk   ! dummy loop indices
      INTEGER           ::   inum     ! local units
      CHARACTER(len=21) ::   clnam    ! filename (mesh and mask informations)
      REAL(wp), DIMENSION(jpi,jpj) ::   z2d   ! workspace
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'cfg_write : create the domain configuration file (', TRIM(cn_domcfg_out),'.nc)'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~'
      !
      !                       ! ============================= !
      !                       !  create 'domcfg_out.nc' file  !
      !                       ! ============================= !
      !
      clnam = cn_domcfg_out  ! filename (configuration information)
      CALL iom_open( TRIM(clnam), inum, ldwrt = .TRUE. )
      !
      !                             !==  Configuration specificities  ==!
      !
      CALL iom_putatt( inum,  'CfgName', TRIM(cn_cfg) )
      CALL iom_putatt( inum, 'CfgIndex',      nn_cfg  )
      !
      !                             !==  domain characteristics  ==!
      !
      !                                   ! lateral boundary of the global domain
      CALL iom_putatt( inum, 'Iperio', COUNT( (/l_Iperio/) ) )
      CALL iom_putatt( inum, 'Jperio', COUNT( (/l_Jperio/) ) )
      CALL iom_putatt( inum,  'NFold', COUNT( (/l_NFold /) ) )
      CALL iom_putatt( inum, 'NFtype',          c_NFtype     )

      !                                   ! type of vertical coordinate
      IF(ln_zco)   CALL iom_putatt( inum, 'VertCoord', 'zco' )
      IF(ln_zps)   CALL iom_putatt( inum, 'VertCoord', 'zps' )
      IF(ln_sco)   CALL iom_putatt( inum, 'VertCoord', 'sco' )
      
      !                                   ! ocean cavities under iceshelves
      CALL iom_putatt( inum, 'IsfCav', COUNT( (/ln_isfcav/) ) )
      !
      !                             !==  horizontal mesh  !
      !
      CALL iom_rstput( 0, 0, inum, 'glamt', glamt, ktype = jp_r8 )   ! latitude
      CALL iom_rstput( 0, 0, inum, 'glamu', glamu, ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'glamv', glamv, ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'glamf', glamf, ktype = jp_r8 )
      !
      CALL iom_rstput( 0, 0, inum, 'gphit', gphit, ktype = jp_r8 )   ! longitude
      CALL iom_rstput( 0, 0, inum, 'gphiu', gphiu, ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'gphiv', gphiv, ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'gphif', gphif, ktype = jp_r8 )
      !
      CALL iom_rstput( 0, 0, inum, 'e1t'  , e1t  , ktype = jp_r8 )   ! i-scale factors (e1.)
      CALL iom_rstput( 0, 0, inum, 'e1u'  , e1u  , ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'e1v'  , e1v  , ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'e1f'  , e1f  , ktype = jp_r8 )
      !
      CALL iom_rstput( 0, 0, inum, 'e2t'  , e2t  , ktype = jp_r8 )   ! j-scale factors (e2.)
      CALL iom_rstput( 0, 0, inum, 'e2u'  , e2u  , ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'e2v'  , e2v  , ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'e2f'  , e2f  , ktype = jp_r8 )
      !
      CALL iom_rstput( 0, 0, inum, 'ff_f' , ff_f , ktype = jp_r8 )   ! coriolis factor
      CALL iom_rstput( 0, 0, inum, 'ff_t' , ff_t , ktype = jp_r8 )
      !
      !                             !==  vertical mesh  ==!
      !
      CALL iom_rstput( 0, 0, inum, 'e3t_1d'  , e3t_1d , ktype = jp_r8 )   ! reference 1D-coordinate
      CALL iom_rstput( 0, 0, inum, 'e3w_1d'  , e3w_1d , ktype = jp_r8 )
      !
      CALL iom_rstput( 0, 0, inum, 'e3t_0'   , e3t_0  , ktype = jp_r8 )   ! vertical scale factors
      CALL iom_rstput( 0, 0, inum, 'e3u_0'   , e3u_0  , ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'e3v_0'   , e3v_0  , ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'e3f_0'   , e3f_0  , ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'e3w_0'   , e3w_0  , ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'e3uw_0'  , e3uw_0 , ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'e3vw_0'  , e3vw_0 , ktype = jp_r8 )
      !
      !                             !==  wet top and bottom level  ==!   (caution: multiplied by ssmask)
      !
      CALL iom_rstput( 0, 0, inum, 'top_level'    , REAL( mikt, wp )*ssmask , ktype = jp_i4 )   ! nb of ocean T-points (ISF)
      CALL iom_rstput( 0, 0, inum, 'bottom_level' , REAL( mbkt, wp )*ssmask , ktype = jp_i4 )   ! nb of ocean T-points
      !
      IF( ln_sco ) THEN             ! s-coordinate: store grid stiffness ratio  (Not required anyway)
         CALL dom_stiff( z2d )
         CALL iom_rstput( 0, 0, inum, 'stiffness', z2d )        !    ! Max. grid stiffness ratio
      ENDIF
      !
      IF( ll_wd ) THEN              ! wetting and drying domain
         CALL iom_rstput( 0, 0, inum, 'ht_0'   , ht_0   , ktype = jp_r8 )
      ENDIF
      !                       ! ============================ !
      !                       !        close the files
      !                       ! ============================ !
      CALL iom_close( inum )
      !
   END SUBROUTINE cfg_write

   !!======================================================================
END MODULE domain
