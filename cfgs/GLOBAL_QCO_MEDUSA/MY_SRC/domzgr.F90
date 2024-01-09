MODULE domzgr
   !!==============================================================================
   !!                       ***  MODULE domzgr   ***
   !! Ocean domain : definition of the vertical coordinate system
   !!==============================================================================
   !! History :  OPA  ! 1995-12  (G. Madec)  Original code : s vertical coordinate
   !!                 ! 1997-07  (G. Madec)  lbc_lnk call
   !!                 ! 1997-04  (J.-O. Beismann) 
   !!            8.5  ! 2002-09  (A. Bozec, G. Madec)  F90: Free form and module
   !!             -   ! 2002-09  (A. de Miranda)  rigid-lid + islands
   !!  NEMO      1.0  ! 2003-08  (G. Madec)  F90: Free form and module
   !!             -   ! 2005-10  (A. Beckmann)  modifications for hybrid s-ccordinates & new stretching function
   !!            2.0  ! 2006-04  (R. Benshila, G. Madec)  add zgr_zco
   !!            3.0  ! 2008-06  (G. Madec)  insertion of domzgr_zps.h90 & conding style
   !!            3.2  ! 2009-07  (R. Benshila) Suppression of rigid-lid option
   !!            3.3  ! 2010-11  (G. Madec) add mbk. arrays associated to the deepest ocean level
   !!            3.4  ! 2012-08  (J. Siddorn) added Siddorn and Furner stretching function
   !!            3.4  ! 2012-12  (R. Bourdalle-Badie and G. Reffray)  modify C1D case  
   !!            3.6  ! 2014-11  (P. Mathiot and C. Harris) add ice shelf capabilitye  
   !!            3.?  ! 2015-11  (H. Liu) Modifications for Wetting/Drying
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dom_zgr       : read or set the ocean vertical coordinate system
   !!   zgr_read      : read the vertical information in the domain configuration file
   !!   zgr_top_bot   : ocean top and bottom level for t-, u, and v-points with 1 as minimum value
   !!---------------------------------------------------------------------
   USE oce            ! ocean variables
   USE dom_oce        ! ocean domain
   USE usrdef_zgr     ! user defined vertical coordinate system
   USE closea         ! closed seas
   USE depth_e3       ! depth <=> e3
   USE wet_dry,   ONLY: ll_wd, ssh_ref  ! Wetting and drying
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! distributed memory computing library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dom_zgr        ! called by dom_init.F90

  !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: domzgr.F90 15556 2021-11-29 15:23:06Z jchanut $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS       

   SUBROUTINE dom_zgr( k_top, k_bot )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE dom_zgr  ***
      !!                   
      !! ** Purpose :   set the depth of model levels and the resulting 
      !!              vertical scale factors.
      !!
      !! ** Method  : - reference 1D vertical coordinate (gdep._1d, e3._1d)
      !!              - read/set ocean depth and ocean levels (bathy, mbathy)
      !!              - vertical coordinate (gdep., e3.) depending on the 
      !!                coordinate chosen :
      !!                   ln_zco=T   z-coordinate   
      !!                   ln_zps=T   z-coordinate with partial steps
      !!                   ln_zco=T   s-coordinate 
      !!
      !! ** Action  :   define gdep., e3., mbathy and bathy
      !!----------------------------------------------------------------------
      INTEGER, DIMENSION(:,:), INTENT(out) ::   k_top, k_bot   ! ocean first and last level indices
      !
      INTEGER  ::   ji,jj,jk            ! dummy loop index
      INTEGER  ::   ikt, ikb            ! top/bot index
      INTEGER  ::   ioptio, ibat, ios   ! local integer
      INTEGER  ::   is_mbkuvf           ! ==0 if mbku, mbkv, mbkf to be computed
      REAL(wp) ::   zrefdep             ! depth of the reference level (~10m)
      REAL(wp), DIMENSION(jpi,jpj  ) ::   zmsk
      REAL(wp), DIMENSION(jpi,jpj,2) ::   ztopbot
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN                     ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'dom_zgr : vertical coordinate'
         WRITE(numout,*) '~~~~~~~'
      ENDIF

      IF( ln_linssh .AND. lwp) WRITE(numout,*) '   linear free surface: the vertical mesh does not change in time'


      IF( ln_read_cfg ) THEN        !==  read in mesh_mask.nc file  ==!
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   ==>>>   Read vertical mesh in ', TRIM( cn_domcfg ), ' file'
         !
         CALL zgr_read   ( ln_zco  , ln_zps  , ln_sco, ln_isfcav,   & 
            &              gdept_1d, gdepw_1d, e3t_1d, e3w_1d   ,   &    ! 1D gridpoints depth
            &              gdept_0 , gdepw_0                    ,   &    ! gridpoints depth 
            &              e3t_0   , e3u_0   , e3v_0 , e3f_0    ,   &    ! vertical scale factors
            &              e3w_0   , e3uw_0  , e3vw_0           ,   &    ! vertical scale factors
            &              k_top   , k_bot                      ,   &    ! 1st & last ocean level
            &              is_mbkuvf, mbku, mbkv, mbkf )                 ! U/V/F points bottom levels
            !
      ELSE                          !==  User defined configuration  ==!
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '          User defined vertical mesh (usr_def_zgr)'
         is_mbkuvf = 0
         !
         CALL usr_def_zgr( ln_zco  , ln_zps  , ln_sco, ln_isfcav,   & 
            &              gdept_1d, gdepw_1d, e3t_1d, e3w_1d   ,   &    ! 1D gridpoints depth
            &              gdept_0 , gdepw_0                    ,   &    ! gridpoints depth 
            &              e3t_0   , e3u_0   , e3v_0 , e3f_0    ,   &    ! vertical scale factors
            &              e3w_0   , e3uw_0  , e3vw_0           ,   &    ! vertical scale factors
            &              k_top   , k_bot            )                  ! 1st & last ocean level
         !
         ! make sure that periodicities are properly applied 
         CALL lbc_lnk( 'dom_zgr', gdept_0, 'T', 1._wp, gdepw_0, 'W', 1._wp,                     &
            &                       e3u_0, 'U', 1._wp,  e3v_0, 'V', 1._wp, e3f_0, 'F', 1._wp,   &
            &                       e3w_0, 'W', 1._wp,  e3uw_0, 'U', 1._wp, e3vw_0, 'V', 1._wp,   &   
            &                     kfillmode = jpfillcopy )   ! do not put 0 over closed boundaries
         CALL lbc_lnk( 'dom_zgr', e3t_0, 'T', 1._dp, kfillmode = jpfillcopy )
         ztopbot(:,:,1) = REAL(k_top, wp)
         ztopbot(:,:,2) = REAL(k_bot, wp)
         CALL lbc_lnk( 'dom_zgr', ztopbot, 'T', 1._wp, kfillmode = jpfillcopy )   ! do not put 0 over closed boundaries
         k_top(:,:) = NINT(ztopbot(:,:,1))
         k_bot(:,:) = NINT(ztopbot(:,:,2))
         !
      ENDIF
      !
      ! the following is mandatory
      ! make sure that closed boundaries are correctly defined in k_top that will be used to compute all mask arrays
      !
      zmsk(:,:) = 1._wp                                       ! default: no closed boundaries
      IF( .NOT. l_Iperio ) THEN                                    ! E-W closed:
         zmsk(  mi0(     1+nn_hls):mi1(     1+nn_hls),:) = 0._wp   ! first column of inner global domain at 0
         zmsk(  mi0(jpiglo-nn_hls):mi1(jpiglo-nn_hls),:) = 0._wp   ! last  column of inner global domain at 0 
      ENDIF
      IF( .NOT. l_Jperio ) THEN                                    ! S closed:
         zmsk(:,mj0(     1+nn_hls):mj1(     1+nn_hls)  ) = 0._wp   ! first   line of inner global domain at 0
      ENDIF
      IF( .NOT. ( l_Jperio .OR. l_NFold ) ) THEN                   ! N closed:
         zmsk(:,mj0(jpjglo-nn_hls):mj1(jpjglo-nn_hls)  ) = 0._wp   ! last    line of inner global domain at 0
      ENDIF
      CALL lbc_lnk( 'usrdef_zgr', zmsk, 'T', 1._wp )             ! set halos
      k_top(:,:) = k_top(:,:) * NINT( zmsk(:,:) )
      !
!!gm to be remove when removing the OLD definition of e3 scale factors so that gde3w disappears
      ! Compute gde3w_0 (vertical sum of e3w)
      gde3w_0(:,:,1) = 0.5_wp * e3w_0(:,:,1)
      DO jk = 2, jpk
         gde3w_0(:,:,jk) = gde3w_0(:,:,jk-1) + e3w_0(:,:,jk)
      END DO
      !
      ! Any closed seas (defined by closea_mask > 0 in domain_cfg file) to be filled 
      ! in at runtime if ln_closea=.false.
      IF( ln_closea ) THEN
         IF ( ln_maskcs ) THEN
            ! mask all the closed sea
            CALL clo_msk( k_top, k_bot, mask_opnsea, 'mask_opensea' )
         ELSE IF ( ln_mask_csundef ) THEN
            ! defined closed sea are kept
            ! mask all the undefined closed sea
            CALL clo_msk( k_top, k_bot, mask_csundef, 'mask_csundef' )
         END IF
      END IF
      !
      IF(lwp) THEN                     ! Control print
         WRITE(numout,*)
         WRITE(numout,*) '   Type of vertical coordinate (read in ', TRIM( cn_domcfg ), ' file or set in userdef_nam) :'
         WRITE(numout,*) '      z-coordinate - full steps      ln_zco    = ', ln_zco
         WRITE(numout,*) '      z-coordinate - partial steps   ln_zps    = ', ln_zps
         WRITE(numout,*) '      s- or hybrid z-s-coordinate    ln_sco    = ', ln_sco
         WRITE(numout,*) '      ice shelf cavities             ln_isfcav = ', ln_isfcav
      ENDIF

      ioptio = 0                       ! Check Vertical coordinate options
      IF( ln_zco      )   ioptio = ioptio + 1
      IF( ln_zps      )   ioptio = ioptio + 1
      IF( ln_sco      )   ioptio = ioptio + 1
      IF( ioptio /= 1 )   CALL ctl_stop( ' none or several vertical coordinate options used' )


      !                                ! top/bottom ocean level indices for t-, u- and v-points (f-point also for top)
      CALL zgr_top_bot( k_top, k_bot, is_mbkuvf )      ! with a minimum value set to 1
      !
      !                                ! ice shelf draft and bathymetry
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
         ikt = mikt(ji,jj)
         ikb = mbkt(ji,jj)
         bathy  (ji,jj) = gdepw_0(ji,jj,ikb+1)
         risfdep(ji,jj) = gdepw_0(ji,jj,ikt  )
      END_2D
      !
      !                                ! deepest/shallowest W level Above/Below ~10m
!!gm BUG in s-coordinate this does not work!
      zrefdep = 10._wp - 0.1_wp * MINVAL( e3w_1d )                   ! ref. depth with tolerance (10% of minimum layer thickness)
      nlb10 = MINLOC( gdepw_1d, mask = gdepw_1d > zrefdep, dim = 1 ) ! shallowest W level Below ~10m
      nla10 = nlb10 - 1                                              ! deepest    W level Above ~10m
!!gm end bug
      !
      IF( lwp )   THEN
         WRITE(numout,*) ' MIN val k_top   ', MINVAL(   k_top(:,:) ), ' MAX ', MAXVAL( k_top(:,:) )
         WRITE(numout,*) ' MIN val k_bot   ', MINVAL(   k_bot(:,:) ), ' MAX ', MAXVAL( k_bot(:,:) )
         WRITE(numout,*) ' MIN val depth t ', MINVAL( gdept_0(:,:,:) ),   &
            &                          ' w ', MINVAL( gdepw_0(:,:,:) ), '3w ', MINVAL( gde3w_0(:,:,:) )
         WRITE(numout,*) ' MIN val e3    t ', MINVAL(   e3t_0(:,:,:) ), ' f ', MINVAL(   e3f_0(:,:,:) ),  &
            &                          ' u ', MINVAL(   e3u_0(:,:,:) ), ' u ', MINVAL(   e3v_0(:,:,:) ),  &
            &                          ' uw', MINVAL(  e3uw_0(:,:,:) ), ' vw', MINVAL(  e3vw_0(:,:,:)),   &
            &                          ' w ', MINVAL(   e3w_0(:,:,:) )

         WRITE(numout,*) ' MAX val depth t ', MAXVAL( gdept_0(:,:,:) ),   &
            &                          ' w ', MAXVAL( gdepw_0(:,:,:) ), '3w ', MAXVAL( gde3w_0(:,:,:) )
         WRITE(numout,*) ' MAX val e3    t ', MAXVAL(   e3t_0(:,:,:) ), ' f ', MAXVAL(   e3f_0(:,:,:) ),  &
            &                          ' u ', MAXVAL(   e3u_0(:,:,:) ), ' u ', MAXVAL(   e3v_0(:,:,:) ),  &
            &                          ' uw', MAXVAL(  e3uw_0(:,:,:) ), ' vw', MAXVAL(  e3vw_0(:,:,:) ),  &
            &                          ' w ', MAXVAL(   e3w_0(:,:,:) )
      ENDIF
      !
   END SUBROUTINE dom_zgr


   SUBROUTINE zgr_read( ld_zco  , ld_zps  , ld_sco  , ld_isfcav,   &   ! type of vertical coordinate
      &                 pdept_1d, pdepw_1d, pe3t_1d , pe3w_1d  ,   &   ! 1D reference vertical coordinate
      &                 pdept , pdepw ,                            &   ! 3D t & w-points depth
      &                 pe3t  , pe3u  , pe3v   , pe3f ,            &   ! vertical scale factors
      &                 pe3w  , pe3uw , pe3vw         ,            &   !     -      -      -
      &                 k_top  , k_bot  ,                          &   ! top & bottom ocean level
      &                 k_mbkuvf  , k_bot_u  , k_bot_v  , k_bot_f  )   ! U/V/F points bottom levels
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE zgr_read  ***
      !!
      !! ** Purpose :   Read the vertical information in the domain configuration file
      !!
      !!----------------------------------------------------------------------
      LOGICAL                   , INTENT(out) ::   ld_zco, ld_zps, ld_sco      ! vertical coordinate flags
      LOGICAL                   , INTENT(out) ::   ld_isfcav                   ! under iceshelf cavity flag
      REAL(wp), DIMENSION(:)    , INTENT(out) ::   pdept_1d, pdepw_1d          ! 1D grid-point depth       [m]
      REAL(wp), DIMENSION(:)    , INTENT(out) ::   pe3t_1d , pe3w_1d           ! 1D vertical scale factors [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(out) ::   pdept, pdepw                ! grid-point depth          [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(out)  :: pe3u, pe3v, pe3f! vertical scale factors    [m]
      REAL(dp), DIMENSION(:,:,:), INTENT(out)  :: pe3t! vertical scale factors    [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(out) ::   pe3w , pe3uw, pe3vw         !    -       -      -
      INTEGER , DIMENSION(:,:)  , INTENT(out) ::   k_top , k_bot               ! first & last ocean level
      INTEGER                   , INTENT(out) ::   k_mbkuvf                    ! ==1 if mbku, mbkv, mbkf are in file
      INTEGER , DIMENSION(:,:)  , INTENT(out) ::   k_bot_u , k_bot_v, k_bot_f  ! bottom levels at U/V/F points
      !
      INTEGER  ::   ji,jj,jk     ! dummy loop index
      INTEGER  ::   inum, iatt
      REAL(WP) ::   z_zco, z_zps, z_sco, z_cav
      REAL(wp), DIMENSION(jpi,jpj) ::   z2d   ! 2D workspace
      CHARACTER(len=7) ::   catt   ! 'zco', 'zps, 'sco' or 'UNKNOWN'
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) '   zgr_read : read the vertical coordinates in ', TRIM( cn_domcfg ), ' file'
         WRITE(numout,*) '   ~~~~~~~~'
      ENDIF
      !
      CALL iom_open( cn_domcfg, inum )
      !
      !                          !* type of vertical coordinate
      CALL iom_getatt( inum, 'VertCoord', catt )   ! returns 'UNKNOWN' if not found
      ld_zco = catt == 'zco'          ! default = .false.
      ld_zps = catt == 'zps'          ! default = .false.
      ld_sco = catt == 'sco'          ! default = .false.
      !                          !* ocean cavities under iceshelves
      CALL iom_getatt( inum,    'IsfCav', iatt )   ! returns -999 if not found
      ld_isfcav = iatt == 1           ! default = .false.
      !
      ! ------- keep compatibility with OLD VERSION... start -------
      IF( catt == 'UNKNOWN' ) THEN
         CALL iom_get( inum,    'ln_zco', z_zco )   ;   ld_zco = z_zco /= 0._wp
         CALL iom_get( inum,    'ln_zps', z_zps )   ;   ld_zps = z_zps /= 0._wp
         CALL iom_get( inum,    'ln_sco', z_sco )   ;   ld_sco = z_sco /= 0._wp
      ENDIF
      IF( iatt == -999 ) THEN
         CALL iom_get( inum, 'ln_isfcav', z_cav )   ;   ld_isfcav = z_cav /= 0._wp
      ENDIF
      ! ------- keep compatibility with OLD VERSION... end -------
      !
      !                          !* ocean top and bottom level
      CALL iom_get( inum, jpdom_global, 'top_level'    , z2d   )   ! 1st wet T-points (ISF)
      k_top(:,:) = NINT( z2d(:,:) )
      CALL iom_get( inum, jpdom_global, 'bottom_level' , z2d   )   ! last wet T-points
      k_bot(:,:) = NINT( z2d(:,:) )
      !
      !                          !* vertical scale factors
      CALL iom_get( inum, jpdom_unknown, 'e3t_1d'  , pe3t_1d  )                     ! 1D reference coordinate
      CALL iom_get( inum, jpdom_unknown, 'e3w_1d'  , pe3w_1d  )
      !
      CALL iom_get( inum, jpdom_global, 'e3t_0'  , pe3t , cd_type = 'T', psgn = 1._dp, kfill = jpfillcopy )    ! 3D coordinate
      CALL iom_get( inum, jpdom_global, 'e3u_0'  , pe3u , cd_type = 'U', psgn = 1._wp, kfill = jpfillcopy )
      CALL iom_get( inum, jpdom_global, 'e3v_0'  , pe3v , cd_type = 'V', psgn = 1._wp, kfill = jpfillcopy )
      CALL iom_get( inum, jpdom_global, 'e3f_0'  , pe3f , cd_type = 'F', psgn = 1._wp, kfill = jpfillcopy )
      CALL iom_get( inum, jpdom_global, 'e3w_0'  , pe3w , cd_type = 'W', psgn = 1._wp, kfill = jpfillcopy )
      CALL iom_get( inum, jpdom_global, 'e3uw_0' , pe3uw, cd_type = 'U', psgn = 1._wp, kfill = jpfillcopy )
      CALL iom_get( inum, jpdom_global, 'e3vw_0' , pe3vw, cd_type = 'V', psgn = 1._wp, kfill = jpfillcopy )
      !
      !                          !* depths
      !                                   !- old depth definition (obsolescent feature)
      IF(  iom_varid( inum, 'gdept_1d', ldstop = .FALSE. ) > 0  .AND.  &
         & iom_varid( inum, 'gdepw_1d', ldstop = .FALSE. ) > 0  .AND.  &
         & iom_varid( inum, 'gdept_0' , ldstop = .FALSE. ) > 0  .AND.  &
         & iom_varid( inum, 'gdepw_0' , ldstop = .FALSE. ) > 0    ) THEN
         CALL ctl_warn( 'zgr_read : old definition of depths and scale factors used ', & 
            &           '           depths at t- and w-points read in the domain configuration file')
         CALL iom_get( inum, jpdom_unknown, 'gdept_1d', pdept_1d )   
         CALL iom_get( inum, jpdom_unknown, 'gdepw_1d', pdepw_1d )
         CALL iom_get( inum, jpdom_global , 'gdept_0' , pdept, kfill = jpfillcopy )
         CALL iom_get( inum, jpdom_global , 'gdepw_0' , pdepw, kfill = jpfillcopy )
         !
      ELSE                                !- depths computed from e3. scale factors
         CALL e3_to_depth( pe3t_1d, pe3w_1d, pdept_1d, pdepw_1d )    ! 1D reference depth
         CALL e3_to_depth( pe3t   , pe3w   , pdept   , pdepw    )    ! 3D depths
#if defined key_qco && key_isf
         DO_3D( nn_hls, nn_hls, nn_hls, nn_hls, 2, jpk )        ! vertical sum at partial cell xxxx other level  
            IF( jk == k_top(ji,jj) ) THEN                               ! first ocean point : partial cell
               pdept(ji,jj,jk) = pdepw(ji,jj,jk  ) + 0.5_wp * pe3w(ji,jj,jk)   ! = risfdep + 1/2 e3w_0(mikt)
            ELSE                                                        !  other levels
               pdept(ji,jj,jk) = pdept(ji,jj,jk-1) +          pe3w(ji,jj,jk) 
            ENDIF
         END_3D
#endif
         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) '              Reference 1D z-coordinate depth and scale factors:'
            WRITE(numout, "(9x,' level  gdept_1d  gdepw_1d  e3t_1d   e3w_1d  ')" )
            WRITE(numout, "(10x, i4, 4f9.2)" ) ( jk, pdept_1d(jk), pdepw_1d(jk), pe3t_1d(jk), pe3w_1d(jk), jk = 1, jpk )
         ENDIF
      ENDIF
      !
      IF( iom_varid( inum, 'mbku', ldstop = .FALSE. ) > 0 ) THEN
         IF(lwp) WRITE(numout,*) '          mbku, mbkv & mbkf read in ', TRIM(cn_domcfg), ' file'
         CALL iom_get( inum, jpdom_global, 'mbku', z2d, cd_type = 'U', psgn = 1._wp, kfill = jpfillcopy )
         k_bot_u(:,:) = NINT( z2d(:,:) )
         CALL iom_get( inum, jpdom_global, 'mbkv', z2d, cd_type = 'V', psgn = 1._wp, kfill = jpfillcopy )
         k_bot_v(:,:) = NINT( z2d(:,:) )
         CALL iom_get( inum, jpdom_global, 'mbkf', z2d, cd_type = 'F', psgn = 1._wp, kfill = jpfillcopy )
         k_bot_f(:,:) = NINT( z2d(:,:) )
         k_mbkuvf = 1
      ELSE
         k_mbkuvf = 0
      ENDIF
      !
      ! reference depth for negative bathy (wetting and drying only)
      IF( ll_wd )  CALL iom_get( inum,  'rn_wd_ref_depth' , ssh_ref   )
      !
      CALL iom_close( inum )
      !
   END SUBROUTINE zgr_read


   SUBROUTINE zgr_top_bot( k_top, k_bot, k_mbkuvf )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE zgr_top_bot  ***
      !!
      !! ** Purpose :   defines the vertical index of ocean bottom (mbk. arrays)
      !!
      !! ** Method  :   computes from k_top and k_bot with a minimum value of 1 over land
      !!
      !! ** Action  :   mikt, miku, mikv :   vertical indices of the shallowest 
      !!                                     ocean level at t-, u- & v-points
      !!                                     (min value = 1)
      !! ** Action  :   mbkt, mbku, mbkv :   vertical indices of the deeptest 
      !!                mbkf                 ocean level at t-, u-, v- & f-points
      !!                                     (min value = 1 over land)
      !!----------------------------------------------------------------------
      INTEGER , DIMENSION(:,:), INTENT(in) ::   k_top, k_bot   ! top & bottom ocean level indices
      INTEGER                 , INTENT(in) ::   k_mbkuvf       ! flag to recompute mbku, mbkv, mbkf
      !
      INTEGER ::   ji, jj   ! dummy loop indices
      REAL(wp), DIMENSION(jpi,jpj) ::   zk   ! workspace
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '    zgr_top_bot : ocean top and bottom k-index of T-, U-, V- and W-levels '
      IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~'
      !
      mikt(:,:) = MAX( k_top(:,:) , 1 )    ! top    ocean k-index of T-level (=1 over land)
      !
      mbkt(:,:) = MAX( k_bot(:,:) , 1 )    ! bottom ocean k-index of T-level (=1 over land)
 
      !                                    ! N.B.  top     k-index of W-level = mikt
      !                                    !       bottom  k-index of W-level = mbkt+1
      DO_2D( 0, 0, 0, 0 )
         miku(ji,jj) = MAX(  mikt(ji+1,jj  ) , mikt(ji,jj)  )
         mikv(ji,jj) = MAX(  mikt(ji  ,jj+1) , mikt(ji,jj)  )
         mikf(ji,jj) = MAX(  mikt(ji  ,jj+1) , mikt(ji,jj), mikt(ji+1,jj  ), mikt(ji+1,jj+1)  )
      END_2D

      IF ( k_mbkuvf==0 ) THEN
         IF(lwp) WRITE(numout,*) '         mbku, mbkv, mbkf computed from mbkt'
         DO_2D( 0, 0, 0, 0 )
            mbku(ji,jj) = MIN(  mbkt(ji+1,jj  ) , mbkt(ji,jj)  )
            mbkv(ji,jj) = MIN(  mbkt(ji  ,jj+1) , mbkt(ji,jj)  )
            mbkf(ji,jj) = MIN(  mbkt(ji  ,jj+1) , mbkt(ji,jj), mbkt(ji+1,jj  ), mbkt(ji+1,jj+1)  )
         END_2D
      ELSE
         IF(lwp) WRITE(numout,*) '         mbku, mbkv, mbkf read from file'
         ! Use mbku, mbkv, mbkf from file
         ! Ensure these are lower than expected bottom level deduced from mbkt
         DO_2D( 0, 0, 0, 0 )
            mbku(ji,jj) = MIN(  mbku(ji,jj), mbkt(ji+1,jj  ) , mbkt(ji,jj)  )
            mbkv(ji,jj) = MIN(  mbkv(ji,jj), mbkt(ji  ,jj+1) , mbkt(ji,jj)  )
            mbkf(ji,jj) = MIN(  mbkf(ji,jj), mbkt(ji  ,jj+1) , mbkt(ji,jj), mbkt(ji+1,jj  ), mbkt(ji+1,jj+1)  )
         END_2D
      ENDIF
      ! convert into REAL to use lbc_lnk ; impose a min value of 1 as a zero can be set in lbclnk 
      DO_2D( 0, 0, 0, 0 )
         zk(ji,jj) = REAL( miku(ji,jj), wp )
      END_2D
      CALL lbc_lnk( 'domzgr', zk, 'U', 1.0_wp )
      miku(:,:) = MAX( NINT( zk(:,:) ), 1 )

      DO_2D( 0, 0, 0, 0 )
         zk(ji,jj) = REAL( mikv(ji,jj), wp )
      END_2D
      CALL lbc_lnk( 'domzgr', zk, 'V', 1.0_wp )
      mikv(:,:) = MAX( NINT( zk(:,:) ), 1 )
      
      DO_2D( 0, 0, 0, 0 )
         zk(ji,jj) = REAL( mikf(ji,jj), wp )
      END_2D
      CALL lbc_lnk( 'domzgr', zk, 'F', 1.0_wp )
      mikf(:,:) = MAX( NINT( zk(:,:) ), 1 )
      !
      DO_2D( 0, 0, 0, 0 )
         zk(ji,jj) = REAL( mbku(ji,jj), wp )
      END_2D
      CALL lbc_lnk( 'domzgr', zk, 'U', 1.0_wp )
      mbku(:,:) = MAX( NINT( zk(:,:) ), 1 )
      
      DO_2D( 0, 0, 0, 0 )
         zk(ji,jj) = REAL( mbkv(ji,jj), wp )
      END_2D
      CALL lbc_lnk( 'domzgr', zk, 'V', 1.0_wp )
      mbkv(:,:) = MAX( NINT( zk(:,:) ), 1 )

      DO_2D( 0, 0, 0, 0 )
         zk(ji,jj) = REAL( mbkf(ji,jj), wp )
      END_2D
      CALL lbc_lnk( 'domzgr', zk, 'F', 1.0_wp )
      mbkf(:,:) = MAX( NINT( zk(:,:) ), 1 )
      !
   END SUBROUTINE zgr_top_bot

   !!======================================================================
END MODULE domzgr
