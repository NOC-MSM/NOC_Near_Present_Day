MODULE dommsk
   !!======================================================================
   !!                       ***  MODULE dommsk   ***
   !! Ocean initialization : domain land/sea mask 
   !!======================================================================
   !! History :  OPA  ! 1987-07  (G. Madec)  Original code
   !!            6.0  ! 1993-03  (M. Guyon)  symetrical conditions (M. Guyon)
   !!            7.0  ! 1996-01  (G. Madec)  suppression of common work arrays
   !!             -   ! 1996-05  (G. Madec)  mask computed from tmask
   !!            8.0  ! 1997-02  (G. Madec)  mesh information put in domhgr.F
   !!            8.1  ! 1997-07  (G. Madec)  modification of kbat and fmask
   !!             -   ! 1998-05  (G. Roullet)  free surface
   !!            8.2  ! 2000-03  (G. Madec)  no slip accurate
   !!             -   ! 2001-09  (J.-M. Molines)  Open boundaries
   !!   NEMO     1.0  ! 2002-08  (G. Madec)  F90: Free form and module
   !!             -   ! 2005-11  (V. Garnier) Surface pressure gradient organization
   !!            3.2  ! 2009-07  (R. Benshila) Suppression of rigid-lid option
   !!            3.6  ! 2015-05  (P. Mathiot) ISF: add wmask,wumask and wvmask
   !!            4.0  ! 2016-06  (G. Madec, S. Flavoni)  domain configuration / user defined interface
   !!            4.x  ! 2020-02  (G. Madec, S. Techene) introduce ssh to h0 ratio
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dom_msk       : compute land/ocean mask
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE domutl         ! 
   USE usrdef_fmask   ! user defined fmask
   USE bdy_oce        ! open boundary
   !
   USE in_out_manager ! I/O manager
   USE iom            ! IOM library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! Massively Parallel Processing library
   USE iom             ! For shlat2d
   USE fldread         ! for sn_shlat2d


   IMPLICIT NONE
   PRIVATE

   PUBLIC   dom_msk    ! routine called by inidom.F90

   !                            !!* Namelist namlbc : lateral boundary condition *
   REAL(wp)        :: rn_shlat   ! type of lateral boundary condition on velocity
   LOGICAL, PUBLIC :: ln_vorlat  !  consistency of vorticity boundary condition 
   !                                            with analytical eqs.

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: dommsk.F90 15556 2021-11-29 15:23:06Z jchanut $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dom_msk( k_top, k_bot )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE dom_msk  ***
      !!
      !! ** Purpose :   Compute land/ocean mask arrays at tracer points, hori-
      !!      zontal velocity points (u & v), vorticity points (f) points.
      !!
      !! ** Method  :   The ocean/land mask  at t-point is deduced from ko_top 
      !!      and ko_bot, the indices of the fist and last ocean t-levels which 
      !!      are either defined in usrdef_zgr or read in zgr_read.
      !!                The velocity masks (umask, vmask, wmask, wumask, wvmask) 
      !!      are deduced from a product of the two neighboring tmask.
      !!                The vorticity mask (fmask) is deduced from tmask taking
      !!      into account the choice of lateral boundary condition (rn_shlat) :
      !!         rn_shlat = 0, free slip  (no shear along the coast)
      !!         rn_shlat = 2, no slip  (specified zero velocity at the coast)
      !!         0 < rn_shlat < 2, partial slip   | non-linear velocity profile
      !!         2 < rn_shlat, strong slip        | in the lateral boundary layer
      !!
      !! ** Action :   tmask, umask, vmask, wmask, wumask, wvmask : land/ocean mask 
      !!                         at t-, u-, v- w, wu-, and wv-points (=0. or 1.)
      !!               fmask   : land/ocean mask at f-point (=0., or =1., or 
      !!                         =rn_shlat along lateral boundaries)
      !!               ssmask , ssumask, ssvmask, ssfmask : 2D ocean mask, i.e. at least 1 wet cell in the vertical
      !!               tmask_i : ssmask * ( excludes halo+duplicated points (NP folding) )
      !!----------------------------------------------------------------------
      INTEGER, DIMENSION(:,:), INTENT(in) ::   k_top, k_bot   ! first and last ocean level
      !
      INTEGER  ::   ji, jj, jk     ! dummy loop indices
      INTEGER  ::   iktop, ikbot   !   -       -
      INTEGER  ::   ios, inum
      !!
      REAL(wp) :: zshlat           !: locally modified shlat for some strait
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zshlat2d
      LOGICAL                         :: ln_shlat2d
      CHARACTER(len = 256)            :: cn_shlat2d_file, cn_shlat2d_var  
      NAMELIST/namlbc/ rn_shlat, ln_vorlat, ln_shlat2d, cn_shlat2d_file, cn_shlat2d_var
      NAMELIST/nambdy/ ln_bdy ,nb_bdy, ln_coords_file, cn_coords_file,         &
         &             ln_mask_file, cn_mask_file, cn_dyn2d, nn_dyn2d_dta,     &
         &             cn_dyn3d, nn_dyn3d_dta, cn_tra, nn_tra_dta,             &
         &             ln_tra_dmp, ln_dyn3d_dmp, rn_time_dmp, rn_time_dmp_out, &
         &             cn_ice, nn_ice_dta,                                     &
         &             ln_vol, nn_volctl, nn_rimwidth
      !!---------------------------------------------------------------------
      !
      READ  ( numnam_ref, namlbc, IOSTAT = ios, ERR = 901 )
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namlbc in reference namelist' )
      READ  ( numnam_cfg, namlbc, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namlbc in configuration namelist' )
      IF(lwm) WRITE ( numond, namlbc )
      
      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'dommsk : ocean mask '
         WRITE(numout,*) '~~~~~~'
         WRITE(numout,*) '   Namelist namlbc'
         WRITE(numout,*) '      lateral momentum boundary cond.    rn_shlat  = ',rn_shlat
         WRITE(numout,*) '      consistency with analytical form   ln_vorlat = ',ln_vorlat 
      ENDIF
      !
      IF(lwp) WRITE(numout,*)
      IF ( ln_shlat2d ) THEN
         IF(lwp) WRITE(numout,*) '         READ shlat as a 2D coefficient in a file '
         ALLOCATE( zshlat2d(jpi,jpj) )
         CALL iom_open(TRIM(cn_shlat2d_file), inum)
         CALL iom_get (inum, jpdom_global, TRIM(cn_shlat2d_var), zshlat2d, 1) !
         CALL iom_close(inum)
       ELSE
         IF     (      rn_shlat == 0.               ) THEN   ;   IF(lwp) WRITE(numout,*) '   ==>>>   ocean lateral  free-slip'
         ELSEIF (      rn_shlat == 2.               ) THEN   ;   IF(lwp) WRITE(numout,*) '   ==>>>   ocean lateral  no-slip'
         ELSEIF ( 0. < rn_shlat .AND. rn_shlat < 2. ) THEN   ;   IF(lwp) WRITE(numout,*) '   ==>>>   ocean lateral  partial-slip'
         ELSEIF ( 2. < rn_shlat                     ) THEN   ;   IF(lwp) WRITE(numout,*) '   ==>>>   ocean lateral  strong-slip'
         ELSE
            CALL ctl_stop( 'dom_msk: wrong value for rn_shlat (i.e. a negalive value). We stop.' )
         ENDIF
      ENDIF

      !  Ocean/land mask at t-point  (computed from ko_top and ko_bot)
      ! ----------------------------
      !
      tmask(:,:,:) = 0._wp
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
         iktop = k_top(ji,jj)
         ikbot = k_bot(ji,jj)
         IF( iktop /= 0 ) THEN       ! water in the column
            tmask(ji,jj,iktop:ikbot) = 1._wp
         ENDIF
      END_2D
      !
      ! Mask corrections for bdy (read in mppini2)
      READ  ( numnam_ref, nambdy, IOSTAT = ios, ERR = 903)
903   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nambdy in reference namelist' )
      READ  ( numnam_cfg, nambdy, IOSTAT = ios, ERR = 904 )
904   IF( ios >  0 )   CALL ctl_nam ( ios , 'nambdy in configuration namelist' )
      ! ------------------------
      IF ( ln_bdy .AND. ln_mask_file ) THEN
         CALL iom_open( cn_mask_file, inum )
         CALL iom_get ( inum, jpdom_global, 'bdy_msk', bdytmask(:,:) )
         CALL iom_close( inum )
         DO_3D( 1, 1, 1, 1, 1, jpkm1 )
            tmask(ji,jj,jk) = tmask(ji,jj,jk) * bdytmask(ji,jj)
         END_3D
      ENDIF
         
      !  Ocean/land mask at u-, v-, and f-points   (computed from tmask)
      ! ----------------------------------------
      ! NB: at this point, fmask is designed for free slip lateral boundary condition
      DO_3D( 0, 0, 0, 0, 1, jpk )
         umask(ji,jj,jk) = tmask(ji,jj  ,jk) * tmask(ji+1,jj  ,jk)
         vmask(ji,jj,jk) = tmask(ji,jj  ,jk) * tmask(ji  ,jj+1,jk)
         fmask(ji,jj,jk) = tmask(ji,jj  ,jk) * tmask(ji+1,jj  ,jk)   &
            &            * tmask(ji,jj+1,jk) * tmask(ji+1,jj+1,jk)
      END_3D
      !
      ! In case of a coarsened grid, account her for possibly aditionnal  
      ! masked points; these have been read in the mesh file and stored in mbku, mbkv, mbkf
      DO_2D( 0, 0, 0, 0 )
         ! Ugly patch to accomodate batropic case (jpk=2)
         ! but with possible masked faces in the coarsening case
         ! A better way would be to keep track of mbku/v/f=0 until here
         ! but these have been assigned a minimum of 1. TBC
         IF (jpk>2) THEN   
            IF (mbku(ji,jj)==1) umask(ji,jj,:) = 0._wp
            IF (mbkv(ji,jj)==1) vmask(ji,jj,:) = 0._wp
            IF (mbkf(ji,jj)==1) fmask(ji,jj,:) = 0._wp
         ENDIF
         IF ( MAXVAL(umask(ji,jj,:))/=0._wp )  umask(ji,jj,mbku(ji,jj)+1:jpk) = 0._wp
         IF ( MAXVAL(vmask(ji,jj,:))/=0._wp )  vmask(ji,jj,mbkv(ji,jj)+1:jpk) = 0._wp
         IF ( MAXVAL(fmask(ji,jj,:))/=0._wp )  fmask(ji,jj,mbkf(ji,jj)+1:jpk) = 0._wp
      END_2D
      CALL lbc_lnk( 'dommsk', umask, 'U', 1.0_wp, vmask, 'V', 1.0_wp, fmask, 'F', 1.0_wp )      ! Lateral boundary conditions
 
      ! Ocean/land mask at wu-, wv- and w points    (computed from tmask)
      !-----------------------------------------
      wmask (:,:,1) = tmask(:,:,1)     ! surface
      wumask(:,:,1) = umask(:,:,1)
      wvmask(:,:,1) = vmask(:,:,1)
      DO jk = 2, jpk                   ! interior values
         wmask (:,:,jk) = tmask(:,:,jk) * tmask(:,:,jk-1)
         wumask(:,:,jk) = umask(:,:,jk) * umask(:,:,jk-1)   
         wvmask(:,:,jk) = vmask(:,:,jk) * vmask(:,:,jk-1)
      END DO

      ! Ocean/land column mask at t-, u-, and v-points   (i.e. at least 1 wet cell in the vertical)
      ! ----------------------------------------------
      ssmask (:,:) = MAXVAL( tmask(:,:,:), DIM=3 )
      ssumask(:,:) = MAXVAL( umask(:,:,:), DIM=3 )
      ssvmask(:,:) = MAXVAL( vmask(:,:,:), DIM=3 )
      ssfmask(:,:) = MAXVAL( fmask(:,:,:), DIM=3 )
      IF( lk_SWE ) THEN      ! Shallow Water Eq. case : redefine ssfmask
         DO_2D( 0, 0, 0, 0 )
            ssfmask(ji,jj) = MAX(  ssmask(ji,jj+1), ssmask(ji+1,jj+1),  & 
               &                   ssmask(ji,jj  ), ssmask(ji+1,jj  )   )
         END_2D
         CALL lbc_lnk( 'dommsk', ssfmask, 'F', 1.0_wp )
      ENDIF
      fe3mask(:,:,:) = fmask(:,:,:)

      ! Interior domain mask  (used for global sum) : 2D ocean mask x (halo+duplicated points) mask 
      ! --------------------
      !
      CALL dom_uniq( tmask_i, 'T' )
      tmask_i(:,:) = ssmask(:,:) * tmask_i(:,:)

      ! Lateral boundary conditions on velocity (modify fmask)
      ! ---------------------------------------  
      IF( rn_shlat /= 0._wp .or. ln_shlat2d ) THEN      ! Not free-slip lateral boundary condition
         !
         IF (  ln_shlat2d ) THEN
            DO_3D( 0, 0, 0, 0, 1, jpk )
               IF( fmask(ji,jj,jk) == 0._wp ) THEN
                  fmask(ji,jj,jk) = zshlat2d(ji,jj) * MIN( 1._wp , MAX( umask(ji,jj,jk), umask(ji,jj+1,jk), &
                     &                                           vmask(ji,jj,jk), vmask(ji+1,jj,jk) ) )
               ENDIF
            END_3D
         ELSE
            DO_3D( 0, 0, 0, 0, 1, jpk )
               IF( fmask(ji,jj,jk) == 0._wp ) THEN
                  fmask(ji,jj,jk) = rn_shlat * MIN( 1._wp , MAX( umask(ji,jj,jk), umask(ji,jj+1,jk), &
                     &                                           vmask(ji,jj,jk), vmask(ji+1,jj,jk) ) )
               ENDIF
            END_3D
         END IF
         !
         IF( ln_shlat2d ) DEALLOCATE( zshlat2d )
         !
         CALL lbc_lnk( 'dommsk', fmask, 'F', 1._wp )      ! Lateral boundary conditions on fmask
         !
         ! CAUTION : The fmask may be further modified in dyn_vor_init ( dynvor.F90 ) depending on ln_vorlat
         !
      ENDIF
      
      ! User defined alteration of fmask (use to reduce ocean transport in specified straits)
      ! -------------------------------- 
      !
      IF ( .not. ln_shlat2d ) THEN      
         CALL usr_def_fmask( cn_cfg, nn_cfg, fmask )
      ENDIF
      !
#if defined key_agrif
      ! Reset masks defining updated points over parent grids
      !  = 1 : updated point from child(s)
      !  = 0 : point not updated
      ! 
      tmask_upd(:,:) = 0._wp
      umask_upd(:,:) = 0._wp
      vmask_upd(:,:) = 0._wp
      !
      ! Reset mask defining actual computationnal domain
      ! e.g. excluding ghosts and updated cells.
      tmask_agrif(:,:) = 1._wp
#endif     
      !
   END SUBROUTINE dom_msk
   
   !!======================================================================
END MODULE dommsk
