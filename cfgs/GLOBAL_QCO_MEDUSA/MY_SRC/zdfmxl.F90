MODULE zdfmxl
   !!======================================================================
   !!                       ***  MODULE  zdfmxl  ***
   !! Ocean physics: mixed layer depth 
   !!======================================================================
   !! History :  1.0  ! 2003-08  (G. Madec)  original code
   !!            3.2  ! 2009-07  (S. Masson, G. Madec)  IOM + merge of DO-loop
   !!            3.7  ! 2012-03  (G. Madec)  make public the density criteria for trdmxl 
   !!             -   ! 2014-02  (F. Roquet)  mixed layer depth calculated using N2 instead of rhop 
   !!----------------------------------------------------------------------
   !!   zdf_mxl      : Compute the turbocline and mixed layer depths.
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE isf_oce        ! ice shelf
   USE dom_oce        ! ocean space and time domain variables
   USE trc_oce  , ONLY: l_offline         ! ocean space and time domain variables
   USE zdf_oce        ! ocean vertical physics
   USE eosbn2         ! for zdf_mxl_zint
   !
   USE in_out_manager ! I/O manager
   USE prtctl         ! Print control
   USE phycst         ! physical constants
   USE iom            ! I/O library
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   zdf_mxl, zdf_mxl_turb, zdf_mxl_alloc   ! called by zdfphy.F90
   PUBLIC   zdf_mxl_zint                           ! called by diahth.F90

   INTEGER , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   nmln    !: number of level in the mixed layer (used by LDF, ZDF, TRD, TOP)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hmld    !: mixing layer depth (turbocline)      [m]   (used by TOP)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hmlp    !: mixed layer depth  (rho=rho0+zdcrit) [m]   (used by LDF)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hmlpt   !: depth of the last T-point inside the mixed layer [m] (used by LDF)
   REAL(wp), PUBLIC, ALLOCATABLE,       DIMENSION(:,:) ::   hmld_zint  !: vertically-interpolated mixed layer depth   [m]
   REAL(wp), PUBLIC, ALLOCATABLE,       DIMENSION(:,:) ::   htc_mld    ! Heat content of hmld_zint
   LOGICAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)    :: ll_found   ! Is T_b to be found by interpolation ?
   LOGICAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  :: ll_belowml ! Flag points below mixed layer when ll_found=F

   REAL(wp), PUBLIC ::   rho_c = 0.01_wp    !: density criterion for mixed layer depth
   REAL(wp), PUBLIC ::   avt_c = 5.e-4_wp   ! Kz criterion for the turbocline depth

   TYPE, PUBLIC :: MXL_ZINT   !: Structure for MLD defs
      INTEGER   :: mld_type   ! mixed layer type     
      REAL(wp)  :: zref       ! depth of initial T_ref
      REAL(wp)  :: dT_crit    ! Critical temp diff
      REAL(wp)  :: iso_frac   ! Fraction of rn_dT_crit 
   END TYPE MXL_ZINT

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "single_precision_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: zdfmxl.F90 15249 2021-09-13 09:59:09Z hadcv $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION zdf_mxl_alloc()
      !!----------------------------------------------------------------------
      !!               ***  FUNCTION zdf_mxl_alloc  ***
      !!----------------------------------------------------------------------
      zdf_mxl_alloc = 0      ! set to zero if no array to be allocated
      IF( .NOT. ALLOCATED( nmln ) ) THEN
         ALLOCATE( nmln(jpi,jpj), hmld(jpi,jpj), hmlp(jpi,jpj), hmlpt(jpi,jpj), hmld_zint(jpi,jpj),     &
   &          htc_mld(jpi,jpj), ll_found(jpi,jpj), ll_belowml(jpi,jpj,jpk), STAT= zdf_mxl_alloc )
         !
         CALL mpp_sum ( 'zdfmxl', zdf_mxl_alloc )
         IF( zdf_mxl_alloc /= 0 )   CALL ctl_stop( 'STOP', 'zdf_mxl_alloc: failed to allocate arrays.' )
         !
      ENDIF
   END FUNCTION zdf_mxl_alloc


   SUBROUTINE zdf_mxl( kt, Kmm )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdfmxl  ***
      !!                   
      !! ** Purpose :   Compute the mixed layer depth with density criteria.
      !!
      !! ** Method  :   The mixed layer depth is the shallowest W depth with 
      !!      the density of the corresponding T point (just bellow) bellow a
      !!      given value defined locally as rho(10m) + rho_c
      !!
      !! ** Action  :   nmln, hmlp, hmlpt
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      INTEGER, INTENT(in) ::   Kmm  ! ocean time level index
      !
      INTEGER  ::   ji, jj, jk      ! dummy loop indices
      INTEGER  ::   iik, ikt        ! local integer
      REAL(wp) ::   zN2_c           ! local scalar
      !!----------------------------------------------------------------------
      !
      IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                       ! Do only on the first tile
         IF( kt == nit000 ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'zdf_mxl : mixed layer depth'
            IF(lwp) WRITE(numout,*) '~~~~~~~ '
            !                             ! allocate zdfmxl arrays
            IF( zdf_mxl_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'zdf_mxl : unable to allocate arrays' )
         ENDIF
      ENDIF
      !
      ! w-level of the mixing and mixed layers
      DO_2D_OVR( nn_hls, nn_hls, nn_hls, nn_hls )
         nmln(ji,jj)  = nlb10                  ! Initialization to the number of w ocean point
         hmlp(ji,jj)  = 0._wp                  ! here hmlp used as a dummy variable, integrating vertically N^2
      END_2D
      zN2_c = grav * rho_c * r1_rho0      ! convert density criteria into N^2 criteria
      DO_3D_OVR( nn_hls, nn_hls, nn_hls, nn_hls, nlb10, jpkm1 )   ! Mixed layer level: w-level
         ikt = mbkt(ji,jj)
         hmlp(ji,jj) =   &
            & hmlp(ji,jj) + MAX( rn2b(ji,jj,jk) , 0._wp ) * e3w(ji,jj,jk,Kmm)
         IF( hmlp(ji,jj) < zN2_c )   nmln(ji,jj) = MIN( jk , ikt ) + 1   ! Mixed layer level
      END_3D
      ! depth of the mixed layer
      DO_2D_OVR( nn_hls, nn_hls, nn_hls, nn_hls )
         iik = nmln(ji,jj)
         hmlp (ji,jj) = gdepw(ji,jj,iik  ,Kmm) * ssmask(ji,jj)    ! Mixed layer depth
         hmlpt(ji,jj) = gdept(ji,jj,iik-1,Kmm) * ssmask(ji,jj)    ! depth of the last T-point inside the mixed layer
      END_2D
      !
      IF( .NOT.l_offline .AND. iom_use("mldr10_1") ) THEN
         IF( .NOT. l_istiled .OR. ntile == nijtile ) THEN         ! Do only on the last tile
            IF( ln_isfcav ) THEN  ;  CALL iom_put( "mldr10_1", hmlp - risfdep)   ! mixed layer thickness
            ELSE                  ;  CALL iom_put( "mldr10_1", hmlp )            ! mixed layer depth
            END IF
         ENDIF
      ENDIF
      !
      IF(sn_cfctl%l_prtctl)   CALL prt_ctl( tab2d_1=REAL(nmln,dp), clinfo1=' nmln : ', tab2d_2=CASTDP(hmlp), clinfo2=' hmlp : ' )
      !
   END SUBROUTINE zdf_mxl

   SUBROUTINE zdf_mxl_zint_mld( sf , Kmm) 
      !!---------------------------------------------------------------------------------- 
      !!                    ***  ROUTINE zdf_mxl_zint_mld  *** 
      !                                                                        
      !   Calculate vertically-interpolated mixed layer depth diagnostic. 
      !            
      !   This routine can calculate the mixed layer depth diagnostic suggested by
      !   Kara et al, 2000, JGR, 105, 16803, but is more general and can calculate
      !   vertically-interpolated mixed-layer depth diagnostics with other parameter
      !   settings set in the namzdf_mldzint namelist.  
      ! 
      !   If mld_type=1 the mixed layer depth is calculated as the depth at which the  
      !   density has increased by an amount equivalent to a temperature difference of  
      !   0.8C at the surface. 
      ! 
      !   For other values of mld_type the mixed layer is calculated as the depth at  
      !   which the temperature differs by 0.8C from the surface temperature.  
      !                                                                        
      !   David Acreman, Daley Calvert                                      
      ! 
      !!----------------------------------------------------------------------------------- 

      TYPE(MXL_ZINT), INTENT(in)  :: sf
      INTEGER, INTENT(in) ::   Kmm  ! ocean time level index

      ! Diagnostic criteria
      INTEGER   :: nn_mld_type   ! mixed layer type     
      REAL(wp)  :: rn_zref       ! depth of initial T_ref
      REAL(wp)  :: rn_dT_crit    ! Critical temp diff
      REAL(wp)  :: rn_iso_frac   ! Fraction of rn_dT_crit used

      ! Local variables
      REAL(wp), PARAMETER :: zepsilon = 1.e-30          ! local small value
      INTEGER, DIMENSION(jpi,jpj) :: ikmt          ! number of active tracer levels 
      INTEGER, DIMENSION(jpi,jpj) :: ik_ref        ! index of reference level 
      INTEGER, DIMENSION(jpi,jpj) :: ik_iso        ! index of last uniform temp level 
      REAL, DIMENSION(jpi,jpj,jpk)  :: zT            ! Temperature or density 
      REAL, DIMENSION(jpi,jpj)    :: ppzdep        ! depth for use in calculating d(rho) 
      REAL, DIMENSION(jpi,jpj)    :: zT_ref        ! reference temperature 
      REAL    :: zT_b                                   ! base temperature 
      REAL, DIMENSION(jpi,jpj,jpk)  :: zdTdz         ! gradient of zT 
      REAL, DIMENSION(jpi,jpj,jpk)  :: zmoddT        ! Absolute temperature difference 
      REAL    :: zdz                                    ! depth difference 
      REAL    :: zdT                                    ! temperature difference 
      REAL, DIMENSION(jpi,jpj)    :: zdelta_T      ! difference critereon 
      REAL, DIMENSION(jpi,jpj)    :: zRHO1, zRHO2  ! Densities 
      INTEGER :: ji, jj, jk                             ! loop counter 

      !!------------------------------------------------------------------------------------- 
      !  
      ! Unpack structure
      nn_mld_type = sf%mld_type
      rn_zref     = sf%zref
      rn_dT_crit  = sf%dT_crit
      rn_iso_frac = sf%iso_frac

      ! Set the mixed layer depth criterion at each grid point 
      IF( nn_mld_type == 0 ) THEN
         zdelta_T(:,:) = rn_dT_crit
         zT(:,:,:) = rhop(:,:,:)
      ELSE IF( nn_mld_type == 1 ) THEN
         ppzdep(:,:)=0.0 
         call eos ( ts(:,:,1,:,Kmm), ppzdep(:,:), zRHO1(:,:) ) 
! Use zT temporarily as a copy of tsn with rn_dT_crit added to SST 
! [assumes number of tracers less than number of vertical levels] 
         zT(:,:,1:jpts)=ts(:,:,1,1:jpts,Kmm) 
         zT(:,:,jp_tem)=zT(:,:,1)+rn_dT_crit 
         CALL eos( zT(:,:,1:jpts), ppzdep(:,:), zRHO2(:,:) ) 
         zdelta_T(:,:) = abs( zRHO1(:,:) - zRHO2(:,:) ) * rho0 
         ! RHO from eos (2d version) doesn't calculate north or east halo: 
         CALL lbc_lnk( 'zdfmxl', zdelta_T, 'T', 1. ) 
         zT(:,:,:) = rhop(:,:,:) 
      ELSE 
         zdelta_T(:,:) = rn_dT_crit                      
         zT(:,:,:) = ts(:,:,:,jp_tem,Kmm)                           
      END IF 

      ! Calculate the gradient of zT and absolute difference for use later 
      DO jk = 1 ,jpk-2 
         zdTdz(:,:,jk)  =    ( zT(:,:,jk+1) - zT(:,:,jk) ) / e3w(:,:,jk+1,Kmm) 
         zmoddT(:,:,jk) = abs( zT(:,:,jk+1) - zT(:,:,jk) ) 
      END DO 

      ! Find density/temperature at the reference level (Kara et al use 10m).          
      ! ik_ref is the index of the box centre immediately above or at the reference level 
      ! Find rn_zref in the array of model level depths and find the ref    
      ! density/temperature by linear interpolation.                                   
      DO jk = jpkm1, 2, -1 
         WHERE ( gdept(:,:,jk,Kmm) > rn_zref ) 
           ik_ref(:,:) = jk - 1 
           zT_ref(:,:) = zT(:,:,jk-1) + zdTdz(:,:,jk-1) * ( rn_zref - gdept(:,:,jk-1,Kmm) ) 
         END WHERE 
      END DO 

      ! If the first grid box centre is below the reference level then use the 
      ! top model level to get zT_ref 
      WHERE ( gdept(:,:,1,Kmm) > rn_zref )  
         zT_ref = zT(:,:,1) 
         ik_ref = 1 
      END WHERE 

      ! The number of active tracer levels is 1 less than the number of active w levels 
      ikmt(:,:) = mbkt(:,:) - 1 

      ! Initialize / reset
      ll_found(:,:) = .false.

      IF ( rn_iso_frac - zepsilon > 0. ) THEN
         ! Search for a uniform density/temperature region where adjacent levels          
         ! differ by less than rn_iso_frac * deltaT.                                      
         ! ik_iso is the index of the last level in the uniform layer  
         ! ll_found indicates whether the mixed layer depth can be found by interpolation 
         ik_iso(:,:)   = ik_ref(:,:) 
         DO jj = 1, jpj      ! Changed from nlcj
            DO ji = 1, jpi   ! Changed from nlci
!CDIR NOVECTOR 
               DO jk = ik_ref(ji,jj), ikmt(ji,jj)-1 
                  IF ( zmoddT(ji,jj,jk) > ( rn_iso_frac * zdelta_T(ji,jj) ) ) THEN 
                     ik_iso(ji,jj)   = jk 
                     ll_found(ji,jj) = ( zmoddT(ji,jj,jk) > zdelta_T(ji,jj) ) 
                     EXIT 
                  END IF 
               END DO 
            END DO 
         END DO 

         ! Use linear interpolation to find depth of mixed layer base where possible 
         hmld_zint(:,:) = rn_zref 
         DO jj = 1, jpj 
            DO ji = 1, jpi 
               IF (ll_found(ji,jj) .and. tmask(ji,jj,1) == 1.0) THEN 
                  zdz =  abs( zdelta_T(ji,jj) / zdTdz(ji,jj,ik_iso(ji,jj)) ) 
                  hmld_zint(ji,jj) = gdept(ji,jj,ik_iso(ji,jj),Kmm) + zdz 
               END IF 
            END DO 
         END DO 
      END IF

      ! If ll_found = .false. then calculate MLD using difference of zdelta_T    
      ! from the reference density/temperature 
 
! Prevent this section from working on land points 
      WHERE ( tmask(:,:,1) /= 1.0 ) 
         ll_found = .true. 
      END WHERE 
 
      DO jk=1, jpk 
         ll_belowml(:,:,jk) = abs( zT(:,:,jk) - zT_ref(:,:) ) >= zdelta_T(:,:)  
      END DO 
 
! Set default value where interpolation cannot be used (ll_found=false)  
      DO jj = 1, jpj 
         DO ji = 1, jpi 
            IF ( .not. ll_found(ji,jj) )  hmld_zint(ji,jj) = gdept(ji,jj,ikmt(ji,jj),Kmm) 
         END DO 
      END DO 

      DO jj = 1, jpj 
         DO ji = 1, jpi 
!CDIR NOVECTOR 
            DO jk = ik_ref(ji,jj)+1, ikmt(ji,jj) 
               IF ( ll_found(ji,jj) ) EXIT 
               IF ( ll_belowml(ji,jj,jk) ) THEN                
                  zT_b = zT_ref(ji,jj) + zdelta_T(ji,jj) * SIGN(1.0, zdTdz(ji,jj,jk-1) ) 
                  zdT  = zT_b - zT(ji,jj,jk-1)                                      
                  zdz  = zdT / zdTdz(ji,jj,jk-1)                                       
                  hmld_zint(ji,jj) = gdept(ji,jj,jk-1,Kmm) + zdz 
                  EXIT                                                   
               END IF 
            END DO 
         END DO 
      END DO 

      hmld_zint(:,:) = hmld_zint(:,:)*tmask(:,:,1) 
      !  
   END SUBROUTINE zdf_mxl_zint_mld

   SUBROUTINE zdf_mxl_zint_htc( kt , Kmm)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdf_mxl_zint_htc  ***
      !! 
      !! ** Purpose :   
      !!
      !! ** Method  :   
      !!----------------------------------------------------------------------

      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      INTEGER, INTENT(in) ::   Kmm  ! ocean time level index

      INTEGER :: ji, jj, jk
      INTEGER :: ikmax
      REAL(wp) :: zc, zcoef
      !
      INTEGER,  ALLOCATABLE, DIMENSION(:,:) ::   ilevel
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   zthick_0, zthick

      !!----------------------------------------------------------------------

      IF( .NOT. ALLOCATED(ilevel) ) THEN
         ALLOCATE( ilevel(jpi,jpj), zthick_0(jpi,jpj), &
         &         zthick(jpi,jpj), STAT=ji )
         IF( lk_mpp  )   CALL mpp_sum( 'zdfmxl', ji )
         IF( ji /= 0 )   CALL ctl_stop( 'STOP', 'zdf_mxl_zint_htc : unable to allocate arrays' )
      ENDIF

      ! Find last whole model T level above the MLD
      ilevel(:,:)   = 0
      zthick_0(:,:) = 0._wp

      DO jk = 1, jpkm1  
         DO jj = 1, jpj
            DO ji = 1, jpi                    
               zthick_0(ji,jj) = zthick_0(ji,jj) + e3t(ji,jj,jk,Kmm)
               IF( zthick_0(ji,jj) < hmld_zint(ji,jj) )   ilevel(ji,jj) = jk
            END DO
         END DO
         WRITE(numout,*) 'zthick_0(jk =',jk,') =',zthick_0(2,2)
         WRITE(numout,*) 'gdepw(jk+1 =',jk+1,') =',gdepw(2,2,jk+1,Kmm)
      END DO

      ! Surface boundary condition
      IF( ln_linssh ) THEN  ;   zthick(:,:) = ssh(:,:,Kmm)   ;   htc_mld(:,:) = ts(:,:,1,jp_tem,Kmm) * ssh(:,:,Kmm) * tmask(:,:,1)   
      ELSE                  ;   zthick(:,:) = 0._wp       ;   htc_mld(:,:) = 0._wp                                   
      ENDIF

      ! Deepest whole T level above the MLD
      ikmax = MIN( MAXVAL( ilevel(:,:) ), jpkm1 )

      ! Integration down to last whole model T level
      DO jk = 1, ikmax
         DO jj = 1, jpj
            DO ji = 1, jpi
               zc = e3t(ji,jj,jk,Kmm) * REAL( MIN( MAX( 0, ilevel(ji,jj) - jk + 1 ) , 1  )  )    ! 0 below ilevel
               zthick(ji,jj) = zthick(ji,jj) + zc
               htc_mld(ji,jj) = htc_mld(ji,jj) + zc * ts(ji,jj,jk,jp_tem,Kmm) * tmask(ji,jj,jk)
            END DO
         END DO
      END DO

      ! Subsequent partial T level
      zthick(:,:) = hmld_zint(:,:) - zthick(:,:)   !   remaining thickness to reach MLD

      DO jj = 1, jpj
         DO ji = 1, jpi
            htc_mld(ji,jj) = htc_mld(ji,jj) + ts(ji,jj,ilevel(ji,jj)+1,jp_tem,Kmm)  & 
      &                      * MIN( e3t(ji,jj,ilevel(ji,jj)+1,Kmm), zthick(ji,jj) ) * tmask(ji,jj,ilevel(ji,jj)+1)
         END DO
      END DO

      WRITE(numout,*) 'htc_mld(after) =',htc_mld(2,2)

      ! Convert to heat content
      zcoef = rho0 * rcp
      htc_mld(:,:) = zcoef * htc_mld(:,:)

   END SUBROUTINE zdf_mxl_zint_htc

   SUBROUTINE zdf_mxl_zint( kt , Kmm)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdf_mxl_zint  ***
      !! 
      !! ** Purpose :   
      !!
      !! ** Method  :   
      !!----------------------------------------------------------------------

      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      INTEGER, INTENT(in) ::   Kmm  ! ocean time level index

      INTEGER :: ios
      INTEGER :: jn

      INTEGER :: nn_mld_diag = 0    ! number of diagnostics

      CHARACTER(len=1) :: cmld

      LOGICAL, SAVE, DIMENSION(5) :: l_mld, l_htc

      TYPE(MXL_ZINT) :: sn_mld1, sn_mld2, sn_mld3, sn_mld4, sn_mld5
      TYPE(MXL_ZINT), SAVE, DIMENSION(5) ::   mld_diags

      NAMELIST/namzdf_mldzint/ nn_mld_diag, sn_mld1, sn_mld2, sn_mld3, sn_mld4, sn_mld5

      !!----------------------------------------------------------------------
      
      IF( kt == nit000 ) THEN
         READ  ( numnam_ref, namzdf_mldzint, IOSTAT = ios, ERR = 901)
901      IF( ios /= 0 ) CALL ctl_nam ( ios , 'namzdf_mldzint in reference namelist' )

         READ  ( numnam_cfg, namzdf_mldzint, IOSTAT = ios, ERR = 902 )
902      IF( ios /= 0 ) CALL ctl_nam ( ios , 'namzdf_mldzint in configuration namelist' )
         IF(lwm) WRITE ( numond, namzdf_mldzint )

         IF( nn_mld_diag > 5 )   CALL ctl_stop( 'STOP', 'zdf_mxl_ini: Specify no more than 5 MLD definitions' )

         l_mld(:) = .FALSE. ; l_htc(:) = .FALSE.

         mld_diags(1) = sn_mld1
         mld_diags(2) = sn_mld2
         mld_diags(3) = sn_mld3
         mld_diags(4) = sn_mld4
         mld_diags(5) = sn_mld5

         IF( nn_mld_diag > 0 ) THEN
            IF( lwp ) THEN
               WRITE(numout,*) '=============== Vertically-interpolated mixed layer ================'
               WRITE(numout,*) '(Diagnostic number, nn_mld_type, rn_zref, rn_dT_crit, rn_iso_frac)'
            ENDIF

            DO jn = 1, nn_mld_diag
               ! Whether the diagnostic is requested
               WRITE(cmld,'(I1)') jn
               IF( iom_use( "mldzint_"//cmld ) ) l_mld(jn) = .TRUE.
               IF( iom_use( "mldhtc_"//cmld  ) ) l_htc(jn) = .TRUE.

               IF( lwp ) THEN
                  WRITE(numout,*) 'MLD criterion',jn,':'
                  WRITE(numout,*) '    nn_mld_type =', mld_diags(jn)%mld_type
                  WRITE(numout,*) '    rn_zref ='    , mld_diags(jn)%zref
                  WRITE(numout,*) '    rn_dT_crit =' , mld_diags(jn)%dT_crit
                  WRITE(numout,*) '    rn_iso_frac =', mld_diags(jn)%iso_frac
               ENDIF
            END DO
            WRITE(numout,*) '===================================================================='
         ENDIF
      ENDIF

      IF( nn_mld_diag > 0 ) THEN
         DO jn = 1, nn_mld_diag
            WRITE(cmld,'(I1)') jn

            IF( l_mld(jn) .OR. l_htc(jn) ) THEN
               CALL zdf_mxl_zint_mld( mld_diags(jn), Kmm)

               IF( l_mld(jn) ) CALL iom_put( "mldzint_"//cmld, hmld_zint(:,:) )

               IF( l_htc(jn) ) THEN
                  CALL zdf_mxl_zint_htc( kt, Kmm )
                  CALL iom_put( "mldhtc_"//cmld, htc_mld(:,:) )
               ENDIF
            ENDIF
         END DO
      ENDIF

   END SUBROUTINE zdf_mxl_zint


   SUBROUTINE zdf_mxl_turb( kt, Kmm )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdf_mxl_turb  ***
      !!
      !! ** Purpose :   Compute the turbocline depth.
      !!
      !! ** Method  :   The turbocline depth is the depth at which the vertical
      !!      eddy diffusivity coefficient (resulting from the vertical physics
      !!      alone, not the isopycnal part, see trazdf.F) fall below a given
      !!      value defined locally (avt_c here taken equal to 5 cm/s2 by default)
      !!
      !! ** Action  :   hmld
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      INTEGER, INTENT(in) ::   Kmm  ! ocean time level index
      !
      INTEGER  ::   ji, jj, jk      ! dummy loop indices
      INTEGER  ::   iik             ! local integer
      INTEGER, DIMENSION(A2D(nn_hls)) ::   imld   ! 2D workspace
      !!----------------------------------------------------------------------
      !
      ! w-level of the turbocline and mixing layer (iom_use)
      imld(:,:) = mbkt(A2D(nn_hls)) + 1                ! Initialization to the number of w ocean point
      DO_3DS( 1, 1, 1, 1, jpkm1, nlb10, -1 )   ! from the bottom to nlb10
         IF( avt (ji,jj,jk) < avt_c * wmask(ji,jj,jk) )   imld(ji,jj) = jk      ! Turbocline
      END_3D
      ! depth of the mixing layer
      DO_2D_OVR( 1, 1, 1, 1 )
         iik = imld(ji,jj)
         hmld (ji,jj) = gdepw(ji,jj,iik  ,Kmm) * ssmask(ji,jj)    ! Turbocline depth
      END_2D
      !
      IF( .NOT.l_offline .AND. iom_use("mldkz5") ) THEN
         IF( .NOT. l_istiled .OR. ntile == nijtile ) THEN         ! Do only on the last tile
            IF( ln_isfcav ) THEN  ;  CALL iom_put( "mldkz5"  , hmld - risfdep )   ! turbocline thickness
            ELSE                  ;  CALL iom_put( "mldkz5"  , hmld )             ! turbocline depth
            END IF
         ENDIF
      ENDIF
      !
   END SUBROUTINE zdf_mxl_turb
   !!======================================================================
END MODULE zdfmxl
