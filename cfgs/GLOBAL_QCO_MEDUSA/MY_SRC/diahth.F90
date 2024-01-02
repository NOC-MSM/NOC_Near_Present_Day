MODULE diahth
   !!======================================================================
   !!                       ***  MODULE  diahth  ***
   !! Ocean diagnostics: thermocline and 20 degree depth
   !!======================================================================
   !! History :  OPA  !  1994-09  (J.-P. Boulanger)  Original code
   !!                 !  1996-11  (E. Guilyardi)  OPA8 
   !!                 !  1997-08  (G. Madec)  optimization
   !!                 !  1999-07  (E. Guilyardi)  hd28 + heat content 
   !!   NEMO     1.0  !  2002-06  (G. Madec)  F90: Free form and module
   !!            3.2  !  2009-07  (S. Masson) hc300 bugfix + cleaning + add new diag
   !!----------------------------------------------------------------------
   !!   dia_hth      : Compute varius diagnostics associated with the mixed layer
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE zdfmxl, ONLY: zdf_mxl_zint
   !
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! MPP library
   USE iom             ! I/O library
   USE timing          ! preformance summary

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dia_hth       ! routine called by step.F90
   PUBLIC   dia_hth_alloc ! routine called by nemogcm.F90

   LOGICAL, SAVE  ::   l_hth     !: thermocline-20d depths flag
   
   ! note: following variables should move to local variables once iom_put is always used 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hth    !: depth of the max vertical temperature gradient [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hd20   !: depth of 20 C isotherm                         [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hd26   !: depth of 26 C isotherm                         [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hd28   !: depth of 28 C isotherm                         [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   htc3   !: heat content of first 300 m                    [W]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   htc7   !: heat content of first 700 m                    [W]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   htc20  !: heat content of first 2000 m                   [W]


   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: diahth.F90 15234 2021-09-08 14:07:02Z clem $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   FUNCTION dia_hth_alloc()
      !!---------------------------------------------------------------------
      INTEGER :: dia_hth_alloc
      !!---------------------------------------------------------------------
      !
      ALLOCATE( hth(jpi,jpj), hd20(jpi,jpj), hd26(jpi,jpj), hd28(jpi,jpj), &
         &      htc3(jpi,jpj), htc7(jpi,jpj), htc20(jpi,jpj), STAT=dia_hth_alloc )
      !
      CALL mpp_sum ( 'diahth', dia_hth_alloc )
      IF(dia_hth_alloc /= 0)   CALL ctl_stop( 'STOP', 'dia_hth_alloc: failed to allocate arrays.' )
      !
   END FUNCTION dia_hth_alloc


   SUBROUTINE dia_hth( kt, Kmm )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dia_hth  ***
      !!
      !! ** Purpose : Computes
      !!      the mixing layer depth (turbocline): avt = 5.e-4
      !!      the depth of strongest vertical temperature gradient
      !!      the mixed layer depth with density     criteria: rho = rho(10m or surf) + 0.03(or 0.01)
      !!      the mixed layer depth with temperature criteria: abs( tn - tn(10m) ) = 0.2       
      !!      the top of the thermochine: tn = tn(10m) - ztem2 
      !!      the pycnocline depth with density criteria equivalent to a temperature variation 
      !!                rho = rho10m + (dr/dT)(T,S,10m)*(-0.2 degC) 
      !!      the barrier layer thickness
      !!      the maximal verical inversion of temperature and its depth max( 0, max of tn - tn(10m) )
      !!      the depth of the 20 degree isotherm (linear interpolation)
      !!      the depth of the 28 degree isotherm (linear interpolation)
      !!      the heat content of first 300 m
      !!
      !! ** Method : 
      !!-------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      INTEGER, INTENT( in ) ::   Kmm     ! ocean time level index
      !!
      INTEGER                      ::   ji, jj, jk            ! dummy loop arguments
      REAL(wp)                     ::   zrho3 = 0.03_wp       ! density     criterion for mixed layer depth
      REAL(wp)                     ::   zrho1 = 0.01_wp       ! density     criterion for mixed layer depth
      REAL(wp)                     ::   ztem2 = 0.2_wp        ! temperature criterion for mixed layer depth
      REAL(wp)                     ::   zztmp, zzdep          ! temporary scalars inside do loop
      REAL(wp)                     ::   zu, zv, zw, zut, zvt  ! temporary workspace
      REAL(wp), DIMENSION(jpi,jpj) ::   zabs2      ! MLD: abs( tn - tn(10m) ) = ztem2 
      REAL(wp), DIMENSION(jpi,jpj) ::   ztm2       ! Top of thermocline: tn = tn(10m) - ztem2     
      REAL(wp), DIMENSION(jpi,jpj) ::   zrho10_3   ! MLD: rho = rho10m + zrho3      
      REAL(wp), DIMENSION(jpi,jpj) ::   zpycn      ! pycnocline: rho = rho10m + (dr/dT)(T,S,10m)*(-0.2 degC)
      REAL(wp), DIMENSION(jpi,jpj) ::   ztinv      ! max of temperature inversion
      REAL(wp), DIMENSION(jpi,jpj) ::   zdepinv    ! depth of temperature inversion
      REAL(wp), DIMENSION(jpi,jpj) ::   zrho0_3    ! MLD rho = rho(surf) = 0.03
      REAL(wp), DIMENSION(jpi,jpj) ::   zrho0_1    ! MLD rho = rho(surf) = 0.01
      REAL(wp), DIMENSION(jpi,jpj) ::   zmaxdzT    ! max of dT/dz
      REAL(wp), DIMENSION(jpi,jpj) ::   zdelr      ! delta rho equivalent to deltaT = 0.2
      !!----------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('dia_hth')

      IF( kt == nit000 ) THEN
         !
         l_hth = iom_use( 'mlddzt'   ) .OR. iom_use( 'mldr0_3'  ) .OR. iom_use( 'mldr0_1'  )    .OR.  & 
            &    iom_use( 'mld_dt02' ) .OR. iom_use( 'topthdep' ) .OR. iom_use( 'mldr10_3' )    .OR.  &    
            &    iom_use( '20d'      ) .OR. iom_use( '26d'      ) .OR. iom_use( '28d'      )    .OR.  &    
            &    iom_use( 'hc300'    ) .OR. iom_use( 'hc700'    ) .OR. iom_use( 'hc2000'   )    .OR.  &    
            &    iom_use( 'pycndep'  ) .OR. iom_use( 'tinv'     ) .OR. iom_use( 'depti'    )
         !
         !                                      ! allocate dia_hth array
         IF( l_hth ) THEN 
            IF( dia_hth_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'dia_hth : unable to allocate standard arrays' )
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'dia_hth : diagnostics of the thermocline depth'
            IF(lwp) WRITE(numout,*) '~~~~~~~ '
            IF(lwp) WRITE(numout,*)
         ENDIF
      ENDIF

      IF( l_hth ) THEN
         !
         ! initialization
         IF( iom_use( 'tinv'   ) )   ztinv  (:,:) = 0._wp  
         IF( iom_use( 'depti'  ) )   zdepinv(:,:) = 0._wp  
         IF( iom_use( 'mlddzt' ) )   zmaxdzT(:,:) = 0._wp  
         IF( iom_use( 'mlddzt' ) .OR. iom_use( 'mld_dt02' ) .OR. iom_use( 'topthdep' )   &
            &                    .OR. iom_use( 'mldr10_3' ) .OR. iom_use( 'pycndep'  ) ) THEN
            DO_2D( 1, 1, 1, 1 )
               zztmp = gdepw(ji,jj,mbkt(ji,jj)+1,Kmm) 
               hth     (ji,jj) = zztmp
               zabs2   (ji,jj) = zztmp
               ztm2    (ji,jj) = zztmp
               zrho10_3(ji,jj) = zztmp
               zpycn   (ji,jj) = zztmp
            END_2D
         ENDIF
         IF( iom_use( 'mldr0_3' ) .OR. iom_use( 'mldr0_1' ) ) THEN
            IF( nla10 > 1 ) THEN 
               DO_2D( 1, 1, 1, 1 )
                  zztmp = gdepw(ji,jj,mbkt(ji,jj)+1,Kmm) 
                  zrho0_3(ji,jj) = zztmp
                  zrho0_1(ji,jj) = zztmp
               END_2D
            ENDIF
         ENDIF
     
         IF( iom_use( 'mlddzt' ) .OR. iom_use( 'mldr0_3' ) .OR. iom_use( 'mldr0_1' ) ) THEN
            ! ------------------------------------------------------------- !
            ! thermocline depth: strongest vertical gradient of temperature !
            ! turbocline depth (mixing layer depth): avt = zavt5            !
            ! MLD: rho = rho(1) + zrho3                                     !
            ! MLD: rho = rho(1) + zrho1                                     !
            ! ------------------------------------------------------------- !
            DO_3DS( 1, 1, 1, 1, jpkm1, 2, -1 )   ! loop from bottom to 2
               !
               zzdep = gdepw(ji,jj,jk,Kmm)
               zztmp = ( ts(ji,jj,jk-1,jp_tem,Kmm) - ts(ji,jj,jk,jp_tem,Kmm) ) &
                      & / zzdep * tmask(ji,jj,jk)   ! vertical gradient of temperature (dT/dz)
               zzdep = zzdep * tmask(ji,jj,1)

               IF( zztmp > zmaxdzT(ji,jj) ) THEN                        
                   zmaxdzT(ji,jj) = zztmp   
                   hth    (ji,jj) = zzdep                ! max and depth of dT/dz
               ENDIF
         
               IF( nla10 > 1 ) THEN 
                  zztmp = rhop(ji,jj,jk) - rhop(ji,jj,1)                       ! delta rho(1)
                  IF( zztmp > zrho3 )   zrho0_3(ji,jj) = zzdep                ! > 0.03
                  IF( zztmp > zrho1 )   zrho0_1(ji,jj) = zzdep                ! > 0.01
               ENDIF
            END_3D
         
            CALL iom_put( 'mlddzt', hth )            ! depth of the thermocline
            IF( nla10 > 1 ) THEN 
               CALL iom_put( 'mldr0_3', zrho0_3 )   ! MLD delta rho(surf) = 0.03
               CALL iom_put( 'mldr0_1', zrho0_1 )   ! MLD delta rho(surf) = 0.01
            ENDIF
            !
         ENDIF
         !
         IF(  iom_use( 'mld_dt02' ) .OR. iom_use( 'topthdep' ) .OR. iom_use( 'mldr10_3' ) .OR.  &    
            &  iom_use( 'pycndep' ) .OR. iom_use( 'tinv'     ) .OR. iom_use( 'depti'    )  ) THEN
            !
            ! Preliminary computation
            ! computation of zdelr = (dr/dT)(T,S,10m)*(-0.2 degC)
            DO_2D( 1, 1, 1, 1 )
               IF( tmask(ji,jj,nla10) == 1. ) THEN
                  zu  =  1779.50 + 11.250 * ts(ji,jj,nla10,jp_tem,Kmm) - 3.80   * ts(ji,jj,nla10,jp_sal,Kmm)  &
                     &           - 0.0745 * ts(ji,jj,nla10,jp_tem,Kmm) * ts(ji,jj,nla10,jp_tem,Kmm)   &
                     &           - 0.0100 * ts(ji,jj,nla10,jp_tem,Kmm) * ts(ji,jj,nla10,jp_sal,Kmm)
                  zv  =  5891.00 + 38.000 * ts(ji,jj,nla10,jp_tem,Kmm) + 3.00   * ts(ji,jj,nla10,jp_sal,Kmm)  &
                     &           - 0.3750 * ts(ji,jj,nla10,jp_tem,Kmm) * ts(ji,jj,nla10,jp_tem,Kmm)
                  zut =    11.25 -  0.149 * ts(ji,jj,nla10,jp_tem,Kmm) - 0.01   * ts(ji,jj,nla10,jp_sal,Kmm)
                  zvt =    38.00 -  0.750 * ts(ji,jj,nla10,jp_tem,Kmm)
                  zw  = (zu + 0.698*zv) * (zu + 0.698*zv)
                  zdelr(ji,jj) = ztem2 * (1000.*(zut*zv - zvt*zu)/zw)
               ELSE
                  zdelr(ji,jj) = 0._wp
               ENDIF
            END_2D
            !
            ! ------------------------------------------------------------- !
            ! MLD: abs( tn - tn(10m) ) = ztem2                              !
            ! Top of thermocline: tn = tn(10m) - ztem2                      !
            ! MLD: rho = rho10m + zrho3                                     !
            ! pycnocline: rho = rho10m + (dr/dT)(T,S,10m)*(-0.2 degC)       !
            ! temperature inversion: max( 0, max of tn - tn(10m) )          !
            ! depth of temperature inversion                                !
            ! ------------------------------------------------------------- !
            DO_3DS( 1, 1, 1, 1, jpkm1, nlb10, -1 )   ! loop from bottom to nlb10
               !
               zzdep = gdepw(ji,jj,jk,Kmm) * tmask(ji,jj,1)
               !
               zztmp = ts(ji,jj,nla10,jp_tem,Kmm) - ts(ji,jj,jk,jp_tem,Kmm)  ! - delta T(10m)
               IF( ABS(zztmp) > ztem2 )      zabs2   (ji,jj) = zzdep   ! abs > 0.2
               IF(     zztmp  > ztem2 )      ztm2    (ji,jj) = zzdep   ! > 0.2
               zztmp = -zztmp                                          ! delta T(10m)
               IF( zztmp >  ztinv(ji,jj) ) THEN                        ! temperature inversion
                  ztinv(ji,jj) = zztmp   
                  zdepinv (ji,jj) = zzdep   ! max value and depth
               ENDIF

               zztmp = rhop(ji,jj,jk) - rhop(ji,jj,nla10)              ! delta rho(10m)
               IF( zztmp > zrho3        )    zrho10_3(ji,jj) = zzdep   ! > 0.03
               IF( zztmp > zdelr(ji,jj) )    zpycn   (ji,jj) = zzdep   ! > equi. delta T(10m) - 0.2
               !
            END_3D

            CALL iom_put( 'mld_dt02', zabs2    )   ! MLD abs(delta t) - 0.2
            CALL iom_put( 'topthdep', ztm2     )   ! T(10) - 0.2
            CALL iom_put( 'mldr10_3', zrho10_3 )   ! MLD delta rho(10m) = 0.03
            CALL iom_put( 'pycndep' , zpycn    )   ! MLD delta rho equi. delta T(10m) = 0.2
            CALL iom_put( 'tinv'    , ztinv    )   ! max. temp. inv. (t10 ref) 
            CALL iom_put( 'depti'   , zdepinv  )   ! depth of max. temp. inv. (t10 ref) 
            !
         ENDIF
 
         ! ------------------------------- !
         !  Depth of 20C/26C/28C isotherm  !
         ! ------------------------------- !
         IF( iom_use ('20d') ) THEN  ! depth of the 20 isotherm
            ztem2 = 20.
            CALL dia_hth_dep( Kmm, ztem2, hd20 )  
            CALL iom_put( '20d', hd20 )    
         ENDIF
         !
         IF( iom_use ('26d') ) THEN  ! depth of the 26 isotherm
            ztem2 = 26.
            CALL dia_hth_dep( Kmm, ztem2, hd26 )  
            CALL iom_put( '26d', hd26 )    
         ENDIF
         !
         IF( iom_use ('28d') ) THEN  ! depth of the 28 isotherm
            ztem2 = 28.
            CALL dia_hth_dep( Kmm, ztem2, hd28 )  
            CALL iom_put( '28d', hd28 )    
         ENDIF
        
         ! ----------------------------- !
         !  Heat content of first 300 m  !
         ! ----------------------------- !
         IF( iom_use ('hc300') ) THEN  
            zzdep = 300.
            CALL  dia_hth_htc( Kmm, zzdep, ts(:,:,:,jp_tem,Kmm), htc3 )
            CALL iom_put( 'hc300', rho0_rcp * htc3 )  ! vertically integrated heat content (J/m2)
         ENDIF
         !
         ! ----------------------------- !
         !  Heat content of first 700 m  !
         ! ----------------------------- !
         IF( iom_use ('hc700') ) THEN  
            zzdep = 700.
            CALL  dia_hth_htc( Kmm, zzdep, ts(:,:,:,jp_tem,Kmm), htc7 )
            CALL iom_put( 'hc700', rho0_rcp * htc7 )  ! vertically integrated heat content (J/m2)
  
         ENDIF
         !
         ! ----------------------------- !
         !  Heat content of first 2000 m  !
         ! ----------------------------- !
         IF( iom_use ('hc2000') ) THEN  
            zzdep = 2000.
            CALL  dia_hth_htc( Kmm, zzdep, ts(:,:,:,jp_tem,Kmm), htc20 )
            CALL iom_put( 'hc2000', rho0_rcp * htc20 )  ! vertically integrated heat content (J/m2)  
         ENDIF
         !
      ENDIF

      ! Vertically-interpolated mixed-layer depth diagnostic
      CALL zdf_mxl_zint( kt, Kmm )

      !
      IF( ln_timing )   CALL timing_stop('dia_hth')
      !
   END SUBROUTINE dia_hth

   SUBROUTINE dia_hth_dep( Kmm, ptem, pdept )
      !
      INTEGER , INTENT(in) :: Kmm      ! ocean time level index
      REAL(wp), INTENT(in) :: ptem
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) :: pdept     
      !
      INTEGER  :: ji, jj, jk, iid
      REAL(wp) :: zztmp, zzdep
      INTEGER, DIMENSION(jpi,jpj) :: iktem
      
      ! --------------------------------------- !
      ! search deepest level above ptem         !
      ! --------------------------------------- !
      iktem(:,:) = 1
      DO_3D( 1, 1, 1, 1, 1, jpkm1 )   ! beware temperature is not always decreasing with depth => loop from top to bottom
         zztmp = ts(ji,jj,jk,jp_tem,Kmm)
         IF( zztmp >= ptem )   iktem(ji,jj) = jk
      END_3D

      ! ------------------------------- !
      !  Depth of ptem isotherm         !
      ! ------------------------------- !
      DO_2D( 1, 1, 1, 1 )
         !
         zzdep = gdepw(ji,jj,mbkt(ji,jj)+1,Kmm)       ! depth of the ocean bottom
         !
         iid = iktem(ji,jj)
         IF( iid /= 1 ) THEN 
             zztmp =     gdept(ji,jj,iid  ,Kmm)   &                     ! linear interpolation
               &  + (    gdept(ji,jj,iid+1,Kmm) - gdept(ji,jj,iid,Kmm)                            )   &
               &  * ( ptem * tmask(ji,jj,iid+1)  - ts(ji,jj,iid,jp_tem,Kmm)                       )   &
               &  / ( ts(ji,jj,iid+1,jp_tem,Kmm) - ts(ji,jj,iid,jp_tem,Kmm) + (1.-tmask(ji,jj,1)) )
            pdept(ji,jj) = MIN( zztmp , zzdep) * tmask(ji,jj,1)       ! bound by the ocean depth
         ELSE 
            pdept(ji,jj) = 0._wp
         ENDIF
      END_2D
      !
   END SUBROUTINE dia_hth_dep


   SUBROUTINE dia_hth_htc( Kmm, pdep, pt, phtc )
      !
      INTEGER , INTENT(in) ::   Kmm      ! ocean time level index
      REAL(wp), INTENT(in) ::   pdep     ! depth over the heat content
      REAL(dp), DIMENSION(jpi,jpj,jpk), INTENT(in)    ::   pt
      REAL(wp), DIMENSION(jpi,jpj),     INTENT(inout) ::   phtc  
      !
      INTEGER  ::   ji, jj, jk, ik
      REAL(wp), DIMENSION(jpi,jpj) ::   zthick
      INTEGER , DIMENSION(jpi,jpj) ::   ilevel


      ! surface boundary condition
      
      IF( .NOT. ln_linssh ) THEN   ;   zthick(:,:) = 0._wp          ;   phtc(:,:) = 0._wp                                   
      ELSE                         ;   zthick(:,:) = ssh(:,:,Kmm)   ;   phtc(:,:) = pt(:,:,1) * ssh(:,:,Kmm) * tmask(:,:,1)   
      ENDIF
      !
      ilevel(:,:) = 1
      DO_3D( 1, 1, 1, 1, 1, jpkm1 )
         IF( ( gdepw(ji,jj,jk+1,Kmm) < pdep ) .AND. ( tmask(ji,jj,jk) == 1 ) ) THEN
             ilevel(ji,jj) = jk+1
             zthick(ji,jj) = zthick(ji,jj) + e3t(ji,jj,jk,Kmm)
             phtc  (ji,jj) = phtc  (ji,jj) + e3t(ji,jj,jk,Kmm) * pt(ji,jj,jk)
         ENDIF
      END_3D
      !
      DO_2D( 1, 1, 1, 1 )
         ik = ilevel(ji,jj)
         IF( tmask(ji,jj,ik) == 1 ) THEN
            zthick(ji,jj) = MIN ( gdepw(ji,jj,ik+1,Kmm), pdep ) - zthick(ji,jj)   ! remaining thickness to reach dephw pdep
            phtc(ji,jj)   = phtc(ji,jj) + pt(ji,jj,ik) * zthick(ji,jj)
         ENDIF
      END_2D
      !
   END SUBROUTINE dia_hth_htc

   !!======================================================================
END MODULE diahth
