MODULE zdftke
   !!======================================================================
   !!                       ***  MODULE  zdftke  ***
   !! Ocean physics:  vertical mixing coefficient computed from the tke
   !!                 turbulent closure parameterization
   !!=====================================================================
   !! History :  OPA  !  1991-03  (b. blanke)  Original code
   !!            7.0  !  1991-11  (G. Madec)   bug fix
   !!            7.1  !  1992-10  (G. Madec)   new mixing length and eav
   !!            7.2  !  1993-03  (M. Guyon)   symetrical conditions
   !!            7.3  !  1994-08  (G. Madec, M. Imbard)  nn_pdl flag
   !!            7.5  !  1996-01  (G. Madec)   s-coordinates
   !!            8.0  !  1997-07  (G. Madec)   lbc
   !!            8.1  !  1999-01  (E. Stretta) new option for the mixing length
   !!  NEMO      1.0  !  2002-06  (G. Madec) add tke_init routine
   !!             -   !  2004-10  (C. Ethe )  1D configuration
   !!            2.0  !  2006-07  (S. Masson)  distributed restart using iom
   !!            3.0  !  2008-05  (C. Ethe,  G.Madec) : update TKE physics:
   !!                 !           - tke penetration (wind steering)
   !!                 !           - suface condition for tke & mixing length
   !!                 !           - Langmuir cells
   !!             -   !  2008-05  (J.-M. Molines, G. Madec)  2D form of avtb
   !!             -   !  2008-06  (G. Madec)  style + DOCTOR name for namelist parameters
   !!             -   !  2008-12  (G. Reffray) stable discretization of the production term
   !!            3.2  !  2009-06  (G. Madec, S. Masson) TKE restart compatible with key_cpl
   !!                 !                                + cleaning of the parameters + bugs correction
   !!            3.3  !  2010-10  (C. Ethe, G. Madec) reorganisation of initialisation phase
   !!            3.6  !  2014-11  (P. Mathiot) add ice shelf capability
   !!            4.0  !  2017-04  (G. Madec)  remove CPP ddm key & avm at t-point only
   !!             -   !  2017-05  (G. Madec)  add top/bottom friction as boundary condition
   !!            4.2  !  2020-12  (G. Madec, E. Clementi) add wave coupling
   !                  !           following Couvelard et al., 2019
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   zdf_tke       : update momentum and tracer Kz from a tke scheme
   !!   tke_tke       : tke time stepping: update tke at now time step (en)
   !!   tke_avn       : compute mixing length scale and deduce avm and avt
   !!   zdf_tke_init  : initialization, namelist read, and parameters control
   !!   tke_rst       : read/write tke restart in ocean restart file
   !!----------------------------------------------------------------------
   USE oce            ! ocean: dynamics and active tracers variables
   USE phycst         ! physical constants
   USE dom_oce        ! domain: ocean
   USE domvvl         ! domain: variable volume layer
   USE sbc_oce        ! surface boundary condition: ocean
   USE zdfdrg         ! vertical physics: top/bottom drag coef.
   USE zdfmxl         ! vertical physics: mixed layer
#if defined key_si3
   USE ice, ONLY: hm_i, h_i
#endif
#if defined key_cice
   USE sbc_ice, ONLY: h_i
#endif
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE prtctl         ! Print control
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)
   USE sbcwave        ! Surface boundary waves

   IMPLICIT NONE
   PRIVATE

   PUBLIC   zdf_tke        ! routine called in step module
   PUBLIC   zdf_tke_init   ! routine called in opa module
   PUBLIC   tke_rst        ! routine called in step module

   !                      !!** Namelist  namzdf_tke  **
   LOGICAL  ::   ln_mxl0   ! mixing length scale surface value as function of wind stress or not
   LOGICAL  ::   ln_mxhsw  ! mixing length scale surface value as a fonction of wave height
   INTEGER  ::   nn_mxlice ! type of scaling under sea-ice (=0/1/2/3)
   REAL(wp) ::   rn_mxlice ! ice thickness value when scaling under sea-ice
   INTEGER  ::   nn_mxl    ! type of mixing length (=0/1/2/3)
   REAL(wp) ::   rn_mxl0   ! surface  min value of mixing length (kappa*z_o=0.4*0.1 m)  [m]
   INTEGER  ::   nn_pdl    ! Prandtl number or not (ratio avt/avm) (=0/1)
   REAL(wp) ::   rn_ediff  ! coefficient for avt: avt=rn_ediff*mxl*sqrt(e)
   REAL(wp) ::   rn_ediss  ! coefficient of the Kolmogoroff dissipation
   REAL(wp) ::   rn_ebb    ! coefficient of the surface input of tke
   REAL(wp) ::   rn_emin   ! minimum value of tke           [m2/s2]
   REAL(wp) ::   rn_emin0  ! surface minimum value of tke   [m2/s2]
   REAL(wp) ::   rn_bshear ! background shear (>0) currently a numerical threshold (do not change it)
   INTEGER  ::   nn_etau   ! type of depth penetration of surface tke (=0/1/2/3)
   INTEGER  ::      nn_htau   ! type of tke profile of penetration (=0/1)
   REAL(wp) ::      rn_htau_scaling   ! a acaling factor to apply to the penetration of TKE
   INTEGER  ::   nn_bc_surf! surface condition (0/1=Dir/Neum) ! Only applicable for wave coupling
   INTEGER  ::   nn_bc_bot ! surface condition (0/1=Dir/Neum) ! Only applicable for wave coupling
   REAL(wp) ::      rn_efr    ! fraction of TKE surface value which penetrates in the ocean
   LOGICAL  ::   ln_lc     ! Langmuir cells (LC) as a source term of TKE or not
   REAL(wp) ::      rn_lc     ! coef to compute vertical velocity of Langmuir cells
   INTEGER  ::   nn_eice   ! attenutaion of langmuir & surface wave breaking under ice (=0/1/2/3)

   REAL(wp) ::   ri_cri    ! critic Richardson number (deduced from rn_ediff and rn_ediss values)
   REAL(wp) ::   rmxl_min  ! minimum mixing length value (deduced from rn_ediff and rn_emin values)  [m]
   REAL(wp) ::   rhftau_add = 1.e-3_wp     ! add offset   applied to HF part of taum  (nn_etau=3)
   REAL(wp) ::   rhftau_scl = 1.0_wp       ! scale factor applied to HF part of taum  (nn_etau=3)

   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:)   ::   htau    ! depth of tke penetration (nn_htau)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   dissl   ! now mixing lenght of dissipation
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   apdlr   ! now mixing lenght of dissipation

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "single_precision_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: zdftke.F90 15071 2021-07-02 13:12:08Z clem $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION zdf_tke_alloc()
      !!----------------------------------------------------------------------
      !!                ***  FUNCTION zdf_tke_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( htau(jpi,jpj) , dissl(jpi,jpj,jpk) , apdlr(jpi,jpj,jpk) ,   STAT= zdf_tke_alloc )
      !
      CALL mpp_sum ( 'zdftke', zdf_tke_alloc )
      IF( zdf_tke_alloc /= 0 )   CALL ctl_stop( 'STOP', 'zdf_tke_alloc: failed to allocate arrays' )
      !
   END FUNCTION zdf_tke_alloc


   SUBROUTINE zdf_tke( kt, Kbb, Kmm, p_sh2, p_avm, p_avt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE zdf_tke  ***
      !!
      !! ** Purpose :   Compute the vertical eddy viscosity and diffusivity
      !!              coefficients using a turbulent closure scheme (TKE).
      !!
      !! ** Method  :   The time evolution of the turbulent kinetic energy (tke)
      !!              is computed from a prognostic equation :
      !!         d(en)/dt = avm (d(u)/dz)**2             ! shear production
      !!                  + d( avm d(en)/dz )/dz         ! diffusion of tke
      !!                  + avt N^2                      ! stratif. destruc.
      !!                  - rn_ediss / emxl en**(2/3)    ! Kolmogoroff dissipation
      !!      with the boundary conditions:
      !!         surface: en = max( rn_emin0, rn_ebb * taum )
      !!         bottom : en = rn_emin
      !!      The associated critical Richardson number is: ri_cri = 2/(2+rn_ediss/rn_ediff)
      !!
      !!        The now Turbulent kinetic energy is computed using the following
      !!      time stepping: implicit for vertical diffusion term, linearized semi
      !!      implicit for kolmogoroff dissipation term, and explicit forward for
      !!      both buoyancy and shear production terms. Therefore a tridiagonal
      !!      linear system is solved. Note that buoyancy and shear terms are
      !!      discretized in a energy conserving form (Bruchard 2002).
      !!
      !!        The dissipative and mixing length scale are computed from en and
      !!      the stratification (see tke_avn)
      !!
      !!        The now vertical eddy vicosity and diffusivity coefficients are
      !!      given by:
      !!              avm = max( avtb, rn_ediff * zmxlm * en^1/2 )
      !!              avt = max( avmb, pdl * avm                 )
      !!              eav = max( avmb, avm )
      !!      where pdl, the inverse of the Prandtl number is 1 if nn_pdl=0 and
      !!      given by an empirical funtion of the localRichardson number if nn_pdl=1
      !!
      !! ** Action  :   compute en (now turbulent kinetic energy)
      !!                update avt, avm (before vertical eddy coef.)
      !!
      !! References : Gaspar et al., JGR, 1990,
      !!              Blanke and Delecluse, JPO, 1991
      !!              Mellor and Blumberg, JPO 2004
      !!              Axell, JGR, 2002
      !!              Bruchard OM 2002
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT(in   ) ::   kt             ! ocean time step
      INTEGER                             , INTENT(in   ) ::   Kbb, Kmm       ! ocean time level indices
      REAL(wp), DIMENSION(A2D(nn_hls),jpk), INTENT(in   ) ::   p_sh2          ! shear production term
      REAL(wp), DIMENSION(:,:,:)          , INTENT(inout) ::   p_avm, p_avt   !  momentum and tracer Kz (w-points)
      !!----------------------------------------------------------------------
      !
      CALL tke_tke( Kbb, Kmm, p_sh2, p_avm, p_avt )   ! now tke (en)
      !
      CALL tke_avn( Kbb, Kmm,        p_avm, p_avt )   ! now avt, avm, dissl
      !
  END SUBROUTINE zdf_tke


   SUBROUTINE tke_tke( Kbb, Kmm, p_sh2, p_avm, p_avt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE tke_tke  ***
      !!
      !! ** Purpose :   Compute the now Turbulente Kinetic Energy (TKE)
      !!
      !! ** Method  : - TKE surface boundary condition
      !!              - source term due to Langmuir cells (Axell JGR 2002) (ln_lc=T)
      !!              - source term due to shear (= Kz dz[Ub] * dz[Un] )
      !!              - Now TKE : resolution of the TKE equation by inverting
      !!                a tridiagonal linear system by a "methode de chasse"
      !!              - increase TKE due to surface and internal wave breaking
      !!             NB: when sea-ice is present, both LC parameterization
      !!                 and TKE penetration are turned off when the ice fraction
      !!                 is smaller than 0.25
      !!
      !! ** Action  : - en : now turbulent kinetic energy)
      !! ---------------------------------------------------------------------
      USE zdf_oce , ONLY : en   ! ocean vertical physics
      !!
      INTEGER                              , INTENT(in   ) ::   Kbb, Kmm       ! ocean time level indices
      REAL(wp), DIMENSION(A2D(nn_hls),jpk) , INTENT(in   ) ::   p_sh2          ! shear production term
      REAL(wp), DIMENSION(:,:,:)           , INTENT(in   ) ::   p_avm, p_avt   ! vertical eddy viscosity & diffusivity (w-points)
      !
      INTEGER ::   ji, jj, jk                  ! dummy loop arguments
      REAL(wp) ::   zetop, zebot, zmsku, zmskv ! local scalars
      REAL(wp) ::   zrhoa  = 1.22              ! Air density kg/m3
      REAL(wp) ::   zcdrag = 1.5e-3            ! drag coefficient
      REAL(wp) ::   zbbrau, zbbirau, zri       ! local scalars
      REAL(wp) ::   zfact1, zfact2, zfact3     !   -      -
      REAL(wp) ::   ztx2  , zty2  , zcof       !   -      -
      REAL(wp) ::   ztau  , zdif               !   -      -
      REAL(wp) ::   zus   , zwlc  , zind       !   -      -
      REAL(wp) ::   zzd_up, zzd_lw             !   -      -
      REAL(wp) ::   ztaui, ztauj, z1_norm
      INTEGER , DIMENSION(A2D(nn_hls))     ::   imlc
      REAL(wp), DIMENSION(A2D(nn_hls))     ::   zice_fra, zhlc, zus3, zWlc2
      REAL(wp), DIMENSION(A2D(nn_hls),jpk) ::   zpelc, zdiag, zd_up, zd_lw
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, SAVE ::   ztmp ! for diags
      REAL(wp) :: zdiv
      !!--------------------------------------------------------------------
      !
      zbbrau  = rn_ebb / rho0       ! Local constant initialisation
      zbbirau = 3.75_wp / rho0
      zfact1  = -.5_wp * rn_Dt
      zfact2  = 1.5_wp * rn_Dt * rn_ediss
      zfact3  = 0.5_wp         * rn_ediss
      !
      zpelc(:,:,:) = 0._wp ! need to be initialised in case ln_lc is not used
      !
      ! ice fraction considered for attenuation of langmuir & wave breaking
      SELECT CASE ( nn_eice )
      CASE( 0 )   ;   zice_fra(:,:) = 0._wp
      CASE( 1 )   ;   zice_fra(:,:) =        TANH( fr_i(A2D(nn_hls)) * 10._wp )
      CASE( 2 )   ;   zice_fra(:,:) =              fr_i(A2D(nn_hls))
      CASE( 3 )   ;   zice_fra(:,:) = MIN( 4._wp * fr_i(A2D(nn_hls)) , 1._wp )
      END SELECT
      !
      !                     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !                     !  Surface/top/bottom boundary condition on tke
      !                     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !
      DO_2D_OVR( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1 )
         en(ji,jj,1) = MAX( rn_emin0, zbbrau * taum(ji,jj) )
         zdiag(ji,jj,1) = 1._wp/en(ji,jj,1)
         zd_lw(ji,jj,1) = 1._wp
         zd_up(ji,jj,1) = 0._wp
      END_2D
      !
      !                     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !                     !  Bottom boundary condition on tke
      !                     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !
      !   en(bot)   = (ebb0/rho0)*0.5*sqrt(u_botfr^2+v_botfr^2) (min value rn_emin)
      ! where ebb0 does not includes surface wave enhancement (i.e. ebb0=3.75)
      ! Note that stress averaged is done using an wet-only calculation of u and v at t-point like in zdfsh2
      !
      IF( .NOT.ln_drg_OFF ) THEN    !== friction used as top/bottom boundary condition on TKE
         !
         DO_2D_OVR( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1 )        ! bottom friction
            zmsku = ( 2. - umask(ji-1,jj,mbkt(ji,jj)) * umask(ji,jj,mbkt(ji,jj)) )
            zmskv = ( 2. - vmask(ji,jj-1,mbkt(ji,jj)) * vmask(ji,jj,mbkt(ji,jj)) )
            !                       ! where 0.001875 = (rn_ebb0/rho0) * 0.5 = 3.75*0.5/1000. (CAUTION CdU<0)
            zebot = - 0.001875_wp * rCdU_bot(ji,jj) * SQRT(  ( zmsku*( uu(ji,jj,mbkt(ji,jj),Kbb)+uu(ji-1,jj,mbkt(ji,jj),Kbb) ) )**2  &
               &                                           + ( zmskv*( vv(ji,jj,mbkt(ji,jj),Kbb)+vv(ji,jj-1,mbkt(ji,jj),Kbb) ) )**2  )
            en(ji,jj,mbkt(ji,jj)+1) = MAX( zebot, rn_emin ) * ssmask(ji,jj)
         END_2D
         IF( ln_isfcav ) THEN
            DO_2D_OVR( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1 )     ! top friction
               zmsku = ( 2. - umask(ji-1,jj,mikt(ji,jj)) * umask(ji,jj,mikt(ji,jj)) )
               zmskv = ( 2. - vmask(ji,jj-1,mikt(ji,jj)) * vmask(ji,jj,mikt(ji,jj)) )
               !                             ! where 0.001875 = (rn_ebb0/rho0) * 0.5 = 3.75*0.5/1000.  (CAUTION CdU<0)
               zetop = - 0.001875_wp * rCdU_top(ji,jj) * SQRT(  ( zmsku*( uu(ji,jj,mikt(ji,jj),Kbb)+uu(ji-1,jj,mikt(ji,jj),Kbb) ) )**2  &
                  &                                           + ( zmskv*( vv(ji,jj,mikt(ji,jj),Kbb)+vv(ji,jj-1,mikt(ji,jj),Kbb) ) )**2  )
               ! (1._wp - tmask(ji,jj,1)) * ssmask(ji,jj) = 1 where ice shelves are present
               en(ji,jj,mikt(ji,jj)) = en(ji,jj,1)           * tmask(ji,jj,1) &
                  &                  + MAX( zetop, rn_emin ) * (1._wp - tmask(ji,jj,1)) * ssmask(ji,jj)
            END_2D
         ENDIF
         !
      ENDIF
      !
      !                     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF( ln_lc ) THEN      !  Langmuir circulation source term added to tke (Axell JGR 2002)
         !                  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
         !
         !                       !* Langmuir velocity scale
         !
         IF ( cpl_sdrftx )  THEN       ! Surface Stokes Drift available
            !                                ! Craik-Leibovich velocity scale Wlc = ( u* u_s )^1/2    with u* = (taum/rho0)^1/2
            !                                ! associated kinetic energy : 1/2 (Wlc)^2 = u* u_s
            !                                ! more precisely, it is the dot product that must be used :
            !                                !     1/2  (W_lc)^2 = MAX( u* u_s + v* v_s , 0 )   only the positive part
!!gm  ! PS: currently we don't have neither the 2 stress components at t-point !nor the angle between u* and u_s
!!gm  ! so we will overestimate the LC velocity....   !!gm I will do the work if !LC have an effect !
!
!  Projection of Stokes drift in the wind stress direction
!
            DO_2D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1 )
                  ztaui   = 0.5_wp * ( utau(ji,jj) + utau(ji-1,jj) )
                  ztauj   = 0.5_wp * ( vtau(ji,jj) + vtau(ji,jj-1) )
                  z1_norm = 1._wp / MAX( ztaui*ztaui+ztauj*ztauj, 1.e-12 ) * tmask(ji,jj,1)
                  zWlc2(ji,jj) = 0.5_wp * z1_norm * ( MAX( ut0sd(ji,jj)*ztaui + vt0sd(ji,jj)*ztauj, 0._wp ) )**2
            END_2D
         ELSE                          ! Surface Stokes drift deduced from surface stress
            !                                ! Wlc = u_s   with u_s = 0.016*U_10m, the surface stokes drift  (Axell 2002, Eq.44)
            !                                ! using |tau| = rho_air Cd |U_10m|^2 , it comes:
            !                                ! Wlc = 0.016 * [|tau|/(rho_air Cdrag) ]^1/2   and thus:
            !                                ! 1/2 Wlc^2 = 0.5 * 0.016 * 0.016 |tau| /( rho_air Cdrag )
            zcof = 0.5 * 0.016 * 0.016 / ( zrhoa * zcdrag )      ! to convert stress in 10m wind using a constant drag
            DO_2D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1 )
               zWlc2(ji,jj) = zcof * taum(ji,jj)
            END_2D
            !
         ENDIF
         !
         !                       !* Depth of the LC circulation  (Axell 2002, Eq.47)
         !                             !- LHS of Eq.47
         DO_2D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1 )
            zpelc(ji,jj,1) =  MAX( rn2b(ji,jj,1), 0._wp ) * gdepw(ji,jj,1,Kmm) * e3w(ji,jj,1,Kmm)
         END_2D
         DO_3D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, 2, jpk )
            zpelc(ji,jj,jk)  = zpelc(ji,jj,jk-1) +   &
               &          MAX( rn2b(ji,jj,jk), 0._wp ) * gdepw(ji,jj,jk,Kmm) * e3w(ji,jj,jk,Kmm)
         END_3D
         !
         !                             !- compare LHS to RHS of Eq.47
         imlc(:,:) = mbkt(A2D(nn_hls)) + 1       ! Initialization to the number of w ocean point (=2 over land)
         DO_3DS( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, jpkm1, 2, -1 )
            IF( zpelc(ji,jj,jk) > zWlc2(ji,jj) )   imlc(ji,jj) = jk
         END_3D
         !                               ! finite LC depth
         DO_2D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1 )
            zhlc(ji,jj) = gdepw(ji,jj,imlc(ji,jj),Kmm)
         END_2D
         !
         zcof = 0.016 / SQRT( zrhoa * zcdrag )
         DO_2D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1 )
            zus = SQRT( 2. * zWlc2(ji,jj) )             ! Stokes drift
            zus3(ji,jj) = MAX( 0._wp, 1._wp - zice_fra(ji,jj) ) * zus * zus * zus * tmask(ji,jj,1) ! zus > 0. ok
         END_2D
         DO_3D_OVR( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, 2, jpkm1 )                  !* TKE Langmuir circulation source term added to en
            IF ( zus3(ji,jj) /= 0._wp ) THEN
               IF ( gdepw(ji,jj,jk,Kmm) - zhlc(ji,jj) < 0 .AND. wmask(ji,jj,jk) /= 0. ) THEN
                  !                                           ! vertical velocity due to LC
                  zwlc = rn_lc * SIN( rpi * gdepw(ji,jj,jk,Kmm) / zhlc(ji,jj) )
                  !                                           ! TKE Langmuir circulation source term
                  en(ji,jj,jk) = en(ji,jj,jk) + rn_Dt * zus3(ji,jj) * ( zwlc * zwlc * zwlc ) / zhlc(ji,jj)
               ENDIF
            ENDIF
         END_3D
         !
      ENDIF
      !
      !                     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !                     !  Now Turbulent kinetic energy (output in en)
      !                     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !                     ! Resolution of a tridiagonal linear system by a "methode de chasse"
      !                     ! computation from level 2 to jpkm1  (e(1) already computed and e(jpk)=0 ).
      !                     ! zdiag : diagonal zd_up : upper diagonal zd_lw : lower diagonal
      !
      IF( nn_pdl == 1 ) THEN          !* Prandtl number = F( Ri )
         DO_3D_OVR( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, 2, jpkm1 )
            !                             ! local Richardson number
            IF (rn2b(ji,jj,jk) <= 0.0_wp) then
                zri = 0.0_wp
            ELSE
                ! This logic is to avoid divide-by-zero errors which can occur for single-precision
                ! The actual value you choose for the denominator instead of zero doesn't really
                ! matter, as long as it is very small and so triggers the same logic below with the
                ! inverse Prandtl number
                zdiv = p_sh2(ji,jj,jk) + rn_bshear
                IF (zdiv == 0.0_wp) THEN
                   zri = rn2b(ji,jj,jk) * p_avm(ji,jj,jk) / rn_bshear
                ELSE
                   zri = rn2b(ji,jj,jk) * p_avm(ji,jj,jk) / zdiv
                ENDIF
            ENDIF
            !                             ! inverse of Prandtl number
            apdlr(ji,jj,jk) = MAX(  0.1_wp,  ri_cri / MAX( ri_cri , zri )  )
         END_3D
      ENDIF
      !
      DO_3D_OVR( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, 2, jpkm1 )   !* Matrix and right hand side in en
         zcof   = zfact1 * tmask(ji,jj,jk)
         !                                   ! A minimum of 2.e-5 m2/s is imposed on TKE vertical
         !                                   ! eddy coefficient (ensure numerical stability)
         zzd_up = zcof * MAX(  p_avm(ji,jj,jk+1) + p_avm(ji,jj,jk  ) , 2.e-5_wp  )   &  ! upper diagonal
            &          /    (  e3t(ji,jj,jk  ,Kmm) * e3w(ji,jj,jk  ,Kmm)  )
         zzd_lw = zcof * MAX(  p_avm(ji,jj,jk  ) + p_avm(ji,jj,jk-1) , 2.e-5_wp  )   &  ! lower diagonal
            &          /    (  e3t(ji,jj,jk-1,Kmm) * e3w(ji,jj,jk  ,Kmm)  )
         !
         zd_up(ji,jj,jk) = zzd_up            ! Matrix (zdiag, zd_up, zd_lw)
         zd_lw(ji,jj,jk) = zzd_lw
         zdiag(ji,jj,jk) = 1._wp - zzd_lw - zzd_up + zfact2 * dissl(ji,jj,jk) * wmask(ji,jj,jk)
         !
         !                                   ! right hand side in en
         en(ji,jj,jk) = en(ji,jj,jk) + rn_Dt * (  p_sh2(ji,jj,jk)                        &   ! shear
            &                                 - p_avt(ji,jj,jk) * rn2(ji,jj,jk)          &   ! stratification
            &                                 + zfact3 * dissl(ji,jj,jk) * en(ji,jj,jk)  &   ! dissipation
            &                                ) * wmask(ji,jj,jk)
      END_3D
      !
      !                     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !                     !  Surface boundary condition on tke if
      !                     !  coupling with waves
      !                     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !
      IF ( cpl_phioc .and. ln_phioc )  THEN
         SELECT CASE (nn_bc_surf) ! Boundary Condition using surface TKE flux from waves

         CASE ( 0 ) ! Dirichlet BC
            DO_2D_OVR( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1 )    ! en(1)   = rn_ebb taum / rho0  (min value rn_emin0)
               IF ( phioc(ji,jj) < 0 )  phioc(ji,jj) = 0._wp
               en(ji,jj,1) = MAX( rn_emin0, .5 * ( 15.8 * phioc(ji,jj) / rho0 )**(2./3.) )  * tmask(ji,jj,1)
               zdiag(ji,jj,1) = 1._wp/en(ji,jj,1)  ! choose to keep coherence with former estimation of
            END_2D

         CASE ( 1 ) ! Neumann BC
            DO_2D_OVR( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1 )
               IF ( phioc(ji,jj) < 0 )  phioc(ji,jj) = 0._wp
               en(ji,jj,2)    = en(ji,jj,2) + ( rn_Dt * phioc(ji,jj) / rho0 ) /e3w(ji,jj,2,Kmm)
               en(ji,jj,1)    = en(ji,jj,2) + (2 * e3t(ji,jj,1,Kmm) * phioc(ji,jj)/rho0) / ( p_avm(ji,jj,1) + p_avm(ji,jj,2) )
               zdiag(ji,jj,2) = zdiag(ji,jj,2) + zd_lw(ji,jj,2)
               zdiag(ji,jj,1) = 1._wp
               zd_lw(ji,jj,2) = 0._wp
            END_2D

         END SELECT

      ENDIF
      !
      !                          !* Matrix inversion from level 2 (tke prescribed at level 1)
      DO_3D_OVR( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, 2, jpkm1 )                ! First recurrence : Dk = Dk - Lk * Uk-1 / Dk-1
         zdiag(ji,jj,jk) = zdiag(ji,jj,jk) - zd_lw(ji,jj,jk) * zd_up(ji,jj,jk-1) / zdiag(ji,jj,jk-1)
      END_3D
!XC : commented to allow for neumann boundary condition
!      DO_2D( 0, 0, 0, 0 )
!         zd_lw(ji,jj,2) = en(ji,jj,2) - zd_lw(ji,jj,2) * en(ji,jj,1)    ! Surface boudary conditions on tke
!      END_2D
      DO_3D_OVR( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, 2, jpkm1 )
         zd_lw(ji,jj,jk) = en(ji,jj,jk) - zd_lw(ji,jj,jk) / zdiag(ji,jj,jk-1) *zd_lw(ji,jj,jk-1)
      END_3D
      DO_2D_OVR( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1 )                          ! thrid recurrence : Ek = ( Lk - Uk * Ek+1 ) / Dk
         en(ji,jj,jpkm1) = zd_lw(ji,jj,jpkm1) / zdiag(ji,jj,jpkm1)
      END_2D
      DO_3DS_OVR( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, jpk-2, 2, -1 )
         en(ji,jj,jk) = ( zd_lw(ji,jj,jk) - zd_up(ji,jj,jk) * en(ji,jj,jk+1) ) / zdiag(ji,jj,jk)
      END_3D
      DO_3D_OVR( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, 2, jpkm1 )                ! set the minimum value of tke
         en(ji,jj,jk) = MAX( en(ji,jj,jk), rn_emin ) * wmask(ji,jj,jk)
      END_3D
      !
      ! Kolmogorov energy of dissipation (W/kg)
      !    ediss = Ce*sqrt(en)/L*en
      !    dissl = sqrt(en)/L
      IF( iom_use('ediss_k') ) THEN
         IF( .NOT. l_istiled .OR. ntile == 1 ) THEN
            ALLOCATE( ztmp(jpi,jpj,jpk) )
            ztmp(:,:,:) = 0._wp
         ENDIF
         DO_3D_OVR( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, 1, jpkm1 )
            ztmp(ji,jj,jk) = zfact3 * dissl(ji,jj,jk) * en(ji,jj,jk) * wmask(ji,jj,jk)
         END_3D
         IF( .NOT. l_istiled .OR. ntile == nijtile ) THEN
            CALL iom_put( 'ediss_k', ztmp )
            DEALLOCATE( ztmp )
         ENDIF
      ENDIF
      !
      !                            !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !                            !  TKE due to surface and internal wave breaking
      !                            !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!!gm BUG : in the exp  remove the depth of ssh !!!
!!gm       i.e. use gde3w in argument (gdepw(:,:,:,Kmm))
      !
      ! penetration is partly switched off below sea-ice if nn_eice/=0
      !
      IF( nn_etau == 1 ) THEN           !* penetration below the mixed layer (rn_efr fraction)
         DO_3D_OVR( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, 2, jpkm1 )
            en(ji,jj,jk) = en(ji,jj,jk) + rn_efr * en(ji,jj,1) * EXP( -gdepw(ji,jj,jk,Kmm) / htau(ji,jj) )   &
               &                                 * MAX( 0._wp, 1._wp - zice_fra(ji,jj) ) * wmask(ji,jj,jk) * tmask(ji,jj,1)
         END_3D
      ELSEIF( nn_etau == 2 ) THEN       !* act only at the base of the mixed layer (jk=nmln)  (rn_efr fraction)
         DO_2D_OVR( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1 )
            jk = nmln(ji,jj)
            en(ji,jj,jk) = en(ji,jj,jk) + rn_efr * en(ji,jj,1) * EXP( -gdepw(ji,jj,jk,Kmm) / htau(ji,jj) )   &
               &                                 * MAX( 0._wp, 1._wp - zice_fra(ji,jj) ) * wmask(ji,jj,jk) * tmask(ji,jj,1)
         END_2D
      ELSEIF( nn_etau == 3 ) THEN       !* penetration belox the mixed layer (HF variability)
         DO_3D_OVR( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, 2, jpkm1 )
            ztx2 = utau(ji-1,jj  ) + utau(ji,jj)
            zty2 = vtau(ji  ,jj-1) + vtau(ji,jj)
            ztau = 0.5_wp * SQRT( ztx2 * ztx2 + zty2 * zty2 ) * tmask(ji,jj,1)    ! module of the mean stress
            zdif = taum(ji,jj) - ztau                            ! mean of modulus - modulus of the mean
            zdif = rhftau_scl * MAX( 0._wp, zdif + rhftau_add )  ! apply some modifications...
            en(ji,jj,jk) = en(ji,jj,jk) + zbbrau * zdif * EXP( -gdepw(ji,jj,jk,Kmm) / htau(ji,jj) )   &
               &                        * MAX( 0._wp, 1._wp - zice_fra(ji,jj) ) * wmask(ji,jj,jk) * tmask(ji,jj,1)
         END_3D
      ENDIF
      !
   END SUBROUTINE tke_tke


   SUBROUTINE tke_avn( Kbb, Kmm, p_avm, p_avt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE tke_avn  ***
      !!
      !! ** Purpose :   Compute the vertical eddy viscosity and diffusivity
      !!
      !! ** Method  :   At this stage, en, the now TKE, is known (computed in
      !!              the tke_tke routine). First, the now mixing lenth is
      !!      computed from en and the strafification (N^2), then the mixings
      !!      coefficients are computed.
      !!              - Mixing length : a first evaluation of the mixing lengh
      !!      scales is:
      !!                      mxl = sqrt(2*en) / N
      !!      where N is the brunt-vaisala frequency, with a minimum value set
      !!      to rmxl_min (rn_mxl0) in the interior (surface) ocean.
      !!        The mixing and dissipative length scale are bound as follow :
      !!         nn_mxl=0 : mxl bounded by the distance to surface and bottom.
      !!                        zmxld = zmxlm = mxl
      !!         nn_mxl=1 : mxl bounded by the e3w and zmxld = zmxlm = mxl
      !!         nn_mxl=2 : mxl bounded such that the vertical derivative of mxl is
      !!                    less than 1 (|d/dz(mxl)|<1) and zmxld = zmxlm = mxl
      !!         nn_mxl=3 : mxl is bounded from the surface to the bottom usings
      !!                    |d/dz(xml)|<1 to obtain lup, and from the bottom to
      !!                    the surface to obtain ldown. the resulting length
      !!                    scales are:
      !!                        zmxld = sqrt( lup * ldown )
      !!                        zmxlm = min ( lup , ldown )
      !!              - Vertical eddy viscosity and diffusivity:
      !!                      avm = max( avtb, rn_ediff * zmxlm * en^1/2 )
      !!                      avt = max( avmb, pdlr * avm )
      !!      with pdlr=1 if nn_pdl=0, pdlr=1/pdl=F(Ri) otherwise.
      !!
      !! ** Action  : - avt, avm : now vertical eddy diffusivity and viscosity (w-point)
      !!----------------------------------------------------------------------
      USE zdf_oce , ONLY : en, avtb, avmb, avtb_2d   ! ocean vertical physics
      !!
      INTEGER                   , INTENT(in   ) ::   Kbb, Kmm       ! ocean time level indices
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   p_avm, p_avt   ! vertical eddy viscosity & diffusivity (w-points)
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zrn2, zraug, zcoef, zav   ! local scalars
      REAL(wp) ::   zdku,   zdkv, zsqen       !   -      -
      REAL(wp) ::   zemxl, zemlm, zemlp, zmaxice       !   -      -
      REAL(wp), DIMENSION(A2D(nn_hls),jpk) ::   zmxlm, zmxld   ! 3D workspace
      !!--------------------------------------------------------------------
      !
      !                     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !                     !  Mixing length
      !                     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !
      !                     !* Buoyancy length scale: l=sqrt(2*e/n**2)
      !
      ! initialisation of interior minimum value (avoid a 2d loop with mikt)
      zmxlm(:,:,:)  = rmxl_min
      zmxld(:,:,:)  = rmxl_min
      !
      IF(ln_sdw .AND. ln_mxhsw) THEN
         ! From Terray et al 1999 and Mellor and Blumberg 2004 it should be 0.85 and not 1.6
         zcoef       = vkarmn * ( (rn_ediff*rn_ediss)**0.25 ) / rn_ediff
         DO_2D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1 )
!            zmxlm(ji,jj,1)= vkarmn * MAX ( 1.6 * hsw(ji,jj) , 0.02 )        ! surface mixing length = F(wave height)
            zmxlm(ji,jj,1)= zcoef * MAX ( 1.6 * hsw(ji,jj) , 0.02 )        ! surface mixing length = F(wave height)
         END_2D
      ELSE
      !
         IF( ln_mxl0 ) THEN            ! surface mixing length = F(stress) : l=vkarmn*2.e5*taum/(rho0*g)
         !
            zraug = vkarmn * 2.e5_wp / ( rho0 * grav )
#if ! defined key_si3 && ! defined key_cice
            DO_2D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1 )                  ! No sea-ice
               zmxlm(ji,jj,1) =  zraug * taum(ji,jj) * tmask(ji,jj,1)
            END_2D
#else
            SELECT CASE( nn_mxlice )             ! Type of scaling under sea-ice
            !
            CASE( 0 )                      ! No scaling under sea-ice
               DO_2D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1 )
                  zmxlm(ji,jj,1) = zraug * taum(ji,jj) * tmask(ji,jj,1)
               END_2D
               !
            CASE( 1 )                      ! scaling with constant sea-ice thickness
               DO_2D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1 )
                  zmxlm(ji,jj,1) =  ( ( 1._wp - fr_i(ji,jj) ) * zraug * taum(ji,jj) + &
                     &                          fr_i(ji,jj)   * rn_mxlice           ) * tmask(ji,jj,1)
               END_2D
               !
            CASE( 2 )                      ! scaling with mean sea-ice thickness
               DO_2D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1 )
#if defined key_si3
                  zmxlm(ji,jj,1) = ( ( 1._wp - fr_i(ji,jj) ) * zraug * taum(ji,jj) + &
                     &                         fr_i(ji,jj)   * hm_i(ji,jj) * 2._wp ) * tmask(ji,jj,1)
#elif defined key_cice
                  zmaxice = MAXVAL( h_i(ji,jj,:) )
                  zmxlm(ji,jj,1) = ( ( 1._wp - fr_i(ji,jj) ) * zraug * taum(ji,jj) + &
                     &                         fr_i(ji,jj)   * zmaxice             ) * tmask(ji,jj,1)
#endif
               END_2D
               !
            CASE( 3 )                      ! scaling with max sea-ice thickness
               DO_2D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1 )
                  zmaxice = MAXVAL( h_i(ji,jj,:) )
                  zmxlm(ji,jj,1) = ( ( 1._wp - fr_i(ji,jj) ) * zraug * taum(ji,jj) + &
                     &                         fr_i(ji,jj)   * zmaxice             ) * tmask(ji,jj,1)
               END_2D
               !
            END SELECT
#endif
            !
            DO_2D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1 )
               zmxlm(ji,jj,1) = MAX( rn_mxl0, zmxlm(ji,jj,1) )
            END_2D
            !
         ELSE
            zmxlm(:,:,1) = rn_mxl0
         ENDIF
      ENDIF
      !
      DO_3D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, 2, jpkm1 )
         zrn2 = MAX( rn2(ji,jj,jk), rsmall )
         zmxlm(ji,jj,jk) = MAX(  rmxl_min,  SQRT( 2._wp * en(ji,jj,jk) / zrn2 )  )
      END_3D
      !
      !                     !* Physical limits for the mixing length
      !
      zmxld(:,:, 1 ) = zmxlm(:,:,1)   ! surface set to the minimum value
      zmxld(:,:,jpk) = rmxl_min       ! last level  set to the minimum value
      !
      SELECT CASE ( nn_mxl )
      !
 !!gm Not sure of that coding for ISF....
      ! where wmask = 0 set zmxlm == e3w(:,:,:,Kmm)
      CASE ( 0 )           ! bounded by the distance to surface and bottom
         DO_3D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, 2, jpkm1 )
            zemxl = MIN( gdepw(ji,jj,jk,Kmm) - gdepw(ji,jj,mikt(ji,jj),Kmm), zmxlm(ji,jj,jk),   &
            &            gdepw(ji,jj,mbkt(ji,jj)+1,Kmm) - gdepw(ji,jj,jk,Kmm) )
            ! wmask prevent zmxlm = 0 if jk = mikt(ji,jj)
            zmxlm(ji,jj,jk) = zemxl * wmask(ji,jj,jk)   &
               &            + MIN( zmxlm(ji,jj,jk) , e3w(ji,jj,jk,Kmm) ) * (1 - wmask(ji,jj,jk))
            zmxld(ji,jj,jk) = zemxl * wmask(ji,jj,jk)   &
               &            + MIN( zmxlm(ji,jj,jk) , e3w(ji,jj,jk,Kmm) ) * (1 - wmask(ji,jj,jk))
         END_3D
         !
      CASE ( 1 )           ! bounded by the vertical scale factor
         DO_3D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, 2, jpkm1 )
            zemxl = MIN( e3w(ji,jj,jk,Kmm), zmxlm(ji,jj,jk) )
            zmxlm(ji,jj,jk) = zemxl
            zmxld(ji,jj,jk) = zemxl
         END_3D
         !
      CASE ( 2 )           ! |dk[xml]| bounded by e3t :
         DO_3D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, 2, jpkm1 )       ! from the surface to the bottom :
            zmxlm(ji,jj,jk) =   &
               &    MIN( zmxlm(ji,jj,jk-1) + e3t(ji,jj,jk-1,Kmm), zmxlm(ji,jj,jk) )
         END_3D
         DO_3DS( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, jpkm1, 2, -1 )   ! from the bottom to the surface :
            zemxl = MIN( zmxlm(ji,jj,jk+1) + e3t(ji,jj,jk+1,Kmm), zmxlm(ji,jj,jk) )
            zmxlm(ji,jj,jk) = zemxl
            zmxld(ji,jj,jk) = zemxl
         END_3D
         !
      CASE ( 3 )           ! lup and ldown, |dk[xml]| bounded by e3t :
         DO_3D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, 2, jpkm1 )        ! from the surface to the bottom : lup
            zmxld(ji,jj,jk) =    &
               &    MIN( zmxld(ji,jj,jk-1) + e3t(ji,jj,jk-1,Kmm), zmxlm(ji,jj,jk) )
         END_3D
         DO_3DS( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, jpkm1, 2, -1 )   ! from the bottom to the surface : ldown
            zmxlm(ji,jj,jk) =   &
               &    MIN( zmxlm(ji,jj,jk+1) + e3t(ji,jj,jk+1,Kmm), zmxlm(ji,jj,jk) )
         END_3D
         DO_3D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, 2, jpkm1 )
            zemlm = MIN ( zmxld(ji,jj,jk),  zmxlm(ji,jj,jk) )
            zemlp = SQRT( zmxld(ji,jj,jk) * zmxlm(ji,jj,jk) )
            zmxlm(ji,jj,jk) = zemlm
            zmxld(ji,jj,jk) = zemlp
         END_3D
         !
      END SELECT
      !
      !                     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !                     !  Vertical eddy viscosity and diffusivity  (avm and avt)
      !                     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      DO_3D_OVR( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, 1, jpkm1 )   !* vertical eddy viscosity & diffivity at w-points
         zsqen = SQRT( en(ji,jj,jk) )
         zav   = rn_ediff * zmxlm(ji,jj,jk) * zsqen
         p_avm(ji,jj,jk) = MAX( zav,                  avmb(jk) ) * wmask(ji,jj,jk)
         p_avt(ji,jj,jk) = MAX( zav, avtb_2d(ji,jj) * avtb(jk) ) * wmask(ji,jj,jk)
         dissl(ji,jj,jk) = zsqen / zmxld(ji,jj,jk)
      END_3D
      !
      !
      IF( nn_pdl == 1 ) THEN          !* Prandtl number case: update avt
         DO_3D_OVR( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, 2, jpkm1 )
            p_avt(ji,jj,jk)   = MAX( apdlr(ji,jj,jk) * p_avt(ji,jj,jk), avtb_2d(ji,jj) * avtb(jk) ) * wmask(ji,jj,jk)
         END_3D
      ENDIF
      !
      IF(sn_cfctl%l_prtctl) THEN
         CALL prt_ctl( tab3d_1=CASTDP(en) , clinfo1=' tke  - e: ', tab3d_2=CASTDP(p_avt), clinfo2=' t: ' )
         CALL prt_ctl( tab3d_1=CASTDP(p_avm), clinfo1=' tke  - m: ' )
      ENDIF
      !
   END SUBROUTINE tke_avn


   SUBROUTINE zdf_tke_init( Kmm )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdf_tke_init  ***
      !!
      !! ** Purpose :   Initialization of the vertical eddy diffivity and
      !!              viscosity when using a tke turbulent closure scheme
      !!
      !! ** Method  :   Read the namzdf_tke namelist and check the parameters
      !!              called at the first timestep (nit000)
      !!
      !! ** input   :   Namlist namzdf_tke
      !!
      !! ** Action  :   Increase by 1 the nstop flag is setting problem encounter
      !!----------------------------------------------------------------------
      USE zdf_oce , ONLY : ln_zdfiwm   ! Internal Wave Mixing flag
      !!
      INTEGER, INTENT(in) ::   Kmm          ! time level index
      INTEGER             ::   ji, jj, jk   ! dummy loop indices
      INTEGER             ::   ios
      !!
      NAMELIST/namzdf_tke/ rn_ediff, rn_ediss , rn_ebb   , rn_emin  ,  &
         &                 rn_emin0, rn_bshear, nn_mxl   , ln_mxl0  ,  &
         &                 rn_mxl0 , nn_mxlice, rn_mxlice,             &
         &                 nn_pdl  , ln_lc    , rn_lc    ,             &
         &                 rn_htau_scaling    ,                        &
         &                 nn_etau , nn_htau  , rn_efr   , nn_eice  ,  &
         &                 nn_bc_surf, nn_bc_bot, ln_mxhsw
      !!----------------------------------------------------------------------
      !
      READ  ( numnam_ref, namzdf_tke, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namzdf_tke in reference namelist' )

      READ  ( numnam_cfg, namzdf_tke, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 ) CALL ctl_nam ( ios , 'namzdf_tke in configuration namelist' )
      IF(lwm) WRITE ( numond, namzdf_tke )
      !
      ri_cri   = 2._wp    / ( 2._wp + rn_ediss / rn_ediff )   ! resulting critical Richardson number
      !
      IF(lwp) THEN                    !* Control print
         WRITE(numout,*)
         WRITE(numout,*) 'zdf_tke_init : tke turbulent closure scheme - initialisation'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namzdf_tke : set tke mixing parameters'
         WRITE(numout,*) '      coef. to compute avt                        rn_ediff  = ', rn_ediff
         WRITE(numout,*) '      Kolmogoroff dissipation coef.               rn_ediss  = ', rn_ediss
         WRITE(numout,*) '      tke surface input coef.                     rn_ebb    = ', rn_ebb
         WRITE(numout,*) '      minimum value of tke                        rn_emin   = ', rn_emin
         WRITE(numout,*) '      surface minimum value of tke                rn_emin0  = ', rn_emin0
         WRITE(numout,*) '      prandl number flag                          nn_pdl    = ', nn_pdl
         WRITE(numout,*) '      background shear (>0)                       rn_bshear = ', rn_bshear
         WRITE(numout,*) '      mixing length type                          nn_mxl    = ', nn_mxl
         WRITE(numout,*) '         surface mixing length = F(stress) or not    ln_mxl0   = ', ln_mxl0
         WRITE(numout,*) '         surface  mixing length minimum value        rn_mxl0   = ', rn_mxl0
         IF( ln_mxl0 ) THEN
            WRITE(numout,*) '      type of scaling under sea-ice               nn_mxlice = ', nn_mxlice
            IF( nn_mxlice == 1 ) &
            WRITE(numout,*) '      ice thickness when scaling under sea-ice    rn_mxlice = ', rn_mxlice
            SELECT CASE( nn_mxlice )             ! Type of scaling under sea-ice
            CASE( 0 )   ;   WRITE(numout,*) '   ==>>>   No scaling under sea-ice'
            CASE( 1 )   ;   WRITE(numout,*) '   ==>>>   scaling with constant sea-ice thickness'
            CASE( 2 )   ;   WRITE(numout,*) '   ==>>>   scaling with mean sea-ice thickness'
            CASE( 3 )   ;   WRITE(numout,*) '   ==>>>   scaling with max sea-ice thickness'
            CASE DEFAULT
               CALL ctl_stop( 'zdf_tke_init: wrong value for nn_mxlice, should be 0,1,2,3 or 4')
            END SELECT
            IF     ( (nn_mxlice>0).AND.(nn_ice==0) ) THEN
               CALL ctl_stop( 'zdf_tke_init: with no ice at all, nn_mxlice must be 0 ') 
            ELSEIF ( (nn_mxlice>1).AND.(nn_ice==1) ) THEN
               CALL ctl_stop( 'zdf_tke_init: with no ice model, nn_mxlice must be 0 or 1')
            ENDIF
         ENDIF
         WRITE(numout,*) '      Langmuir cells parametrization              ln_lc     = ', ln_lc
         WRITE(numout,*) '         coef to compute vertical velocity of LC     rn_lc  = ', rn_lc
         IF ( cpl_phioc .and. ln_phioc )  THEN
            SELECT CASE( nn_bc_surf)             ! Type of scaling under sea-ice
            CASE( 0 )   ;   WRITE(numout,*) '  nn_bc_surf=0 ==>>> DIRICHLET SBC using surface TKE flux from waves'
            CASE( 1 )   ;   WRITE(numout,*) '  nn_bc_surf=1 ==>>> NEUMANN SBC using surface TKE flux from waves'
            END SELECT
         ENDIF
         WRITE(numout,*) '      test param. to add tke induced by wind      nn_etau   = ', nn_etau
         WRITE(numout,*) '          type of tke penetration profile            nn_htau   = ', nn_htau
         WRITE(numout,*) '          scaling factor for tke penetration depth   rn_htau_scaling   = ', rn_htau_scaling
         WRITE(numout,*) '          fraction of TKE that penetrates            rn_efr    = ', rn_efr
         WRITE(numout,*) '      langmuir & surface wave breaking under ice  nn_eice = ', nn_eice
         SELECT CASE( nn_eice )
         CASE( 0 )   ;   WRITE(numout,*) '   ==>>>   no impact of ice cover on langmuir & surface wave breaking'
         CASE( 1 )   ;   WRITE(numout,*) '   ==>>>   weigthed by 1-TANH( fr_i(:,:) * 10 )'
         CASE( 2 )   ;   WRITE(numout,*) '   ==>>>   weighted by 1-fr_i(:,:)'
         CASE( 3 )   ;   WRITE(numout,*) '   ==>>>   weighted by 1-MIN( 1, 4 * fr_i(:,:) )'
         CASE DEFAULT
            CALL ctl_stop( 'zdf_tke_init: wrong value for nn_eice, should be 0,1,2, or 3')
         END SELECT
         WRITE(numout,*)
         WRITE(numout,*) '   ==>>>   critical Richardson nb with your parameters  ri_cri = ', ri_cri
         WRITE(numout,*)
      ENDIF
      !
      IF( ln_zdfiwm ) THEN          ! Internal wave-driven mixing
         rn_emin  = 1.e-10_wp             ! specific values of rn_emin & rmxl_min are used
         rmxl_min = 1.e-03_wp             ! associated avt minimum = molecular salt diffusivity (10^-9 m2/s)
         IF(lwp) WRITE(numout,*) '   ==>>>   Internal wave-driven mixing case:   force   rn_emin = 1.e-10 and rmxl_min = 1.e-3'
      ELSE                          ! standard case : associated avt minimum = molecular viscosity (10^-6 m2/s)
         rmxl_min = 1.e-6_wp / ( rn_ediff * SQRT( rn_emin ) )    ! resulting minimum length to recover molecular viscosity
         IF(lwp) WRITE(numout,*) '   ==>>>   minimum mixing length with your parameters rmxl_min = ', rmxl_min
      ENDIF
      !
      !                              ! allocate tke arrays
      IF( zdf_tke_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'zdf_tke_init : unable to allocate arrays' )
      !
      !                               !* Check of some namelist values
      IF( nn_mxl  < 0   .OR.  nn_mxl  > 3 )   CALL ctl_stop( 'bad flag: nn_mxl is  0, 1, 2 or 3' )
      IF( nn_pdl  < 0   .OR.  nn_pdl  > 1 )   CALL ctl_stop( 'bad flag: nn_pdl is  0 or 1' )
      IF( ( nn_htau < 0   .OR.  nn_htau > 1 ) .AND. nn_htau .NE. 4 .AND. nn_htau .NE. 5 )   CALL ctl_stop( 'bad flag: nn_htau is 0, 1 , 4 or 5 ' )
      IF( nn_etau == 3 .AND. .NOT. ln_cpl )   CALL ctl_stop( 'nn_etau == 3 : HF taum only known in coupled mode' )
      !
      IF( ln_mxl0 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   ==>>>   use a surface mixing length = F(stress) :   set rn_mxl0 = rmxl_min'
         rn_mxl0 = rmxl_min
      ENDIF
      !                               !* depth of penetration of surface tke
      IF( nn_etau /= 0 ) THEN
         SELECT CASE( nn_htau )             ! Choice of the depth of penetration
         CASE( 0 )                                 ! constant depth penetration (here 10 meters)
            htau(:,:) = rn_htau_scaling*10._wp 
         CASE( 1 )                                 ! F(latitude) : 0.5m to 30m poleward of 40 degrees
            htau(:,:) = MAX(  0.5_wp, MIN( 30._wp, 45._wp* ABS( SIN( rpi/180._wp * gphit(:,:) ) ) )   )
         CASE( 4 )                                 ! F(latitude) : 0.5m to 10m/30m poleward of 13/40 degrees north/south
            DO_2D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1 )
                  IF( gphit(ji,jj) <= 0._wp ) THEN
                     htau(ji,jj) = MAX(  0.5_wp, MIN( 30._wp, 45._wp* rn_htau_scaling*ABS( SIN( rpi/180._wp * gphit(ji,jj) ) ) )   )
                  ELSE
                     htau(ji,jj) = MAX(  0.5_wp, MIN( 10._wp, 45._wp* rn_htau_scaling*ABS( SIN( rpi/180._wp * gphit(ji,jj) ) ) )   )
                  ENDIF
            END_2D
         CASE( 5 )                                 ! Variation on case 4 with a steeper ramp further south in Southern Hemisphere
            DO_2D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1 )
                  htau(ji,jj) = MAX(  0.5_wp, MIN( 10._wp, 45._wp* rn_htau_scaling*ABS( SIN( rpi/180._wp * gphit(ji,jj) ) ) )   )
                  IF( gphit(ji,jj) <= -40._wp ) THEN
                      htau(ji,jj) = htau(ji,jj) + MIN( 20._wp, 135._wp * ABS( SIN( rpi/180._wp * (gphit(ji,jj) + 40.0) ) ) )
                  ENDIF
            END_2D
         END SELECT
      ENDIF
      !                                !* read or initialize all required files
      CALL tke_rst( nit000, 'READ' )      ! (en, avt_k, avm_k, dissl)
      !
   END SUBROUTINE zdf_tke_init


   SUBROUTINE tke_rst( kt, cdrw )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE tke_rst  ***
      !!
      !! ** Purpose :   Read or write TKE file (en) in restart file
      !!
      !! ** Method  :   use of IOM library
      !!                if the restart does not contain TKE, en is either
      !!                set to rn_emin or recomputed
      !!----------------------------------------------------------------------
      USE zdf_oce , ONLY : en, avt_k, avm_k   ! ocean vertical physics
      !!
      INTEGER         , INTENT(in) ::   kt     ! ocean time-step
      CHARACTER(len=*), INTENT(in) ::   cdrw   ! "READ"/"WRITE" flag
      !
      INTEGER ::   jit, jk              ! dummy loop indices
      INTEGER ::   id1, id2, id3, id4   ! local integers
      !!----------------------------------------------------------------------
      !
      IF( TRIM(cdrw) == 'READ' ) THEN        ! Read/initialise
         !                                   ! ---------------
         IF( ln_rstart ) THEN                   !* Read the restart file
            id1 = iom_varid( numror, 'en'   , ldstop = .FALSE. )
            id2 = iom_varid( numror, 'avt_k', ldstop = .FALSE. )
            id3 = iom_varid( numror, 'avm_k', ldstop = .FALSE. )
            id4 = iom_varid( numror, 'dissl', ldstop = .FALSE. )
            !
            IF( MIN( id1, id2, id3, id4 ) > 0 ) THEN      ! fields exist
               CALL iom_get( numror, jpdom_auto, 'en'   , en    )
               CALL iom_get( numror, jpdom_auto, 'avt_k', avt_k )
               CALL iom_get( numror, jpdom_auto, 'avm_k', avm_k )
               CALL iom_get( numror, jpdom_auto, 'dissl', dissl )
            ELSE                                          ! start TKE from rest
               IF(lwp) WRITE(numout,*)
               IF(lwp) WRITE(numout,*) '   ==>>>   previous run without TKE scheme, set en to background values'
               en   (:,:,:) = rn_emin * wmask(:,:,:)
               dissl(:,:,:) = 1.e-12_wp
               ! avt_k, avm_k already set to the background value in zdf_phy_init
            ENDIF
         ELSE                                   !* Start from rest
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) '   ==>>>   start from rest: set en to the background value'
            en   (:,:,:) = rn_emin * wmask(:,:,:)
            dissl(:,:,:) = 1.e-12_wp
            ! avt_k, avm_k already set to the background value in zdf_phy_init
         ENDIF
         !
      ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN   ! Create restart file
         !                                   ! -------------------
         IF(lwp) WRITE(numout,*) '---- tke_rst ----'
         CALL iom_rstput( kt, nitrst, numrow, 'en'   , en    )
         CALL iom_rstput( kt, nitrst, numrow, 'avt_k', avt_k )
         CALL iom_rstput( kt, nitrst, numrow, 'avm_k', avm_k )
         CALL iom_rstput( kt, nitrst, numrow, 'dissl', dissl )
         !
      ENDIF
      !
   END SUBROUTINE tke_rst

   !!======================================================================
END MODULE zdftke
