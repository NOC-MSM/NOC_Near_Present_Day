MODULE in_out_manager
   !!======================================================================
   !!                       ***  MODULE  in_out_manager  ***
   !! I/O manager utilities : Defines run parameters together with logical units
   !!=====================================================================
   !! History :   1.0  !  2002-06  (G. Madec)   original code
   !!             2.0  !  2006-07  (S. Masson)  iom, add ctl_stop, ctl_warn
   !!             3.0  !  2008-06  (G. Madec)   add ctmp4 to ctmp10
   !!             3.2  !  2009-08  (S. MAsson)  add new ctl_opn
   !!             3.3  !  2010-10  (A. Coward)  add NetCDF4 usage
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   USE par_oce       ! ocean parameter
   USE nc4interface  ! NetCDF4 interface

   IMPLICIT NONE
   PUBLIC

   !!----------------------------------------------------------------------
   !!                   namrun namelist parameters
   !!----------------------------------------------------------------------
   CHARACTER(lc) ::   cn_exp           !: experiment name used for output filename
   CHARACTER(lc) ::   cn_ocerst_in     !: suffix of ocean restart name (input)
   CHARACTER(lc) ::   cn_ocerst_indir  !: restart input directory
   CHARACTER(lc) ::   cn_ocerst_out    !: suffix of ocean restart name (output)
   CHARACTER(lc) ::   cn_ocerst_outdir !: restart output directory
   LOGICAL       ::   ln_rstart        !: start from (F) rest or (T) a restart file
   LOGICAL       ::   ln_reset_ts      !: use TS from initial condition file? (F) no or (T) yes
   LOGICAL       ::   ln_rst_list      !: output restarts at list of times (T) or by frequency (F)
   LOGICAL       ::   ln_rst_eos       !: check equation of state used for the restart is consistent with model
   INTEGER       ::   nn_rstctl        !: control of the time step (0, 1 or 2)
   INTEGER       ::   nn_rstssh   = 0  !: hand made initilization of ssh or not (1/0)
   INTEGER       ::   nn_it000         !: index of the first time step
   INTEGER       ::   nn_itend         !: index of the last time step
   INTEGER       ::   nn_date0         !: initial calendar date aammjj
   INTEGER       ::   nn_time0         !: initial time of day in hhmm
   INTEGER       ::   nn_leapy         !: Leap year calendar flag (0/1 or 30)
   INTEGER       ::   nn_istate        !: initial state output flag (0/1)
   INTEGER       ::   nn_write         !: model standard output frequency
   INTEGER       ::   nn_stock         !: restart file frequency
   INTEGER, DIMENSION(10) :: nn_stocklist  !: restart dump times
   LOGICAL       ::   ln_mskland       !: mask land points in NetCDF outputs (costly: + ~15%)
   LOGICAL       ::   ln_rstdate       !: T=> stamp output restart files with date instead of timestep
   LOGICAL       ::   ln_cfmeta        !: output additional data to netCDF files required for compliance with the CF metadata standard
   LOGICAL       ::   ln_clobber       !: clobber (overwrite) an existing file
   INTEGER       ::   nn_chunksz       !: chunksize (bytes) for NetCDF file (works only with iom_nf90 routines)
   LOGICAL       ::   ln_xios_read     !: use xios to read single file restart
   INTEGER       ::   nn_wxios         !: write resart using xios 0 - no, 1 - single, 2 - multiple file output
   INTEGER       ::   nn_no            !: Assimilation cycle

#if defined key_netcdf4
   !!----------------------------------------------------------------------
   !!                   namnc4 namelist parameters                         (key_netcdf4)
   !!----------------------------------------------------------------------
   ! The following four values determine the partitioning of the output fields
   ! into netcdf4 chunks. They are unrelated to the nn_chunk_sz setting which is
   ! for runtime optimisation. The individual netcdf4 chunks can be optionally
   ! gzipped (recommended) leading to significant reductions in I/O volumes
   !                         !!!**  variables only used with iom_nf90 routines and key_netcdf4 **
   INTEGER ::   nn_nchunks_i   !: number of chunks required in the i-dimension
   INTEGER ::   nn_nchunks_j   !: number of chunks required in the j-dimension
   INTEGER ::   nn_nchunks_k   !: number of chunks required in the k-dimension
   INTEGER ::   nn_nchunks_t   !: number of chunks required in the t-dimension
   LOGICAL ::   ln_nc4zip      !: netcdf4 usage: (T) chunk and compress output using the HDF5 sublayers of netcdf4
   !                           !                 (F) ignore chunking request and use the netcdf4 library
   !                           !                     to produce netcdf3-compatible files
#endif

!$AGRIF_DO_NOT_TREAT
   TYPE(snc4_ctl)     :: snc4set        !: netcdf4 chunking control structure (always needed for decision making)
!$AGRIF_END_DO_NOT_TREAT


   !! conversion of DOCTOR norm namelist name into model name
   !! (this should disappear in a near futur)

   CHARACTER(lc) ::   cexper                      !: experiment name used for output filename
   INTEGER       ::   nrstdt                      !: control of the time step (0, 1 or 2)
   INTEGER       ::   nit000                      !: index of the first time step
   INTEGER       ::   nitend                      !: index of the last time step
   INTEGER       ::   ndate0                      !: initial calendar date aammjj
   INTEGER       ::   nleapy                      !: Leap year calendar flag (0/1 or 30)
   INTEGER       ::   ninist                      !: initial state output flag (0/1)

   !!----------------------------------------------------------------------
   !! was in restart but moved here because of the OFF line... better solution should be found...
   !!----------------------------------------------------------------------
   INTEGER ::   nitrst                !: time step at which restart file should be written
   LOGICAL ::   lrst_oce              !: logical to control the oce restart write
   LOGICAL ::   lrst_ice              !: logical to control the ice restart write
   LOGICAL ::   lrst_abl              !: logical to control the abl restart write
   INTEGER ::   numror = 0            !: logical unit for ocean restart (read). Init to 0 is needed for SAS (in daymod.F90)
   INTEGER ::   numrir = 0            !: logical unit for ice   restart (read)
   INTEGER ::   numrar = 0            !: logical unit for abl   restart (read)
   INTEGER ::   numrow = 0            !: logical unit for ocean restart (write)
   INTEGER ::   numriw = 0            !: logical unit for ice   restart (write)
   INTEGER ::   numraw = 0            !: logical unit for abl   restart (write)
   INTEGER ::   numrtr = 0            !: trc restart (read )
   INTEGER ::   numrtw = 0            !: trc restart (write )
   INTEGER ::   numrsr = 0            !: logical unit for sed restart (read)
   INTEGER ::   numrsw = 0            !: logical unit for sed restart (write)

   INTEGER ::   nrst_lst              !: number of restart to output next

   !!----------------------------------------------------------------------
   !!                    output monitoring
   !!----------------------------------------------------------------------
   TYPE :: sn_ctl                !: structure for control over output selection
      LOGICAL :: l_runstat = .FALSE.  !: Produce/do not produce run.stat file (T/F)
      LOGICAL :: l_trcstat = .FALSE.  !: Produce/do not produce tracer.stat file (T/F)
      LOGICAL :: l_oceout  = .FALSE.  !: Produce all ocean.outputs    (T) or just one (F)
      LOGICAL :: l_layout  = .FALSE.  !: Produce all layout.dat files (T) or just one (F)
      LOGICAL :: l_prtctl  = .FALSE.  !: Produce/do not produce mpp.output_XXXX files (T/F)
      LOGICAL :: l_prttrc  = .FALSE.  !: Produce/do not produce mpp.top.output_XXXX files (T/F)
      LOGICAL :: l_oasout  = .FALSE.  !: Produce/do not write oasis setup info to ocean.output (T/F)
                                      !  Optional subsetting of processor report files
                                      !  Default settings of 0/1000000/1 should ensure all areas report.
                                      !  Set to a more restrictive range to select specific areas
      INTEGER :: procmin   = 0        !: Minimum narea to output
      INTEGER :: procmax   = 1000000  !: Maximum narea to output
      INTEGER :: procincr  = 1        !: narea increment to output
      INTEGER :: ptimincr  = 1        !: timestep increment to output (time.step and run.stat)
   END TYPE
   TYPE(sn_ctl), SAVE :: sn_cfctl     !: run control structure for selective output, must have SAVE for default init. of sn_ctl
   LOGICAL ::   ln_timing        !: run control for timing
   LOGICAL ::   ln_diacfl        !: flag whether to create CFL diagnostics
   INTEGER ::   nn_ictls         !: Start i indice for the SUM control
   INTEGER ::   nn_ictle         !: End   i indice for the SUM control
   INTEGER ::   nn_jctls         !: Start j indice for the SUM control
   INTEGER ::   nn_jctle         !: End   j indice for the SUM control
   INTEGER ::   nn_isplt         !: number of processors following i
   INTEGER ::   nn_jsplt         !: number of processors following j

   !!----------------------------------------------------------------------
   !!                        logical units
   !!----------------------------------------------------------------------
   INTEGER ::   numstp          =   -1      !: logical unit for time step
   INTEGER ::   numtime         =   -1      !: logical unit for timing
   INTEGER ::   numout          =    6      !: logical unit for output print; Set to stdout to ensure any
   INTEGER ::   numnul          =   -1      !: logical unit for /dev/null
      !                                     !  early output can be collected; do not change
   INTEGER ::   numond          =   -1      !: logical unit for Output Namelist Dynamics
   INTEGER ::   numoni          =   -1      !: logical unit for Output Namelist Ice
   INTEGER ::   numevo_ice      =   -1      !: logical unit for ice variables (temp. evolution)
   INTEGER ::   numrun          =   -1      !: logical unit for run statistics
   INTEGER ::   numdct_in       =   -1      !: logical unit for transports computing
   INTEGER ::   numdct_vol      =   -1      !: logical unit for volume transports output
   INTEGER ::   numdct_heat     =   -1      !: logical unit for heat   transports output
   INTEGER ::   numdct_salt     =   -1      !: logical unit for salt   transports output
   INTEGER ::   numfl           =   -1      !: logical unit for floats ascii output
   INTEGER ::   numflo          =   -1      !: logical unit for floats ascii output
      !
   CHARACTER(LEN=:), ALLOCATABLE :: numnam_ref      !: character buffer for reference namelist
   CHARACTER(LEN=:), ALLOCATABLE :: numnam_cfg      !: character buffer for configuration specific namelist
   CHARACTER(LEN=:), ALLOCATABLE :: numnam_ice_ref  !: character buffer for ice reference namelist
   CHARACTER(LEN=:), ALLOCATABLE :: numnam_ice_cfg  !: character buffer for ice configuration specific namelist

   !!----------------------------------------------------------------------
   !!                          Run control
   !!----------------------------------------------------------------------
   INTEGER       ::   no_print = 0          !: optional argument of fld_fill (if present, suppress some control print)
   INTEGER       ::   nstop = 0             !: error flag (=number of reason for a premature stop run)
!$AGRIF_DO_NOT_TREAT
   INTEGER       ::   ngrdstop = -1         !: grid number having nstop > 1
!$AGRIF_END_DO_NOT_TREAT
   INTEGER       ::   nwarn = 0             !: warning flag (=number of warning found during the run)
   CHARACTER(lc) ::   ctmp1, ctmp2, ctmp3   !: temporary characters 1 to 3
   CHARACTER(lc) ::   ctmp4, ctmp5, ctmp6   !: temporary characters 4 to 6
   CHARACTER(lc) ::   ctmp7, ctmp8, ctmp9   !: temporary characters 7 to 9
   CHARACTER(lc) ::   ctmp10                !: temporary character 10
   LOGICAL       ::   lwm      = .FALSE.    !: boolean : true on the 1st processor only (always)
   LOGICAL       ::   lwp      = .FALSE.    !: boolean : true on the 1st processor only .OR. sn_cfctl%l_oceout=T
   LOGICAL       ::   lsp_area = .TRUE.     !: to make a control print over a specific area
   CHARACTER(LEN=lc) ::   cxios_context     !: context name used in xios
   CHARACTER(LEN=lc) ::   cr_ocerst_cxt     !: context name used in xios to read OCE restart
   CHARACTER(LEN=lc) ::   cw_ocerst_cxt     !: context name used in xios to write OCE restart file
   CHARACTER(LEN=lc) ::   cr_icerst_cxt     !: context name used in xios to read SI3 restart
   CHARACTER(LEN=lc) ::   cw_icerst_cxt     !: context name used in xios to write SI3 restart file
   CHARACTER(LEN=lc) ::   cr_ablrst_cxt     !: context name used in xios to read ABL restart
   CHARACTER(LEN=lc) ::   cw_ablrst_cxt     !: context name used in xios to write ABL restart file
   CHARACTER(LEN=lc) ::   cr_toprst_cxt     !: context name used in xios to read TOP restart
   CHARACTER(LEN=lc) ::   cw_toprst_cxt     !: context name used in xios to write TOP restart file
   CHARACTER(LEN=lc) ::   cr_sedrst_cxt     !: context name used in xios to read SEDIMENT restart
   CHARACTER(LEN=lc) ::   cw_sedrst_cxt     !: context name used in xios to write SEDIMENT restart file




   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: in_out_manager.F90 14553 2021-02-26 17:01:43Z gsamson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!=====================================================================
END MODULE in_out_manager
