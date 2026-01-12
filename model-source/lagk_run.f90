! !     ==============================================================================
! !     WRC MODIFICATION NOTES - RESTART CAPABILITY ADDITIONS (January 2026)
! !     ==============================================================================
! !     The following modifications enable Lag-K routing restart capability
! !     for operational forecasting workflows. Changes allow:
! !     - Warm start from previous C array states (carryover values)
! !     - Continuation of routing operations across simulation boundaries
! !     - Support for ensemble forecasting from same initial conditions
! !     
! !     Key additions:
! !     1. New input arguments: c_array_in (100x n_hrus), use_c_array_restart (flag)
! !     2. New output argument: c_array_out (100x n_hrus) for state persistence
! !     3. Conditional C array initialization logic
! !     4. Final state capture for each HRU/tributary
! !     
! !     Technical details:
! !     - The C array contains carryover values used by flag7 and fka7 subroutines
! !     - C array stores: initial conditions, inflow, outflow, storage, carryover
! !     - pin7 always initializes P and C arrays from parameters
! !     - If restarting, C array is overridden with saved state after pin7
! !     - c_cpy is used during routing to preserve original C state
! !     ==============================================================================
! !     WRC: RESTART WORKFLOW TECHNICAL NOTES
! !     =================================================================
! !     Cold Start Workflow (use_c_array_restart = 0):
! !       1. pin7 initializes P and C arrays from parameters
! !       2. C array uses ico, iinfl, ioutfl, istor as initial conditions
! !       3. flag7 and fka7 perform routing with default initialization
! !       4. c_array_out captures final state for potential future restart
! !
! !     Warm Start Workflow (use_c_array_restart = 1):
! !       1. pin7 initializes P array from parameters (required for routing tables)
! !       2. pin7 also initializes C array, but this is immediately overridden
! !       3. C array loaded from c_array_in (previous simulation's c_array_out)
! !       4. flag7 and fka7 perform routing with restored carryover state
! !       5. c_array_out captures new final state for next restart
! !
! !     Continuous Simulation Chain:
! !       Run 1: lagk(..., c_array_dummy, 0, ..., c_array_1)
! !       Run 2: lagk(..., c_array_1, 1, ..., c_array_2)
! !       Run 3: lagk(..., c_array_2, 1, ..., c_array_3)
! !       etc.
! !
! !     Why pin7 is always called:
! !       - P array contains routing tables (lag/Q, K/Q, storage curves)
! !       - P array structure depends on jlag, jk, and timestep parameters
! !       - P array must be reconstructed each run even for warm starts
! !       - Only C array (carryover values) persists between runs
! !
! !     Array sizing note:
! !       - C array is fixed at 100 elements (oversized for safety)
! !       - Actual used length stored in c(1) by pin7
! !       - Entire 100-element array saved/restored to preserve all state
! !     =================================================================
! !     WRC: USAGE EXAMPLES
! !     =================================================================
! !     Example 1 - Cold Start (Original Behavior):
! !       
! !       double precision, dimension(100, n_hrus) :: c_dummy, c_out
! !       c_dummy = 0.0d0  ! Not used, but required argument
! !       
! !       call lagk(n_hrus, ita, itb, &
! !           lagtbl_a, lagtbl_b, lagtbl_c, lagtbl_d, &
! !           ktbl_a, ktbl_b, ktbl_c, ktbl_d, &
! !           lagk_lagmax, lagk_kmax, lagk_qmax, &
! !           lagk_lagmin, lagk_kmin, lagk_qmin, &
! !           ico, iinfl, ioutfl, istor, &
! !           c_dummy, 0, &  ! Cold start: use_c_array_restart = 0
! !           qa, sim_length, return_states, &
! !           lagk_out, co_st_out, inflow_st_out, storage_st_out, &
! !           c_out)  ! Save final state
! !
! !     Example 2 - Warm Start (Restart from Previous Run):
! !       
! !       double precision, dimension(100, n_hrus) :: c_previous, c_new
! !       ! c_previous contains c_array_out from previous simulation
! !       
! !       call lagk(n_hrus, ita, itb, &
! !           lagtbl_a, lagtbl_b, lagtbl_c, lagtbl_d, &
! !           ktbl_a, ktbl_b, ktbl_c, ktbl_d, &
! !           lagk_lagmax, lagk_kmax, lagk_qmax, &
! !           lagk_lagmin, lagk_kmin, lagk_qmin, &
! !           ico, iinfl, ioutfl, istor, &
! !           c_previous, 1, &  ! Warm start: use_c_array_restart = 1
! !           qa, sim_length, return_states, &
! !           lagk_out, co_st_out, inflow_st_out, storage_st_out, &
! !           c_new)  ! Save new final state
! !
! !     Example 3 - Forecast Ensemble from Same Initial Conditions:
! !       
! !       double precision, dimension(100, n_hrus) :: c_initial
! !       ! Run historical simulation to get c_initial
! !       
! !       do ensemble_member = 1, n_ensemble
! !         ! Each ensemble member starts from same c_initial state
! !         call lagk(..., c_initial, 1, ..., c_temp)
! !         ! Store results for this ensemble member
! !       end do
! !     =================================================================
! !     WRC: COMPARISON - ORIGINAL vs. RESTART-ENABLED VERSION
! !     =================================================================
! !     Feature                    | Original        | Restart-Enabled
! !     ---------------------------+-----------------+------------------
! !     Subroutine arguments       | 13              | 16 (+3)
! !     C array initialization     | pin7 only       | pin7 or c_array_in
! !     State persistence          | None            | c_array_out
! !     Warm start capability      | No              | Yes
! !     Forecast workflows         | Limited         | Full support
! !     Backward compatibility     | N/A             | Yes (flag=0)
! !     
! !     New Arguments:
! !       - c_array_in (input, 100 x n_hrus)
! !       - use_c_array_restart (input, integer flag)
! !       - c_array_out (output, 100 x n_hrus)
! !     
! !     Modified Logic:
! !       - Conditional C array override after pin7
! !       - Final C array state capture before end of HRU loop
! !     =================================================================

subroutine lagk(n_hrus, ita, itb, &
    lagtbl_a_in, lagtbl_b_in, lagtbl_c_in, lagtbl_d_in,&
    ktbl_a_in, ktbl_b_in, ktbl_c_in, ktbl_d_in, &
    lagk_lagmax_in, lagk_kmax_in, lagk_qmax_in, &
    lagk_lagmin_in, lagk_kmin_in, lagk_qmin_in, &
    ico_in, iinfl_in, ioutfl_in, istor_in, &
    c_array_in, use_c_array_restart, & ! WRC: NEW - restart state inputs
    qa_in, sim_length, &
    return_states, &
    lagk_out, co_st_out, &
    inflow_st_out, storage_st_out, &
    c_array_out) ! WRC: NEW - restart state output

    ! !There are three subroutines to execute:  pin7, flag7, fka7
    ! !subroutines should be ran in the order presented
    ! !
    ! !pin7:  Formats inputs (via a P and C array) run LagK operations
    ! !
    ! !         pin7(p,c,ita,itb,jlag,jk,meteng,lagtbl,ktbl,ico,iinfl,ioutfl,istor)
    ! !          p:  output from pin7 subroutine.  Contains lag/q, k/q, and 
    ! !              2*S/(DT/4)+O, O tables.  Also specifies the timestep of
    ! !              Input output
    ! !          c:   output from pin7 subroutine. Contains initial, inflow,
    ! !               utflow, storage, and carryover values
    ! !         ita:  Data time interval of the inflow time series in HR (integer)
    ! !         itb:  Data time interval of the outflow time series HR (integer)
    ! !         jlag: If > 0 - number of pairs of Lag and Q values used to define the
    ! !               variable Lag vs Q curve.  If = 0 - constant Lag will be used
    ! !         jk:   If > 0 - number of pairs of K and Q values used to define the
    ! !               variable K vs Q curve.  If = 0 - constant K will be used
    ! !         meteng:  Code specifying whether units of q, lag, k parameters and initial
    ! !                 values are English or metric:  'ENGL' = enter flow in CFS and 
    ! !                 volume in CFSD,'METR' = enter flow in CMS and volume in CMSD, 
    ! !                 Default is metric.  Note output:  P and C output is ALWAYS converted
    ! !                 to metric
    ! !         lagtbl:  If jlag=0, constant Lag value.  If jlag>0, lag and q pairs, in
    ! !                   that order, in a single column array.  Must be in ascending order
    ! !                   of q (example: [6, 0, 4, 10000, 3.5, 20000]).
    ! !         ktbl:  If jk=0, constant K value.  If jk>0, K and q pairs, in
    ! !                   that order, in a single column array.  Must be in ascending order
    ! !                   of q (example: [1, 100, 1, 40000, 3, 100000]).
    ! !         ico:  Initial carry over 
    ! !         iinfl:  Initial inflow
    ! !         ioutfl:  Initial outflow
    ! !         istor:  Initial storage
    ! !
    ! !  !!Notes: 
    ! !           1) This routine handles all the variables which would be optimized:
    ! !              lagtbl, ktble, ico, iinfl, ioutfl, istor
    ! !           2) %EDIT% This wrapper Fortran code uses parameters:  a, b c, d to develop
    ! !              lagtbl and ktbl.  The equation is lag/k_table_entry=a*(Q-d)**2+b*Q+c
    ! !              Q is the flow table entry
    ! !           3) The pin7.f subroutine was edited to start with a empty c/p array
    ! !              far larger than which should be needed [p(500),c(100)].  For anyone
    ! !              interestd below is python code I used to chop the unused lines after
    ! !              executing the subroutine. This is not necessary for flagk, and fka
    ! !              subroutines to properly run.  I used this document, pg 1-3, as reference:
    ! !              https://www.weather.gov/media/owp/oh/hrl/docs/833lagk.pdf
    ! !               k_start=int(p[17])
    ! !               k_len=int(p[k_start-1])
    ! !               pina7_len=int(p[k_start+2*k_len])
    ! !               p_end=k_start+2*(k_len+pina7_len)+1
    ! !               p=p[:p_end]
    ! !               c=c[:int(c[0])]
    ! !
    ! !flag7:  Controls the Lag Operation
    ! !             flag7(p,c,qa,qb,CO_ST,[ndt])
    ! !
    ! !            qb: downstream streamflow values (single column array) with only
    ! !                lag applied, time step is assumed to correspond to itb.
    ! !            p:  output from pin7 subroutine.  Contains lag/q, k/q, and
    ! !                2*S/(DT/4)+O, O tables.  Also specifies the timestep of
    ! !                Input output
    ! !            c:  output from pin7 subroutine. Contains initial, inflow, 
    ! !                outflow, storage, and carryover values.  !!NOTE!! use a copy of
    ! !                 the original c array as it get edited during the subroutine
    ! !            qa: Upstream streamflow values (single column array), time step is
    ! !                assumed to correspond to ita
    ! !            ntd:  Optional variable, total number to time steps to process.
    ! !                  if less than full qa array is desired 
    ! !           CO_ST:  Lagk state used for a warm start.  This is a timeseries of lag
    ! !                   time for flow input
    ! !
    ! !fka7:    Perform the attenuation (K) computations 
    ! !             flag7(p,c,qb,qc,STOR_ST,[ndt])
    ! !
    ! !            qc: downstream streamflow values (single column array) with both
    ! !                lag and attenuation applied, time step is assumed to correspond
    ! !                to itb.  
    ! !            p:  output from pin7 subroutine.  Contains lag/q, k/q, and
    ! !                2*S/(DT/4)+O, O tables.  Also specifies the timestep of
    ! !                Input output
    ! !            c:  output from pin7 subroutine. Contains initial, inflow,
    ! !                outflow, storage, and carryover values
    ! !            qb: downstream streamflow values (single column array) with only
    ! !                lag applied, time step is assumed to correspond to itb
    ! !            ntd: Optional variable, total number to time steps to process.
    ! !                  if less than full qa array is desired 
    ! !           STOR_ST:  Lagk state used for a warm start.  This is a timeseries of attenuation
    ! !                   storage

    ! !Wrapper varible:
    ! !            return_states:  Binary option to return co_st_out, inflow_st_out,storage_st_out.
    ! !                            1:  Return states, 0: Return only routed flows
    ! !            n_hrus:  The number of time steps of the simulation (integer)

    ! !            UNITS CONVERSION
    ! !             1 CFS to 0.0283168 CMS
    ! !             1 CMS to 35.3147 CFS
    ! !             1 CFD to 0.0283168 CMD

    ! !     WRC: RESTART CAPABILITY ARGUMENTS (Added January 2026)
    ! !     =================================================================
    ! !     INPUTS:
    ! !       c_array_in: (100 x n_hrus double array) WRC: Saved C array states 
    ! !                   from previous simulation run. Contains carryover values
    ! !                   for lag and attenuation operations. Only used if 
    ! !                   use_c_array_restart = 1.
    ! !
    ! !       use_c_array_restart: (integer flag) WRC: Controls restart behavior
    ! !                            1 = Warm start - use c_array_in to initialize C
    ! !                            0 = Cold start - use pin7 defaults (original behavior)
    ! !
    ! !     OUTPUTS:
    ! !       c_array_out: (100 x n_hrus double array) WRC: Final C array states
    ! !                    at end of simulation. Save these values to enable
    ! !                    warm start for subsequent simulations. Contains updated
    ! !                    carryover values after routing operations complete.
    ! !
    ! !     C ARRAY STRUCTURE (from pin7 subroutine):
    ! !       c(1)   = Total length of C array in use
    ! !       c(2)   = Initial carryover value
    ! !       c(3)   = Initial inflow value
    ! !       c(4)   = Initial outflow value
    ! !       c(5)   = Initial storage value
    ! !       c(6+)  = Additional carryover/state values (structure depends on 
    ! !                lag/K table configuration and timestep settings)
    ! !     =================================================================


  implicit none

  ! ! inputs
  integer, intent(in):: n_hrus, ita, itb, sim_length
  integer, intent(in):: use_c_array_restart   ! WRC: NEW - flag to use restart C array (1=warm start, 0=cold start)
  double precision, dimension(100, n_hrus), intent(in):: c_array_in ! WRC: NEW - saved C array states from previous run
  character(len = 4), parameter:: meteng = 'METR'
  logical:: return_states
  double precision, dimension(n_hrus), intent(in):: ico_in, iinfl_in, ioutfl_in, istor_in
  double precision, dimension(n_hrus), intent(in):: lagtbl_a_in, lagtbl_b_in, lagtbl_c_in, lagtbl_d_in  
  double precision, dimension(n_hrus), intent(in):: ktbl_a_in, ktbl_b_in, ktbl_c_in, ktbl_d_in
  double precision, dimension(n_hrus), intent(in):: lagk_kmax_in, lagk_lagmax_in, lagk_qmax_in
  double precision, dimension(n_hrus), intent(in):: lagk_kmin_in, lagk_lagmin_in, lagk_qmin_in
  double precision, dimension(sim_length, n_hrus), intent(in):: qa_in
  
  ! ! local variables
  real, dimension(n_hrus):: ico, iinfl, ioutfl, istor
  real, dimension(n_hrus):: lagtbl_a, lagtbl_b, lagtbl_c, lagtbl_d    
  real, dimension(n_hrus):: ktbl_a, ktbl_b, ktbl_c, ktbl_d
  real, dimension(n_hrus):: lagk_kmax, lagk_lagmax, lagk_qmax
  real, dimension(n_hrus):: lagk_kmin, lagk_lagmin, lagk_qmin
  real, dimension(22, n_hrus):: lagtbl, ktbl
  real, dimension(sim_length, n_hrus):: qa 
  real, dimension(500,n_hrus):: p
  real, dimension(100,n_hrus):: c
  real, dimension(100):: c_cpy
  real, dimension(sim_length, n_hrus):: qb, qc
  real, dimension(sim_length, n_hrus):: storage_st, co_st
  integer, dimension(n_hrus):: jlag, jk
  integer:: nh, i
  real:: ndq, lag_entry, k_entry
  
  ! ! output 
  double precision, dimension(sim_length, n_hrus), intent(out):: lagk_out
  double precision, dimension(sim_length, n_hrus), intent(out):: inflow_st_out, storage_st_out, co_st_out
  double precision, dimension(100, n_hrus), intent(out):: c_array_out ! WRC: NEW - final C array states for next restart

  ! Convert double precision to single precision
  ico=real(ico_in)*0.0283168e0
  iinfl=real(iinfl_in)*0.0283168e0
  ioutfl=real(ioutfl_in)*0.0283168e0
  istor=real(istor_in)*0.0283168e0
  
  lagtbl_a=real(lagtbl_a_in)
  lagtbl_b=real(lagtbl_b_in)
  lagtbl_c=real(lagtbl_c_in)
  lagtbl_d=real(lagtbl_d_in) 
  ktbl_a=real(ktbl_a_in) 
  ktbl_b=real(ktbl_b_in) 
  ktbl_c=real(ktbl_c_in) 
  ktbl_d=real(ktbl_d_in)
  lagk_kmax=real(lagk_kmax_in)
  lagk_lagmax=real(lagk_lagmax_in)
  lagk_qmax=real(lagk_qmax_in)*0.0283168e0
  lagk_kmin=real(lagk_kmin_in)
  lagk_lagmin=real(lagk_lagmin_in) 
  lagk_qmin=real(lagk_qmin_in)*0.0283168e0
  
  qa=real(qa_in)*0.0283168e0

  lagk_out = 0
  lagtbl = 0 
  ktbl = 0
  lag_entry = 0
  k_entry = 0
  p = 0
  c = 0
  c_cpy = 0
  qb = 0
  qc = 0
  if(return_states)then
    storage_st = 0
    co_st = 0
  end if
  
  ! Populate Lag and K tables  
  
  do nh=1,n_hrus  
   ndq=0
   do i=1,11
    lagtbl(i*2,nh)=ndq*(lagk_qmax(nh)-lagk_qmin(nh))+lagk_qmin(nh)
    ktbl(i*2,nh)=ndq*(lagk_qmax(nh)-lagk_qmin(nh))+lagk_qmin(nh)
   
    lag_entry=lagtbl_a(nh)*(ndq-lagtbl_d(nh))**2+lagtbl_b(nh)*ndq+lagtbl_c(nh)
    k_entry=ktbl_a(nh)*(ndq-ktbl_d(nh))**2+ktbl_b(nh)*ndq+ktbl_c(nh)
   
    if (lag_entry > 0 .AND. lag_entry < 1) then
     lagtbl(i*2-1,nh)=lag_entry*(lagk_lagmax(nh)-lagk_lagmin(nh))+lagk_lagmin(nh)
    else if (lag_entry >= 1) then
     lagtbl(i*2-1,nh)=lagk_lagmax(nh)
    else
     lagtbl(i*2-1,nh)=lagk_lagmin(nh)
    end if
    
    if (k_entry > 0 .AND. k_entry < 1) then
     ktbl(i*2-1,nh)=k_entry*(lagk_kmax(nh)-lagk_kmin(nh))+lagk_kmin(nh)
    else if (k_entry >= 1) then
     ktbl(i*2-1,nh)=lagk_kmax(nh)
    else
     ktbl(i*2-1,nh)=lagk_kmin(nh)
    end if
   
    ndq=ndq+.1
   end do
  end do

  ! Get length of K and Lag Table
  do nh=1,n_hrus  
    if (MAXVAL(lagtbl(::2,nh))==MINVAL(lagtbl(::2,nh))) then
      jlag(nh)=0
    else
      jlag(nh)=size(lagtbl,1)/2
    end if
    if (MAXVAL(ktbl(::2,nh))==MINVAL(ktbl(::2,nh))) then
      jk(nh)=0
    else
      jk(nh)=size(ktbl,1)/2
    end if
  end do

  ! Loop through each reach and calculate lag 
  do nh=1,n_hrus
    
    ! WRC: Always call pin7 to initialize P array and default C array
    ! WRC: pin7 sets up routing tables and initial conditions from parameters
    call pin7(p(:,nh),c(:,nh),int(ita,4),int(itb,4),jlag(nh),jk(nh),meteng,lagtbl(:,nh), &
       ktbl(:,nh),ico(nh),iinfl(nh),ioutfl(nh),istor(nh))

    ! WRC: ===== CONDITIONAL C ARRAY RESTART INITIALIZATION (BEGIN) =====
    ! WRC: If restarting, override the C array with saved state from previous run
    ! WRC: This must happen AFTER pin7 (to ensure P array is set) but BEFORE routing
    if (use_c_array_restart .eq. 1) then
      c(:,nh) = real(c_array_in(:,nh))
    end if
    ! WRC: ===== CONDITIONAL C ARRAY RESTART INITIALIZATION (END) =====
 
    c_cpy=c(:,nh) ! WRC: Create working copy - flag7 and fka7 modify C in place

    call flag7(p(:,nh),c_cpy,qa(:,nh),qb(:,nh),int(sim_length,4), &
       co_st(:,nh))
    
    call fka7(p(:,nh),c_cpy,qb(:,nh),qc(:,nh),int(sim_length,4), &
       storage_st(:,nh))
       
    ! WRC: ===== SAVE FINAL C ARRAY STATE FOR RESTART (BEGIN) =====
    ! WRC: Capture the final C array state for this tributary/HRU
    ! WRC: c_cpy has been modified by flag7 and fka7 to contain updated carryover values
    ! WRC: Convert from single to double precision for output
    
    c_array_out(:,nh) = dble(c_cpy)
    
    ! WRC: ===== SAVE FINAL C ARRAY STATE FOR RESTART (END) =====

  end do
  
  lagk_out=dble(qc)*35.3147d0
  if(return_states)then
    inflow_st_out=dble(qb)*35.3147d0
    storage_st_out=dble(storage_st)*35.3147d0
    co_st_out=dble(co_st)
  end if
  
end subroutine