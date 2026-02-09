! !     ==============================================================================
! !     WRC MODIFICATION NOTES - RESTART CAPABILITY ADDITIONS (January 2026)
! !     ==============================================================================
! !     The following modifications enable Lag-K routing restart capability
! !     for operational forecasting workflows. Changes allow:
! !     - Warm start from previous C array states (carryover values)
! !     - Continuation of routing operations across simulation boundaries
! !     - Support for ensemble forecasting from same initial conditions
! !     - NEW: Optional C-array time series output for exact snapshot states
! !     
! !     Key additions:
! !     1. New input arguments: c_array_in (100x n_hrus), use_c_array_restart (flag)
! !     2. New output argument: c_array_out (100x n_hrus) for state persistence
! !     3. NEW: c_array_ts (100, sim_length, n_hrus) for snapshot capability
! !        Only populated when save_states=1 to minimize memory overhead
! !     4. Conditional C array initialization logic
! !     5. Final state capture for each HRU/tributary
! !     
! !     FIX (Feb 2026): flag7 now takes QT as a passed-in argument to avoid
! !     ALLOCATABLE ABI conflicts between system gfortran (which compiles
! !     flag7.f) and conda gfortran (which compiles this file and links).
! !     QT is allocated HERE (by conda gfortran) and passed to flag7.
! !     ==============================================================================

subroutine lagk(n_hrus, ita, itb, &
    lagtbl_a_in, lagtbl_b_in, lagtbl_c_in, lagtbl_d_in,&
    ktbl_a_in, ktbl_b_in, ktbl_c_in, ktbl_d_in, &
    lagk_lagmax_in, lagk_kmax_in, lagk_qmax_in, &
    lagk_lagmin_in, lagk_kmin_in, lagk_qmin_in, &
    ico_in, iinfl_in, ioutfl_in, istor_in, &
    c_array_in, use_c_array_restart, &
    qa_in, sim_length, &
    return_states, save_states, &
    lagk_out, co_st_out, &
    inflow_st_out, storage_st_out, &
    c_array_out, c_array_ts)

  implicit none

  ! inputs
  integer, intent(in):: n_hrus, ita, itb, sim_length
  integer, intent(in):: use_c_array_restart
  integer, intent(in):: save_states  ! flag to save c_array at every timestep
  double precision, dimension(100, n_hrus), intent(in):: c_array_in
  character(len = 4), parameter:: meteng = 'METR'
  logical:: return_states
  double precision, dimension(n_hrus), intent(in):: ico_in, iinfl_in, ioutfl_in, istor_in
  double precision, dimension(n_hrus), intent(in):: lagtbl_a_in, lagtbl_b_in, lagtbl_c_in, lagtbl_d_in  
  double precision, dimension(n_hrus), intent(in):: ktbl_a_in, ktbl_b_in, ktbl_c_in, ktbl_d_in
  double precision, dimension(n_hrus), intent(in):: lagk_kmax_in, lagk_lagmax_in, lagk_qmax_in
  double precision, dimension(n_hrus), intent(in):: lagk_kmin_in, lagk_lagmin_in, lagk_qmin_in
  double precision, dimension(sim_length, n_hrus), intent(in):: qa_in
  
  ! local variables
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
  integer:: nh, i, nqtsz
  real:: ndq, lag_entry, k_entry

  ! Local single-element arrays for single-timestep flag7/fka7 calls
  real :: qa_1(1), qb_1(1), co_1(1), qc_1(1), st_1(1)

  ! QT work array for flag7, allocated here to avoid mixed-compiler
  ! ALLOCATABLE ABI issues (this file = conda gfortran, flag7.f = system gfortran).
  ! Size: max(sim_length*3, 500) covers both full-timeseries and NDT=1 modes.
  real, allocatable :: qt_work(:)
  
  ! output 
  double precision, dimension(sim_length, n_hrus), intent(out):: lagk_out
  double precision, dimension(sim_length, n_hrus), intent(out):: inflow_st_out, storage_st_out, co_st_out
  double precision, dimension(100, n_hrus), intent(out):: c_array_out
  ! c_array time series for snapshot capability
  double precision, dimension(100, sim_length, n_hrus), intent(out):: c_array_ts

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
  c_array_out = 0
  c_array_ts = 0
  
  if(return_states)then
    storage_st = 0
    co_st = 0
  end if

  ! Allocate QT work array for flag7.
  ! Must be large enough for the full time series case (sim_length*3)
  ! AND for NDT=1 carryover pairs (need at least 500).
  nqtsz = sim_length * 3
  if (nqtsz < 500) nqtsz = 500
  allocate(qt_work(nqtsz))
  
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
    
    ! Always call pin7 to initialize P array and default C array
    call pin7(p(:,nh),c(:,nh),int(ita,4),int(itb,4),jlag(nh),jk(nh),meteng,lagtbl(:,nh), &
       ktbl(:,nh),ico(nh),iinfl(nh),ioutfl(nh),istor(nh))

    ! If restarting, override the C array with saved state from previous run
    if (use_c_array_restart .eq. 1) then
      c(:,nh) = real(c_array_in(:,nh))
    end if
 
    c_cpy=c(:,nh)

    ! ==================================================================
    ! When save_states=1, process one timestep at a time so we can
    ! capture the exact C array state after each (flag7 + fka7) pair.
    !
    ! flag7 now takes QT as an argument (qt_work, nqtsz) instead of
    ! allocating it internally. This avoids ABI conflicts between the
    ! system gfortran (which compiles flag7.f) and conda gfortran
    ! (which compiles this file). QT is allocated here once and reused
    ! for all 62,820 calls per reach.
    ! ==================================================================
    
    if (save_states == 1) then
      do i = 1, sim_length
        qa_1(1) = qa(i, nh)
        qb_1(1) = 0.0
        co_1(1) = 0.0
        qc_1(1) = 0.0
        st_1(1) = 0.0
        
        ! Lag operation (single timestep) - pass qt_work as work array
        call flag7(p(:,nh), c_cpy, qa_1, qb_1, 1, co_1, qt_work, nqtsz)
        
        ! K attenuation (single timestep)
        call fka7(p(:,nh), c_cpy, qb_1, qc_1, 1, st_1)
        
        ! Copy results to output arrays
        qb(i,nh) = qb_1(1)
        qc(i,nh) = qc_1(1)
        co_st(i,nh) = co_1(1)
        storage_st(i,nh) = st_1(1)
        
        ! Capture C array state at this timestep
        c_array_ts(:,i,nh) = dble(c_cpy)
      end do
    else
      ! Original fast approach - call on entire time series at once
      ! Pass qt_work as work array to flag7
      call flag7(p(:,nh), c_cpy, qa(:,nh), qb(:,nh), sim_length, &
                 co_st(:,nh), qt_work, nqtsz)
      call fka7(p(:,nh), c_cpy, qb(:,nh), qc(:,nh), sim_length, storage_st(:,nh))
    end if
       
    ! Save final C array state for restart
    c_array_out(:,nh) = dble(c_cpy)

  end do

  deallocate(qt_work)
  
  lagk_out=dble(qc)*35.3147d0
  if(return_states)then
    inflow_st_out=dble(qb)*35.3147d0
    storage_st_out=dble(storage_st)*35.3147d0
    co_st_out=dble(co_st)
  end if
  
end subroutine