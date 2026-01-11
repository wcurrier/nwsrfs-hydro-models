subroutine lagk(n_hrus, ita, itb, &
    lagtbl_a_in, lagtbl_b_in, lagtbl_c_in, lagtbl_d_in,&
    ktbl_a_in, ktbl_b_in, ktbl_c_in, ktbl_d_in, &
    lagk_lagmax_in, lagk_kmax_in, lagk_qmax_in, &
    lagk_lagmin_in, lagk_kmin_in, lagk_qmin_in, &
    ico_in, iinfl_in, ioutfl_in, istor_in, &
    c_array_in, use_c_array_restart, &
    qa_in, sim_length, &
    return_states, &
    lagk_out, co_st_out, &
    inflow_st_out, storage_st_out, &
    c_array_out)

  implicit none

  ! ! inputs
  integer, intent(in):: n_hrus, ita, itb, sim_length
  integer, intent(in):: use_c_array_restart  ! ADD THIS
  double precision, dimension(100, n_hrus), intent(in):: c_array_in
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
  double precision, dimension(100, n_hrus), intent(out):: c_array_out

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
    
    ! Always call pin7 to initialize P array
    call pin7(p(:,nh),c(:,nh),int(ita,4),int(itb,4),jlag(nh),jk(nh),meteng,lagtbl(:,nh), &
       ktbl(:,nh),ico(nh),iinfl(nh),ioutfl(nh),istor(nh))
    
    ! If restarting, override the C array with saved state
    if (use_c_array_restart .eq. 1) then
      c(:,nh) = real(c_array_in(:,nh))
    end if
    
    c_cpy=c(:,nh)

    call flag7(p(:,nh),c_cpy,qa(:,nh),qb(:,nh),int(sim_length,4), &
       co_st(:,nh))
    
    call fka7(p(:,nh),c_cpy,qb(:,nh),qc(:,nh),int(sim_length,4), &
       storage_st(:,nh))
    
    ! Save the final C array state for this tributary
    c_array_out(:,nh) = dble(c_cpy)
    
  end do
  
  lagk_out=dble(qc)*35.3147d0
  if(return_states)then
    inflow_st_out=dble(qb)*35.3147d0
    storage_st_out=dble(storage_st)*35.3147d0
    co_st_out=dble(co_st)
  end if
  
end subroutine