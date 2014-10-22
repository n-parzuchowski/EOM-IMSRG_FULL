program main_IMSRG
  use basic_IMSRG
  use HF_mod
  use IMSRG_ODE
  use IMSRG_MAGNUS
  ! ground state IMSRG calculation for nuclear system 
  implicit none
  
  type(spd) :: jbasis
  type(sq_op) :: HS,ETA,DH,w1,w2
  type(cross_coupled_31_mat) :: CCHS,CCETA,WCC
  character(200) :: inputs_from_command
  integer :: i,j,T,P,JT,a,b,c,d,g,q,ham_type,j3
  integer :: np,nh,nb,k,l,m,n
  real(8) :: hw ,sm,omp_get_wtime,t1,t2,bet_off,d6ji,gx
  logical :: hartree_fock,magnus_exp 
  external :: dHds_white_gs,dHds_TDA_shell

!============================================================
! READ INPUTS SET UP STORAGE STRUCTURE
!============================================================
  call getarg(1,inputs_from_command) 
  call read_main_input_file(inputs_from_command,HS,ham_type,&
       hartree_fock,magnus_exp,hw)
 
  HS%herm = 1

  call read_sp_basis(jbasis,HS%Aprot,HS%Aneut) 
  call allocate_blocks(jbasis,HS)   
  call read_interaction(HS,jbasis,ham_type,hw)
  
!============================================================
! BUILD BASIS
!============================================================
  call calculate_h0_harm_osc(hw,jbasis,HS,ham_type) 
  
  if (hartree_fock) then 
     call calc_HF(HS,jbasis) 
    ! calc_HF normal orders the hamiltonian
  else 
     call normal_order(HS,jbasis) 
  end if
    
!============================================================
! IM-SRG CALCULATION 
!============================================================ 

  if (magnus_exp) then 
     call magnus_decouple(HS,jbasis)
  else
     call decouple_hamiltonian(HS,jbasis,dHds_white_gs) 
  end if 

!  call TDA_decouple(HS,jbasis,dHds_TDA_shell) 
  
end program main_IMSRG




