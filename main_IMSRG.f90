program main_IMSRG
  use basic_IMSRG
  use HF_mod
  use generators
  use commutators 
  ! ground state IMSRG calculation for nuclear system 
  implicit none
  
  type(spd) :: jbasis 
  type(sq_op) :: HS,ETA,DH
  character(200) :: sp_input_file,interaction_file
  character(200) :: inputs_from_command
  integer :: i,j,T,P,JT,a,b,c,d,ham_type,j3
  real(8) :: hw ,sm,omp_get_wtime,t1,t2
  logical :: hartree_fock 

!============================================================
! READ INPUTS SET UP STORAGE STRUCTURE
!============================================================
  call getarg(1,inputs_from_command) 
  call read_main_input_file(inputs_from_command,HS,ham_type,&
       hartree_fock,hw,sp_input_file,interaction_file)
 
  HS%herm = 1

  call read_sp_basis(jbasis,sp_input_file,HS%Aprot,HS%Aneut) 
  call allocate_blocks(jbasis,HS)   
  call read_interaction(HS,interaction_file,jbasis,ham_type,hw)
 
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
! PLAYGROUND
!============================================================ 

  print*, HS%E0

  ! for testing commutators
  do i = 1,HS%Nsp-HS%belowEF
     do j = 1,HS%belowEF
        
        call random_number(HS%fph(i,j))
       
     end do 
  end do 
  
  call duplicate_sq_op(HS,ETA) 
  call duplicate_sq_op(HS,DH) 
  call build_gs_white(HS,ETA,jbasis) 

  t1 = omp_get_wtime()
  call commutator_111(ETA,HS,DH,jbasis)
  t2 = omp_get_wtime()
  print*, t2-t1
  call print_matrix(DH%fph(1:6,:))
  call print_matrix(DH%fhh)
  call print_matrix(DH%fpp(1:10,1:10))
  
  
end program main_IMSRG

