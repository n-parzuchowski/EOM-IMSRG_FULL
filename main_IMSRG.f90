program main_IMSRG
  use basic_IMSRG
  use HF_mod
  ! ground state IMSRG calculation for nuclear system 
  implicit none
  
  type(spd) :: jbasis 
  type(sq_op) :: HS 
  character(200) :: sp_input_file,interaction_file
  character(200) :: inputs_from_command
  integer :: i,j,T,P,JT,a,b,c,d,ham_type,j3
  real(8) :: hw ,sm
  logical :: hartree_fock 
  
  call getarg(1,inputs_from_command) 
  call read_main_input_file(inputs_from_command,HS,ham_type,&
       hartree_fock,hw,sp_input_file,interaction_file)
 
  HS%herm = 1
  
  call read_sp_basis(jbasis,sp_input_file,HS%Aprot,HS%Aneut) 

  call allocate_blocks(jbasis,HS)
   
  call read_interaction(HS,interaction_file,jbasis,ham_type,hw)
 
  call calculate_h0_harm_osc(hw,jbasis,HS,ham_type) 
 
  IF (hartree_fock) then 
     call calc_HF(HS,jbasis) 
     ! calc_HF normal orders the hamiltonian
  else 
     call normal_order(HS,jbasis) 
  END IF
  
  print*, HS%E0
     
     
end program main_IMSRG

