program main_IMSRG
  use basic_IMSRG
  use HF_mod
  ! ground state IMSRG calculation for nuclear system 
  implicit none
  
  type(spd) :: jbasis 
  type(sq_op) :: HS 
  character(50) :: sp_input_file,interaction_file
  integer :: i,j,T,P,JT,a,b,c,d,ham_type,j3
  real(8) :: hw ,sm
  
  HS%Aprot = 20
  HS%Aneut = 20
  HS%herm = 1  
  ham_type = 1
  hw = 28.0
  sp_input_file ='nl4.sps'
  interaction_file = 'vsrg.int' 
  
  call read_sp_basis(jbasis,sp_input_file,HS%Aprot,HS%Aneut) 

  call allocate_blocks(jbasis,HS)
   
  call read_interaction(HS,interaction_file,jbasis,ham_type,hw)
 
  call calculate_h0_harm_osc(hw,jbasis,HS,ham_type) 
 
  call calc_HF(HS,jbasis) 
  
end program main_IMSRG

