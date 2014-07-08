program main_IMSRG
  use basic_IMSRG
  ! ground state IMSRG calculation for nuclear system 
  implicit none
  
  type(spd) :: jbasis 
  type(sq_op) :: HS 
  character(50) :: sp_input_file,interaction_file
  integer :: i
  
  HS%nbody = 2
  HS%herm = 1  
  sp_input_file ='nl4.sps'
  interaction_file = 'vsrg.int' 
  
  call read_sp_basis(jbasis,sp_input_file,HS%nbody) 
   
  call allocate_blocks(jbasis,HS) 
  
  call read_interaction(HS,interaction_file,jbasis) 
  
end program
