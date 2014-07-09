program main_IMSRG
  use basic_IMSRG
  ! ground state IMSRG calculation for nuclear system 
  implicit none
  
  type(spd) :: jbasis 
  type(sq_op) :: HS 
  character(50) :: sp_input_file,interaction_file
  integer :: i,T,P,J,a,b,c,d
  real(8) :: hw 
  
  HS%Abody = 2
  HS%herm = 1  
  hw = 28.0
  sp_input_file ='nl4.sps'
  interaction_file = 'vsrg.int' 
  
  call read_sp_basis(jbasis,sp_input_file,HS%Abody) 

  call allocate_blocks(jbasis,HS)
   
  call read_interaction(HS,interaction_file,jbasis)
 
  call calculate_h0_harm_osc(hw,jbasis,HS,1) 
  
  print*, v_elem(4,14,14,6,6,1,0,HS,jbasis) 
 
end program

