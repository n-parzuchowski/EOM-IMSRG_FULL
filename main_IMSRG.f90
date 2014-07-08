program main_IMSRG
  use basic_IMSRG
  ! ground state IMSRG calculation for nuclear system 
  implicit none
  
  type(single_particle_descript) :: jbasis 
  character(50) :: sp_input_file 
  integer :: i 
  
  sp_input_file ='nl4.sps'
  call read_sp_basis(jbasis,sp_input_file) 
  
  do i = 1, jbasis%total_orbits
     print*, jbasis%nn(i),jbasis%ll(i),jbasis%jj(i),jbasis%itzp(i),jbasis%e(i)
  end do 
end program
