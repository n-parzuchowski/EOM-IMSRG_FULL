program add_sp_energies
  !for HO basis
  !calculates sp energies for a given hw and vls 
  implicit none
  
  character(50) :: filename
  real(8) :: hw,vls,e
  integer :: ist,lbl,n,l,j,tz
  
  print*, 'name of file in sp directory: ' 
  read*, filename
  
  print*, 'hbar-omega: '
  read*, hw
  
  print*, 'spin_orbit potential: '
  read*, vls
  
  open(unit=39,file = '../../sp_inputs/'//trim(adjustl(filename))) 
  
  open(unit=40,file = '../../sp_inputs/temporary.sps') 
  
  ! read in each state, calculate it's energy, write it to a temp file
  do 
     read(39,*,iostat=ist) lbl,n,l,j,tz,e 
     if (ist > 0) stop 'file not formatted correctly'
     if (ist < 0) exit
     
     e = hw *(2*n+l+1.5) + &
          0.5*vls*(0.25 * j * (j+1) - l * (l+1) -0.75)  
     
     write(40,'(5(I5),e17.7)')  lbl,n,l,j,tz,e
     
  end do 
  
  close(39) 
  close(40) 
  
  ! replace real file with temp file
  call system('cp ../../sp_inputs/temporary.sps '&
       //'../../sp_inputs/'//trim(adjustl(filename))) 
  
  call system('rm ../../sp_inputs/temporary.sps') 
  
end program
     
     
