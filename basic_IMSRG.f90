module basic_IMSRG 
  ! contains basic type declarations, constants, and adminstrative routines
  
  TYPE :: single_particle_descript
     INTEGER :: total_orbits
     INTEGER, ALLOCATABLE,DIMENSION(:) :: nn, ll, jj, itzp, nshell, mvalue
     ! for clarity:  nn, ll, nshell are all the true value
     ! jj is j+1/2 (so it's an integer) 
     ! likewise itzp is 2*tz 
     CHARACTER (LEN=10), ALLOCATABLE,DIMENSION(:) :: orbit_status, model_space,basis
     REAL(8), ALLOCATABLE,DIMENSION(:) :: e, evalence, e_original
  END TYPE single_particle_descript

contains 
!====================================================
!====================================================
subroutine read_sp_basis(jbas,sp_input_file)
  ! fills jscheme_basis array with single particle data from file
  ! file format must be: 5 integers and one real 
  ! for each state we have:   | label |  n  | l | 2 * j | 2*tz | E_sp |  
  implicit none 
  
  type(single_particle_descript) :: jbas
  character(50) :: sp_input_file
  integer :: ist,i,label,ni,li,ji,tzi,ix
  real(8) :: e
  
  open(unit=39,file='../sp_inputs/'//trim(adjustl(sp_input_file)))
  
  ix = 0 
  
  ! count the number of states in the file. 
  do 
  
     read(39,*,iostat=ist) label,ni,li,ji,tzi,e
     
     if (ist>0) stop 'input error in .sps file. See &
          subroutine read_sp_basis in basic_IMSRG.f90' 
     if (ist<0) exit
  
     ix = ix + 1 
  
  end do
  
  ! build the jscheme descriptor
  jbas%total_orbits=ix
  allocate(jbas%nn(ix)) 
  allocate(jbas%ll(ix))
  allocate(jbas%jj(ix))
  allocate(jbas%itzp(ix))
  allocate(jbas%nshell(ix))
  allocate(jbas%e(ix)) 
  
  ! go back to the start and read them in. 
  rewind(39) 
  
  do i = 1,ix
     
     read(39,*) label,ni,li,ji,tzi,e
     
     jbas%nn(i) = ni 
     jbas%ll(i) = li
     jbas%jj(i) = ji
     jbas%itzp(i) = tzi 
     jbas%nshell(i) = 2*ni + li 
     jbas%e(i) =  e 

  end do 
  
  close(39) 
  
end subroutine 
!==============================================
!==============================================  
end module
       
  
