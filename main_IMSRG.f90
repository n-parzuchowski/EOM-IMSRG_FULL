program main_IMSRG
  use basic_IMSRG
  use HF_mod
  use generators
  use commutators 
  ! ground state IMSRG calculation for nuclear system 
  implicit none
  
  type(spd) :: jbasis 
  type(sq_op) :: HS,ETA,DH,w1,w2
  character(200) :: sp_input_file,interaction_file
  character(200) :: inputs_from_command
  integer :: i,j,T,P,JT,a,b,c,d,g,q,ham_type,j3
  integer :: np,nh,nb
  real(8) :: hw ,sm,omp_get_wtime,t1,t2,bet_off
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

  np = HS%mat(1)%npp
     nh = HS%mat(1)%nhh
     nb = HS%mat(1)%nph
  ! for testing commutators
 
  call duplicate_sq_op(HS,ETA) !generator
  call duplicate_sq_op(HS,DH) !derivative
  call duplicate_sq_op(HS,w1) !workspace
  call duplicate_sq_op(HS,w2) !workspace
  call build_gs_white(HS,ETA,jbasis) 

 
   do i = 1,np
     do j = 1,nb
        
        call random_number(ETA%mat(1)%gam(2)%X(i,j))
       
     end do 
  end do 
   
  DH%fph = 0.
  DH%fpp = 0.
  DH%fhh = 0. 

  t1 = omp_get_wtime()
  call commutator_222_pp_hh(ETA,HS,DH,w1,w2,jbasis)
  call commutator_221(ETA,HS,DH,w1,w2,jbasis) 
  t2 = omp_get_wtime()
  print*, t2-t1

  call print_matrix(DH%mat(1)%gam(4)%X(1:4,1:4))
  
  call print_matrix(DH%fhh)
  call print_matrix(DH%fpp(1:10,1:10))
  call print_matrix(DH%fph(1:6,1:6)) 
  
  do q = 1,DH%nblocks
     do g = 1,6
        DH%mat(q)%gam(g)%X = 0.d0 
     end do 
  end do 
  
  

end program main_IMSRG

