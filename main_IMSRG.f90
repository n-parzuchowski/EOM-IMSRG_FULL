program main_IMSRG
  use basic_IMSRG
  use HF_mod
  use generators
  use commutators 
  ! ground state IMSRG calculation for nuclear system 
  implicit none
  
  type(spd) :: jbasis
  type(sq_op) :: HS,ETA,DH,w1,w2
  type(cross_coupled_31_mat) :: CCHS,CCETA,WCC
  character(200) :: sp_input_file,interaction_file
  character(200) :: inputs_from_command
  integer :: i,j,T,P,JT,a,b,c,d,g,q,ham_type,j3
  integer :: np,nh,nb,k,l,m,n
  real(8) :: hw ,sm,omp_get_wtime,t1,t2,bet_off,d6ji
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

!  print*, HS%E0

     np = HS%mat(1)%npp
     nh = HS%mat(1)%nhh
     nb = HS%mat(1)%nph
  ! for testing commutators
        
  call duplicate_sq_op(HS,ETA) !generator
  call duplicate_sq_op(HS,DH) !derivative
  call duplicate_sq_op(HS,w1) !workspace
  call duplicate_sq_op(HS,w2) !workspace
  call allocate_CCMAT(HS,CCHS,jbasis) !cross coupled ME
  call duplicate_CCMAT(CCHS,CCETA) !cross coupled ME
  call allocate_CC_wkspc(CCHS,WCC)
  call build_gs_white(HS,ETA,jbasis) 
 
!  t1 = omp_get_wtime()
  call calculate_cross_coupled(HS,CCHS,jbasis,.true.)
  call calculate_cross_coupled(ETA,CCETA,jbasis,.false.) 
!  t2 = omp_get_wtime()
!  print*, t2-t1

  DH%E0 = commutator_110(ETA,HS,jbasis) + commutator_220(ETA,HS,jbasis)
!  t2 = omp_get_wtime()
!  print*, t2-t1
  call commutator_111(ETA,HS,DH,jbasis) 
!  t2 = omp_get_wtime()
!  print*, t2-t1
  call commutator_121(ETA,HS,DH,jbasis)
!  t2 = omp_get_wtime()
!  print*, t2-t1
  call commutator_122(ETA,HS,DH,jbasis)
!  t2 = omp_get_wtime()
!  print*, t2-t1
  call commutator_222_pp_hh(ETA,HS,DH,w1,w2,jbasis)
!  t2 = omp_get_wtime()
!  print*, t2-t1
  call commutator_221(ETA,HS,DH,w1,w2,jbasis)
!  t2 = omp_get_wtime()
!  print*, t2-t1
  call commutator_222_ph(CCETA,CCHS,DH,WCC,jbasis)
!  t2 = omp_get_wtime()
!  print*, t2-t1

 ! call print_matrix(DH%mat(1)%gam(2)%X)
end program main_IMSRG

