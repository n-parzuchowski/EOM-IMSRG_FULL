program main_IMSRG
  use basic_IMSRG
  use HF_mod
  use IMSRG_ODE
  use IMSRG_MAGNUS
  ! ground state IMSRG calculation for nuclear system 
  implicit none
  
  type(spd) :: jbasis
  type(sq_op) :: HS,ETA,DH,w1,w2,Hcm
  type(cross_coupled_31_mat) :: CCHS,CCETA,WCC
  character(200) :: inputs_from_command
  integer :: i,j,T,P,JT,a,b,c,d,g,q,ham_type,j3
  integer :: np,nh,nb,k,l,m,n
  real(8) :: hw ,sm,omp_get_wtime,t1,t2,bet_off,d6ji,gx
  logical :: hartree_fock,magnus_exp,tda_calculation,COM_calc
  external :: dHds_white_gs,dHds_TDA_shell

!============================================================
! READ INPUTS SET UP STORAGE STRUCTURE
!============================================================
  call getarg(1,inputs_from_command) 
  call read_main_input_file(inputs_from_command,HS,ham_type,&
       hartree_fock,magnus_exp,tda_calculation,COM_calc,hw)
  
  HS%herm = 1
  HS%hospace = hw

  call read_sp_basis(jbasis,HS%Aprot,HS%Aneut) 
  call allocate_blocks(jbasis,HS)   
  
  ! for calculating COM expectation value
  if (COM_calc) then  
     call duplicate_sq_op(HS,Hcm)   
     call read_interaction(HS,jbasis,ham_type,hw,Hcm,hw)
     call calculate_h0_harm_osc(hw,jbasis,Hcm,2) 
  else 
     call read_interaction(HS,jbasis,ham_type,hw)
  end if 
  
!============================================================
! BUILD BASIS
!============================================================
  call calculate_h0_harm_osc(hw,jbasis,HS,ham_type) 

  if (hartree_fock) then 
     if (COM_calc) then 
        call calc_HF(HS,jbasis,Hcm)
        call normal_order(Hcm,jbasis)
        Hcm%E0 = Hcm%E0 - 1.5d0*hw
     else
        call calc_HF(HS,jbasis)
     end if 
    ! calc_HF normal orders the hamiltonian
  else 
     call normal_order(HS,jbasis) 
  end if

!============================================================
! IM-SRG CALCULATION 
!============================================================ 
  if (magnus_exp) then 
    
     if (COM_calc) then 
        call magnus_decouple(HS,jbasis,Hcm) 
     else 
        call magnus_decouple(HS,jbasis) 
     end if 
     
  else
     
     if (COM_calc) then 
        call decouple_hamiltonian(HS,jbasis,dHds_white_gs) 
     else 
        call decouple_hamiltonian(HS,jbasis,dHds_white_gs) 
     end if
     
  end if
  
  if (tda_calculation) then 
     call TDA_decouple(HS,jbasis,dHds_TDA_shell) 
  end if 
  
end program main_IMSRG




