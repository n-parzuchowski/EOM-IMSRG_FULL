program main_IMSRG
  use basic_IMSRG
  use HF_mod
  use IMSRG_ODE
  use IMSRG_MAGNUS
  use IMSRG_CANONICAL
  use operators
  ! ground state IMSRG calculation for nuclear system 
  implicit none
  
  type(spd) :: jbasis
  type(sq_op) :: HS,ETA,DH,w1,w2,Hcm,rirj,pipj,r2_rms
  type(cross_coupled_31_mat) :: CCHS,CCETA,WCC
  type(full_sp_block_mat) :: coefs,TDA,ppTDA,rrTDA
  character(200) :: inputs_from_command
  integer :: i,j,T,P,JT,a,b,c,d,g,q,ham_type,j3
  integer :: np,nh,nb,k,l,m,n,method_int
  real(8) :: hw ,sm,omp_get_wtime,t1,t2,bet_off,d6ji,gx
  logical :: hartree_fock,tda_calculation,COM_calc,r2rms_calc
  external :: dHds_white_gs,dHds_TDA_shell,dHds_TDA_shell_w_1op
  external :: dHds_white_gs_with_1op,dHds_white_gs_with_2op
  external :: dHds_TDA_shell_w_2op

!============================================================
! READ INPUTS SET UP STORAGE STRUCTURE
!============================================================
  call getarg(1,inputs_from_command) 
  call read_main_input_file(inputs_from_command,HS,ham_type,&
       hartree_fock,method_int,tda_calculation,COM_calc,r2rms_calc,hw)
  
  HS%herm = 1
  HS%hospace = hw

  call read_sp_basis(jbasis,HS%Aprot,HS%Aneut) 
  call allocate_blocks(jbasis,HS)   
  call print_header
  ! for calculating COM expectation value
  if (COM_calc) then  
     call duplicate_sq_op(HS,rirj)
     call duplicate_sq_op(HS,pipj)
     call read_interaction(HS,jbasis,ham_type,hw,rr=rirj,pp=pipj)
     call calculate_h0_harm_osc(hw,jbasis,pipj,4)
     call calculate_h0_harm_osc(hw,jbasis,rirj,5)
  else if (r2rms_calc) then
     call duplicate_sq_op(HS,rirj)
     call duplicate_sq_op(HS,r2_rms) 
     call read_interaction(HS,jbasis,ham_type,hw,rr=rirj)
     call initialize_CM_radius(r2_rms,rirj,jbasis) 
  else    
     call read_interaction(HS,jbasis,ham_type,hw)
  end if 
  
  
!============================================================
! BUILD BASIS
!============================================================
  call calculate_h0_harm_osc(hw,jbasis,HS,ham_type) 

  if (hartree_fock) then 
     
     if (COM_calc) then 
        call calc_HF(HS,jbasis,coefs,pipj,rirj)
        call normal_order(pipj,jbasis)
        call normal_order(rirj,jbasis)
     else if (r2rms_calc) then 
        call calc_HF(HS,jbasis,coefs,r2_rms)
        call normal_order(r2_rms,jbasis)
     else
        call calc_HF(HS,jbasis,coefs)
     end if 
    ! calc_HF normal orders the hamiltonian
  
  else 
     call normal_order(HS,jbasis) 
  end if
  
!============================================================
! IM-SRG CALCULATION 
!============================================================ 
 
! just a large series of IF statements deciding exactly which type of
! calculation to run. (The call statements run the full calculation) 
  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ground state decoupling
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  select case (method_int) 
  
  case (1) ! magnus 
    
     if (COM_calc) then 
        call magnus_decouple(HS,jbasis,pipj,rirj)
        call calculate_CM_energy(pipj,rirj,hw) 
     else if (r2rms_calc) then
        call magnus_decouple(HS,jbasis,r2_rms)
        call write_tilde_from_Rcm(r2_rms)
     else 
        call magnus_decouple(HS,jbasis) 
     end if 
     
  case (2) ! traditional
     if (COM_calc) then 
        call decouple_hamiltonian(HS,jbasis,dHds_white_gs_with_2op,pipj,rirj)
        call calculate_CM_energy(pipj,rirj,hw)  ! this writes to file
     else if (r2rms_calc) then 
        call decouple_hamiltonian(HS,jbasis,dHds_white_gs_with_1op,r2_rms)
        call write_tilde_from_Rcm(r2_rms)
     else 
        call decouple_hamiltonian(HS,jbasis,dHds_white_gs) 
    !    call discrete_decouple(HS,jbasis) 
     end if
     
  case (3) 
  
     if (COM_calc) then 
        call discrete_decouple(HS,jbasis,pipj,rirj)
        call calculate_CM_energy(pipj,rirj,hw)  ! this writes to file
     else if (r2rms_calc) then 
        call discrete_decouple(HS,jbasis,r2_rms)
        call write_tilde_from_Rcm(r2_rms)     
     else 
        call discrete_decouple(HS,jbasis) 
     end if

  end select
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! excited state decoupling
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if (tda_calculation) then 
     
     call initialize_TDA(TDA,jbasis,HS%Jtarg,HS%Ptarg,HS%valcut)
     allocate(HS%exlabels(TDA%map(1),2))
     HS%exlabels=TDA%blkM(1)%labels
    
     select case (method_int) 
        case(1) !magnus
    
        if (COM_calc) then 
           call magnus_TDA(HS,TDA,jbasis,pipj,ppTDA,rirj,rrTDA) 
           call calculate_CM_energy_TDA(TDA,ppTDA,rrTDA,hw) 
        else if (r2rms_calc) then
           call magnus_TDA(HS,TDA,jbasis,r2_rms,rrTDA)
        else 
           call magnus_TDA(HS,TDA,jbasis) 
        end if
        
        case(2) !traditional
     
        if (COM_calc) then 
           call TDA_decouple(HS,TDA,jbasis,dHds_TDA_shell_w_2op, &
                pipj,ppTDA,rirj,rrTDA) 
           call calculate_CM_energy_TDA(TDA,ppTDA,rrTDA,hw) 
        else if (r2rms_calc) then
           call TDA_decouple(HS,TDA,jbasis,dHds_TDA_shell_w_1op,&
                r2_rms,rrTDA)
        else 
           call TDA_decouple(HS,TDA,jbasis,dHds_TDA_shell) 
        end if 
        
        case(3) !discrete
      
         
        if (COM_calc) then 
           call discrete_TDA(HS,TDA,jbasis,pipj,ppTDA,rirj,rrTDA) 
           call calculate_CM_energy_TDA(TDA,ppTDA,rrTDA,hw) 
        else if (r2rms_calc) then
           call discrete_TDA(HS,TDA,jbasis,r2_rms,rrTDA)
        else 
           call discrete_TDA(HS,TDA,jbasis) 
        end if 
      
     end select
     
     
  end if 
  
end program main_IMSRG
!=========================================================================
subroutine print_header
  implicit none 
  
  print* 
  print*, 'Constructing Basis...' 
  print*
  print*, '=============================================================='
  print*, '  iter      s             E0        E0+MBPT(2)     |MBPT(2)|  '  
  print*, '=============================================================='

end subroutine   


