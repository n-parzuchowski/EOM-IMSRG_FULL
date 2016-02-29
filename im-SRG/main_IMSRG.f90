program main_IMSRG
  use basic_IMSRG
  use HF_mod
  use IMSRG_ODE
  use IMSRG_MAGNUS
  use IMSRG_CANONICAL
  use operators
  use me2j_format
  use EOM_IMSRG
  use brute_force_testing
  use three_body_routines
  ! ground state IMSRG calculation for nuclear system 
  implicit none
  
  type(spd) :: jbasis
  type(sq_op) :: HS,ETA,DH,w1,w2,rirj,pipj,r2_rms,Otrans,exp_omega,num
  type(sq_op),allocatable,dimension(:) :: ladder_ops 
  type(cc_mat) :: CCHS,CCETA,WCC
  type(full_sp_block_mat) :: coefs,TDA,ppTDA,rrTDA
  type(three_body_force) :: threebod
  character(200) :: inputs_from_command
  character(1) :: quads,trips,trans_type
  integer :: i,j,T,JTot,a,b,c,d,g,q,ham_type,j3,ix,jx,kx,lx,PAR,Tz,trans_rank
  integer :: np,nh,nb,k,l,m,n,method_int,mi,mj,ma,mb,j_min,ex_Calc_int
  integer :: na,la,lb
  real(8) :: hw ,sm,omp_get_wtime,t1,t2,bet_off,d6ji,gx,dcgi,dcgi00,pre,x
  logical :: hartree_fock,COM_calc,r2rms_calc,me2j,me2b,trans_calc
  logical :: skip_setup,skip_gs,writing,TEST_commutators,mortbin,write_omega
  external :: build_gs_white,build_specific_space,build_hartree_fock_gen
  integer :: heiko(30)
!============================================================
! READ INPUTS SET UP STORAGE STRUCTURE
!============================================================
  
  t1 = omp_get_wtime()
  heiko = (/1,2,5,6,3,4,11,12,9,10,7,8,19,20,17,18,15,16,&
       13,14,29,30,27,28,25,26,23,24,21,22/)   

  call getarg(1,inputs_from_command) 
  
  if (trim(inputs_from_command) == 'X') then 
     test_commutators = .true.
     inputs_from_command = ''
  else
     test_commutators = .false.
  end if

  call read_main_input_file(inputs_from_command,HS,ham_type,&
       hartree_fock,method_int,ex_calc_int,COM_calc,r2rms_calc,me2j,&
       me2b,mortbin,hw,skip_setup,skip_gs,quads,trips,&
       trans_type,trans_rank,threebod%e3max)

  call read_sp_basis(jbasis,HS%Aprot,HS%Aneut,method_int)

  if (TEST_COMMUTATORS) then 
     ! run this by typing ' X' after the input file in the command line
     ! This takes forever, you might want to comment some of this out. 
     call test
     stop
  end if 
  
  call allocate_blocks(jbasis,HS)
   
  HS%herm = 1
  HS%hospace = hw

  call initialize_transition_operator&
       (trans_type,trans_rank,Otrans,HS,jbasis,trans_calc)  
  
  print*, f_tensor_elem(1,27,Otrans,jbasis),f_tensor_elem(27,1,Otrans,jbasis)
  ! for calculating COM expectation value
  if (COM_calc) then  
     
     call duplicate_sq_op(HS,rirj)
     call duplicate_sq_op(HS,pipj)
     
     if (me2j) then  ! heiko format or scott format? 
        call read_me2j_interaction(HS,jbasis,ham_type,hw,rr=rirj,pp=pipj) 
     else if (mortbin) then 
        call read_gz(HS,jbasis,ham_type,hw,rr=rirj,pp=pipj)
     else
        call read_interaction(HS,jbasis,ham_type,hw,rr=rirj,pp=pipj)
     end if
     
     call calculate_h0_harm_osc(hw,jbasis,pipj,4)
     call calculate_h0_harm_osc(hw,jbasis,rirj,5)
 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  else if (r2rms_calc) then ! radius of some sort 
  
     call duplicate_sq_op(HS,rirj)
     call duplicate_sq_op(HS,r2_rms) 
     
     if (me2j) then 
        call read_me2j_interaction(HS,jbasis,ham_type,hw,rr=rirj)
     else if (mortbin) then 
        call read_gz(HS,jbasis,ham_type,hw,rr=rirj)  
     else
        call read_interaction(HS,jbasis,ham_type,hw,rr=rirj)
     end if
     
     call initialize_rms_radius(r2_rms,rirj,jbasis) 
 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  else    ! normal boring
  
     if (me2j) then 
        call read_me2j_interaction(HS,jbasis,ham_type,hw) 
     else if (me2b) then
        ! pre normal ordered interaction with three body included at No2b       
        call read_me2b_interaction(HS,jbasis,ham_type,hw) 
        goto 12 ! skip the normal ordering. 
        ! it's already done.  line 128 or search "bare" 
     else if (mortbin) then 
        call read_gz(HS,jbasis,ham_type,hw) 
     else
        call read_interaction(HS,jbasis,ham_type,hw)
     end if
     
  end if
    
!============================================================
! BUILD BASIS
!============================================================
  
  call calculate_h0_harm_osc(hw,jbasis,HS,ham_type) 
  
  writing = .false. 
  if (threebod%e3Max.ne.0) then 
     writing = read_twobody_operator( HS,'bare' )
     if ( writing.or.COM_calc.or.r2rms_calc.or.trans_calc ) then 
        print*, 'Reading Three Body Force'
        call allocate_three_body_storage(jbasis,threebod)
        call read_me3j(threebod,jbasis)
     else 
        goto 77 
     end if
  end if 

  
  if (hartree_fock) then 
  
    if (COM_calc) then 
       call calc_HF(HS,threebod,jbasis,coefs,pipj,rirj)
       call normal_order(pipj,jbasis)
       call normal_order(rirj,jbasis)
     else if (r2rms_calc) then 
        call calc_HF(HS,threebod,jbasis,coefs,r2_rms)
        call normal_order(r2_rms,jbasis)
     else if (trans_calc) then 
        call calc_HF(HS,threebod,jbasis,coefs,Otrans)
     else
        call calc_HF(HS,threebod,jbasis,coefs)
     end if 
    ! calc_HF normal orders the hamiltonian
  
  else 
     if (COM_calc) then
        call normal_order(pipj,jbasis)
        call normal_order(rirj,jbasis)
     else if (r2rms_calc) then 
        call normal_order(r2_rms,jbasis)
     end if 

     call normal_order(HS,jbasis) 
  end if

  call deallocate_3b(threebod)

  ! lawson 0b term
12 HS%E0 = HS%E0 - HS%lawson_beta * 1.5d0* HS%com_hw
77 if (writing) call write_twobody_operator(HS,'bare')
!============================================================
! IM-SRG CALCULATION 
!============================================================ 
 
! just a large series of IF statements deciding exactly which type of
! calculation to run. (The call statements run the full calculation) 
  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ground state decoupling
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  call print_header
  select case (method_int) 
  
  case (1) ! magnus 
     
     call duplicate_sq_op(HS,exp_omega)
     exp_omega%herm = -1
     write_omega = read_twobody_operator( exp_omega ,'omega' )     
!     write_omega=.true.
     if ( write_omega ) then 
        call magnus_decouple(HS,exp_omega,jbasis,quads,trips,build_gs_white)    
        call write_twobody_operator(exp_omega,'omega')
     else
        print*, 'READ TRANSFORMATION FROM FILE, SKIPPING IMSRG...' 
        call transform_observable_BCH(HS,exp_omega,jbasis,quads) 
     end if 
     
     if (COM_calc) then 
        print*, 'TRANSFORMING Hcm'
        call transform_observable_BCH(pipj,exp_omega,jbasis,quads)
        call transform_observable_BCH(rirj,exp_omega,jbasis,quads)
        print*, 'CALCULTING Ecm' 
        call calculate_CM_energy(pipj,rirj,hw) 
     end if 
     
     if (r2rms_calc) then
        print*, 'TRANSFORMING RADIUS'
        call transform_observable_BCH(r2_rms,exp_omega,jbasis,quads)
        call write_tilde_from_Rcm(r2_rms) 
     end if
    
     if (trans_calc) then
        print*, 'TRANSFORMING TRANSITION OPERATOR'
        call transform_observable_BCH(Otrans,exp_omega,jbasis,quads)
     end if
    
  case (2) ! traditional
     if (COM_calc) then 
        call decouple_hamiltonian(HS,jbasis,build_gs_white,pipj,rirj)
        call calculate_CM_energy(pipj,rirj,hw)  ! this writes to file
     else if (r2rms_calc) then 
        call decouple_hamiltonian(HS,jbasis,build_gs_white,r2_rms)
!        call write_tilde_from_Rcm(r2_rms)
        print*, sqrt(r2_rms%E0)
    else 
        call decouple_hamiltonian(HS,jbasis,build_gs_white) 
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
!============================================================
! store hamiltonian in easiest format for quick reading
!============================================================
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  equation of motion calculation 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  t2 = omp_get_wtime() 
  write(*,'(A5,f12.7)') 'TIME:', t2-t1
  
  if (ex_calc_int==1) then 
     call calculate_excited_states(HS%Jtarg,HS%Ptarg,10,HS,jbasis,Otrans) 
     t2 = omp_get_wtime() 
     write(*,'(A5,f12.7)') 'TIME:', t2-t1
  end if
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! excited state decoupling
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (ex_calc_int==2) then 

     call initialize_TDA(TDA,jbasis,HS%Jtarg,HS%Ptarg,HS%valcut)

     deallocate(HS%exlabels) 

     allocate(HS%exlabels(TDA%map(1),2))
     HS%exlabels=TDA%blkM(1)%labels

     select case (method_int) 
        case(1) !magnus
           
           call magnus_TDA(HS,TDA,exp_omega,jbasis,quads,build_specific_space)     
       
           if (COM_calc) then 
              call calculate_CM_energy_TDA(TDA,rirj,pipj,ppTDA,rrTDA,hw) 
           else if (r2rms_calc) then
              print*, 'fuck'
           else 
              print*, 'superfuck'
        end if
        
        case(2) !traditional
     
        if (COM_calc) then 
           call TDA_decouple(HS,TDA,jbasis,build_gs_white, &
                pipj,ppTDA,rirj,rrTDA) 
           call calculate_CM_energy_TDA(TDA,rirj,pipj,ppTDA,rrTDA,hw) 
        else if (r2rms_calc) then
           call TDA_decouple(HS,TDA,jbasis,build_gs_white,&
                r2_rms,rrTDA)
           print*, rrTDA%blkM(1)%eigval
        else 
           call TDA_decouple(HS,TDA,jbasis,build_gs_white) 
        end if 
        
        case(3) !discrete
      
         
        if (COM_calc) then 
           call discrete_TDA(HS,TDA,jbasis,pipj,ppTDA,rirj,rrTDA) 
           call calculate_CM_energy_TDA(TDA,rirj,pipj,ppTDA,rrTDA,hw) 
        else if (r2rms_calc) then
           call discrete_TDA(HS,TDA,jbasis,r2_rms,rrTDA)
        else 
           call discrete_TDA(HS,TDA,jbasis) 
        end if 
      
     end select
     
!     t2 = omp_get_wtime() 
     write(*,'(A5,f12.7)') 'TIME:', t2-t1
  
  end if 

contains

subroutine test

  call test_scalar_scalar_commutator(jbasis,-1,1) 
  deallocate(jbasis%xmap)
  call test_EOM_scalar_scalar_commutator(jbasis,1,1)
  deallocate(jbasis%xmap)
  call test_EOM_scalar_tensor_commutator(jbasis,1,1,4,0)  
  deallocate(jbasis%xmap,jbasis%xmap_tensor,phase_hh,phase_pp)
  deallocate(half6j%tp_mat)
  call test_scalar_tensor_commutator(jbasis,1,1,4,0) 

end subroutine test
end program main_IMSRG
!=========================================================================
subroutine print_header
  implicit none 
  
  print* 
  print*, 'Constructing Basis...' 
  print*
  print*, '================================'//&
       '==================================='
  print*, '  iter        s            E0      '//&
       '    E0+MBPT(2)      |MBPT(2)|  '  
  print*, '================================'//&
       '==================================='

end subroutine   


