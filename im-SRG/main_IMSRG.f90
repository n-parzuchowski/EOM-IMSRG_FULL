program main_IMSRG
  use isospin_operators
  use HF_mod
  use IMSRG_ODE
  use IMSRG_MAGNUS
  use IMSRG_CANONICAL
  use operators
  use interaction_IO
  use EOM_IMSRG
  use brute_force_testing
  use three_body_routines
  ! ground state IMSRG calculation for nuclear system 
  implicit none
  
  type(spd) :: jbas,jbx
  type(sq_op) :: HS,ETA,DH,w1,w2,rirj,pipj,r2_rms,Otrans,exp_omega,num,cr,H0,Hcm
  type(sq_op),allocatable,dimension(:) :: ladder_ops 
  type(iso_ladder),allocatable,dimension(:) :: isoladder_ops
  type(cc_mat) :: CCHS,CCETA,WCC
  type(full_sp_block_mat) :: coefs,TDA,ppTDA,rrTDA
  type(three_body_force) :: threebod
  type(obsv_mgr) :: trans,moments
  type(eom_mgr) :: eom_states
  character(200) :: inputs_from_command
  character(1) :: quads,trips,trans_type
  integer :: i,j,T,JTot,a,b,c,d,g,q,ham_type,j3,ix,jx,kx,lx,PAR,Tz,trans_rank
  integer :: np,nh,nb,k,l,m,n,method_int,mi,mj,ma,mb,j_min,ex_Calc_int
  integer :: na,la,lb,totstates,numstates,oldnum,qx,dTZ,oldnum_dTz,numstates_dTz
  real(8) :: hw ,sm,omp_get_wtime,t1,t2,bet_off,d6ji,gx,dcgi,dcgi00,pre,x,corr,de_trips
  logical :: hartree_fock,COM_calc,r2rms_calc,me2j,me2b,trans_calc
  logical :: skip_setup,skip_gs,do_HF,TEST_commutators,mortbin,decouple
  external :: build_gs_white,build_specific_space,build_gs_atan,build_gs_w2
  external :: build_ex_imtime,build_sd_shellmodel
!============================================================
! READ INPUTS SET UP STORAGE STRUCTURE
!============================================================
  t1 = omp_get_wtime()

  call getarg(1,inputs_from_command) 
  call getarg(2,resubmitter) 
  
  if (trim(inputs_from_command) == 'X') then 
     test_commutators = .true.
     inputs_from_command = ''
  else
     test_commutators = .false.
  end if


  call read_main_input_file(inputs_from_command,HS,ham_type,&
       hartree_fock,method_int,ex_calc_int,COM_calc,r2rms_calc,me2j,&
       me2b,mortbin,hw,skip_setup,skip_gs,quads,trips,threebod%e3max)

  call read_sp_basis(jbas,HS%Aprot,HS%Aneut,HS%eMax,HS%lmax,trips,jbx)

  call print_system(jbas) 
  
  if (TEST_COMMUTATORS) then 
     ! run this by typing ' X' after the input file in the command line
     ! This takes forever, you might want to comment some of this out. 
     call test
     stop
  end if 
  
  call allocate_blocks(jbas,HS)

  HS%herm = 1
  HS%hospace = hw

!=================================================================
! SET UP OPERATORS
!=================================================================
  if (COM_calc) then 
     call duplicate_sq_op(HS,rirj)
     call duplicate_sq_op(HS,pipj)
     call calculate_pipj(pipj,jbas)
     call calculate_rirj(rirj,jbas)
  end if
  
  if (r2rms_calc) then 
     call duplicate_sq_op(HS,rirj)
     call duplicate_sq_op(HS,r2_rms) 
     call initialize_rms_radius(r2_rms,rirj,jbas) 
  end if 
  
!============================================================
!  CAN WE SKIP STUFF?  
!============================================================
  
  do_hf = .true. 
  IF (reading_decoupled) then 
     do_hf = read_twobody_operator(HS,'decoupled') 
     if (.not. do_hf) goto 91 
  end if  
  
  if (reading_bare) then 
     do_hf = read_twobody_operator(HS,'bare')

     if (.not. do_hf) then
        call read_umat(coefs,jbas)
        goto 90
     end if
  end if 
  ! yes, goto 
!=============================================================
! READ INTERACTION 
!=============================================================
  print*, 'reading 2-body interaction'
  ! for calculating COM expectation value
  
  if (me2j) then 
     call read_me2j_interaction(HS,jbas,jbx,ham_type) 
  else if (me2b) then
     ! pre normal ordered interaction with three body included at No2b       
     print*, 'READING PRE-NORMAL ORDERED INTERACTION FROM HEIKO' 
     call read_me2b_interaction(HS,jbas,ham_type,hw) 
     goto 12 ! skip the normal ordering. 
     ! it's already done.  line 128 or search "bare" 
  else if (mortbin) then 
     call read_gz(HS,jbas,ham_type) 
  else
     call read_interaction(HS,jbas,ham_type)
  end if
    
!============================================================
! BUILD BASIS
!============================================================
  
  call calculate_h0_harm_osc(hw,jbas,HS,ham_type) 
  
  if (threebod%e3Max.ne.0) then 
     print*, 'Reading Three Body Force From file'
     call allocate_three_body_storage(jbas,jbx,threebod,HS%eMax,HS%lmax)
     call read_me3j(threebod,jbas,jbx,HS%eMax,HS%lmax)
  end if 
    
  if (hartree_fock) then 
     call calc_HF(HS,threebod,jbas,coefs)   
  else 
     call normal_order(HS,jbas) 
  end if

  call deallocate_3b(threebod)

  ! lawson 0b term
12 HS%E0 = HS%E0 - HS%lawson_beta * 1.5d0* HS%com_hw
  if (writing_bare) then 
     call write_twobody_operator(HS,'bare')   
     call write_umat(coefs)
  end if

  ! Normal Order Observables 
90 if (hartree_fock) then    
     call observable_to_HF(pipj,coefs,jbas)
     call observable_to_HF(rirj,coefs,jbas)
     call observable_to_HF(r2_rms,coefs,jbas)
  else 
     call normal_order(pipj,jbas) 
     call normal_order(rirj,jbas)
     call normal_order(r2_rms,jbas)
  end if

print*, 'BASIS SETUP COMPLETE' 
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
     decouple=.true.

     if ( reading_omega ) then
        decouple = read_twobody_operator( exp_omega ,'omega' )
     end if 

     if (trips .ne. 'n') then 
        call duplicate_Sq_op(HS,H0)
        call copy_sq_op(HS,H0)
     end if

     if ( decouple ) then 

        call magnus_decouple(HS,exp_omega,jbas,quads,trips,build_gs_w2)    
        if (writing_omega) call write_twobody_operator(exp_omega,'omega')
     else
        print*, 'READ TRANSFORMATION FROM FILE, SKIPPING IMSRG...' 
        call transform_observable_BCH(HS,exp_omega,jbas,quads) 
     end if 
     
     if (COM_calc) then 
        print*, 'TRANSFORMING Hcm'
        call transform_observable_BCH(pipj,exp_omega,jbas,quads)
        call transform_observable_BCH(rirj,exp_omega,jbas,quads)
        print*, 'CALCULTING Ecm' 
        call calculate_CM_energy(pipj,rirj,hw) 
     end if 
     
     if (r2rms_calc) then
        print*, 'TRANSFORMING RADIUS'
        call transform_observable_BCH(r2_rms,exp_omega,jbas,quads)
        call write_tilde_from_Rcm(r2_rms) 
     end if
        
  case (2) ! traditional
     if (COM_calc) then 
        call decouple_hamiltonian(HS,jbas,build_gs_white,pipj,rirj)
        call calculate_CM_energy(pipj,rirj,hw)  ! this writes to file
     else if (r2rms_calc) then 
        call decouple_hamiltonian(HS,jbas,build_gs_white,r2_rms)
!        call write_tilde_from_Rcm(r2_rms)
        print*, sqrt(r2_rms%E0)
    else 
        call decouple_hamiltonian(HS,jbas,build_gs_white) 
    !    call discrete_decouple(HS,jbas) 
     end if
     
  case (3) 
     if (COM_calc) then 
        call discrete_decouple(HS,jbas,pipj,rirj)
        call calculate_CM_energy(pipj,rirj,hw)  ! this writes to file
     else if (r2rms_calc) then 
        call discrete_decouple(HS,jbas,r2_rms)
        call write_tilde_from_Rcm(r2_rms)     
     else 
        call discrete_decouple(HS,jbas) 
     end if

  end select
!============================================================
! store hamiltonian in easiest format for quick reading
!============================================================
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  equation of motion calculation 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (writing_decoupled) then 
     call write_twobody_operator(HS,'decoupled')
  end if


!============TRIPLES MAGNUS=============================================
  if (trips == 'y') then 
!     call enumerate_three_body(threebas,jbas)
     print*, 'computing triples'
     t1 = omp_get_wtime()
     corr =  restore_triples(H0,exp_omega,jbas) 
     t2 = omp_get_wtime()
     print*, 'FINAL ENERGY:', corr + HS%E0,t2-t1
     open(unit=39,file=trim(OUTPUT_DIR)//&
       trim(adjustl(prefix))//'_magnus_triples.dat')
     write(39,'(2(e15.7))') HS%E0,HS%E0+corr
     close(39)
     
  else if (trips == 'C') then
     t1 = omp_get_wtime()
     call duplicate_sq_op(H0,CR)         
     ! completely renormalized bit.
     call CR_EXPAND(CR,exp_omega,H0,jbas,quads) 
     corr =  restore_triples(CR,exp_omega,jbas)
     t2 = omp_get_wtime()
     print*, 'FINAL ENERGY:', corr + HS%E0,t2-t1
     open(unit=39,file=trim(OUTPUT_DIR)//&
       trim(adjustl(prefix))//'_magnus_cr_triples.dat')
     write(39,'(2(e15.7))') HS%E0,HS%E0+corr
     close(39)
  end if
!=======================================================================

91 t2 = omp_get_wtime() 
  write(*,'(A5,f12.7)') 'TIME:', t2-t1

  
  if (ex_calc_int==1) then

     totstates=read_eom_file(trans,moments,eom_states,jbas)! total number of states
     
     allocate(ladder_ops(totstates-eom_states%total_dTz))
     allocate(isoladder_ops(eom_states%total_dTz))     

     oldnum = 0
     oldnum_dTz = 0
     numstates = 0
     numstates_dTz = 0
     do q = 1,eom_states%num

        if (eom_states%dTz(q) == 0 ) then 
           oldnum = oldnum + Numstates
           Numstates = eom_states%number_requested(q)        
           ladder_ops(1+oldnum:Numstates+oldnum)%xindx = q
           call calculate_excited_states(eom_states%ang_mom(q),eom_states%par(q),numstates,HS,&
                jbas,ladder_ops(1+oldnum:Numstates+oldnum))
        else
           oldnum_dTz = oldnum_dTz + Numstates_dTz
           Numstates_dTz = eom_states%number_requested(q)        
           isoladder_ops(1+oldnum_dTz:Numstates_dTz+oldnum_dTz)%xindx = q
           call calculate_isospin_states(eom_states%ang_mom(q),eom_states%par(q),eom_states%dTz(q),&
                numstates_dTZ,HS,jbas,isoladder_ops(1+oldnum_dTz:Numstates_dTz+oldnum_dTz))
        
        end if        
        
        
        
       ! print*
       ! print*, '================================================'
       ! print*, '  J^Pi          E            E+dE       time    '
       ! print*, '================================================'
       ! do qx = 1+oldnum,Numstates+oldnum
       !    t1= omp_get_wtime()
       !    dE_trips=EOM_triples(HS,ladder_ops(qx),jbas)  
       !    t2= omp_get_wtime()
       !    write(*,'(A2,3(f20.10))') eom_states%name(q),&
        !         ladder_ops(qx)%E0,ladder_ops(qx)%E0 + dE_trips,t2-t1
        ! end do
        
     end do

     t2 = omp_get_wtime() 
     write(*,'(A5,f12.7)') 'TIME:', t2-t1

     Otrans%xindx = eom_states%num+1
     trans_type = trans%oper(1:1)

     read(trans%oper(2:2),'(I1)') trans_rank

     call initialize_transition_operator(trans_type,trans_rank,Otrans,HS,jbas,trans_calc)
     
     if (trans_calc) then 
        if (hartree_fock) then 
           call observable_to_HF(Otrans,coefs,jbas)
        end if
        print*, 'Transforming transition operator...' 

        if ( trans%num + moments%num > 0 ) call transform_observable_BCH(Otrans,exp_omega,jbas,quads)

        if (com_calc) then 
           Hcm%rank = 0
           Hcm%dpar = 0
           Hcm%xindx = Otrans%xindx + 1 
           call build_Hcm(pipj,rirj,Hcm,jbas)
        end if
        
        call EOM_observables( ladder_ops, isoladder_ops, Otrans, HS, Hcm,trans, moments,eom_states,jbas)
        
     end if
     
  end if
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! excited state decoupling
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (ex_calc_int==2) then 

     call initialize_TDA(TDA,jbas,HS%Jtarg,HS%Ptarg,HS%valcut)

     deallocate(HS%exlabels) 

     allocate(HS%exlabels(TDA%map(1),2))
     HS%exlabels=TDA%blkM(1)%labels

     select case (method_int) 
        case(1) !magnus
           
           call magnus_TDA(HS,TDA,exp_omega,jbas,quads,build_specific_space)     
       
           if (COM_calc) then 
              call calculate_CM_energy_TDA(TDA,rirj,pipj,ppTDA,rrTDA,hw) 
           else if (r2rms_calc) then
              print*, 'fail'
           else 
              print*, 'superfail'
        end if
        
        case(2) !traditional
     
        if (COM_calc) then 
           call TDA_decouple(HS,TDA,jbas,build_specific_space, &
                pipj,ppTDA,rirj,rrTDA) 
           call calculate_CM_energy_TDA(TDA,rirj,pipj,ppTDA,rrTDA,hw) 
        else if (r2rms_calc) then
           call TDA_decouple(HS,TDA,jbas,build_specific_space,&
                r2_rms,rrTDA)
           print*, rrTDA%blkM(1)%eigval
        else 
           call TDA_decouple(HS,TDA,jbas,build_specific_space) 
!           call calculate_excited_states(HS%Jtarg,HS%Ptarg,10,HS,jbas,Otrans) 
        end if 
        
        case(3) !discrete
      
         
        if (COM_calc) then 
           call discrete_TDA(HS,TDA,jbas,pipj,ppTDA,rirj,rrTDA) 
           call calculate_CM_energy_TDA(TDA,rirj,pipj,ppTDA,rrTDA,hw) 
        else if (r2rms_calc) then
           call discrete_TDA(HS,TDA,jbas,r2_rms,rrTDA)
        else 
           call discrete_TDA(HS,TDA,jbas) 
        end if 
      
     end select
     
!     t2 = omp_get_wtime() 
     write(*,'(A5,f12.7)') 'TIME:', t2-t1
  
  end if 

  ! free it up brah
  deallocate(isoladder_ops,ladder_ops)
  
contains

subroutine test
  ! call compare_tensor_scalar_commutator(jbas,-1,1) 
  ! stop
  ! deallocate(jbas%xmap) 
  ! call test_scalar_scalar_commutator(jbas,-1,1) 
  ! deallocate(jbas%xmap)
  ! call test_EOM_scalar_scalar_commutator(jbas,1,1)
  ! deallocate(jbas%xmap)
!   call test_EOM_scalar_tensor_commutator(jbas,1,1,6,2)  
!   deallocate(jbas%xmap,jbas%xmap_tensor,phase_hh,phase_pp)
!   deallocate(half6j%tp_mat)
!  call test_scalar_tensor_commutator(jbas,-1,1,6,2) 
 !  call test_tensor_product(jbas,1,1,2,4,6,2,2,0) 
!  call test_EOM_iso_commutator(jbas,1,1,4,0,0)
!  call test_scalar_iso_commutator(jbas,-1,1,6,2,1) !butt
  call test_tensor_dTZ_product(jbas,1,1,4,4,4,2,0,2,1) 
  
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


