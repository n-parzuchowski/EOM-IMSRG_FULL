module IMSRG_ODE 
  use commutators
  use adams_ode 
  use operators
  use generators
  use basic_IMSRG
  implicit none 
  
contains
!==============================================================
subroutine decouple_hamiltonian( H , jbas, deriv_calculator ) 
  ! runs IMSRG using the specified derivative calculator
  implicit none 

  ! SRG convergence / failsafe / error tolerances
  integer,parameter :: max_steps = 10000
  real(8),parameter :: conv_crit = 1.d-6
  real(8),parameter :: relerr = 1.d-6, abserr = 1.d-6

  type(spd) :: jbas
  type(sq_op) :: H ,HOD
  type(cross_coupled_31_mat) :: HCC
  type(full_sp_block_mat) :: TDA
  integer,dimension(5) :: iwork
  real(8),allocatable,dimension(:) :: cur_vec,work
  integer :: neq,iflag,Atot,Ntot,nh,np,nb,q,steps  
  real(8) :: ds,s,E_old,crit
  character(200) :: spfile,intfile,prefix
  logical :: com_calc 
  external :: deriv_calculator 
  common /files/ spfile,intfile,prefix
  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! first figure out how many equations there are:
  Atot = H%belowEF
  Ntot = H%Nsp
  

  neq = 1 + Atot*Atot + Atot*(Ntot-Atot) + (Ntot - Atot)**2 
  
  do q = 1, H%nblocks
     
     nh = H%mat(q)%nhh
     np = H%mat(q)%npp
     nb = H%mat(q)%nph 
     
     neq = neq + (nh*nh+nh +  nb*nb+nb + np*np+np)/2 + nb*np + nh*np + nh*nb 
  end do 

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  H%neq = neq
  
  allocate(cur_vec(neq)) ! carries the system into SG solver
  allocate(work(100+21*neq))  ! memory eater
     
  ! parameters for solver
  work = 0.d0
  iwork = 0 
  iflag = 1
  
  ! flow equation variables
  ds = .1d0
  s = 0.d0 
  
  steps = 0 
  
  open(unit=36,file='../../output/'//&
       trim(adjustl(prefix))//'_0bflow.dat')
  
  write(36,'(I6,3(e14.6))') steps,s,H%E0,crit

! main loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
  do while (steps < max_steps) 
    
     E_old = H%E0 
     ! send info to SG solver
     call vectorize(H,cur_vec)
    
     call ode(deriv_calculator,neq,cur_vec,H,jbas,&
          s,s+ds,relerr,abserr,iflag,work,iwork) 
          
     call repackage(H,cur_vec) 
        
     steps = steps + 1
  
     ! weak convergence criteria, but it works
     crit = abs(H%E0 - E_old) 
     
     write(36,'(I6,3(e14.6))') steps,s,H%E0,crit     
     write(*,'(I6,3(e14.6))') steps,s,H%E0,crit
     if (crit < conv_crit) exit

  end do 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  close(36)
end subroutine   
!================================================
!==============================================================
subroutine TDA_decouple( H , jbas, deriv_calculator ) 
  ! runs IMSRG using the specified derivative calculator
  implicit none 

  ! SRG convergence / failsafe / error tolerances
  integer,parameter :: max_steps = 2000
  real(8),parameter :: conv_crit = 1.d-5
  real(8),parameter :: relerr = 1.d-5, abserr = 1.d-5

  type(spd) :: jbas
  type(sq_op) :: H ,HOD
  type(cross_coupled_31_mat) :: HCC
  type(full_sp_block_mat) :: TDA
  integer,dimension(5) :: iwork
  real(8),allocatable,dimension(:) :: cur_vec,work,E_old
  integer :: neq,iflag,Atot,Ntot,nh,np,nb,q,steps ,i 
  real(8) :: ds,s,crit,min_crit
  character(200) :: spfile,intfile,prefix
  external :: deriv_calculator 
  common /files/ spfile,intfile,prefix

 ! E_old = H%E0 ! for convergence check
  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! first figure out how many equations there are:
  Atot = H%belowEF
  Ntot = H%Nsp
  
  neq = 1 + Atot*Atot + Atot*(Ntot-Atot) + (Ntot - Atot)**2 
  
  do q = 1, H%nblocks
     
     nh = H%mat(q)%nhh
     np = H%mat(q)%npp
     nb = H%mat(q)%nph 
     
     neq = neq + (nh*nh+nh +  nb*nb+nb + np*np+np)/2 + nb*np + nh*np + nh*nb 
  end do 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  H%neq = neq
  
  allocate(cur_vec(neq)) ! carries the system into SG solver
  allocate(work(100+21*neq))  ! memory eater
     
  ! parameters for solver
  work = 0.d0
  iwork = 0 
  iflag = 1
  
  ! flow equation variables
  ds = .01d0
  s = 0.d0 
  
  steps = 0 

  call allocate_CCMAT(H,HCC,jbas) 
  call initialize_TDA(TDA,jbas,H%Jtarg,H%Ptarg,H%valcut)
  allocate(H%exlabels(TDA%map(1),2))
  H%exlabels=TDA%blkM(1)%labels
  allocate(E_old(TDA%map(1)))
  call calculate_cross_coupled(H,HCC,jbas,.true.) 
  call calc_TDA(TDA,H,HCC,jbas) 
  call diagonalize_blocks(TDA)
  
     
  E_old = TDA%blkM(1)%eigval
    
  open(unit=37,file='../../output/'//&
       trim(adjustl(prefix))//'_excited.dat')
  
  call write_excited_states(steps,s,TDA,H%E0,37) 
  
!===========================================================
!===========================================================
! main loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
  min_crit = 10000.d0
  do while (steps < max_steps) 
    

     ! send info to SG solver
     call vectorize(H,cur_vec)
     call ode(deriv_calculator,neq,cur_vec,H,jbas,&
          s,s+ds,relerr,abserr,iflag,work,iwork) 
     call repackage(H,cur_vec) 
        
     steps = steps + 1
     
     call calculate_cross_coupled(H,HCC,jbas,.true.) 
     call calc_TDA(TDA,H,HCC,jbas) 
     call diagonalize_blocks(TDA)
  
     call write_excited_states(steps,s,TDA,H%E0,37) 
     
     ! convergence criteria
     crit = sum(abs(E_old-TDA%blkM(1)%eigval))/TDA%map(1)
     write(*,'(7(e14.5))') crit,E_old(1:6)
     E_old = TDA%blkM(1)%eigval
 

     if (crit > 100*min_crit) then
        print*, 'convergence failed' 
        exit
     end if 
     min_crit = min(min_crit,crit) 
     if (crit < conv_crit) exit
  end do 
  
  open(unit=38,file='../../output/'//&
       trim(adjustl(prefix))//'_final_excited.dat')
  
  do q = 1, TDA%blocks
     write(38,'(2(I4))') TDA%blkM(q)%lmda(1:2) 
     do i = 1, TDA%map(q)
        write(38,'(e14.6)') TDA%blkM(q)%eigval(i) + H%E0 
     end do 
  end do 
     
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  close(37)
  close(38)
end subroutine   
!================================================
!================================================
end module
!================================================
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!================================================
subroutine dHds_white_gs(t,yp,HS,jbas) 
  ! calculates the derivatives inside the solver
  ! use in shampine and gordon, modified to include HS,jbas
  use basic_IMSRG
  use commutators
  use generators
  implicit none 
  
  real(8) :: t,ex,ex2
  integer :: i,j,neq,bytes,n,m,p,q,r,s,px,qx,rx,sx,k,l
  real(8) :: yp(*)
  type(spd) :: jbas
  type(sq_op) :: HS,ETA,DH,w1,w2
  type(cross_coupled_31_mat) :: WCC,ETACC,HSCC 

!!! we need the sq_op structure to compute the derivatives at max speed
!!! so we allocate a bunch of those to work in 

! ALLOCATE A BUNCH OF WORKSPACE
  call duplicate_sq_op(HS,ETA) !generator
  call duplicate_sq_op(HS,DH) !derivative
  call duplicate_sq_op(HS,w1) !workspace
  call duplicate_sq_op(HS,w2) !workspace
  call allocate_CCMAT(HS,HSCC,jbas) ! cross coupled ME
  call duplicate_CCMAT(HSCC,ETACC) !cross coupled ME
  call allocate_CC_wkspc(HSCC,WCC) ! workspace for CCME

  call build_gs_white(HS,ETA,jbas) ! constructs generator
  
  call calculate_cross_coupled(HS,HSCC,jbas,.true.)
  call calculate_cross_coupled(ETA,ETACC,jbas,.false.) 
   
  DH%E0 = commutator_110(ETA,HS,jbas) + commutator_220(ETA,HS,jbas)

  call commutator_111(ETA,HS,DH,jbas) 
  call commutator_121(ETA,HS,DH,jbas)
  call commutator_122(ETA,HS,DH,jbas)
  
  call commutator_222_pp_hh(ETA,HS,DH,w1,w2,jbas)
  
  call commutator_221(ETA,HS,DH,w1,w2,jbas)
  call commutator_222_ph(ETACC,HSCC,DH,WCC,jbas)

  ! rewrite in a form that shampine and gordon are comfortable with.
  call vectorize(DH,yp)
  
end subroutine 
!================================================
subroutine dHds_TDA_shell(t,yp,HS,jbas) 
  ! calculates the derivatives inside the solver
  ! use in shampine and gordon, modified to include HS,jbas
  use basic_IMSRG
  use commutators
  use generators
  implicit none 
  
  real(8) :: t,ex,ex2
  integer :: i,j,neq,bytes,n,m,p,q,r,s,px,qx,rx,sx,k,l
  real(8) :: yp(*)
  type(spd) :: jbas
  type(sq_op) :: HS,ETA,DH,w1,w2
  type(cross_coupled_31_mat) :: WCC,ETACC,HSCC 

!!! we need the sq_op structure to compute the derivatives at max speed
!!! so we allocate a bunch of those to work in 

! ALLOCATE A BUNCH OF WORKSPACE
  call duplicate_sq_op(HS,ETA) !generator
  call duplicate_sq_op(HS,DH) !derivative
  call duplicate_sq_op(HS,w1) !workspace
  call duplicate_sq_op(HS,w2) !workspace
  call allocate_CCMAT(HS,HSCC,jbas) ! cross coupled ME
  call duplicate_CCMAT(HSCC,ETACC) !cross coupled ME
  call allocate_CC_wkspc(HSCC,WCC) ! workspace for CCME

  !call build_valence_decouple(HS,ETA,jbas) ! constructs generator
  ! call build_ex_imtime(HS,ETA,jbas)
  call build_specific_space(HS,ETA,jbas)
  
  call calculate_cross_coupled(HS,HSCC,jbas,.true.)
  call calculate_cross_coupled(ETA,ETACC,jbas,.false.) 
   
  DH%E0 = commutator_110(ETA,HS,jbas) + commutator_220(ETA,HS,jbas)
  
  call commutator_111(ETA,HS,DH,jbas) 
  call commutator_121(ETA,HS,DH,jbas)
  call commutator_122(ETA,HS,DH,jbas)
  
  call commutator_222_pp_hh(ETA,HS,DH,w1,w2,jbas)
  
  call commutator_221(ETA,HS,DH,w1,w2,jbas)
  call commutator_222_ph(ETACC,HSCC,DH,WCC,jbas)

  ! rewrite in a form that shampine and gordon are comfortable with.
  call vectorize(DH,yp)
  
end subroutine 
!==============================================================
!==============================================================
