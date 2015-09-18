module IMSRG_ODE 
  use commutators
  use adams_ode 
  use operators
  use generators
  use basic_IMSRG
  implicit none 
  
contains
!==============================================================
subroutine decouple_hamiltonian( H , jbas, deriv_calculator,O1,O2) 
  ! runs IMSRG using the specified derivative calculator
  implicit none 

  ! SRG convergence / failsafe / error tolerances
  integer,parameter :: max_steps = 10000
  real(8),parameter :: conv_crit = 1.d-6
  real(8),parameter :: relerr = 1.d-6, abserr = 1.d-6

  type(spd) :: jbas
  type(sq_op) :: H ,HOD
  type(sq_op),optional :: O1,O2
  type(cross_coupled_31_mat) :: HCC
  type(full_sp_block_mat) :: TDA
  integer,dimension(5) :: iwork
  real(8),allocatable,dimension(:) :: cur_vec,work
  integer :: neq,iflag,Atot,Ntot,nh,np,nb,q,steps  
  real(8) :: ds,s,E_old,crit,E_mbpt2
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

  ! add more room if there are other operators
  if (present(O1)) then 
     O1%neq=H%neq
     if (present(O2)) then
        O2%neq = H%neq
        neq = 3*neq
     else 
        neq = 2*neq
     end if 
  end if 
  
  allocate(cur_vec(neq)) ! carries the system into SG solver
  allocate(work(100+21*neq))  ! memory eater
     
  ! parameters for solver
  work = 0.d0
  iwork = 0 
  iflag = 1
  
  ! flow equation variables
  ds = 0.1d0
  s = 0.d0 
  
  steps = 0 
  
  open(unit=36,file='../../output/'//&
       trim(adjustl(prefix))//'_0bflow.dat')
  
  E_mbpt2 = mbpt2(H,jbas) 
  crit = abs(E_mbpt2)
 
  write(36,'(I6,4(e15.7))') steps,s,H%E0,H%E0+E_mbpt2,crit
  write(*,'(I6,4(e15.7))') steps,s,H%E0,H%E0+E_mbpt2,crit
  
  if (present(O1)) then 
     if (present(O2)) then 
        
! main loop   (two operators) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
  do while (steps < max_steps) 
    
     E_old = H%E0 
     ! send info to SG solver
     call vectorize(H,cur_vec(1:H%neq))
     call vectorize(O1,cur_vec(H%neq+1:2*H%neq))
     call vectorize(O2,cur_vec(2*H%neq+1:3*H%neq)) 

     call ode(deriv_calculator,neq,cur_vec,H,jbas,&
          s,s+ds,relerr,abserr,iflag,work,iwork) 

     call repackage(H,cur_vec(1:H%neq)) 
     call repackage(O1,cur_vec(H%neq+1:2*H%neq))
     call repackage(O2,cur_vec(2*H%neq+1:3*H%neq)) 

     steps = steps + 1
  
     ! weak convergence criteria, but it works

     E_mbpt2 = mbpt2(H,jbas) 
     crit = abs(E_mbpt2)
  
     write(36,'(I6,4(e15.7))') steps,s,H%E0,H%E0+E_mbpt2,crit
     write(*,'(I6,4(e15.7))') steps,s,H%E0,H%E0+E_mbpt2,crit
    
     if (crit < conv_crit) exit

  end do 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  else 

  ! main loop (one operator) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
  do while (steps < max_steps) 
    
     E_old = H%E0 
     ! send info to SG solver
     call vectorize(H,cur_vec(1:H%neq))
     call vectorize(O1,cur_vec(H%neq+1:2*H%neq))
    
     call ode(deriv_calculator,neq,cur_vec,H,jbas,&
          s,s+ds,relerr,abserr,iflag,work,iwork) 

     call repackage(H,cur_vec(1:H%neq)) 
     call repackage(O1,cur_vec(H%neq+1:2*H%neq))   
     steps = steps + 1
  
     ! weak convergence criteria, but it works

     E_mbpt2 = mbpt2(H,jbas) 
     crit = abs(E_mbpt2)
  
     write(36,'(I6,5(e15.7))') steps,s,H%E0,H%E0+E_mbpt2,crit
     write(*,'(I6,5(e15.7))') steps,s,H%E0,H%E0+E_mbpt2,crit
    
     if (crit < conv_crit) exit

  end do 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  end if 
else
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

     E_mbpt2 = mbpt2(H,jbas) 
     crit = abs(E_mbpt2)
  
     write(36,'(I6,4(e15.7))') steps,s,H%E0,H%E0+E_mbpt2,crit
     write(*,'(I6,4(e15.7))') steps,s,H%E0,H%E0+E_mbpt2,crit
    
     if (crit < conv_crit) exit

  end do 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end if 
 
  close(36)
end subroutine   
!================================================
!==============================================================
subroutine TDA_decouple( H , TDA, jbas, deriv_calculator,O1,O1TDA,O2,O2TDA ) 
  ! runs IMSRG using the specified derivative calculator
  implicit none 

  ! SRG convergence / failsafe / error tolerances
  integer,parameter :: max_steps = 50
  real(8),parameter :: conv_crit = 1.d-6
  real(8),parameter :: relerr = 1.d-6, abserr = 1.d-6

  type(spd) :: jbas
  type(sq_op) :: H ,HOD
  type(sq_op),optional :: O1,O2 
  type(cross_coupled_31_mat) :: HCC,OeCC
  type(full_sp_block_mat) :: TDA
  type(full_sp_block_mat),optional :: O1TDA,O2TDA
  integer,dimension(5) :: iwork
  real(8),allocatable,dimension(:) :: cur_vec,work,E_old
  integer :: neq,iflag,Atot,Ntot,nh,np,nb,q,steps ,i 
  real(8) :: ds,s,crit,min_crit
  character(200) :: spfile,intfile,prefix
  character(1) :: Jlabel,Plabel
  integer :: Jsing
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
  
  ! add more room if there are other operators
  if (present(O1)) then 
     if (present(O2)) then 
        neq = 3*neq
     else 
        neq = 2*neq
     end if 
  end if 

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
  allocate(E_old(TDA%map(1)))
  call calculate_cross_coupled(H,HCC,jbas,.true.) 
  call calc_TDA(TDA,H,HCC,jbas) 
 
  call diagonalize_blocks(TDA)
  
  E_old = TDA%blkM(1)%eigval
    
  Jsing = H%Jtarg/2
  write(Jlabel,'(I1)') Jsing
  if (H%Ptarg == 0) then 
     Plabel = '+'
  else
     Plabel = '-'
  end if 
  open(unit=37,file='../../output/'//&
       trim(adjustl(prefix))//'_'//Jlabel//Plabel//'_excited.dat')
  
  call write_excited_states(steps,s,TDA,H%E0,37) 
  
!===========================================================
!===========================================================
min_crit = 10000.d0
if (present(O1)) then 
   if (present(O2)) then 
      
! main loop   (two operators) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
  do while (steps < max_steps) 
    
    
     ! send info to SG solver
     call vectorize(H,cur_vec(1:H%neq))
     call vectorize(O1,cur_vec(H%neq+1:2*H%neq))
     call vectorize(O2,cur_vec(2*H%neq+1:3*H%neq)) 
     
     call ode(deriv_calculator,neq,cur_vec,H,jbas,&
          s,s+ds,relerr,abserr,iflag,work,iwork) 
          
     call repackage(H,cur_vec(1:H%neq)) 
     call repackage(O1,cur_vec(H%neq+1:2*H%neq))
     call repackage(O2,cur_vec(2*H%neq+1:3*H%neq)) 
   
     steps = steps + 1
  
     call calculate_cross_coupled(H,HCC,jbas,.true.) 
     call calc_TDA(TDA,H,HCC,jbas) 
     call diagonalize_blocks(TDA)
  

     call write_excited_states(steps,s,TDA,H%E0,37) 
     
     ! convergence criteria
     crit = sum(abs(E_old-TDA%blkM(1)%eigval))/TDA%map(1)
     write(*,'(I4,7(e14.5))') steps,crit,E_old(1:6) 
     E_old = TDA%blkM(1)%eigval
 

     if (crit > 100*min_crit) then
        print*, 'convergence failed' 
        exit
     end if 
     min_crit = min(min_crit,crit) 
     if (crit < conv_crit) exit

  end do 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  else 

  ! main loop (one operator) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
  do while (steps < max_steps) 
    
     ! send info to SG solver
     call vectorize(H,cur_vec(1:H%neq))
     call vectorize(O1,cur_vec(H%neq+1:2*H%neq))
    
     call ode(deriv_calculator,neq,cur_vec,H,jbas,&
          s,s+ds,relerr,abserr,iflag,work,iwork) 

     call repackage(H,cur_vec(1:H%neq)) 
     call repackage(O1,cur_vec(H%neq+1:2*H%neq))   
     steps = steps + 1
  
     call calculate_cross_coupled(H,HCC,jbas,.true.) 
     call calc_TDA(TDA,H,HCC,jbas) 
     call diagonalize_blocks(TDA)
  
     call write_excited_states(steps,s,TDA,H%E0,37) 
     
     ! convergence criteria
     crit = sum(abs(E_old-TDA%blkM(1)%eigval))/TDA%map(1)
     write(*,'(I4,7(e14.5))') steps,crit,E_old(1:6)
!     write(*,'(7(e14.5))') crit,E_old(1:6)
     E_old = TDA%blkM(1)%eigval
 

     if (crit > 100*min_crit) then
        print*, 'convergence failed' 
        exit
     end if 
     min_crit = min(min_crit,crit) 
     if (crit < conv_crit) exit

  end do 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  end if 
else
! main loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
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
     write(*,'(I4,2(e14.5))') steps,crit,E_old(1)
     E_old = TDA%blkM(1)%eigval
 

     if (crit > 100*min_crit) then
        print*, 'convergence failed' 
        exit
     end if 
     min_crit = min(min_crit,crit) 
     if (crit < conv_crit) exit
  end do 
end if 
  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  close(37)
! calculate TDA matrices for operators
  if (present(O1)) then 
     call duplicate_CCMAT(HCC,OeCC)

     if (present(O2)) then 
        call duplicate_sp_mat(TDA,O2TDA)
        allocate(O2TDA%blkM(1)%labels(TDA%map(1),2)) 
        O2TDA%blkM(1)%labels = TDA%blkM(1)%labels   
        call calculate_cross_coupled(O2,OeCC,jbas,.true.)
        call calc_TDA(O2TDA,O2,OeCC,jbas)
     end if 

     ! transform observable
        call duplicate_sp_mat(TDA,O1TDA)
        allocate(O1TDA%blkM(1)%labels(TDA%map(1),2)) 
        O1TDA%blkM(1)%labels = TDA%blkM(1)%labels      
        call calculate_cross_coupled(O1,OeCC,jbas,.true.)
        call calc_TDA(O1TDA,O1,OeCC,jbas)
        
  end if 
  
end subroutine   
!================================================
!================================================
end module
!================================================
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!================================================
subroutine dHds_white_gs(t,yy,yp,HS,jbas) 
  ! calculates the derivatives inside the solver
  ! use in shampine and gordon, modified to include HS,jbas
  use basic_IMSRG
  use commutators
  use TS_commutators
  use generators
  use operators
  implicit none 
  
  real(8) :: t,ex,ex2
  integer :: i,j,k,l,neq,bytes,n,m,p,q,r,s,px,qx,rx,sx,IX,JX,a,b,c,d
  integer :: ji,jj,jk,jl,ja,jb,jc,jd,J1,J2
  real(8) :: yp(*),yy(*),sm,x,d6ji
  type(spd) :: jbas
  type(sq_op) :: HS,ETA,DH,w1,w2,w3,w4
  type(cross_coupled_31_mat) :: WCC,ETACC,HSCC 

!!! we need the sq_op structure to compute the derivatives at max speed
!!! so we allocate a bunch of those to work in 

! ALLOCATE A BUNCH OF WORKSPACE
  call duplicate_sq_op(HS,ETA) !generator
  call duplicate_sq_op(HS,DH) !derivative
!  call duplicate_sq_op(HS,w1) !workspace
!  call duplicate_sq_op(HS,w2) !workspace
  call allocate_CCMAT(HS,HSCC,jbas) ! cross coupled ME
  call duplicate_CCMAT(HSCC,ETACC) !cross coupled ME
  call allocate_CC_wkspc(HSCC,WCC) ! workspace for CCME

 ! call build_gs_wegner(HS,ETA,jbas,HSCC,ETACC,WCC,w1,w2) 
  call build_gs_white(HS,ETA,jbas) ! constructs generator
  
! THIS IS WHERE THE TS COMMUTATOR TESTGROUND BEGINS. 
  
  ETA%herm = -1
  w1%herm = 1
  w1%rank =4 
  call allocate_tensor(jbas,w1,HS) 
  call calculate_EX(w1,jbas)
  
  call duplicate_sq_op(w1,w2) 
  call duplicate_sq_op(w1,w3) 
  call duplicate_sq_op(w1,w4) 
  w2%herm = -1*ETA%herm
  
  
 do i = 1, w1%belowEF
     do j = i, w1%belowEF
        
        ji = jbas%jj(i) 
        jj = jbas%jj(j)

        if (ji .ne. jj) cycle        
        if (jbas%itzp(i) .ne. jbas%itzp(j) ) cycle
        if (jbas%ll(i) .ne. jbas%ll(j) ) cycle

        call random_number(x) 
        x = 10.*x-5.
        ETA%fhh(i,j) = x
        ETA%fhh(j,i) = ETA%herm * x
        
        
     end do 
     
     ETA%fhh(i,i) = (ETA%herm+1)*ETA%fhh(i,i)/2.d0
     
  end do 

 do i = 1, w1%nsp-w1%belowEF
     do j = 1, w1%belowEF
        
        ji = jbas%jj(i+6) 
        jj = jbas%jj(j)

        if (ji .ne. jj) cycle        
        if (jbas%itzp(i+6) .ne. jbas%itzp(j) ) cycle
        if (jbas%ll(i+6) .ne. jbas%ll(j) ) cycle

        call random_number(x) 
        x = 10.*x-5.
        ETA%fph(i,j) = x

     end do 
  end do 

 do i = 1,w1%nsp- w1%belowEF
     do j = i,w1%nsp- w1%belowEF
        
        ji = jbas%jj(i+6) 
        jj = jbas%jj(j+6)

        if (ji .ne. jj) cycle        
        if (jbas%itzp(i+6) .ne. jbas%itzp(j+6) ) cycle
        if (jbas%ll(i+6) .ne. jbas%ll(j+6) ) cycle

        call random_number(x) 
        x = 10.*x-5.
        ETA%fpp(i,j) = x 
        ETA%fpp(j,i) = ETA%herm * x
     end do 
     ETA%fpp(i,i) = (ETA%herm+1)*ETA%fpp(i,i)/2.d0
  end do 
  
  do q=1,ETA%nblocks
     do i = 1, 6
        do IX = 1,size(ETA%mat(q)%gam(i)%X(:,1))
           do JX = 1,size(ETA%mat(q)%gam(i)%X(1,:))
              
              call random_number(x) 
              x = 10.*x-5.
              
              ETA%mat(q)%gam(i)%X(IX,JX) =  x * (1 + ETA%herm *kron_del(IX,JX) ) 
              if (sqs(i)) ETA%mat(q)%gam(i)%X(JX,IX) = ETA%herm*x*  (1 + ETA%herm *kron_del(IX,JX) )  
              
              
           end do 
        end do 
    end do 
  end do 
  
   do i = 1, w1%belowEF
     do j = i, w1%belowEF
        
        ji = jbas%jj(i) 
        jj = jbas%jj(j)

        if (.not. triangle(ji,jj,w1%rank)) cycle        
        if (jbas%itzp(i) .ne. jbas%itzp(j) ) cycle
        if (mod(jbas%ll(i),2) .ne. mod(jbas%ll(j),2)) cycle

        call random_number(x) 
        x = 10.*x-5.
        w1%fhh(i,j) = x
        w1%fhh(j,i) = (-1)**((ji-jj)/2) * x
     end do 
  end do 

  do i = 1, w1%nsp- w1%belowEF
     do j = i, w1%nsp-w1%belowEF
        
        ji = jbas%jj(i+6) 
        jj = jbas%jj(j+6)

        if (.not. triangle(ji,jj,w1%rank)) cycle        
        if (jbas%itzp(i+6) .ne. jbas%itzp(j+6) ) cycle
        if (mod(jbas%ll(i+6),2) .ne. mod(jbas%ll(j+6),2)) cycle

        call random_number(x) 
        x = 10.*x-5.
        w1%fpp(i,j) = x
        w1%fpp(j,i) = (-1)**((ji-jj)/2) * x
     end do 
  end do 
  
  do i = 1, w1%nsp- w1%belowEF
     do j = 1, w1%belowEF
        
        ji = jbas%jj(i+6) 
        jj = jbas%jj(j)

        if (.not. triangle(ji,jj,w1%rank)) cycle        
        if (jbas%itzp(i+6) .ne. jbas%itzp(j) ) cycle
        if (mod(jbas%ll(i+6),2) .ne. mod(jbas%ll(j),2)) cycle

        call random_number(x) 
        x = 10.*x-5.
        w1%fph(i,j) = x
     end do 
  end do 
      
  do q = 1, w1%nblocks
     do i = 1, 9
         
        do IX = 1, size(w1%tblck(q)%tgam(i)%X(:,1) )
           do JX = 1, size(w1%tblck(q)%tgam(i)%X(1,:) )
              
              call random_number(x)
              
              x = 10.*x-5.
              w1%tblck(q)%tgam(i)%X(IX,JX) = x
              
           end do 
        end do 
        
       
           if ( (w1%tblck(q)%Jpair(1) == w1%tblck(q)%Jpair(2)).and.(sqs(i)) ) then
              w1%tblck(q)%tgam(i)%X = &
                   0.5*(w1%tblck(q)%tgam(i)%X + Transpose(w1%tblck(q)%tgam(i)%X))
             
           end if 
       
        
     end do 
  if ( w1%tblck(q)%Jpair(1) == w1%tblck(q)%Jpair(2) ) then
     w1%tblck(q)%tgam(7)%X = Transpose(w1%tblck(q)%tgam(3)%X)
     w1%tblck(q)%tgam(8)%X = Transpose(w1%tblck(q)%tgam(2)%X)
     w1%tblck(q)%tgam(9)%X = Transpose(w1%tblck(q)%tgam(6)%X)
  end if 
  end do 
 
 !  !call calculate_cross_coupled(ETA,ETACC,jbas,.false.) 
  
  call TS_commutator_122(ETA,w1,w2,jbas) 
  
  

   print*, w2%tblck(1)%Jpair
   print*
   print*, w2%tblck(1)%tensor_qn(3,1)%Y(:,1)
   print*, w2%tblck(1)%tensor_qn(3,1)%Y(:,2)
   print*
   print*, w2%tblck(1)%tensor_qn(2,2)%Y(:,1)
   print*, w2%tblck(1)%tensor_qn(2,2)%Y(:,2)
  
   print*
   print*, w2%tblck(8)%Jpair
   print*
   print*, w2%tblck(8)%tensor_qn(3,1)%Y(:,1)
   print*, w2%tblck(8)%tensor_qn(3,1)%Y(:,2)
   print*
   print*, w2%tblck(8)%tensor_qn(2,2)%Y(:,1)
   print*, w2%tblck(8)%tensor_qn(2,2)%Y(:,2)

  do 
     read*, i,j,k,l, J1,J2
     ji=jbas%jj(i)
     jj=jbas%jj(j)
     jk=jbas%jj(k)
     jl=jbas%jj(l)
     
     sm = 0.d0 
     do a = 1, 30 
        sm = sm + f_elem(i,a,ETA,jbas)*tensor_elem(a,j,k,l,J1,J2,w1,jbas) &
                + f_elem(j,a,ETA,jbas)*tensor_elem(i,a,k,l,J1,J2,w1,jbas) &
                - f_elem(a,k,ETA,jbas)*tensor_elem(i,j,a,l,J1,J2,w1,jbas) &
                - f_elem(a,l,ETA,jbas)*tensor_elem(i,j,k,a,J1,J2,w1,jbas)

     end do 
     
     print*, tensor_elem(i,j,k,l,J1,J2,w2,jbas),sm
!    print*, f_elem(i,j,w2,jbas),sm
  end do
   stop
  

  call calculate_cross_coupled(HS,HSCC,jbas,.true.)
  call calculate_cross_coupled(ETA,ETACC,jbas,.false.) 
   
  DH%E0 = commutator_110(ETA,HS,jbas) + commutator_220(ETA,HS,jbas)
  !print*, 'dE0/ds =', DH%E0
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
subroutine dHds_TDA_shell(t,yy,yp,HS,jbas) 
  ! calculates the derivatives inside the solver
  ! use in shampine and gordon, modified to include HS,jbas
  use basic_IMSRG
  use commutators
  use generators
  implicit none 
  
  real(8) :: t,ex,ex2
  integer :: i,j,neq,bytes,n,m,p,q,r,s,px,qx,rx,sx,k,l
  real(8) :: yp(*),yy(*)
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
subroutine dHds_white_gs_with_1op(t,yy,yp,HS,jbas) 
  ! calculates the derivatives inside the solver
  ! use in shampine and gordon, modified to include HS,jbas
  use basic_IMSRG
  use commutators
  use generators
  implicit none 
  
  real(8) :: t,ex,ex2
  integer :: i,j,neq,bytes,n,m,p,q,r,s,px,qx,rx,sx,k,l,num_ops 
  real(8) :: yp(*),yy(*)
  type(spd) :: jbas
  type(sq_op) :: HS,ETA,DH,w1,w2,O1,O2
  type(cross_coupled_31_mat) :: WCC,ETACC,HSCC 

!!! we need the sq_op structure to compute the derivatives at max speed
!!! so we allocate a bunch of those to work in 

  neq = 2*HS%neq
 
! ALLOCATE A BUNCH OF WORKSPACE
  call duplicate_sq_op(HS,ETA) !generator
  call duplicate_sq_op(HS,DH) !derivative
  call duplicate_sq_op(HS,w1) !workspace
  call duplicate_sq_op(HS,w2) !workspace
 
  call duplicate_sq_op(HS,O1)   
  call repackage(O1,yy(HS%neq+1:2*HS%neq)) 
  
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
  call vectorize(DH,yp(1:HS%neq))
  
!===================================================================
  
  call calculate_cross_coupled(O1,HSCC,jbas,.true.)
  
  DH%E0 = commutator_110(ETA,O1,jbas) + commutator_220(ETA,O1,jbas)
 
  call commutator_111(ETA,O1,DH,jbas) 
  call commutator_121(ETA,O1,DH,jbas)
  call commutator_122(ETA,O1,DH,jbas)
  
  call commutator_222_pp_hh(ETA,O1,DH,w1,w2,jbas)
  
  call commutator_221(ETA,O1,DH,w1,w2,jbas)
  call commutator_222_ph(ETACC,HSCC,DH,WCC,jbas)
  
  ! rewrite in a form that shampine and gordon are comfortable with.
  call vectorize(DH,yp(HS%neq+1:2*HS%neq))

end subroutine
!==================================================================
!==================================================================
subroutine dHds_white_gs_with_2op(t,yy,yp,HS,jbas) 
  ! calculates the derivatives inside the solver
  ! use in shampine and gordon, modified to include HS,jbas
  use basic_IMSRG
  use commutators
  use generators
  implicit none 
  
  real(8) :: t,ex,ex2
  integer :: i,j,neq,bytes,n,m,p,q,r,s,px,qx,rx,sx,k,l,num_ops 
  real(8) :: yp(*),yy(*)
  type(spd) :: jbas
  type(sq_op) :: HS,ETA,DH,w1,w2,O1,O2
  type(cross_coupled_31_mat) :: WCC,ETACC,HSCC 

!!! we need the sq_op structure to compute the derivatives at max speed
!!! so we allocate a bunch of those to work in 

  neq = 3*HS%neq
 
! ALLOCATE A BUNCH OF WORKSPACE
  call duplicate_sq_op(HS,ETA) !generator
  call duplicate_sq_op(HS,DH) !derivative
  call duplicate_sq_op(HS,w1) !workspace
  call duplicate_sq_op(HS,w2) !workspace
  
  call duplicate_sq_op(HS,O1) 
  call repackage(O1,yy(HS%neq+1:2*HS%neq)) 

  call duplicate_sq_op(HS,O2)
  call repackage(O2,yy(2*HS%neq+1:3*HS%neq))
  
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
  call vectorize(DH,yp(1:HS%neq))
  
!===================================================================
  
  call calculate_cross_coupled(O1,HSCC,jbas,.true.)
  
  DH%E0 = commutator_110(ETA,O1,jbas) + commutator_220(ETA,O1,jbas)
 
  call commutator_111(ETA,O1,DH,jbas) 
  call commutator_121(ETA,O1,DH,jbas)
  call commutator_122(ETA,O1,DH,jbas)
  
  call commutator_222_pp_hh(ETA,O1,DH,w1,w2,jbas)
  
  call commutator_221(ETA,O1,DH,w1,w2,jbas)
  call commutator_222_ph(ETACC,HSCC,DH,WCC,jbas)
  
  ! rewrite in a form that shampine and gordon are comfortable with.
  call vectorize(DH,yp(HS%neq+1:2*HS%neq))

!===================================================================
  call calculate_cross_coupled(O2,HSCC,jbas,.true.)
   
  DH%E0 = commutator_110(ETA,O2,jbas) + commutator_220(ETA,O2,jbas)

  call commutator_111(ETA,O2,DH,jbas) 
  call commutator_121(ETA,O2,DH,jbas)
  call commutator_122(ETA,O2,DH,jbas)
  
  call commutator_222_pp_hh(ETA,O2,DH,w1,w2,jbas)
  
  call commutator_221(ETA,O2,DH,w1,w2,jbas)
  call commutator_222_ph(ETACC,HSCC,DH,WCC,jbas)
  
  ! rewrite in a form that shampine and gordon are comfortable with.
  call vectorize(DH,yp(2*HS%neq+1:3*HS%neq))
end subroutine
!================================================
subroutine dHds_TDA_shell_w_1op(t,yy,yp,HS,jbas) 
  ! calculates the derivatives inside the solver
  ! use in shampine and gordon, modified to include HS,jbas
  use basic_IMSRG
  use commutators
  use generators
  implicit none 
  
  real(8) :: t,ex,ex2
  integer :: i,j,neq,bytes,n,m,p,q,r,s,px,qx,rx,sx,k,l
  real(8) :: yp(*),yy(*)
  type(spd) :: jbas
  type(sq_op) :: HS,ETA,DH,w1,w2,O1
  type(cross_coupled_31_mat) :: WCC,ETACC,HSCC 

!!! we need the sq_op structure to compute the derivatives at max speed
!!! so we allocate a bunch of those to work in 

  neq = 2*HS%neq
  
! ALLOCATE A BUNCH OF WORKSPACE
  call duplicate_sq_op(HS,ETA) !generator
  call duplicate_sq_op(HS,DH) !derivative
  call duplicate_sq_op(HS,w1) !workspace
  call duplicate_sq_op(HS,w2) !workspace
  call allocate_CCMAT(HS,HSCC,jbas) ! cross coupled ME
  call duplicate_CCMAT(HSCC,ETACC) !cross coupled ME
  call allocate_CC_wkspc(HSCC,WCC) ! workspace for CCME

  call duplicate_sq_op(HS,O1) 
  call repackage(O1,yy(HS%neq+1:2*HS%neq)) 

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
  call vectorize(DH,yp(1:HS%neq))
  
!===================================================================
  
  call calculate_cross_coupled(O1,HSCC,jbas,.true.)
  
  DH%E0 = commutator_110(ETA,O1,jbas) + commutator_220(ETA,O1,jbas)
 
  call commutator_111(ETA,O1,DH,jbas) 
  call commutator_121(ETA,O1,DH,jbas)
  call commutator_122(ETA,O1,DH,jbas)
  
  call commutator_222_pp_hh(ETA,O1,DH,w1,w2,jbas)
  
  call commutator_221(ETA,O1,DH,w1,w2,jbas)
  call commutator_222_ph(ETACC,HSCC,DH,WCC,jbas)
  
  ! rewrite in a form that shampine and gordon are comfortable with.
  call vectorize(DH,yp(HS%neq+1:2*HS%neq))
  
end subroutine 
!==============================================================
!==============================================================
!================================================
subroutine dHds_TDA_shell_w_2op(t,yy,yp,HS,jbas) 
  ! calculates the derivatives inside the solver
  ! use in shampine and gordon, modified to include HS,jbas
  use basic_IMSRG
  use commutators
  use generators
  implicit none 
  
  real(8) :: t,ex,ex2
  integer :: i,j,neq,bytes,n,m,p,q,r,s,px,qx,rx,sx,k,l
  real(8) :: yp(*),yy(*)
  type(spd) :: jbas
  type(sq_op) :: HS,ETA,DH,w1,w2,O1,O2
  type(cross_coupled_31_mat) :: WCC,ETACC,HSCC 

!!! we need the sq_op structure to compute the derivatives at max speed
!!! so we allocate a bunch of those to work in 

  neq = 3*HS%neq
! ALLOCATE A BUNCH OF WORKSPACE
  call duplicate_sq_op(HS,ETA) !generator
  call duplicate_sq_op(HS,DH) !derivative
  call duplicate_sq_op(HS,w1) !workspace
  call duplicate_sq_op(HS,w2) !workspace
  call allocate_CCMAT(HS,HSCC,jbas) ! cross coupled ME
  call duplicate_CCMAT(HSCC,ETACC) !cross coupled ME
  call allocate_CC_wkspc(HSCC,WCC) ! workspace for CCME

  call duplicate_sq_op(HS,O1) 
  call repackage(O1,yy(HS%neq+1:2*HS%neq)) 

  call duplicate_sq_op(HS,O2)
  call repackage(O2,yy(2*HS%neq+1:3*HS%neq))

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
  call vectorize(DH,yp(1:HS%neq))
  
!===================================================================
  
  call calculate_cross_coupled(O1,HSCC,jbas,.true.)
  
  DH%E0 = commutator_110(ETA,O1,jbas) + commutator_220(ETA,O1,jbas)
 
  call commutator_111(ETA,O1,DH,jbas) 
  call commutator_121(ETA,O1,DH,jbas)
  call commutator_122(ETA,O1,DH,jbas)
  
  call commutator_222_pp_hh(ETA,O1,DH,w1,w2,jbas)
  
  call commutator_221(ETA,O1,DH,w1,w2,jbas)
  call commutator_222_ph(ETACC,HSCC,DH,WCC,jbas)
  
  ! rewrite in a form that shampine and gordon are comfortable with.
  call vectorize(DH,yp(HS%neq+1:2*HS%neq))

!===================================================================
  call calculate_cross_coupled(O2,HSCC,jbas,.true.)
   
  DH%E0 = commutator_110(ETA,O2,jbas) + commutator_220(ETA,O2,jbas)

  call commutator_111(ETA,O2,DH,jbas) 
  call commutator_121(ETA,O2,DH,jbas)
  call commutator_122(ETA,O2,DH,jbas)
  
  call commutator_222_pp_hh(ETA,O2,DH,w1,w2,jbas)
  
  call commutator_221(ETA,O2,DH,w1,w2,jbas)
  call commutator_222_ph(ETACC,HSCC,DH,WCC,jbas)
  
  ! rewrite in a form that shampine and gordon are comfortable with.
  call vectorize(DH,yp(2*HS%neq+1:3*HS%neq))
  
end subroutine 
!==============================================================
!==============================================================

