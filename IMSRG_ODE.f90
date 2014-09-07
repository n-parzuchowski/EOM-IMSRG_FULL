module IMSRG_ODE 
  use commutators
  use adams_ode 
  use basic_IMSRG
  implicit none 
  
contains

subroutine decouple_hamiltonian( H , jbas, deriv_calculator ) 
  ! runs IMSRG using the specified derivative calculator
  implicit none 
  
  integer :: neq,iflag,Atot,Ntot,nh,np,nb,q,steps
  type(spd) :: jbas
  type(sq_op) :: H 
  external :: deriv_calculator 
  integer,dimension(5) :: iwork
  real(8),allocatable,dimension(:) :: cur_vec,work
  real(8),parameter :: relerr = 1.d-6, abserr = 1.d-6
  real(8) :: ds,s,E_old,crit

  E_old = H%E0 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! first figure out how many equations there are:
  Atot = H%belowEF
  Ntot = H%Nsp
  
  neq = 1 + Atot*Atot + Atot*(Ntot-Atot) * (Ntot - Atot)**2 
  
  do q = 1, H%nblocks
     
     nh = H%mat(q)%nhh
     np = H%mat(q)%npp
     nb = H%mat(q)%nph 
     
     neq = neq + (nh*nh+nh +  nb*nb+nb + np*np+np)/2 + nb*np + nh*np * nh*nb 
  end do 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  H%neq = neq
  print*, neq
  allocate(cur_vec(neq)) 
  allocate(work(100+21*neq)) 
  iflag = 1
  
  ds = .1d0
  s = 0.d0 
  
  steps = 0 
 
  do while (steps < 10000) 
    
     E_old = H%E0 
     call vectorize(H,cur_vec) 
     call ode(deriv_calculator,neq,cur_vec,H,jbas,&
          s,s+ds,relerr,abserr,iflag,work,iwork) 
     call repackage(H,cur_vec) 
  
     steps = steps + 1
     !crit = sqrt(sum(H%mat(1)%gam(3)%X**2))
     crit = abs(H%E0 - E_old)
     print*, steps,s,H%E0,crit
     
   !  crit = abs((H%E0 - E_old)/E_old)
     if (crit < 1e-6) exit
  end do 
  
end subroutine   
end module
!================================================
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!================================================
subroutine dHds_white_gs(t,yp,HS,jbas) 
  ! calculates the derivatives inside the solver
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

!!! we need the full_ham structure to compute the derivatives at max speed
!!! so we allocate a bunch of those to work in 

! ALLOCATE A BUNCH OF WORKSPACE

  call duplicate_sq_op(HS,ETA) !generator
  call duplicate_sq_op(HS,DH) !derivative
  call duplicate_sq_op(HS,w1) !workspace
  call duplicate_sq_op(HS,w2) !workspace
  call allocate_CCMAT(HS,HSCC,jbas) !cross coupled ME
  call duplicate_CCMAT(HSCC,ETACC) !cross coupled ME
  call allocate_CC_wkspc(HSCC,WCC) ! workspace for CCME
    
  call build_gs_white(HS,ETA,jbas) 
 

  call calculate_cross_coupled(HS,HSCC,jbas,.true.)
  call calculate_cross_coupled(ETA,ETACC,jbas,.false.) 
 
  DH%E0 = commutator_110(ETA,HS,jbas) + commutator_220(ETA,HS,jbas)

  call commutator_111(ETA,HS,DH,jbas) 
  call commutator_121(ETA,HS,DH,jbas)
  call commutator_122(ETA,HS,DH,jbas)
  

  call commutator_222_pp_hh(ETA,HS,DH,w1,w2,jbas)
  

  call commutator_221(ETA,HS,DH,w1,w2,jbas)
  call commutator_222_ph(ETACC,HSCC,DH,WCC,jbas)

  ! the problem is with this matrix element. out of phase somehow? 
  ! seems to be the case with most phph, phhh,phpp
  ! figure it out hoe.
  !    do i = 1, HS%nsp
  !       do j = i ,HS%nsp
  !          do k = 1, HS%nsp
  !             do l = k, HS%nsp
              
  !                if (jbas%jj(i) .ne. jbas%jj(j)) cycle 
  !                if (jbas%jj(k) .ne. jbas%jj(l)) cycle 
  !                if (jbas%itzp(i) + jbas%itzp(j) .ne. jbas%itzp(k) + jbas%itzp(l) ) cycle
  !                if (mod(jbas%ll(i)+jbas%ll(j),2) .ne. mod(jbas%ll(k)+jbas%ll(l),2)) cycle
  !                write(72,'(4(I4),e14.6)') i,j,k,l,v_elem(i,j,k,l,0,DH,jbas)
              
  !    end do ; end do ; end do; end do
  ! ! ! ! ! ME is -.141992 should be positive.
  !    do i = 1, HS%nsp
  !       do j = i, HS%nsp
  !          if (jbas%jj(i) .ne. jbas%jj(j) ) cycle
  !          if (jbas%ll(i) .ne. jbas%ll(j) ) cycle
  !          if (jbas%itzp(i) .ne. jbas%itzp(j) ) cycle
  !          write(73,'(2(I4),e14.6)') i,j,f_elem(i,j,DH,jbas) 
  !       end do 
  !    end do 
    
  !    print*, DH%E0
  !    stop
  !print*, 'moo'
  call vectorize(DH,yp)

end subroutine 
