module IMSRG_MAGNUS
  use commutators
  use operators
  use generators
  use basic_IMSRG
  use HF_mod 
  implicit none 
  
contains

subroutine magnus_decouple(HS,jbas,O1,O2) 
  ! runs IMSRG using magnus expansion method
  implicit none 
  
  integer :: Atot,Ntot,nh,np,nb,q,steps,i,j
  type(spd) :: jbas
  type(sq_op),optional :: O1,O2
  type(sq_op) :: H , G ,ETA, HS,INT1,INT2,AD,w1,w2,DG,G0,ETA0,H0,Oevolv
  type(cross_coupled_31_mat) :: GCC,ADCC,WCC 
  real(8) :: ds,s,E_old,E_mbpt2,crit,nrm1,nrm2,wTs(2),Ecm(3)
  character(200) :: spfile,intfile,prefix
  common /files/ spfile,intfile,prefix
  
  call duplicate_sq_op(HS,ETA) !generator
  call duplicate_sq_op(HS,H) !evolved hamiltonian
  call duplicate_sq_op(HS,w1) !workspace
  call duplicate_sq_op(HS,w2) !workspace
  call duplicate_sq_op(HS,INT1) !workspace
  call duplicate_sq_op(HS,INT2) !workspace
  call duplicate_sq_op(HS,AD) !workspace
  call duplicate_sq_op(HS,G) !magnus operator
  call duplicate_sq_op(HS,DG) !magnus operator
  call duplicate_sq_op(HS,G0) ! backup copy of G
  call duplicate_sq_op(HS,ETA0) ! backup copy of G
  call duplicate_sq_op(HS,H0) ! backup copy of G
  G%herm = -1 
  DG%herm = -1
  G0%herm = -1 
  ETA0%herm = -1
  
  call allocate_CCMAT(HS,ADCC,jbas) !cross coupled ME
  call duplicate_CCMAT(ADCC,GCC) !cross coupled ME
  call allocate_CC_wkspc(ADCC,WCC) ! workspace for CCME
  
  !call build_gs_wegner(HS,ETA,jbas,ADCC,GCC,WCC,w1,w2) 
  call build_gs_white(HS,ETA,jbas) 
  call copy_sq_op(HS,H) 
  
  s = 0.d0 
  ds = 1.0d0
  crit = 10.
  steps = 0

  open(unit=36,file='../../output/'//&
       trim(adjustl(prefix))//'_0b_magnus_flow.dat')

  E_mbpt2 = mbpt2(HS,jbas) 
  crit=abs(E_mbpt2)

  write(36,'(I6,4(e15.7))') steps,s,H%E0,HS%E0+E_mbpt2,crit
  write(*,'(I6,4(e15.7))') steps,s,HS%E0,HS%E0+E_mbpt2,crit

  nrm1 = mat_frob_norm(ETA)
  do while (crit > 1e-6) 
     
     call copy_sq_op(G,G0) 
     
     call MAGNUS_EXPAND(DG,G,ETA,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s)
     call euler_step(G,DG,s,ds) 
  
     call copy_sq_op(HS,H0) 
     call BCH_EXPAND(HS,G,H,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s) 

     call copy_sq_op(ETA,ETA0)
 !   call build_gs_wegner(HS,ETA,jbas,ADCC,GCC,WCC,w1,w2)  
     call build_gs_white(HS,ETA,jbas) 

     nrm2 = HS%E0 !mat_frob_norm(ETA)
     nrm2 = mat_frob_norm(ETA)
     if ( nrm1 < nrm2 )  then
        s = s-ds
        call copy_sq_op(G0,G)
        call copy_sq_op(ETA0,ETA)
        call copy_sq_op(H0,HS)
        ds = ds/2.d0 
        cycle 
     end if 

     E_mbpt2 = mbpt2(HS,jbas) 

     crit = abs(E_mbpt2) 
     nrm1 = nrm2 
     steps = steps + 1
     write(36,'(I6,4(e15.7))') steps,s,HS%E0,HS%E0+E_mbpt2,crit
     write(*,'(I6,4(e15.7))') steps,s,HS%E0,HS%E0+E_mbpt2,crit
  end do

! calculate any observables which have been requested =====================

  if (present(O1)) then 
     call duplicate_sq_op(O1,Oevolv)
    
     if (present(O2)) then 

       call BCH_EXPAND(Oevolv,G,O2,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s) 
       
       call copy_sq_op(Oevolv,O2) 

     end if 
 
     call clear_sq_op(Oevolv)
     call BCH_EXPAND(Oevolv,G,O1,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s) 
     
     call copy_sq_op(Oevolv,O1) 
          
  end if

!===========================================================================  
  close(36)
  
end subroutine  
!===========================================================================
!===========================================================================
subroutine magnus_TDA(HS,TDA,jbas,O1,O1TDA,O2,O2TDA) 
  ! runs IMSRG TDA decoupling using magnus expansion method
  implicit none 
  
  integer :: Atot,Ntot,nh,np,nb,q,steps,i
  type(spd) :: jbas
  type(sq_op),optional :: O1,O2
  type(full_sp_block_mat) :: TDA
  type(full_sp_block_mat),optional :: O1TDA,O2TDA
  type(sq_op) :: H , G ,ETA, HS,INT1,INT2,AD,w1,w2,DG,G0,ETA0,H0,Oevolv
  type(cross_coupled_31_mat) :: GCC,ADCC,WCC,HCC,OeCC
  real(8) :: ds,s,crit,nrm1,nrm2,wTs(2),Ecm(3)
  real(8),allocatable,dimension(:) :: E_old,wTvec
  character(200) :: spfile,intfile,prefix
  character(3) :: args
  common /files/ spfile,intfile,prefix
  
  call duplicate_sq_op(HS,ETA) !generator
  call duplicate_sq_op(HS,H) !evolved hamiltonian
  call duplicate_sq_op(HS,w1) !workspace
  call duplicate_sq_op(HS,w2) !workspace
  call duplicate_sq_op(HS,INT1) !workspace
  call duplicate_sq_op(HS,INT2) !workspace
  call duplicate_sq_op(HS,AD) !workspace
  call duplicate_sq_op(HS,G) !magnus operator
  call duplicate_sq_op(HS,DG) !magnus operator
  call duplicate_sq_op(HS,G0) ! backup copy of G
  call duplicate_sq_op(HS,ETA0) ! backup copy of G
  call duplicate_sq_op(HS,H0) ! backup copy of G

  G%herm = -1 
  DG%herm = -1
  G0%herm = -1 
  ETA0%herm = -1
    
  call allocate_CCMAT(HS,ADCC,jbas) !cross coupled ME
  call duplicate_CCMAT(ADCC,GCC) !cross coupled ME
  call duplicate_CCMAT(ADCC,HCC) 
  call allocate_CC_wkspc(ADCC,WCC) ! workspace for CCME
  
  ! TDA stuff
  call calculate_cross_coupled(HS,HCC,jbas,.true.)
  call calc_TDA(TDA,HS,HCC,jbas)
  call diagonalize_blocks(TDA)
  allocate(E_old(TDA%map(1)))

  s = 0.d0 
  ds = 0.001d0
  crit = 10.
  steps = 0

  E_old = TDA%blkM(1)%eigval
  
  open(unit=37,file='../../output/'//&
       trim(adjustl(prefix))//'_excited.dat')
  
  call write_excited_states(steps,s,TDA,HS%E0,37)
  call build_specific_space(HS,ETA,jbas) 
  call copy_sq_op(HS,H) 
  
  nrm1 = mat_frob_norm(ETA) 
 
  do while (crit > 1e-5) 
     
     call copy_sq_op(G,G0) 
     call MAGNUS_EXPAND(DG,G,ETA,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s)
     call euler_step(G,DG,s,ds) 
 
     call copy_sq_op(HS,H0) 
     call BCH_EXPAND(HS,G,H,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s) 
     
     call copy_sq_op(ETA,ETA0)
    
     call build_specific_space(HS,ETA,jbas)
   
     nrm2 = mat_frob_norm(ETA)
   
     if ( ( 10*nrm1 < nrm2 ) .and. (ds > 1e-2))  then
        s = s-ds
        call copy_sq_op(G0,G)
        call copy_sq_op(ETA0,ETA)
        call copy_sq_op(H0,HS)
        ds = ds/2.d0 
        cycle 
     end if 
    
     call calculate_cross_coupled(HS,HCC,jbas,.true.) 
     call calc_TDA(TDA,HS,HCC,jbas) 
     call diagonalize_blocks(TDA)
  
     call write_excited_states(steps,s,TDA,HS%E0,37) 
     
     crit = sum(abs(E_old-TDA%blkM(1)%eigval))/TDA%map(1)
     write(*,'(I6,6(e15.7))') steps,s,TDA%blkM(1)%eigval(1:4),crit
     E_old = TDA%blkM(1)%eigval

     nrm1 = nrm2 
     steps = steps + 1
     
     
  end do
 
! calculate any observables which have been requested =====================
  
  if (present(O1)) then 
     call duplicate_sq_op(O1,Oevolv) 
     call duplicate_CCMAT(HCC,OeCC)

     if (present(O2)) then 
        
        ! transform observable
        call BCH_EXPAND(Oevolv,G,O2,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s) 
        call duplicate_sp_mat(TDA,O2TDA) 
        allocate(O2TDA%blkM(1)%labels(TDA%map(1),2)) 
        O2TDA%blkM(1)%labels = TDA%blkM(1)%labels      
        call calculate_cross_coupled(Oevolv,OeCC,jbas,.true.) 
        call calc_TDA(O2TDA,Oevolv,OeCC,jbas)
        call copy_sq_op(Oevolv,O2)
     end if 
     
     ! transform observable
     call BCH_EXPAND(Oevolv,G,O1,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s) 
     call duplicate_sp_mat(TDA,O1TDA)
     allocate(O1TDA%blkM(1)%labels(TDA%map(1),2)) 
     O1TDA%blkM(1)%labels = TDA%blkM(1)%labels      
     call calculate_cross_coupled(Oevolv,OeCC,jbas,.true.) 
     call calc_TDA(O1TDA,Oevolv,OeCC,jbas) 
     call copy_sq_op(Oevolv,O1) 
  end if 

!===========================================================================  
  close(37)
end subroutine
!=========================================================================
!=========================================================================
subroutine BCH_EXPAND(HS,G,H,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s) 
  implicit none 
  
  real(8), parameter :: conv = 1e-8
  integer :: trunc,i,m,n,q,j,k,l
  type(spd) :: jbas
  type(sq_op) :: H , G, ETA, INT1, INT2, HS, AD,w1,w2
  type(cross_coupled_31_mat) :: WCC,ADCC,GCC
  real(8) ::  cof(15),adnorm,fullnorm,s,advals(15)
  character(3) :: args
  
  advals = 0.d0 
  
  cof = (/1.d0,1.d0,0.5d0,0.166666666666666666d0, &
       0.04166666666666666d0,0.0083333333333333333d0,&
       .001388888888888d0,1.984126984d-4,2.48015873016d-5,&
       2.75573192239d-6,2.75573192239d-7,2.505210839d-8, &
       2.087675698d-9,1.6059043837d-10,1.1470745598d-11/) 

  ! intermediates must be HERMITIAN
  INT2%herm = 1
  INT1%herm = 1 
  AD%herm = 1
  
  !! so here:  H is the current hamiltonian and we 
  !! copy it onto HS.  We copy this onto 
  !! INT2, just to make things easier for everyone.

  call copy_sq_op( H , HS )  !basic_IMSRG
  call copy_sq_op( HS , INT2 )
 
  advals(1) = abs(H%E0)   
 
  do i = 2 ,15

     ! current value of HS is renamed INT1 
     ! INT2 is renamed AD, for the AD parameters in BCH and magnus expansions
     ! AD refers AD_n and INT2 becomes AD_{n+1} . 
     call copy_sq_op( HS , INT1) 
     call copy_sq_op( INT2 , AD ) 
     ! so to start, AD is equal to H
     call clear_sq_op(INT2)    
     !now: INT2 = [ G , AD ]  
        
! zero body commutator
 
     call calculate_cross_coupled(AD,ADCC,jbas,.true.)
     call calculate_cross_coupled(G,GCC,jbas,.false.) 
 
     INT2%E0 = commutator_110(G,AD,jbas) + commutator_220(G,AD,jbas)

     call commutator_111(G,AD,INT2,jbas) 
     call commutator_121(G,AD,INT2,jbas)
     call commutator_122(G,AD,INT2,jbas)    

     call commutator_222_pp_hh(G,AD,INT2,w1,w2,jbas)
  
     call commutator_221(G,AD,INT2,w1,w2,jbas)
     call commutator_222_ph(GCC,ADCC,INT2,WCC,jbas)
     
     ! so now just add INT1 + c_n * INT2 to get current value of HS
     
     call add_sq_op(INT1 , 1.d0 , INT2 , cof(i) , HS )   !basic_IMSRG
    
     advals(i) = abs(INT2%E0*cof(i))
     if (advals(i) < conv) exit
     
  end do 
 
end subroutine 
!===============================================================
!===============================================================
subroutine MAGNUS_EXPAND(DG,G,ETA,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s)
  implicit none 
  
  real(8), parameter :: conv = 1e-8
  integer :: trunc,i,q,j,k,l,ry
  type(spd) :: jbas
  type(sq_op) :: H , G, ETA, INT1, INT2, HS, AD,w1,w2,DG
  type(cross_coupled_31_mat) :: WCC,ADCC,GCC
  real(8) ::  cof(7),adnorm,fullnorm,s,advals(7) 
  character(3) :: args
  
  advals = 0.d0 
  ! Intermediates are ANTI-HERMITIAN 
  INT2%herm = -1
  INT1%herm = -1 
  AD%herm = -1
  
  cof = (/1.d0,-0.5d0,.0833333333d0,0.d0,-0.00138888888d0,0.d0,3.306878d-5/) 
  ! nested derivatives not so important
    
  !! same deal as BCH expansion, which is explained ad nauseam above. 
  call copy_sq_op( ETA, DG )  !ME_general
  call copy_sq_op( DG , INT2 )
  advals(1) = mat_frob_norm(INT2)  
  
  fullnorm = mat_frob_norm(G)   
  if (fullnorm < 1e-9) return
  
  q = 1
 ! return
  do i = 2 , 7

     call copy_sq_op( DG , INT1) 
     call copy_sq_op( INT2 , AD ) 
  
     adnorm = advals(i-q) 
  
     if  (abs(cof(i)) > 1e-6) then  
        q = 1 
     else 
        q = 2
     end if 
     
     if ( abs(adnorm/fullnorm) < conv ) exit
          
     call calculate_cross_coupled(AD,ADCC,jbas,.true.)
     call calculate_cross_coupled(G,GCC,jbas,.false.) 
 
     call commutator_111(G,AD,INT2,jbas) 
     call commutator_121(G,AD,INT2,jbas)
     call commutator_122(G,AD,INT2,jbas)    
  
     call commutator_222_pp_hh(G,AD,INT2,w1,w2,jbas)   
     call commutator_221(G,AD,INT2,w1,w2,jbas)     
     call commutator_222_ph(GCC,ADCC,INT2,WCC,jbas)
                 
     call add_sq_op(INT1 , 1.d0 , INT2 , cof(i) , DG ) !ME_general
     
     advals(i) = mat_frob_norm(INT2)*abs(cof(i))
    
  end do   
     
end subroutine 
!=====================================================
!=====================================================
subroutine euler_step(G,DG,s,stp)
  implicit none 

  integer :: i
  real(8) :: s , stp
  type(sq_op) :: G , DG

  call add_sq_op(G,1.d0,DG,stp,G) !durr probably wrong. 
  s = s + stp 

end subroutine 
end module
!================================================
!================================================
  
  
  
