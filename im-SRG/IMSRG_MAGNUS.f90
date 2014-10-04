module IMSRG_MAGNUS
  use commutators
  use generators
  use basic_IMSRG
  implicit none 
  
contains

subroutine magnus_decouple( HS , jbas)!, deriv_calculator) 
  ! runs IMSRG using magnus expansion method
  implicit none 
  
  integer :: Atot,Ntot,nh,np,nb,q,steps
  type(spd) :: jbas
  type(sq_op) :: H , G ,ETA, HS,INT1,INT2,AD,w1,w2,DG,G0,ETA0,H0
  type(cross_coupled_31_mat) :: GCC,ADCC,WCC 
  real(8) :: ds,s,E_old,crit,nrm1,nrm2
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
  
  call build_gs_white(HS,ETA,jbas) 
  call copy_sq_op(HS,H) 
  
  nrm1 = mat_frob_norm(ETA) 
  s = 0.d0 
  ds = .7d0
  crit = 10.
  steps = 0
  
  open(unit=36,file='../../output/'//&
       trim(adjustl(prefix))//'_magnus0bflow.dat')
  write(36,'(I6,3(e14.6))') steps,s,H%E0,crit
  
  do while (crit > 1e-4) 
     
     call copy_sq_op(G,G0) 
     call MAGNUS_EXPAND(DG,G,ETA,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,.true.)
     call euler_step(G,DG,s,ds) 
     
     call copy_sq_op(HS,H0) 
     call BCH_EXPAND(HS,G,H,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas) 
     
     call copy_sq_op(ETA,ETA0) 
     call build_gs_white(HS,ETA,jbas) 
     nrm2 = mat_frob_norm(ETA)
     
     if ( nrm1 < nrm2 )  then 
        s = s-ds
        call copy_sq_op(G0,G)
        call copy_sq_op(ETA0,ETA)
        call copy_sq_op(H0,HS)
        ds = ds/2.d0 
        cycle 
     end if 
     
     crit = mat_frob_norm(ETA) 
     steps = steps + 1
     
     write(36,'(I6,3(e14.6))') steps,s,HS%E0,crit
  
  end do
  close(36)
end subroutine  
!=========================================================================
!=========================================================================
subroutine BCH_EXPAND(HS,G,H,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas) 
  implicit none 
  
  real(8), parameter :: conv = 1e-4 
  integer :: trunc,i,m,n
  type(spd) :: jbas
  type(sq_op) :: H , G, ETA, INT1, INT2, HS, AD,w1,w2
  type(cross_coupled_31_mat) :: WCC,ADCC,GCC
  real(8) ::  cof(11)
 
  cof = (/1.d0,0.5d0,0.166666666666666666d0, &
       0.04166666666666666d0,0.0083333333333333333d0,&
       .001388888888888d0,1.984126984d-4,2.48015873d-5,&
       2.755731922d-6,2.755731922d-7,2.505210839d-8/) 

  ! intermediates must be HERMITIAN
  INT2%herm = 1
  INT1%herm = 1 
  AD%herm = 1
  
  !! so here:  H is the current hamiltonian and we 
  !! copy it onto HS.  We copy this onto 
  !! INT2, just to make things easier for everyone.

  call copy_sq_op( H , HS )  !basic_IMSRG
  call copy_sq_op( HS , INT2 )
 
  do i = 2 , 12
          
     ! current value of HS is renamed INT1 
     ! INT2 is renamed AD, for the AD parameters in BCH and magnus expansions
     ! AD refers AD_n and INT2 becomes AD_{n+1} . 
     call copy_sq_op( HS , INT1) 
     call copy_sq_op( INT2 , AD ) 
     ! so to start, AD is equal to H
  
     !now: INT2 = [ G , AD ]  

! zero body commutator
!     call set_to_zero(INT2)
 
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
     call add_sq_op(INT1 , 1.d0 , INT2 , cof(i-1) , HS )   !basic_IMSRG
     if ( abs(mat_frob_norm(INT2)*cof(i-1)/mat_frob_norm(INT1)) < conv ) exit
     
  end do 
end subroutine 
!===============================================================
!===============================================================
subroutine MAGNUS_EXPAND(DG,G,ETA,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,ZCOM)
  implicit none 
  
  real(8), parameter :: conv = 1e-4
  integer :: trunc,i
  type(spd) :: jbas
  type(sq_op) :: H , G, ETA, INT1, INT2, HS, AD,w1,w2,DG
  type(cross_coupled_31_mat) :: WCC,ADCC,GCC
  real(8) ::  cof(4)
  logical ::  ZCOM 
  ! if ZCOM is set to true, only uses d(OMEGA) = ETA 
 
  ! Intermediates are ANTI-HERMITIAN 
  INT2%herm = -1
  INT1%herm = -1 
  AD%herm = -1
  
  cof = (/-0.5d0,.0833333333d0,0.d0,-0.00138888888d0/) 
  ! nested derivatives not so important
    
  !! same deal as BCH expansion, which is explained ad nauseam above. 
  call copy_sq_op( ETA, DG )  !ME_general
  call copy_sq_op( DG , INT2 )
  
  if (ZCOM) return 
  do i = 2 , 5
     !exit
     call copy_sq_op( DG , INT1) 
     call copy_sq_op( INT2 , AD ) 
   

     call calculate_cross_coupled(AD,ADCC,jbas,.true.)
     call calculate_cross_coupled(G,GCC,jbas,.false.) 
 
     call commutator_111(G,AD,INT2,jbas) 
     call commutator_121(G,AD,INT2,jbas)
     call commutator_122(G,AD,INT2,jbas)    

     call commutator_222_pp_hh(G,AD,INT2,w1,w2,jbas)
  
     call commutator_221(G,AD,INT2,w1,w2,jbas)
     call commutator_222_ph(GCC,ADCC,INT2,WCC,jbas)

     call add_sq_op(INT1 , 1.d0 , INT2 , cof(i-1) , DG ) !ME_general

     if ( abs(mat_frob_norm(INT2)*cof(i-1) / mat_frob_norm(INT1)) &
          < conv ) exit

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
!================================================
!================================================
end module
