module IMSRG_MAGNUS
  use commutators
  use operators
  use generators
  use basic_IMSRG
  use HF_mod 
  use mscheme
  implicit none 
  
contains

subroutine magnus_decouple(HS,jbas,O1,O2,quads,trips) 
  ! runs IMSRG using magnus expansion method
  implicit none 
  
  integer :: Atot,Ntot,nh,np,nb,q,steps,i,j
  type(spd) :: jbas
  type(tpd),allocatable,dimension(:) :: threebas
  type(mscheme_3body) :: threebd
  type(sq_op),optional :: O1,O2
  type(sq_op) :: H,H2,G1b,G2b,G,ETA,HS,INT1,INT2,AD,w1,w2,DG,G0,ETA0,H0,Oevolv,CR
  type(cross_coupled_31_mat) :: GCC,ADCC,WCC 
  real(8) :: ds,s,E_old,E_mbpt2,crit,nrm1,nrm2,wTs(2),Ecm(3),corr,dcgi00,xxx
  real(8) :: omp_get_wtime,t1,t2
  character(200) :: spfile,intfile,prefix
  character(1),optional :: quads,trips
  logical :: qd_calc,trip_calc,xxCR
  common /files/ spfile,intfile,prefix
  
  qd_calc = .false. 
  if (present(quads)) then 
     qd_calc = .true.
  end if 
  
  trip_calc = .false. 
  if (present(trips)) then 
     trip_calc=.true.
     xxCr = .false. 
     if (trips == 'C') xxCr = .true.
  end if 
     

  HS%neq = 1
!  call duplicate_sq_op(HS,ETA) !generator
  call duplicate_sq_op(HS,H) !evolved hamiltonian
!  call duplicate_sq_op(HS,H2) !evolved hamiltonian
  call duplicate_sq_op(HS,w1) !workspace
  call duplicate_sq_op(HS,w2) !workspace
  call duplicate_sq_op(HS,INT1) !workspace
  call duplicate_sq_op(HS,INT2) !workspace
  call duplicate_sq_op(HS,AD) !workspace
  call duplicate_sq_op(HS,G) !magnus operator
  call duplicate_sq_op(HS,DG) !magnus operator
!  call duplicate_sq_op(HS,G0) ! backup copy of G
!  call duplicate_sq_op(HS,ETA0) ! backup copy of G
!  call duplicate_sq_op(HS,H0) ! backup copy of G
!  call duplicate_sq_op(HS,G1b)
!  call duplicate_sq_op(HS,G2b) 
  G%herm = -1 
  DG%herm = -1
  G0%herm = -1 
  ETA0%herm = -1
  
  call allocate_CCMAT(HS,ADCC,jbas) !cross coupled ME
  call duplicate_CCMAT(ADCC,GCC) !cross coupled ME
  call allocate_CC_wkspc(ADCC,WCC) ! workspace for CCME
  
  !call build_gs_wegner(HS,ETA,jbas,ADCC,GCC,WCC,w1,w2) 
 
  call build_gs_white(HS,DG,jbas)  
  !call build_gs_imtime(HS,DG,jbas) 
  
  call copy_sq_op(HS,H) 
  
  s = 0.d0 
    
  if (HS%lawson_beta < 3.0) then 
     ds = 1.0d0
  else if (HS%lawson_beta < 6.0) then 
     ds = 0.5d0
  else
     ds = 0.1d0
  end if 
  
  crit = 10.
  steps = 0

  open(unit=36,file='../../output/'//&
       trim(adjustl(prefix))//'_0b_magnus_flow.dat')

  E_mbpt2 = mbpt2(HS,jbas) 
  crit=abs(E_mbpt2)

  write(36,'(I6,4(e15.7))') steps,s,H%E0,HS%E0+E_mbpt2,crit
  write(*,'(I6,4(e15.7))') steps,s,HS%E0,HS%E0+E_mbpt2,crit

  nrm1 = mat_frob_norm(DG)
  
  do while (crit > 1e-6) 
     
    ! call copy_sq_op(G,G0) 
    ! call copy_sq_op(DG,ETA0)

     call MAGNUS_EXPAND(DG,G,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s)
     
     call euler_step(G,DG,s,ds) 
     !call copy_sq_op(HS,H0)
     
     if (qd_calc) then 
        ! calculate quadrupoles correction
       ! call split_1b_2b(G,G1b,G2b) 
       ! call BCH_EXPAND(H2,G1b,H,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s) 
       ! call BCH_EXPAND(HS,G2b,H2,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s,'y')   
        call BCH_EXPAND(HS,G,H,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s,'y') 
     else
        ! make two separate operators
        !call split_1b_2b(G,G1b,G2b) 
        !call BCH_EXPAND(H2,G1b,H,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s) 
        !call BCH_EXPAND(HS,G2b,H2,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s) 
        call BCH_EXPAND(HS,G,H,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s) 
     end if
     
     !   call build_gs_wegner(HS,ETA,jbas,ADCC,GCC,WCC,w1,w2)  
  
     call build_gs_white(HS,DG,jbas)   
     
     nrm2 = HS%E0 !mat_frob_norm(ETA)
     nrm2 = mat_frob_norm(DG)
     ! if ( nrm1 < nrm2 )  then
     !    s = s-ds
     !    call copy_sq_op(G0,G)
     !    call copy_sq_op(ETA0,DG)
     !    call copy_sq_op(H0,HS)
     !    ds = ds/2.d0 
     !    cycle 
     ! end if 

     E_mbpt2 = mbpt2(HS,jbas) 

     crit = abs(E_mbpt2) 
     nrm1 = nrm2 
     steps = steps + 1
     write(36,'(I6,4(e15.7))') steps,s,HS%E0,HS%E0+E_mbpt2,crit
     write(*,'(I6,4(e15.7))') steps,s,HS%E0,HS%E0+E_mbpt2,crit

 end do

!  call split_1b_2b(G,G1b,G2b) 
!  call BCH_EXPAND(H2,G1b,H,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s) 
!  call BCH_EXPAND(HS,G2b,H2,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s) 
 
! calculate any observables which have been requested =====================

  ! triples correction
  if (trip_calc) then 
     call enumerate_three_body(threebas,jbas)
     t1 = omp_get_wtime()
     if ( xxCR ) then 
        call duplicate_sq_op(H,CR)         
        call CR_EXPAND(CR,G,H,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s) 
        corr =  restore_triples(CR,G,threebas,jbas)
     else
        corr =  restore_triples(H,G,threebas,jbas) 
     end if 
     t2 = omp_get_wtime()
     print*, 'FINAL ENERGY:', corr + HS%E0,t2-t1
     open(unit=39,file='../../output/'//&
       trim(adjustl(prefix))//'_magnus_triples.dat')
     write(39,'(I6,4(e15.7))') steps,s,HS%E0+corr,HS%E0+E_mbpt2,crit
     close(39)
  end if
  
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
subroutine magnus_TDA(HS,TDA,jbas,O1,O1TDA,O2,O2TDA,quads) 
  ! runs IMSRG TDA decoupling using magnus expansion method
  implicit none 
  
  integer :: Atot,Ntot,nh,np,nb,q,steps,i,Jsing
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
  character(1),optional :: quads 
  logical :: qd_calc 
  character(1) :: Jlabel,Plabel
  common /files/ spfile,intfile,prefix
  
  qd_calc = .false.
  if (present(quads)) qd_calc =.true. 
  
!  call duplicate_sq_op(HS,ETA) !generator
  call duplicate_sq_op(HS,H) !evolved hamiltonian
  call duplicate_sq_op(HS,w1) !workspace
  call duplicate_sq_op(HS,w2) !workspace
  call duplicate_sq_op(HS,INT1) !workspace
  call duplicate_sq_op(HS,INT2) !workspace
  call duplicate_sq_op(HS,AD) !workspace
  call duplicate_sq_op(HS,G) !magnus operator
  call duplicate_sq_op(HS,DG) !magnus operator
!  call duplicate_sq_op(HS,G0) ! backup copy of G
!  call duplicate_sq_op(HS,ETA0) ! backup copy of G
!  call duplicate_sq_op(HS,H0) ! backup copy of G

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
  ds = 0.005d0
  crit = 10.
  steps = 0

  E_old = TDA%blkM(1)%eigval
  Jsing = H%Jtarg/2

  write( Jlabel ,'(I1)') Jsing

  if (H%Ptarg == 0 ) then 
     Plabel ='+'
  else 
     Plabel ='-'
  end if
  open(unit=37,file='../../output/'//&
       trim(adjustl(prefix))//'_'//Jlabel//Plabel//'_excited.dat')
  
  call write_excited_states(steps,s,TDA,HS%E0,37)
  call build_specific_space(HS,DG,jbas) 
  call copy_sq_op(HS,H) 
  
  !nrm1 = mat_frob_norm(DG) 
 
  do while (crit > 1e-5) 
     
   !  call copy_sq_op(G,G0) 
   !  call copy_sq_op(DG,ETA0)
  
     call MAGNUS_EXPAND(DG,G,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s)
     call euler_step(G,DG,s,ds) 
 
   !  call copy_sq_op(HS,H0) 
     if (qd_calc) then 
        call BCH_EXPAND(HS,G,H,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s,'y') 
     else
        call BCH_EXPAND(HS,G,H,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s) 
     end if
    
         
     call build_specific_space(HS,DG,jbas)
   
 !    nrm2 = mat_frob_norm(DG)
   
     ! if ( ( 10*nrm1 < nrm2 ) .and. (ds > 1e-2))  then
     !    s = s-ds
     !    print*, 'damnit'
     !    call copy_sq_op(G0,G)
     !    call copy_sq_op(ETA0,DG)
     !    call copy_sq_op(H0,HS)
     !    ds = ds/2.d0 
     !    cycle 
     ! end if 
    
     call calculate_cross_coupled(HS,HCC,jbas,.true.) 
     call calc_TDA(TDA,HS,HCC,jbas) 
     call diagonalize_blocks(TDA)
  
     call write_excited_states(steps,s,TDA,HS%E0,37) 
     
     crit = sum(abs(E_old-TDA%blkM(1)%eigval))/TDA%map(1)
     write(*,'(I6,7(e15.7))') steps,s,TDA%blkM(1)%eigval(1:5),crit
     E_old = TDA%blkM(1)%eigval

    ! nrm1 = nrm2 
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
subroutine BCH_EXPAND(HS,G,H,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s,quads) 
  implicit none 
  
  real(8), parameter :: conv = 1e-8
  integer :: trunc,i,m,n,q,j,k,l,a,b,c,d,iw
  integer :: ix,jx,kx,lx,ax,cx,bx,dx,jmin,jmax,Jtot
  integer :: mi,mj,mk,ml,ma,mc,mb,md,ja,jb,jj,ji,JT,MT
  type(spd) :: jbas
  type(sq_op) :: H , G, ETA, INT1, INT2, INT3,HS, AD,w1,w2
  type(cross_coupled_31_mat) :: WCC,ADCC,GCC
  real(8) ::  cof(15),adnorm,fullnorm,s,advals(15),sm,sm2,dcgi,dcgi00
  character(3) :: args
  character(1),optional :: quads ! enter some character to restore quadrupoles 
  logical :: qd_calc   

  qd_calc = .false. 
  if (present(quads)) then 
     qd_calc = .true. 
  end if 
  
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
  call duplicate_sq_op(H,INT3)

  call copy_sq_op( H , HS )  !basic_IMSRG
  call copy_sq_op( HS , INT2 )
 
  advals(1) = abs(H%E0)   

  do iw = 2 ,15

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
          
     call add_sq_op(INT3, 1.d0 , INT2, 1.d0, INT2)     
     call add_sq_op(INT1 ,1.d0 , INT2 , cof(iw) , HS )   !basic_IMSRG
    
     if (qd_calc) then 
        call restore_quadrupoles(AD,G,w1,w2,INT3,jbas) 
     end if 
     
     advals(iw) = abs(INT2%E0*cof(iw))
     if (advals(iw) < conv) exit
     
  end do 
 
end subroutine 
!=========================================================================
!=========================================================================
subroutine CR_EXPAND(HS,G,H,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s,quads) 
  implicit none 
  
  real(8), parameter :: conv = 1e-8
  integer :: trunc,i,m,n,q,j,k,l,a,b,c,d,iw
  integer :: ix,jx,kx,lx,ax,cx,bx,dx,jmin,jmax,Jtot
  integer :: mi,mj,mk,ml,ma,mc,mb,md,ja,jb,jj,ji,JT,MT
  type(spd) :: jbas
  type(sq_op) :: H , G, ETA, INT1, INT2, INT3,HS, AD,w1,w2
  type(cross_coupled_31_mat) :: WCC,ADCC,GCC
  real(8) ::  cof(14),adnorm,fullnorm,s,advals(14),sm,sm2,dcgi,dcgi00
  character(3) :: args
  character(1),optional :: quads ! enter some character to restore quadrupoles 
  logical :: qd_calc   

  qd_calc = .false. 
  if (present(quads)) then 
     qd_calc = .true. 
  end if 
  
  advals = 0.d0 
  
  cof = (/1.d0,0.5d0,0.166666666666666666d0, &
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
  call duplicate_sq_op(H,INT3)

  call copy_sq_op( H , HS )  !basic_IMSRG
  call copy_sq_op( HS , INT2 )
 
  advals(1) = abs(H%E0)   

  do iw = 2 ,14

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
          
     call add_sq_op(INT3, 1.d0 , INT2, 1.d0, INT2)     
     call add_sq_op(INT1 ,1.d0 , INT2 , cof(iw) , HS )   !basic_IMSRG
    
     if (qd_calc) then 
        call restore_quadrupoles(AD,G,w1,w2,INT3,jbas) 
     end if 
     
     advals(iw) = abs(INT2%E0*cof(iw))
     if (advals(iw) < conv) exit
     
  end do 
 
end subroutine 
!===============================================================
!===============================================================
subroutine MAGNUS_EXPAND(DG,G,INT1,INT2,AD,w1,w2,ADCC,GCC,WCC,jbas,s)
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
!=====================================================
!=====================================================
real(8) function restore_triples(H,OM,threebas,jbas) 
  implicit none 
  
  type(spd) :: jbas
  type(tpd),dimension(:) :: threebas
  type(sq_op) :: H,OM
  integer :: a,b,c,i,j,k,Jtot,Jab,Jij,g1
  integer :: ja,jb,jc,ji,jj,jk,AAA,q
  integer :: ax,bx,cx,ix,jx,kx,III
  integer :: jab_min,jab_max,jij_min,jij_max
  integer :: J_min, J_max,x,total_threads,thread
  real(8) :: sm,denom,dlow
  
  sm = 0.d0   
  total_threads = size(threebas(1)%direct_omp) - 1

!$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) SHARED(threebas,jbas,H,OM) & 
!$OMP& REDUCTION(+:sm)  
  
  do thread = 1, total_threads
  do q = 1+threebas(1)%direct_omp(thread),&
       threebas(1)%direct_omp(thread+1)
  
     Jtot = threebas(q)%chan(1) 
     do AAA = 1, size(threebas(q)%ppp(:,1)) 
        
        a = threebas(q)%ppp(AAA,1)
        b = threebas(q)%ppp(AAA,2)
        c = threebas(q)%ppp(AAA,3)

        ja = jbas%jj(a)      
        jb = jbas%jj(b) 
        jc = jbas%jj(c) 
        
        jab_min = abs(ja-jb) 
        jab_max = ja+jb
   
        dlow = f_elem(a,a,H,jbas)&
                +f_elem(b,b,H,jbas)+f_elem(c,c,H,jbas)
        
        do III = 1, size(threebas(q)%hhh(:,1)) 
             
           i = threebas(q)%hhh(III,1)
           j = threebas(q)%hhh(III,2)
           k = threebas(q)%hhh(III,3)

           ji = jbas%jj(i)       
           jj = jbas%jj(j) 
           jk = jbas%jj(k)  
     
           jij_min = abs(ji-jj) 
           jij_max = ji+jj
               
           denom = f_elem(i,i,H,jbas)+f_elem(j,j,H,jbas)&
                +f_elem(k,k,H,jbas)-dlow
                
           do jab = jab_min,jab_max,2
               
              if ( .not. (triangle(Jtot,jc,jab))) cycle
               
              do jij = jij_min, jij_max,2
                  
                 if ( .not. (triangle(Jtot,jk,jij))) cycle
                      
                 sm = sm + &
                      commutator_223_single(OM,H,a,b,c,i,j,k,Jtot,jab,jij,jbas)**2 &
                      /denom*(Jtot+1.d0)
           
              end do
           end do

        end do
     end do
  end do
  end do 
 !$OMP END PARALLEL DO 
  restore_triples = sm / 36.d0 

end function
!================================================
!================================================
end module
!================================================
!================================================
subroutine restore_quadrupoles( X , OM, w1,w2, RES,jbas ) 
  use basic_IMSRG
  implicit none
  
  type(sq_op) :: X,OM,INT1,INT2,RES
  type(spd) :: jbas
  type(sq_op) :: w1,w2
  integer :: p1x,p2x,h1x,h2x
  integer :: i,j,b,q,Abody,Ntot,nh,np,nb,a,c,p1,p2,h1,h2,cx
  integer :: AA,II,Jtot,ik,jk,ck,ji,jj,ti,tj,li,lj,jc,JT,pm,jp1,jp2,ja,jb
  integer :: kx,dx,k,d
  real(8) :: sm,sm2
  
  Abody = X%belowEF
  Ntot = X%Nsp
  
  ! I think this is all I need here (these are the one -body insertions)
  allocate(INT1%fpp(Ntot-Abody,Ntot-Abody),INT1%fhh(Abody,Abody)) 
  INT1%fpp = 0.d0
  INT1%fhh = 0.d0 
 
  allocate(INT2%fpp(Ntot-Abody,Ntot-Abody),INT2%fhh(Abody,Abody)) 
  INT2%fpp = 0.d0
  INT2%fhh = 0.d0
 
  
  pm = X%herm*OM%herm

  ! w1=X.OM w2=OM.OM matrices are made 
  call build_intermediates_For_intermediates(X,OM,w1,w2,jbas)

! fhh
  do i = 1 , Abody
     ik = jbas%holes(i) 
     ji = jbas%jj(ik) 
     li = jbas%ll(ik) 
     ti = jbas%itzp(ik) 
     
     do j = 1 , Abody
        
        jk = jbas%holes(j) 
        jj = jbas%jj(jk) 
        if (ji .ne. jj)  cycle
        lj = jbas%ll(jk) 
        if (li .ne. lj)  cycle
        tj = jbas%itzp(jk)
        if (tj .ne. ti) cycle 
                
        sm = 0.d0 
        sm2 = 0.d0 
        do c = 1, Abody
           ck = jbas%holes(c) 
           jc = jbas%jj(ck)
           do JT = abs(jc - ji),jc+ji,2
              sm = sm + v_elem(ck,ik,ck,jk,JT,w1,jbas)*(JT + 1)
              sm2 = sm2 + v_elem(ck,ik,ck,jk,JT,w2,jbas)*(JT + 1) 
              ! w2 is subtracted, because the commutator in this case has a minus sign. 
           end do
        end do 
           
        INT1%fhh(i,j) = INT1%fhh(i,j) + sm / (ji + 1.d0 )
        INT2%fhh(i,j) = INT2%fhh(i,j) + sm2 / (ji + 1.d0 )
       ! nothing is hermitian or anti-hermitian here
     end do 
  end do       

! fpp
  do i = 1 , Ntot - Abody
     ik = jbas%parts(i) 
     ji = jbas%jj(ik) 
     li = jbas%ll(ik) 
     ti = jbas%itzp(ik) 
     
     do j = i , Ntot - Abody
        
        jk = jbas%parts(j) 
        jj = jbas%jj(jk) 
        if (ji .ne. jj)  cycle
        lj = jbas%ll(jk) 
        if (li .ne. lj)  cycle
        tj = jbas%itzp(jk)
        if (tj .ne. ti) cycle 
                
       
        sm = 0.d0 
        sm2 = 0.d0
        do c = 1, Ntot - Abody
           ck = jbas%parts(c) 
           jc = jbas%jj(ck)
           ! w2 matrix results from multiplying the hh channel
           do JT = abs(jc - ji),jc+ji,2
              ! the hermiticity of X and OM are exploited here. 
              sm = sm + X%herm*OM%herm*v_elem(ck,jk,ck,ik,JT,w1,jbas)*(JT + 1) 
              sm2 = sm2 + v_elem(ck,jk,ck,ik,JT,w2,jbas)*(JT + 1)
           end do 
        end do 
     
        INT1%fpp(i,j) = INT1%fpp(i,j) + sm / (ji + 1.d0) 
        INT2%fpp(i,j) = INT2%fpp(i,j) + sm2 / (ji + 1.d0) 
        
     end do 
  end do       

  !!! now add the new stuff to RES
  
  do q = 1, RES%nblocks
     
     JTot = RES%mat(q)%lam(1) 
     
     np = RES%mat(q)%npp
     nh = RES%mat(q)%nhh
     
     do AA = 1, np
        a = RES%mat(q)%qn(1)%Y(AA,1)
        b = RES%mat(q)%qn(1)%Y(AA,2)
        ja = jbas%jj(a)
        jb = jbas%jj(b)
      
        do II = 1, nh    
           i = RES%mat(q)%qn(3)%Y(II,1)
           j = RES%mat(q)%qn(3)%Y(II,2)
           ji = jbas%jj(i)
           jj = jbas%jj(j)
           
           RES%mat(q)%gam(3)%X(AA,II) = 0.d0
           do c = 1, Abody
              
              cx = jbas%holes(c) 
              h1 = jbas%holesb4(i)+1 
              h2 = jbas%holesb4(j)+1
              
              
              ! i've already built the INT2 minus signs into the one-body
              ! insertion
              RES%mat(q)%gam(3)%X(AA,II) = RES%mat(q)%gam(3)%X(AA,II) &
                   - INT1%fhh(h1,c) * v_elem(a,b,cx,j,Jtot,OM,jbas)  &
                   + (-1)**((ji+jj-Jtot)/2)*INT1%fhh(h2,c) * v_elem(a,b,cx,i,Jtot,OM,jbas) &
                   + INT2%fhh(h1,c) * v_elem(a,b,cx,j,Jtot,X,jbas) &
                   - (-1)**((ji+jj-Jtot)/2)*INT2%fhh(h2,c) * v_elem(a,b,cx,i,Jtot,X,jbas)
                   
          end do 
       
          do c = 1, Ntot-Abody
              
              cx = jbas%parts(c) 
              p1 = jbas%partsb4(a)+1 
              p2 = jbas%partsb4(b)+1
              
              
              ! i've already built the INT2 minus signs into the one-body
              ! insertion
               RES%mat(q)%gam(3)%X(AA,II) = RES%mat(q)%gam(3)%X(AA,II) &
                   - INT1%fpp(p1,c) * v_elem(cx,b,i,j,Jtot,OM,jbas) &
                    + (-1)**((ja+jb-Jtot)/2)* INT1%fpp(p2,c) * v_elem(cx,a,i,j,Jtot,OM,jbas) &
                    + INT2%fpp(p1,c) * v_elem(cx,b,i,j,Jtot,X,jbas) &
                    -(-1)**((ja+jb-Jtot)/2)*INT2%fpp(p2,c) * v_elem(cx,a,i,j,Jtot,X,jbas)
                   
          end do 
       end do 
    end do 
 end do 

end subroutine 
!========================================================================
!========================================================================
subroutine build_intermediates_For_intermediates(L,R,w1,w2,jbas) 
  !VERIFIED
  !NEEDS TO BE RUN BEFORE 221, because it sets up the 
  !intermediary matrices for the calculation of 
  ! one body insertions in restore_quadrupoles
  use basic_IMSRG
  implicit none
  
  type(spd) :: jbas
  type(sq_op) ::  L,R,w1,w2
  integer :: q
  integer :: np,nb,nh,pm
  real(8) :: bet_off,al_off
  
  pm = R%herm*L%herm
   !construct temporary matrices
  do q = 1, L%nblocks
     
     nh = L%mat(q)%nhh
     np = L%mat(q)%npp
     nb = L%mat(q)%nph
        
     if (nh + np == 0 ) cycle 
     if (np + nb == 0 ) cycle 
     if (nh + nb == 0 ) cycle
     
  
             
     if (np*nh .ne. 0) then 
     !L_hhpp . R_pphh = W1_hhhh
     call dgemm('T','N',nh,nh,np,al,L%mat(q)%gam(3)%X,np,&
          R%mat(q)%gam(3)%X,np,bet,w1%mat(q)%gam(5)%X,nh)
    
     w1%mat(q)%gam(5)%X = w1%mat(q)%gam(5)%X * L%herm 
     
     !L_pphh . R_hhpp = W1_pppp
     call dgemm('N','T',np,np,nh,al,L%mat(q)%gam(3)%X,np,&
          R%mat(q)%gam(3)%X,np,bet,w1%mat(q)%gam(1)%X,np)
     
     w1%mat(q)%gam(1)%X = w1%mat(q)%gam(1)%X * R%herm
     
     
     !R_hhpp . R_pphh = W2_hhhh
     call dgemm('T','N',nh,nh,np,al,R%mat(q)%gam(3)%X,np,&
          R%mat(q)%gam(3)%X,np,bet,w2%mat(q)%gam(5)%X,nh)
    
     w2%mat(q)%gam(5)%X = w2%mat(q)%gam(5)%X * R%herm 
     
     !R_pphh . R_hhpp = W2_pppp
     call dgemm('N','T',np,np,nh,al,R%mat(q)%gam(3)%X,np,&
          R%mat(q)%gam(3)%X,np,bet,w2%mat(q)%gam(1)%X,np)
     
     w2%mat(q)%gam(1)%X = w2%mat(q)%gam(1)%X * R%herm
     
     end if
    
    
  end do

end subroutine 


