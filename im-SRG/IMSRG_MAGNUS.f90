module IMSRG_MAGNUS
  use commutators
  use TS_commutators
  use operators
  use HF_mod 
  implicit none 
  
contains

subroutine magnus_decouple(HS,G,jbas,quads,trips,build_generator)
  ! runs IMSRG using magnus expansion method
  implicit none 
  
  integer :: Atot,Ntot,nh,np,nb,q,steps,i,j
  type(spd) :: jbas
  type(tpd),allocatable,dimension(:) :: threebas
  type(sq_op) :: H,G,HS,AD
  type(sq_op) :: DG,CR
  type(cc_mat) :: GCC,ADCC,WCC 
  real(8) :: ds,s,sx,Eold,E_mbpt2,crit,nrm1,nrm2,wTs(2),Ecm(3),corr,dcgi00,xxx
  real(8) :: omp_get_wtime,t1,t2
  character(1) :: quads,trips
  logical :: trip_calc,xxCR,chkpoint_restart
  external :: build_generator 
  
  trip_calc = .false. 
  xxCr = .false. 

  if (trips=='y') then 
     ! triples correction
     trip_calc=.true.
  else if (trips == 'C') then 
     ! triples correction with full internal BCH
     trip_calc=.true. 
     xxCr = .true.
  end if 
     
  HS%neq = 1
  call duplicate_sq_op(HS,H) !evolved hamiltonian
  if (.not. allocated(G%mat)) call duplicate_sq_op(HS,G) !magnus operator
  call duplicate_sq_op(HS,AD) !workspace

  call duplicate_sq_op(HS,DG) !magnus operator
  G%herm = -1 
  DG%herm = -1 
   
  call build_generator(HS,DG,jbas)  
  call copy_sq_op(HS,H) 
  
  s = 0.d0 
    
  if (HS%lawson_beta < 3.0) then 
     ds = 0.5d0
  else if (HS%lawson_beta < 6.0) then 
     ds = 0.5d0
  else
     ds = 0.1d0
  end if 

  crit = 10.
  steps = 0
  
  if (checkpointing) then 
     chkpoint_restart = read_omega_checkpoint(G,sx) 
  else
     chkpoint_restart = .true. 
  end if 
  
  if (chkpoint_restart) then 
  

     open(unit=36,file=trim(OUTPUT_DIR)//&
       trim(adjustl(prefix))//'_0b_magnus_flow.dat')

     
     write(36,'(I6,4(e15.7))') steps,s,H%E0,HS%E0+E_mbpt2,crit
     write(*,'(I6,4(e15.7))') steps,s,HS%E0,HS%E0+E_mbpt2,crit
     Eold=0.
     E_mbpt2 = mbpt2(HS,jbas) 
     crit=abs(E_mbpt2)
  
  else 
     
     ! CHECKPOINT RESTART
     open(unit=36,file=trim(OUTPUT_DIR)//&
       trim(adjustl(prefix))//'_0b_magnus_flow.dat',position='append')
     
     call BCH_EXPAND(HS,G,H,jbas,quads) 
     call build_generator(HS,DG,jbas)  
     
     s= sx
  
     steps = nint(sx/ds)

  end if 
  
     
  do while (crit > 1e-6) 
     
     steps = steps + 1
     
     if (checkpointing) then 
        if (mod(steps,4)==0) then 
           call write_omega_checkpoint(G,s)
        end if
     end if 
     
     call MAGNUS_EXPAND(DG,G,AD,jbas)
 
     call euler_step(G,DG,s,ds) 
     
     call BCH_EXPAND(HS,G,H,jbas,quads) 
     
     call build_generator(HS,DG,jbas)   

     E_mbpt2 = mbpt2(HS,jbas) 
     
     crit = abs(E_mbpt2)  
       
     write(36,'(I6,4(e15.7))') steps,s,HS%E0,HS%E0+E_mbpt2,crit
     write(*,'(I6,4(e15.7))') steps,s,HS%E0,HS%E0+E_mbpt2,crit

  end do
 
! triples correction ===========================================================
  if (trip_calc) then 
     call enumerate_three_body(threebas,jbas)
     t1 = 0.!fuckfaces
     if ( xxCR ) then 
        call duplicate_sq_op(H,CR)         
        ! completely renormalized bit.
        call CR_EXPAND(CR,G,H,jbas,quads) 
        corr =  restore_triples(CR,G,threebas,jbas)
     else
        corr =  restore_triples(H,G,threebas,jbas) 
     end if 
     t2 = 0.!fuckfaces
   
     print*, 'FINAL ENERGY:', corr + HS%E0,t2-t1
     
     open(unit=39,file=trim(OUTPUT_DIR)//&
       trim(adjustl(prefix))//'_magnus_triples.dat')
     write(39,'(I6,4(e15.7))') steps,s,HS%E0+corr,HS%E0+E_mbpt2,crit
     close(39)
  end if
  
end subroutine  
!===========================================================================
!===========================================================================
subroutine magnus_TDA(HS,TDA,G,jbas,quads,build_generator) 
  ! runs IMSRG TDA decoupling using magnus expansion method
  implicit none 
  
  integer :: Atot,Ntot,nh,np,nb,q,steps,i,Jsing
  type(spd) :: jbas
  type(full_sp_block_mat) :: TDA
  type(sq_op) :: H , G ,ETA, HS,INT1,INT2,AD,w1,w2,DG,G0,ETA0,H0,Oevolv
  type(cc_mat) :: GCC,ADCC,WCC,HCC,OeCC
  real(8) :: ds,s,crit,nrm1,nrm2,wTs(2),Ecm(3)
  real(8),allocatable,dimension(:) :: E_old,wTvec
  character(3) :: args
  character(1) :: quads 
  character(1) :: Jlabel,Plabel
  external :: build_generator 
  
  call duplicate_sq_op(HS,H) !evolved hamiltonian
  call duplicate_sq_op(HS,AD) !workspace
  call duplicate_sq_op(HS,DG) !magnus operator

  G%herm = -1 
  DG%herm = -1
    
  call init_ph_mat(HS,HCC,jbas) !cross coupled ME
  
  ! TDA stuff
  call calculate_cross_coupled(HS,HCC,jbas)
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
  open(unit=37,file=trim(OUTPUT_DIR)//&
       trim(adjustl(prefix))//'_'//Jlabel//Plabel//'_excited.dat')
  
  call write_excited_states(steps,s,TDA,HS%E0,37)
  call build_generator(HS,DG,jbas) 
  call copy_sq_op(HS,H) 
  
  !nrm1 = mat_frob_norm(DG) 
 
  do while (crit > 1e-5) 
  
     call MAGNUS_EXPAND(DG,G,AD,jbas)
     call euler_step(G,DG,s,ds) 
 
     call BCH_EXPAND(HS,G,H,jbas,quads) 
         
     call build_generator(HS,DG,jbas)
    
     call calculate_cross_coupled(HS,HCC,jbas) 
     call calc_TDA(TDA,HS,HCC,jbas) 
     call diagonalize_blocks(TDA)
  
     call write_excited_states(steps,s,TDA,HS%E0,37) 
     
     crit = sum(abs(E_old-TDA%blkM(1)%eigval))/TDA%map(1)
     write(*,'(I6,7(e15.7))') steps,s,TDA%blkM(1)%eigval(1:5),crit
     E_old = TDA%blkM(1)%eigval

     steps = steps + 1
          
  end do
 
!===========================================================================  
  close(37)
end subroutine
!============================================================================
!============================================================================
subroutine magnus_HF(HS,G,jbas,build_generator)
  ! runs IMSRG using magnus expansion method
  implicit none 
  
  integer :: Atot,Ntot,nh,np,nb,q,steps,i,j
  type(spd) :: jbas
  type(tpd),allocatable,dimension(:) :: threebas
  type(sq_op) :: H,G,HS,AD
  type(sq_op) :: DG,CR
  type(cc_mat) :: GCC,ADCC,WCC 
  real(8) :: ds,s,Eold,E_mbpt2,crit,nrm1,nrm2,wTs(2),Ecm(3),corr,dcgi00,xxx
  real(8) :: omp_get_wtime,t1,t2
  character(1) :: quads,trips
  logical :: trip_calc,xxCR,ecrit
  external :: build_generator 
     
  HS%neq = 1
  call duplicate_sq_op(HS,H) !evolved hamiltonian
  if (.not. allocated(G%mat)) call duplicate_sq_op(HS,G) !magnus operator
  call duplicate_sq_op(HS,AD) !workspace

  call duplicate_sq_op(HS,DG) !magnus operator
  G%herm = -1 
  DG%herm = -1 
   
  call build_generator(HS,DG,jbas)  
  call copy_sq_op(HS,H) 
  
  s = 0.d0 
    
  if (HS%lawson_beta < 3.0) then 
     ds = 1.0d0
  else if (HS%lawson_beta < 6.0) then 
     ds = 0.5d0
  else
     ds = 0.1d0
  end if 

  ds = 1.0
  crit = 10.
  steps = 0
  
  open(unit=36,file=trim(OUTPUT_DIR)//&
       trim(adjustl(prefix))//'_0b_magnus_flow.dat')

  write(36,'(I6,4(e15.7))') steps,s,H%E0,HS%E0+E_mbpt2,crit
  write(*,'(I6,4(e15.7))') steps,s,HS%E0,HS%E0+E_mbpt2,crit
  Eold=0.
  do while (crit > 1e-12) 
     
     call MAGNUS_EXPAND_1b(DG,G,AD,jbas)
     
     call euler_step(G,DG,s,ds) 
     
     call BCH_EXPAND_1b(HS,G,H,jbas,quads) 
    
     call build_generator(HS,DG,jbas)   
     
     crit = abs(HS%E0-Eold) 
     Eold = HS%E0
     
     steps = steps + 1
 !    write(36,'(I6,4(e15.7))') steps,s,HS%E0,HS%E0+E_mbpt2,crit
     write(*,'(I6,4(e15.7))') steps,s,HS%E0,HS%E0+E_mbpt2,crit

  end do
  
end subroutine  
!===========================================================================
!===========================================================================
!=========================================================================
!=========================================================================
subroutine transform_observable_BCH(Op,G,jbas,quads)
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: Op,G,Oevolved
  character(1) :: quads
  
  call duplicate_sq_op(Op,Oevolved) 
  
  if (Op%rank > 0) then      
     call BCH_TENSOR(Oevolved,G,Op,jbas,quads)    
  else
     call BCH_EXPAND(Oevolved,G,Op,jbas,quads)
  end if 
  
  call copy_sq_op(Oevolved,Op) 
 
end subroutine   
!=========================================================================
!=========================================================================
subroutine BCH_EXPAND(HS,G,H,jbas,quads) 
  implicit none 
  
  real(8), parameter :: conv = 1e-6
  integer :: trunc,i,m,n,q,j,k,l,a,b,c,d,iw
  integer :: ix,jx,kx,lx,ax,cx,bx,dx,jmin,jmax,Jtot
  integer :: mi,mj,mk,ml,ma,mc,mb,md,ja,jb,jj,ji,JT,MT
  type(spd) :: jbas
  type(sq_op) :: H , G, ETA, INT1, INT2, INT3,HS, AD,w1,w2
  type(cc_mat) :: WCC,ADCC,GCC
  real(8) :: adnorm,fullnorm,s,advals(15),sm,sm2,coef
  character(3) :: args
  character(1) :: quads ! enter some character to restore quadrupoles 

  call duplicate_sq_op(HS,w1) !workspace
  call duplicate_sq_op(HS,w2) !workspace
  call duplicate_sq_op(HS,INT1) !workspace
  call duplicate_sq_op(HS,INT2) !workspace
  call duplicate_sq_op(HS,AD) !workspace
  INT2%herm = 1
  INT1%herm = 1 
  AD%herm = 1
  call init_ph_mat(AD,ADCC,jbas) !cross coupled ME
  call duplicate_ph_mat(ADCC,GCC) !cross coupled ME
  call init_ph_wkspc(ADCC,WCC) ! workspace for CCME

  advals = 0.d0 
  coef = 1.d0 

  ! intermediates must be HERMITIAN
  INT2%herm = 1
  INT1%herm = 1 
  AD%herm = 1
  
  !! so here:  H is the current hamiltonian and we 
  !! copy it onto HS.  We copy this onto 
  !! INT2, just to make things easier for everyone.
  if (quads=='y') call duplicate_sq_op(H,INT3)

  call copy_sq_op( H , HS )  !basic_IMSRG
  call copy_sq_op( HS , INT2 )
 
  advals(1) = abs(H%E0)   

  do iw = 2 ,30

     coef = coef/(iw-1.d0) 
     ! current value of HS is renamed INT1 
     ! INT2 is renamed AD, for the AD parameters in BCH and magnus expansions
     ! AD refers AD_n and INT2 becomes AD_{n+1} . 
     call copy_sq_op( HS , INT1) 
     call copy_sq_op( INT2 , AD ) 
     ! so to start, AD is equal to H
     call clear_sq_op(INT2)    
     !now: INT2 = [ G , AD ]  
        
! zero body commutator 
     call calculate_cross_coupled(AD,ADCC,jbas)
     call calculate_cross_coupled(G,GCC,jbas) 
  
     INT2%E0 = commutator_110(G,AD,jbas) + commutator_220(G,AD,jbas)
     
     call commutator_111(G,AD,INT2,jbas) 
     call commutator_121(G,AD,INT2,jbas)
     call commutator_122(G,AD,INT2,jbas)    
          
     call commutator_222_pp_hh(G,AD,INT2,w1,w2,jbas)

     call commutator_221(G,AD,INT2,w1,w2,jbas)
     call commutator_222_ph(GCC,ADCC,INT2,WCC,jbas)
     
     ! so now just add INT1 + c_n * INT2 to get current value of HS
     if (quads=='y') then 
        call add_sq_op(INT3, 1.d0 , INT2, 1.d0, INT2)          
        call restore_quadrupoles(AD,G,w1,w2,INT3,jbas) 
     end if 

     call add_sq_op(INT1 ,1.d0 , INT2 , coef , HS )   !basic_IMSRG
          
     advals(iw) = mat_frob_norm(INT2)*coef
     if (advals(iw) < conv) exit
     
  end do 
  
end subroutine BCH_EXPAND
!=========================================================================
!=========================================================================
subroutine BCH_EXPAND_1b(HS,G,H,jbas,quads) 
  ! for 1-body transformations
  implicit none 
  
  real(8), parameter :: conv = 1e-6
  integer :: trunc,i,m,n,q,j,k,l,a,b,c,d,iw
  integer :: ix,jx,kx,lx,ax,cx,bx,dx,jmin,jmax,Jtot
  integer :: mi,mj,mk,ml,ma,mc,mb,md,ja,jb,jj,ji,JT,MT
  type(spd) :: jbas
  type(sq_op) :: H , G, ETA, INT1, INT2, INT3,HS, AD,w1,w2
  type(cc_mat) :: WCC,ADCC,GCC
  real(8) :: adnorm,fullnorm,s,advals(15),sm,sm2,coef
  character(3) :: args
  character(1) :: quads ! enter some character to restore quadrupoles 

  call duplicate_sq_op(HS,INT1) !workspace
  call duplicate_sq_op(HS,INT2) !workspace
  call duplicate_sq_op(HS,AD) !workspace
  INT2%herm = 1
  INT1%herm = 1 
  AD%herm = 1

  advals = 0.d0 
  coef = 1.d0 

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

  do iw = 2 ,30
     coef = coef/(iw-1.d0) 
     ! current value of HS is renamed INT1 
     ! INT2 is renamed AD, for the AD parameters in BCH and magnus expansions
     ! AD refers AD_n and INT2 becomes AD_{n+1} . 
     call copy_sq_op( HS , INT1) 
     call copy_sq_op( INT2 , AD ) 
     ! so to start, AD is equal to H
     call clear_sq_op(INT2)    
     !now: INT2 = [ G , AD ]  
        
! zero body commutator 
     INT2%E0 = commutator_110(G,AD,jbas)
     
     call commutator_111(G,AD,INT2,jbas) 
     call commutator_121(G,AD,INT2,jbas)
     call commutator_122(G,AD,INT2,jbas)    
     
     ! so now just add INT1 + c_n * INT2 to get current value of HS
     call add_sq_op(INT1 ,1.d0 , INT2 , coef , HS )   !basic_IMSRG
          
     advals(iw) = mat_frob_norm(INT2)*coef
     if (advals(iw) < conv) exit
     
  end do 
  
end subroutine BCH_EXPAND_1b
!=========================================================================
!=========================================================================
subroutine BCH_TENSOR(HS,G,H,jbas,quads) 
  implicit none 
  
  real(8), parameter :: conv = 1e-6
  integer :: trunc,i,m,n,q,j,k,l,a,b,c,d,iw
  integer :: ix,jx,kx,lx,ax,cx,bx,dx,jmin,jmax,Jtot
  integer :: mi,mj,mk,ml,ma,mc,mb,md,ja,jb,jj,ji,JT,MT
  type(spd) :: jbas
  type(sq_op) :: H , G, ETA, INT1, INT2, INT3,HS, AD,w1,w2
  type(pandya_mat) :: WCC,ADCC
  type(cc_mat) :: GCC 
  real(8) ::  coef,adnorm,fullnorm,s,advals(15),sm,sm2,dcgi,dcgi00
  character(3) :: args
  character(1) :: quads ! enter some character to restore quadrupoles 
  
  call duplicate_sq_op(HS,w1) !workspace
  call duplicate_sq_op(HS,w2) !workspace
  call duplicate_sq_op(HS,INT1) !workspace
  call duplicate_sq_op(HS,INT2) !workspace
  call duplicate_sq_op(HS,AD) !workspace
  INT2%herm = 1
  INT1%herm = 1 
  AD%herm = 1
  call init_ph_mat(AD,ADCC,jbas) !cross coupled ME
  call init_ph_mat(G,GCC,jbas) !cross coupled ME
  call init_ph_wkspc(ADCC,WCC) ! workspace for CCME
  
  advals = 0.d0 
  coef = 1.d0
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

     
  do iw = 2 ,15
     
     coef = coef/(iw-1.d0)
     ! current value of HS is renamed INT1 
     ! INT2 is renamed AD, for the AD parameters in BCH and magnus expansions
     ! AD refers AD_n and INT2 becomes AD_{n+1} . 
     call copy_sq_op( HS , INT1) 
     call copy_sq_op( INT2 , AD ) 
     ! so to start, AD is equal to H
     call clear_sq_op(INT2)    
     !now: INT2 = [ G , AD ]  
        
! zero body commutator 
     call calculate_generalized_pandya(AD,ADCC,jbas)
     call calculate_cross_coupled(G,GCC,jbas) 
       
     call TS_commutator_111(G,AD,INT2,jbas) 
     call TS_commutator_121(G,AD,INT2,jbas)
     call TS_commutator_211(GCC,AD,INT2,jbas)      
     call TS_commutator_122(G,AD,INT2,jbas)   
     call TS_commutator_212(G,AD,INT2,jbas)
   
     call TS_commutator_222_pp_hh(G,AD,INT2,w1,w2,jbas)

     call TS_commutator_221(w1,w2,G%herm*AD%herm,INT2,jbas)
     call TS_commutator_222_ph(GCC,ADCC,INT2,WCC,jbas)
     
     ! so now just add INT1 + c_n * INT2 to get current value of HS
     call add_sq_op(INT1 ,1.d0 , INT2 , coef , HS )   !basic_IMSRG
         
     advals(iw) = mat_frob_norm(INT2)*coef
    if (advals(iw) < conv) exit
     
  end do 
 
end subroutine BCH_TENSOR
!=========================================================================
!=========================================================================
subroutine CR_EXPAND(HS,G,H,jbas,quads) 
  implicit none 
  
  real(8), parameter :: conv = 1e-6
  integer :: trunc,i,m,n,q,j,k,l,a,b,c,d,iw
  integer :: ix,jx,kx,lx,ax,cx,bx,dx,jmin,jmax,Jtot
  integer :: mi,mj,mk,ml,ma,mc,mb,md,ja,jb,jj,ji,JT,MT
  type(spd) :: jbas
  type(sq_op) :: H , G, ETA, INT1, INT2, INT3,HS, AD,w1,w2
  type(cc_mat) :: WCC,ADCC,GCC
  real(8) ::  coef,adnorm,fullnorm,s,advals(14),sm,sm2,dcgi,dcgi00
  character(3) :: args
  character(1) :: quads ! enter some character to restore quadrupoles 

  call duplicate_sq_op(HS,w1) !workspace
  call duplicate_sq_op(HS,w2) !workspace
  call duplicate_sq_op(HS,INT1) !workspace
  call duplicate_sq_op(HS,INT2) !workspace
  call duplicate_sq_op(HS,AD) !workspace
  INT2%herm = 1
  INT1%herm = 1 
  AD%herm = 1
  call init_ph_mat(AD,ADCC,jbas) !cross coupled ME
  call duplicate_ph_mat(ADCC,GCC) !cross coupled ME
  call init_ph_wkspc(ADCC,WCC) ! workspace for CCME


  advals = 0.d0 
  
  coef = 1.d0
  ! intermediates must be HERMITIAN
  INT2%herm = 1
  INT1%herm = 1 
  AD%herm = 1
  
  !! so here:  H is the current hamiltonian and we 
  !! copy it onto HS.  We copy this onto 
  !! INT2, just to make things easier for everyone.
  if (quads=='y') call duplicate_sq_op(H,INT3)

  call copy_sq_op( H , HS )  !basic_IMSRG
  call copy_sq_op( HS , INT2 )
 
  advals(1) = abs(H%E0)   

  do iw = 2 ,14
     coef= coef/float(iw)
     ! current value of HS is renamed INT1 
     ! INT2 is renamed AD, for the AD parameters in BCH and magnus expansions
     ! AD refers AD_n and INT2 becomes AD_{n+1} . 
     call copy_sq_op( HS , INT1) 
     call copy_sq_op( INT2 , AD ) 
     ! so to start, AD is equal to H
     call clear_sq_op(INT2)    
     !now: INT2 = [ G , AD ]  
        
! zero body commutator
 
     call calculate_cross_coupled(AD,ADCC,jbas)
     call calculate_cross_coupled(G,GCC,jbas) 
 
     INT2%E0 = commutator_110(G,AD,jbas) + commutator_220(G,AD,jbas)

     call commutator_111(G,AD,INT2,jbas) 
     call commutator_121(G,AD,INT2,jbas)
     call commutator_122(G,AD,INT2,jbas)    

     call commutator_222_pp_hh(G,AD,INT2,w1,w2,jbas)
  
     call commutator_221(G,AD,INT2,w1,w2,jbas)
     call commutator_222_ph(GCC,ADCC,INT2,WCC,jbas)

     ! so now just add INT1 + c_n * INT2 to get current value of HS
     if (quads=='y') then
        call add_sq_op(INT3, 1.d0 , INT2, 1.d0, INT2)        
        call restore_quadrupoles(AD,G,w1,w2,INT3,jbas) 
     end if           
     
     call add_sq_op(INT1 ,1.d0 , INT2 , coef , HS )   !basic_IMSRG
         
     advals(iw) = mat_frob_norm(INT2)*coef
     if (advals(iw) < conv) exit
     
  end do 
 
end subroutine 
!===============================================================
!===============================================================
subroutine MAGNUS_EXPAND(DG,G,AD,jbas)
  implicit none 
  
  real(8), parameter :: conv = 1e-6
  integer :: trunc,i,q,j,k,l,ry
  type(spd) :: jbas
  type(sq_op) :: H , G, ETA, INT1, INT2, HS, AD,w1,w2,DG
  type(cc_mat) :: WCC,ADCC,GCC
  real(8) ::  cof(7),adnorm,fullnorm,s,advals(7) 
  character(3) :: args
  
  call duplicate_sq_op(DG,w1) !workspace
  call duplicate_sq_op(DG,w2) !workspace
  call duplicate_sq_op(DG,INT1) !workspace
  call duplicate_sq_op(DG,INT2) !workspace
  INT2%herm = -1
  INT1%herm = -1 
  AD%herm = -1
  
  call init_ph_mat(AD,ADCC,jbas) !cross coupled ME
  call duplicate_ph_mat(ADCC,GCC) !cross coupled ME
  call init_ph_wkspc(ADCC,WCC) ! workspace for CCME

  advals = 0.d0 
  ! Intermediates are ANTI-HERMITIAN 
  
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
          
     call calculate_cross_coupled(AD,ADCC,jbas)
     call calculate_cross_coupled(G,GCC,jbas) 
 
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
subroutine MAGNUS_EXPAND_1b(DG,G,AD,jbas)
  implicit none 
  
  real(8), parameter :: conv = 1e-6
  integer :: trunc,i,q,j,k,l,ry
  type(spd) :: jbas
  type(sq_op) :: H , G, ETA, INT1, INT2, HS, AD,w1,w2,DG
  type(cc_mat) :: WCC,ADCC,GCC
  real(8) ::  cof(7),adnorm,fullnorm,s,advals(7) 
  character(3) :: args
  
  call duplicate_sq_op(DG,INT1) !workspace
  call duplicate_sq_op(DG,INT2) !workspace
  INT2%herm = -1
  INT1%herm = -1 
  AD%herm = -1
  
  advals = 0.d0 
  ! Intermediates are ANTI-HERMITIAN 
  
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
  
     call commutator_111(G,AD,INT2,jbas) 
     call commutator_121(G,AD,INT2,jbas)
     call commutator_122(G,AD,INT2,jbas)    
                 
     call add_sq_op(INT1 , 1.d0 , INT2 , cof(i) , DG ) !ME_general
     
     advals(i) = mat_frob_norm(INT2)*abs(cof(i))
    
  end do   
     
end subroutine MAGNUS_EXPAND_1b
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
  real(8) :: sm,denom,dlow,w
  
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
                 ! w = commutator_223_single(OM,H,a,b,c,i,j,k,Jtot,jab,jij,jbas)
                 ! sm = sm + w*w/denom*(Jtot+1.d0)
                 sm =sm  + commutator_223_single(OM,H,a,b,c,i,j,k,Jtot,jab,jij,jbas)**2&
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
              h1 = hb4(i)+1 
              h2 = hb4(j)+1
              
              
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
              p1 = pb4(a)+1 
              p2 = pb4(b)+1
              
              
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
