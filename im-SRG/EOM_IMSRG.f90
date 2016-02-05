module EOM_IMSRG
  use basic_IMSRG
  use EOM_scalar_commutators
  use EOM_TS_commutators
  use operators
  implicit none
  
contains 

subroutine calculate_excited_states( J, PAR, Numstates, HS , jbas,O1) 
  implicit none
  
  real(8) :: BE2,Mfi 
  type(spd) :: jbas
  type(sq_op) :: HS 
  type(sq_op),optional :: O1
  type(sq_op),allocatable,dimension(:) :: ladder_ops 
  integer :: J,PAR,Numstates,i,q
  character(2) :: Jlabel,Plabel,betalabel  
  
  allocate(ladder_ops(numstates)) 
  ladder_ops%herm = 1
 
  ladder_ops%rank = J 
  ladder_ops%dpar = 2*PAR
  
  ladder_ops%pphh_ph = .true. 
  if ( ladder_ops(1)%rank .ne. 0 ) then 
    
     call allocate_tensor(jbas,ladder_ops(1),HS)   
     
     do q = 1,ladder_ops(1)%nblocks
        ladder_ops(1)%tblck(q)%lam(1) = 1 
     end do
  else 
     call duplicate_sq_op(HS,ladder_ops(1)) 
  end if
  
  do i = 2, Numstates
     call duplicate_sq_op(ladder_ops(1),ladder_ops(i),'y')
  end do
  
  print* 
  write(*,'((A55),(I1),(A3),(I1))') 'EXECUTING EOM CALCULATION'// &
       ' FOR EXCITED STATES: J=',J/2,' P=',PAR  
  print*
  call lanczos_diagonalize(jbas,HS,ladder_ops,Numstates)  
  
  
  
  ! if ( present(O1) ) then 
     
  !    print*
  !    write(*,'((A21),(f16.9))') 'Ground State Energy: ',HS%E0 
  !    print*
  !    print*, 'EXCITED STATE ENERGIES:'
  !    print*, '================================================='
  !    print*, '      dE             E_0 + dE        BE(2)       '
  !    print*, '================================================='
  
  !    do i = 1, Numstates
  !       Mfi = transition_to_ground_ME( O1 , ladder_ops(i),jbas )
  !       BE2 = Mfi**2/(J+1.d0) 
  !       write(*,'(3(f16.9))') ladder_ops(i)%E0 ,ladder_ops(i)%E0+HS%E0,BE2 
  !    end do
  
!  else
     
     print*
     write(*,'((A21),(f16.9))') 'Ground State Energy: ',HS%E0 
     print*
     print*, 'EXCITED STATE ENERGIES:'
     print*, '=================================='
     print*, '      dE             E_0 + dE'
     print*, '=================================='
     do i = 1, Numstates
        write(*,'(2(f16.9))') ladder_ops(i)%E0 ,ladder_ops(i)%E0+HS%E0
     end do
 ! end if 
  
  ! WRITE STUFF TO FILES. 
  write( Jlabel ,'(I2)') HS%Jtarg
  write( betalabel ,'(I2)') nint(HS%lawson_beta)
  if (HS%Ptarg == 0 ) then 
     Plabel ='+'
  else 
     Plabel ='-'
  end if
  
  Jlabel = adjustl(Jlabel)
  Plabel = adjustl(Plabel)
  betalabel = adjustl(betalabel)
  open(unit=72,file=trim(OUTPUT_DIR)//&
       trim(adjustl(prefix))//&
       '_EOM_spec_law'//trim(betalabel)//'.dat')
  
  do i = 1, Numstates
     write(72,'(2(f16.9))') ladder_ops(i)%E0 ,ladder_ops(i)%E0+HS%E0
  end do

  close(72)
  
  open(unit=75,file=trim(OUTPUT_DIR)//&
       trim(adjustl(prefix))//'_lawson_check.dat'&
       ,position='append')
  write(75,*) HS%lawson_beta, HS%E0,ladder_ops(1:3)%E0+HS%E0
  close(75)
  
end subroutine 

subroutine LANCZOS_DIAGONALIZE(jbas,OP,Vecs,nev)    
  
  integer :: N 
  integer,intent(in) :: nev
  type(spd) :: jbas
  type(sq_op) :: op,V1,Q1,Q2,w1,w2
  type(sq_op),dimension(nev) :: Vecs
  type(cc_mat) :: OpCC,QCC,WCC
  type(ex_pandya_mat) :: QPP,WPP
  type(ex_cc_mat) :: OpPP
  real(8),allocatable,dimension(:) :: workl,D,eigs,resid,work,workD
  real(8),allocatable,dimension(:,:) :: V,Z
  integer :: i,j,ix,jx,lwork,info,ido,ncv,ldv,iparam(11),ipntr(11),q,II,JJ
  integer :: ishift,mxiter,nb,nconv,mode,np,lworkl,ldz,p,h,sps,tps,jp,jh
  real(8) ::  x,tol,y,sigma,t1,t2
  character(1) :: BMAT,HOWMNY 
  character(2) :: which
  logical :: rvec
  logical,allocatable,dimension(:) :: selct

  Q1%pphh_ph=.true.
  Q2%pphh_ph=.true.
  w1%pphh_ph=.false.
  w2%pphh_ph=.false.

  call duplicate_sq_op(vecs(1),w1,'w') !workspace
  call duplicate_sq_op(vecs(1),w2,'w') !workspace
  call duplicate_sq_op(vecs(1),Q1,'y') !workspace
  call duplicate_sq_op(vecs(1),Q2,'y') !workspace
 
  if (vecs(1)%rank == 0) then 
     call init_ph_mat(Op,OpCC,jbas) !cross coupled ME
     call duplicate_ph_mat(OpCC,QCC) !cross coupled ME     
     call init_ph_wkspc(QCC,WCC) 
  else 
     call init_ph_mat(Op,OpPP,jbas) !cross coupled ME
     call init_ph_mat(vecs(1),QPP,jbas) !cross coupled ME
     call init_ph_wkspc(QPP,WPP) 
  end if
  
  h = OP%belowEF !holes
  p = OP%Nsp-h  !particles
 
  sps = 0 
  do ix = 1,p
     do jx = 1,h
        
        i = jbas%parts(ix)
        j = jbas%holes(jx) 
  
        if (triangle(jbas%jj(i),jbas%jj(j),vecs(1)%rank)) then  

           if (jbas%itzp(i) .ne. jbas%itzp(j) ) cycle
           if (mod(jbas%ll(i) + jbas%ll(j) + vecs(1)%dpar/2,2) .ne. 0 ) cycle
        
           sps = sps+1
        end if 
     end do
  end do
  
  if (vecs(1)%rank == 0) then 
     ! scalar operator
     tps = 0 
     do q = 1, OP%nblocks
     
        do II = 1,OP%mat(q)%npp
           do JJ = 1, OP%mat(q)%nhh 
           
              if (mod(OP%mat(q)%lam(1)/2,2) == 1) then 
                 if ( OP%mat(q)%qn(1)%Y(II,1) == &
                      OP%mat(q)%qn(1)%Y(II,2) ) cycle
                 if ( OP%mat(q)%qn(3)%Y(JJ,1) == &
                      OP%mat(q)%qn(3)%Y(JJ,2) ) cycle
              end if
           
              tps = tps+ 1
           end do
        end do
     end do
  
  else 
     ! tensor case
     tps = 0 
     do q = 1, vecs(1)%nblocks
     
        do II = 1,size( vecs(1)%tblck(q)%tgam(3)%X(:,1) ) 
           do JJ = 1, size( vecs(1)%tblck(q)%tgam(3)%X(1,:) )  
           
              if (mod(vecs(1)%tblck(q)%Jpair(1)/2,2) == 1) then 
                 
                 if ( vecs(1)%tblck(q)%tensor_qn(1,1)%Y(II,1) == &
                      vecs(1)%tblck(q)%tensor_qn(1,1)%Y(II,2) ) cycle
              end if 

              if (mod(vecs(1)%tblck(q)%Jpair(2)/2,2) == 1) then 
   
                 if ( vecs(1)%tblck(q)%tensor_qn(3,2)%Y(JJ,1) == &
                      vecs(1)%tblck(q)%tensor_qn(3,2)%Y(JJ,2) ) cycle
              end if
           
              tps = tps+ 1
           end do
        end do
        
        if (vecs(1)%tblck(q)%Jpair(1) == vecs(1)%tblck(q)%Jpair(2)) cycle
 
        do II = 1,size( vecs(1)%tblck(q)%tgam(7)%X(:,1) ) 
           do JJ = 1, size( vecs(1)%tblck(q)%tgam(7)%X(1,:) )  
           
              if (mod(vecs(1)%tblck(q)%Jpair(1)/2,2) == 1) then 
                 
                 if ( vecs(1)%tblck(q)%tensor_qn(3,1)%Y(II,1) == &
                      vecs(1)%tblck(q)%tensor_qn(3,1)%Y(II,2) ) cycle
              end if 

              if (mod(vecs(1)%tblck(q)%Jpair(2)/2,2) == 1) then 
   
                 if ( vecs(1)%tblck(q)%tensor_qn(1,2)%Y(JJ,1) == &
                      vecs(1)%tblck(q)%tensor_qn(1,2)%Y(JJ,2) ) cycle
              end if
           
              tps = tps+ 1
           end do
        end do

     end do
  end if    
           
  N = sps + tps ! number of ph and pphh SDs 
  
  ido = 0  ! status integer is 0 at start
  BMAT = 'I' ! standard eigenvalue problem (N for generalized) 
  which = 'SM' ! compute smallest eigenvalues in magnitude ('SA') is algebraic. 
  tol = 0.0E+00 ! error tolerance? (wtf zero?) 
  info = 0
  ncv = 5*nev ! number of lanczos vectors I guess
  lworkl = ncv*(ncv+8) 
  allocate(V(N,NCV),workl(lworkl))
  LDV = N  
  ishift = 1
  mxiter = 500 
  mode = 1
   
  allocate(eigs(N),resid(N),work(10*N),workD(3*N)) 

  iparam(1) = ishift
  iparam(3) = mxiter
  iparam(7) = mode
  i = 0

  do 
     ! so V is the krylov subspace matrix that is being diagonalized
     ! it does not need to be initialized, so long as you have the other 
     ! stuff declare right, the code should know this. 
     call dsaupd ( ido, bmat, N, which, nev, tol, resid, &
      ncv, v, ldv, iparam, ipntr, workd, workl, &
      lworkl, info )
     ! The actual matrix only gets multiplied with the "guess" vector in "matvec_prod" 
     call progress_bar( i )
     i=i+1 
  
    if ( ido /= -1 .and. ido /= 1 ) then
      exit
    end if

    if ( vecs(1)%rank == 0 ) then 
       call matvec_prod(N,OP,Q1,Q2,w1,w2,OpCC,QCC,WCC,jbas, workd(ipntr(1)), workd(ipntr(2)) ) 
    else
       call matvec_nonzeroX_prod(N,OP,Q1,Q2,w1,w2,OpPP,QPP,WPP,jbas, workd(ipntr(1)), workd(ipntr(2)) ) 
    end if
    
    end do 
    write(6,*) 
  ! the ritz values are out of order right now. Need to do post
  ! processing to fix this, and get the eigenvectors
  rvec= .true. 
  howmny = 'A'
  
  allocate(selct(NCV)) 
  allocate(D(NEV)) 
  allocate(Z(N,NEV)) 
  ldz = N  
  call dseupd( rvec, howmny, selct, d, Z, ldv, sigma, &
      bmat, n, which, nev, tol, resid, ncv, v, ldv, &
      iparam, ipntr, workd, workl, lworkl, info )
  
  ! right now Z contains the eigenvectors in the columns
  ! d contains the eigenvalues in the same order. 
  do i = 1, nev
     if (Vecs(i)%rank .ne. 0 ) then
        call unwrap_tensor(Z(:,i),Vecs(i),N,jbas) 
     else
        call unwrap(Z(:,i),Vecs(i),N,jbas)
     end if 
     Vecs(i)%E0 = d(i)
  end do 
      
end subroutine

subroutine matvec_prod(N,OP,Q_op,Qout,w1,w2,OpCC,QCC,WCC,jbas,v,w) 
  implicit none 
  
  integer :: N ,q,a,b,c,d,i,j,k,l,Jtot
  real(8) :: sm
  type(sq_op) :: OP, Q_op ,Qout,w1,w2
 ! type(ph_mat) :: QCC,WCC
  type(cc_mat) :: OpCC,QCC,WCC
  type(spd) :: jbas
  real(8),dimension(N) :: v,w 
  real(8) :: coef9,dfact0

  ! the name says mat-vec product, but really this
  ! is a commutator. 
  
  ! the commutator here is equivalent to the matrix-vector product 
  ! with only connected diagrams retained. 

  
  ! FIRST WE NEED TO CONVERT v TO a (SQ_OP) variable
  
  call unwrap(v,Q_op,N,jbas)

  ! now we have two sq_op operators which can be used with my commutator expressions. Noice. 
    
  call calculate_cross_coupled(Op,OpCC,jbas)
  call EOM_scalar_cross_coupled(Q_op,QCC,jbas) 
 
  !print*, EOM_scalar_commutator_110(Op,Q_op,jbas) + &
   !    EOM_scalar_commutator_220(Op,Q_op,jbas) ,sum(v**2)

  call EOM_scalar_commutator_111(Op,Q_op,Qout,jbas) ! verified   
  call EOM_scalar_commutator_121(Op,Q_op,Qout,jbas) ! verified
  call EOM_scalar_commutator_122(Op,Q_op,Qout,jbas)  ! verified, damnit.  
  
  call EOM_scalar_commutator_222_pp_hh(Op,Q_op,Qout,w1,w2,jbas) ! verified   
  call EOM_scalar_commutator_221(Op,Q_op,Qout,w1,w2,jbas)  ! verified.    
  call EOM_scalar_commutator_222_ph(OpCC,QCC,Qout,WCC,jbas) 
 
  ! Okay, now we have the "matrix vector product" So lets put it back in vector form:
  call rewrap(w,Qout,N,jbas) 
end subroutine
!======================================================================================
!======================================================================================
subroutine matvec_nonzeroX_prod(N,OP,Q_op,Qout,w1,w2,OpCC,QCC,WCC,jbas,v,w) 
  implicit none 
  
  integer :: N ,q,a,b,c,d,i,j,k,l,Jtot
  real(8) :: sm
  type(sq_op) :: OP, Q_op ,Qout,w1,w2
  type(ex_pandya_mat) :: QCC,WCC
  type(ex_cc_mat) :: OpCC
  type(spd) :: jbas
  real(8),dimension(N) :: v,w 
  real(8) :: coef9,dfact0

  ! the name says mat-vec product, but really this
  ! is a commutator. 
  
  ! the commutator here is equivalent to the matrix-vector product 
  ! with only connected diagrams retained. 

  
  ! FIRST WE NEED TO CONVERT v TO a (SQ_OP) variable
  
  call unwrap_tensor(v,Q_op,N,jbas)
  ! now we have two sq_op operators which can be used with my commutator expressions. Noice. 
  
  call EOM_generalized_pandya(Q_op,QCC,jbas)
  call calculate_cross_coupled_pphh(Op,OpCC,jbas) 
  
  call EOM_TS_commutator_111(Op,Q_op,Qout,jbas) 
  call EOM_TS_commutator_121(Op,Q_op,Qout,jbas)
  call EOM_TS_commutator_211(OpCC,Q_op,Qout,jbas) 
  call EOM_TS_commutator_122(Op,Q_op,Qout,jbas)
  call EOM_TS_commutator_212(Op,Q_op,Qout,jbas)
  
  call EOM_TS_commutator_222_pp_hh(Op,Q_op,Qout,w1,w2,jbas)  
  call EOM_TS_commutator_221(w1,w2,Op%herm*Q_op%herm,Qout,jbas)
  call EOM_TS_commutator_222_ph(OpCC,QCC,Qout,WCC,jbas)
  
  ! Okay, now we have the "matrix vector product" So lets put it back in vector form:
  call rewrap_tensor(w,Qout,N,jbas) 
end subroutine
!======================================================================================
!======================================================================================
subroutine unwrap( v, AX ,N ,jbas) 
  implicit none 
  
  type(spd) :: jbas
  integer :: N ,i, II,JJ, parts,holes,q,IX,JX
  real(8),dimension(N) :: v
  type(sq_op) :: AX 
  
  i = 1
  
  holes = AX%belowEF
  parts = AX%Nsp- holes 
  
  ! one body
  do IX = 1,parts
     do JX = 1,holes
        
        II = jbas%parts(ix)
        JJ = jbas%holes(jx) 
        
        ! i have to do this because i don't have the 
        ! one body piece blocked by conserved quantities
        if (jbas%jj(II) .ne. jbas%jj(JJ) ) cycle
        if (jbas%itzp(II) .ne. jbas%itzp(JJ) ) cycle
        if (jbas%ll(II) .ne. jbas%ll(JJ) ) cycle
        
        AX%fph(IX,JX) = v(i) 
        i = i + 1
     end do
  end do

  
  ! two body 

  do q = 1, AX%nblocks
     
     do II = 1,AX%mat(q)%npp
        do JJ = 1, AX%mat(q)%nhh 
           
           if (mod(AX%mat(q)%lam(1)/2,2) == 1) then 
              if ( AX%mat(q)%qn(1)%Y(II,1) == &
                   AX%mat(q)%qn(1)%Y(II,2) ) cycle
             if ( AX%mat(q)%qn(3)%Y(JJ,1) == &
                   AX%mat(q)%qn(3)%Y(JJ,2) ) cycle
           end if
           AX%mat(q)%gam(3)%X(II,JJ) = v(i)
           i = i + 1
        end do 
     end do 
  end do 

end subroutine 
!======================================================================================
!======================================================================================
subroutine unwrap_tensor( v, AX ,N ,jbas) 
  implicit none 
  
  type(spd) :: jbas
  integer :: N ,i, II,JJ, parts,holes,q,IX,JX,qx
  real(8),dimension(N) :: v
  type(sq_op) :: AX 
  
  i = 1
  
  holes = AX%belowEF
  parts = AX%Nsp- holes 
  

  do ix = 1,parts
     do jx = 1,holes
        
        ii = jbas%parts(ix)
        JJ = jbas%holes(jx) 
  
        if (triangle(jbas%jj(II),jbas%jj(JJ),AX%rank)) then  

           if (jbas%itzp(II) .ne. jbas%itzp(JJ) ) cycle
           if (mod(jbas%ll(II) + jbas%ll(JJ) + AX%dpar/2,2) .ne. 0 ) cycle
        
           AX%fph(IX,JX) = v(i) 
           i = i+ 1
        end if 
     end do
  end do
  
  ! tensor case
  
  do q = 1, AX%nblocks
     
     do II = 1,size( AX%tblck(q)%tgam(3)%X(:,1) ) 
        do JJ = 1, size( AX%tblck(q)%tgam(3)%X(1,:) )  
           
           if (mod(AX%tblck(q)%Jpair(1)/2,2) == 1) then 
              
              if ( AX%tblck(q)%tensor_qn(1,1)%Y(II,1) == &
                   AX%tblck(q)%tensor_qn(1,1)%Y(II,2) ) cycle
           end if

           if (mod(AX%tblck(q)%Jpair(2)/2,2) == 1) then 

              if ( AX%tblck(q)%tensor_qn(3,2)%Y(JJ,1) == &
                   AX%tblck(q)%tensor_qn(3,2)%Y(JJ,2) ) cycle
           end if

           AX%tblck(q)%tgam(3)%X(II,JJ) = v(i)
           i = i + 1
        end do
     end do

     ! IF THE Js ARE THE SAME THEN WE NEED TO MAKE SURE THAT IS REFLECTED IN 
     ! THE TRANSPOSE
     if (AX%tblck(q)%Jpair(1) == AX%tblck(q)%Jpair(2)) then 
        if (AX%dpar == 2 ) then 
           qx = tensor_block_index(AX%tblck(q)%Jpair(1),AX%tblck(q)%Jpair(2)&
                ,AX%rank,AX%tblck(q)%lam(3),mod(AX%tblck(q)%lam(2)+1,2))
           AX%tblck(qx)%tgam(7)%X = Transpose(AX%tblck(q)%tgam(3)%X)  
        else       
           AX%tblck(q)%tgam(7)%X = Transpose(AX%tblck(q)%tgam(3)%X)  
        end if
        
        cycle 
     
     end if
     
     ! OTHERWISE BUSINESS AS USUAL.
     do II = 1,size( AX%tblck(q)%tgam(7)%X(:,1) ) 
        do JJ = 1, size( AX%tblck(q)%tgam(7)%X(1,:) )  

           if (mod(AX%tblck(q)%Jpair(1)/2,2) == 1) then 

              if ( AX%tblck(q)%tensor_qn(3,1)%Y(II,1) == &
                   AX%tblck(q)%tensor_qn(3,1)%Y(II,2) ) cycle
           end if

           if (mod(AX%tblck(q)%Jpair(2)/2,2) == 1) then 

              if ( AX%tblck(q)%tensor_qn(1,2)%Y(JJ,1) == &
                   AX%tblck(q)%tensor_qn(1,2)%Y(JJ,2) ) cycle
           end if

           AX%tblck(q)%tgam(7)%X(II,JJ) = v(i)
           i = i + 1
        end do
     end do

  end do

end subroutine 
!======================================================================================
!======================================================================================
subroutine rewrap_tensor( v, AX ,N ,jbas) 
  implicit none 
  
  type(spd) :: jbas
  integer :: N ,i, II,JJ, parts,holes,q,IX,JX,qx
  real(8),dimension(N) :: v
  type(sq_op) :: AX 
  
  i = 1
  
  holes = AX%belowEF
  parts = AX%Nsp- holes 
  
  do ix = 1,parts
     do jx = 1,holes
        
        ii = jbas%parts(ix)
        JJ = jbas%holes(jx) 
  
        if (triangle(jbas%jj(II),jbas%jj(JJ),AX%rank)) then  

           if (jbas%itzp(II) .ne. jbas%itzp(JJ) ) cycle
           if (mod(jbas%ll(II) + jbas%ll(JJ) + AX%dpar/2,2) .ne. 0 ) cycle
        
           v(i) = AX%fph(IX,JX) 
           i = i + 1
        end if 
     end do
  end do
  
  ! tensor case
  
  do q = 1, AX%nblocks
     
     do II = 1,size( AX%tblck(q)%tgam(3)%X(:,1) ) 
        do JJ = 1, size( AX%tblck(q)%tgam(3)%X(1,:) )  
           
           if (mod(AX%tblck(q)%Jpair(1)/2,2) == 1) then 
              
              if ( AX%tblck(q)%tensor_qn(1,1)%Y(II,1) == &
                   AX%tblck(q)%tensor_qn(1,1)%Y(II,2) ) cycle
           end if

           if (mod(AX%tblck(q)%Jpair(2)/2,2) == 1) then 

              if ( AX%tblck(q)%tensor_qn(3,2)%Y(JJ,1) == &
                   AX%tblck(q)%tensor_qn(3,2)%Y(JJ,2) ) cycle
           end if

           v(i) = AX%tblck(q)%tgam(3)%X(II,JJ) 
           i = i + 1
        end do
     end do

     ! IF THE Js ARE THE SAME THEN WE NEED TO MAKE SURE THAT IS REFLECTED IN 
     ! THE TRANSPOSE
     if (AX%tblck(q)%Jpair(1) == AX%tblck(q)%Jpair(2)) then 
        
        cycle 
     
     end if
     
     ! OTHERWISE BUSINESS AS USUAL.
     do II = 1,size( AX%tblck(q)%tgam(7)%X(:,1) ) 
        do JJ = 1, size( AX%tblck(q)%tgam(7)%X(1,:) )  

           if (mod(AX%tblck(q)%Jpair(1)/2,2) == 1) then 

              if ( AX%tblck(q)%tensor_qn(3,1)%Y(II,1) == &
                   AX%tblck(q)%tensor_qn(3,1)%Y(II,2) ) cycle
           end if

           if (mod(AX%tblck(q)%Jpair(2)/2,2) == 1) then 

              if ( AX%tblck(q)%tensor_qn(1,2)%Y(JJ,1) == &
                   AX%tblck(q)%tensor_qn(1,2)%Y(JJ,2) ) cycle
           end if

           v(i) = AX%tblck(q)%tgam(7)%X(II,JJ) 
           i = i + 1
        end do
     end do

  end do

end subroutine 
!============================================================================  
!============================================================================
subroutine rewrap( v, AX ,N ,jbas) 
  implicit none 
  
  type(spd) :: jbas
  integer :: N ,i, II,JJ, parts,holes,q,IX,JX
  real(8),dimension(N) :: v
  type(sq_op) :: AX 
  
  i = 1
  
  holes = AX%belowEF
  parts = AX%Nsp- holes 
  
  ! one body
  do IX = 1,parts
     do JX = 1,holes
        
        II = jbas%parts(ix)
        JJ = jbas%holes(jx) 
        
        ! i have to do this because i don't have the 
        ! one body piece blocked by conserved quantities
        if (jbas%jj(II) .ne. jbas%jj(JJ) ) cycle
        if (jbas%itzp(II) .ne. jbas%itzp(JJ) ) cycle
        if (jbas%ll(II) .ne. jbas%ll(JJ) ) cycle
        
        v(i) = AX%fph(IX,JX) 
        i = i + 1
     end do
  end do
  
  ! two body 

  do q = 1, AX%nblocks
     
     do II = 1,AX%mat(q)%npp
        do JJ = 1, AX%mat(q)%nhh 
           
           if (mod(AX%mat(q)%lam(1)/2,2) == 1) then 
              if ( AX%mat(q)%qn(1)%Y(II,1) == &
                   AX%mat(q)%qn(1)%Y(II,2) ) cycle
             if ( AX%mat(q)%qn(3)%Y(JJ,1) == &
                   AX%mat(q)%qn(3)%Y(JJ,2) ) cycle
           end if
           v(i) = AX%mat(q)%gam(3)%X(II,JJ)
           i = i + 1
        end do 
     end do 
  end do 
           

end subroutine 


subroutine progress_bar( step  )  
  implicit none
  
  integer :: step
 
  if ( step .ne. 0) then 
     ! hold the backspace key down
     write(6,'(A)',advance='no') char(8)//char(8)//char(8)//char(8)//char(8) &
       //char(8)//char(8)//char(8)//char(8)//char(8)//char(8)//char(8)&
       //char(8)//char(8)//char(8)//char(8)//char(8)//char(8)//char(8)//char(8) &
       //char(8)//char(8)//char(8)//char(8)
     flush 6
  end if 

  write(6,'((A19),(I5))',advance='no') 'Lanczos iteration: ',step
  flush 6 

end subroutine 

end module
