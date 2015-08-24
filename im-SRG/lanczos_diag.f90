module lanczos_diag
  use basic_IMSRG
  use EOM_scalar_commutators
  implicit none
  

contains 

subroutine LANCZOS_DIAGONALIZE(jbas,OP,Vecs,nev)    
  
  integer :: N 
  integer,intent(in) :: nev
  type(spd) :: jbas
  type(sq_op) :: op,V1,Q1,Q2,w1,w2
  type(sq_op),dimension(nev) :: Vecs
  type(cross_coupled_31_mat) :: OpCC,QCC,WCC
  real(8),allocatable,dimension(:) :: workl,D,eigs,resid,work,workD
  real(8),allocatable,dimension(:,:) :: V,Z
  integer :: i,j,ix,jx,lwork,info,ido,ncv,ldv,iparam(11),ipntr(11),q,II,JJ
  integer :: ishift,mxiter,nb,nconv,mode,np,lworkl,ldz,p,h,sps,tps,jp,jh
  real(8) ::  x,tol,y,sigma,t1,t2
  character(1) :: BMAT,HOWMNY 
  character(2) :: which
  logical :: rvec
  logical,allocatable,dimension(:) :: selct

  call duplicate_sq_op(Op,w1) !workspace
  call duplicate_sq_op(Op,w2) !workspace
  call duplicate_sq_op(Op,Q1) !workspace
  call duplicate_sq_op(Op,Q2) !workspace

  call allocate_CCMAT(Op,OpCC,jbas) !cross coupled ME
  call duplicate_CCMAT(OpCC,QCC) !cross coupled ME
  call allocate_CC_wkspc(OpCC,WCC) ! workspace for CCME
  
  h = OP%belowEF !holes
  p = OP%Nsp-h  !particles
 
  sps = 0 
  do ix = 1,p
     do jx = 1,h
        
        i = jbas%parts(ix)
        j = jbas%holes(jx) 
        
        if (jbas%jj(i) .ne. jbas%jj(j) ) cycle
        if (jbas%itzp(i) .ne. jbas%itzp(j) ) cycle
        if (jbas%ll(i) .ne. jbas%ll(j) ) cycle
        
        sps = sps+1
     end do
  end do
  
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
     print*, i
     i=i+1 
  
    if ( ido /= -1 .and. ido /= 1 ) then
      exit
    end if

    call matvec_prod(N,OP,Q1,Q2,w1,w2,OpCC,QCC,WCC,jbas, workd(ipntr(1)), workd(ipntr(2)) ) 
  end do 
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
     call unwrap(Z(:,i),Vecs(i),N,jbas) 
     Vecs(i)%E0 = d(i)
  end do 
      
end subroutine

subroutine matvec_prod(N,OP,Q_op,Qout,w1,w2,OpCC,QCC,WCC,jbas,v,w) 
  implicit none 
  
  integer :: N ,q,a,b,c,d,i,j,k,l,Jtot
  real(8) :: sm
  type(sq_op) :: OP, Q_op ,Qout,w1,w2
  type(cross_coupled_31_mat) :: OpCC,QCC,WCC
  
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
    
  call calculate_cross_coupled(Op,OpCC,jbas,.true.)
 
  call EOM_scalar_cross_coupled(Q_op,QCC,jbas,.false.) 
 
  !print*, EOM_scalar_commutator_110(Op,Q_op,jbas) + &
   !    EOM_scalar_commutator_220(Op,Q_op,jbas) ,sum(v**2)

  call EOM_scalar_commutator_111(Op,Q_op,Qout,jbas) ! verified   
  call EOM_scalar_commutator_121(Op,Q_op,Qout,jbas) ! verified
  call EOM_scalar_commutator_122(Op,Q_op,Qout,jbas)  ! verified, damnit.  
  
  call EOM_scalar_commutator_222_pp_hh(Op,Q_op,Qout,w1,w2,jbas)   
  call EOM_scalar_commutator_221(Op,Q_op,Qout,w1,w2,jbas)  ! verified.    
  call EOM_scalar_commutator_222_ph(OpCC,QCC,Qout,WCC,jbas)
  
  

  ! Okay, now we have the "matrix vector product" So lets put it back in vector form:
  call rewrap(w,Qout,N,jbas) 
end subroutine


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

end module
