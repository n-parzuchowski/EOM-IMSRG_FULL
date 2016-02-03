module brute_force_testing
  use cross_coupled
  use commutators
  use TS_commutators
  use EOM_TS_commutators
  use EOM_scalar_commutators
  
  
contains
!============================================================
!============================================================
subroutine construct_random_rank0(OP,HERM,jbas) 
  implicit none 
  
  integer,intent(in):: HERM
  type(sq_op) :: OP 
  type(spd) :: jbas
  integer :: i,j,k,l,ji,jj,jk,jl
  integer :: J1,J2,IX,JX,q,ig,jg 
  real(8) :: x
  
  OP%rank = 0
  OP%herm = HERM 
  
  do ig = 1, OP%belowEF
     do jg = ig, OP%belowEF
        
        i = jbas%holes(ig)
        j = jbas%holes(jg)
 
        ji = jbas%jj(i) 
        jj = jbas%jj(j)

        if (ji .ne. jj) cycle        
        if (jbas%itzp(i) .ne. jbas%itzp(j) ) cycle
        if (jbas%ll(i) .ne. jbas%ll(j) ) cycle

        call random_number(x) 
        x = 10.*x-5.
        OP%fhh(ig,jg) = x
        OP%fhh(jg,ig) = OP%herm * x
        
        
     end do
     
     OP%fhh(ig,ig) = (OP%herm+1)*OP%fhh(ig,ig)/2.d0
     
  end do

  do ig = 1, OP%nsp-OP%belowEF
     do jg = 1, OP%belowEF

        i = jbas%parts(ig)
        j = jbas%holes(jg)
        
        ji = jbas%jj(i) 
        jj = jbas%jj(j)

        if (ji .ne. jj) cycle        
        if (jbas%itzp(i) .ne. jbas%itzp(j) ) cycle
        if (jbas%ll(i) .ne. jbas%ll(j) ) cycle

        call random_number(x) 
        x = 10.*x-5.
        OP%fph(ig,jg) = x

     end do
  end do

  do ig = 1,OP%nsp- OP%belowEF
     do jg = ig,OP%nsp- OP%belowEF
        
        i = jbas%parts(ig)
        j = jbas%parts(jg)

        ji = jbas%jj(i) 
        jj = jbas%jj(j)

        if (ji .ne. jj) cycle        
        if (jbas%itzp(i) .ne. jbas%itzp(j) ) cycle
        if (jbas%ll(i) .ne. jbas%ll(j) ) cycle

        call random_number(x) 
        x = 10.*x-5.
        OP%fpp(ig,jg) = x 
        OP%fpp(jg,ig) = OP%herm * x
     end do
     OP%fpp(ig,ig) = (OP%herm+1)*OP%fpp(ig,ig)/2.d0
  end do
  
  do q=1,OP%nblocks
     do i = 1, 6
        do IX = 1,size(OP%mat(q)%gam(i)%X(:,1))
           do JX = 1,size(OP%mat(q)%gam(i)%X(1,:))
              
              call random_number(x) 
              x = 10.*x-5.

              ! PAULI PRINCIPLE 
              if ( OP%mat(q)%qn(sea1(i))%Y(IX,1) == OP%mat(q)%qn(sea1(i))%Y(IX,2) ) then 
                 if ( mod( OP%mat(q)%lam(1)/2, 2) == 1) x = 0.d0 
              end if 

              if ( OP%mat(q)%qn(sea2(i))%Y(JX,1) == OP%mat(q)%qn(sea2(i))%Y(JX,2) ) then 
                 if ( mod( OP%mat(q)%lam(1)/2, 2) == 1) x = 0.d0 
              end if 
              
              OP%mat(q)%gam(i)%X(IX,JX) =  x * (1 + OP%herm *kron_del(IX,JX) ) 
              if (sqs(i)) OP%mat(q)%gam(i)%X(JX,IX) = OP%herm*x*  (1 + OP%herm *kron_del(IX,JX) )  
              
              
           end do 
        end do 
     end do
  end do 

end subroutine 
!============================================================
!============================================================
subroutine construct_random_rankX(OP,HERM,jbas) 
  implicit none 
  
  integer,intent(in) :: HERM
  type(sq_op) :: OP 
  type(spd) :: jbas
  integer :: i,j,k,l,ji,jj,jk,jl,ig,jg
  integer :: J1,J2,IX,JX,q,qx,rank
  real(8) :: x

  rank = OP%rank
  OP%herm = HERM 

  do ig = 1, OP%belowEF
     do jg = ig, OP%belowEF
        
        i = jbas%holes(ig)
        j = jbas%holes(jg)

        ji = jbas%jj(i) 
        jj = jbas%jj(j)

        if (.not. triangle(ji,jj,OP%rank)) cycle        
        if (jbas%itzp(i) .ne. jbas%itzp(j) ) cycle
        if (mod(jbas%ll(i),2) .ne. mod(jbas%ll(j)+op%dpar/2,2)) cycle

        call random_number(x) 
        x = 10.*x-5.
        OP%fhh(ig,jg) = x
        OP%fhh(jg,ig) = (-1)**((ji-jj)/2) * x* OP%herm
         
     end do
     OP%fhh(ig,ig) = (OP%herm+1)*OP%fhh(ig,ig)/2.d0      
  end do

  do ig = 1, OP%nsp- OP%belowEF
     do jg = ig, OP%nsp-OP%belowEF
        
        i = jbas%parts(ig)
        j = jbas%parts(jg)

        ji = jbas%jj(i) 
        jj = jbas%jj(j)

        if (.not. triangle(ji,jj,OP%rank)) cycle        
        if (jbas%itzp(i) .ne. jbas%itzp(j) ) cycle
        if (mod(jbas%ll(i),2) .ne. mod(jbas%ll(j)+op%dpar/2,2)) cycle

        call random_number(x) 
        x = 10.*x-5.
        OP%fpp(ig,jg) = x
        OP%fpp(jg,ig) = (-1)**((ji-jj)/2) * x* OP%herm 
     end do 
     OP%fpp(ig,ig) = (OP%herm+1)*OP%fpp(ig,ig)/2.d0 
  end do 
  
  do ig = 1, OP%nsp- OP%belowEF
     do jg = 1, OP%belowEF
        
        i = jbas%parts(ig)
        j = jbas%holes(jg)

        ji = jbas%jj(i) 
        jj = jbas%jj(j)

        if (.not. triangle(ji,jj,OP%rank)) cycle        
        if (jbas%itzp(i) .ne. jbas%itzp(j) ) cycle
        if (mod(jbas%ll(i),2) .ne. mod(jbas%ll(j)+op%dpar/2,2)) cycle

        call random_number(x) 
        x = 10.*x-5.
        OP%fph(ig,jg) = x
     end do 
  end do 

  do q = 1, OP%nblocks
     do i = 1, 9
         
        do IX = 1, size(OP%tblck(q)%tgam(i)%X(:,1) )
           do JX = 1, size(OP%tblck(q)%tgam(i)%X(1,:) )
              
              call random_number(x)

              x = 10.*x-5.              
              ! PAULI PRINCIPLE 
              if ( OP%tblck(q)%tensor_qn(sea1(i),1)%Y(IX,1) == OP%tblck(q)%tensor_qn(sea1(i),1)%Y(IX,2) ) then 
                 if ( mod( OP%tblck(q)%Jpair(1)/2, 2) == 1) x = 0.d0 
              end if
              
              if ( OP%tblck(q)%tensor_qn(sea2(i),2)%Y(JX,1) == OP%tblck(q)%tensor_qn(sea2(i),2)%Y(JX,2) ) then 
                 if ( mod( OP%tblck(q)%Jpair(2)/2, 2) == 1) x = 0.d0 
              end if
              
              OP%tblck(q)%tgam(i)%X(IX,JX) = x
              
           end do
        end do
        
        

        if ( (OP%tblck(q)%Jpair(1) == OP%tblck(q)%Jpair(2)).and.(sqs(i)) ) then
           if (mod(OP%dpar/2,2) == 0 ) then
              OP%tblck(q)%tgam(i)%X = &
                   0.5*(OP%tblck(q)%tgam(i)%X + OP%herm * Transpose(OP%tblck(q)%tgam(i)%X))
           
           else
              qx = tensor_block_index(OP%tblck(q)%Jpair(1),OP%tblck(q)%Jpair(2)&
                   ,rank,OP%tblck(q)%lam(3),mod(OP%tblck(q)%lam(2)+1,2))
              
              OP%tblck(qx)%tgam(i)%X = Transpose( OP%tblck(q)%tgam(i)%X ) * OP%herm 
              
           end if 
       end if 
        
     end do
     
     
     if ( OP%tblck(q)%Jpair(1) == OP%tblck(q)%Jpair(2) ) then
        if (mod(OP%dpar/2,2) == 0 ) then 

           OP%tblck(q)%tgam(7)%X = Transpose(OP%tblck(q)%tgam(3)%X) * OP%herm
           OP%tblck(q)%tgam(8)%X = Transpose(OP%tblck(q)%tgam(2)%X) * OP%herm
           OP%tblck(q)%tgam(9)%X = Transpose(OP%tblck(q)%tgam(6)%X) * OP%herm
        else
           

           qx = tensor_block_index(OP%tblck(q)%Jpair(1),OP%tblck(q)%Jpair(2)&
                   ,rank,OP%tblck(q)%lam(3),mod(OP%tblck(q)%lam(2)+1,2))
           
           OP%tblck(qx)%tgam(7)%X = Transpose(OP%tblck(q)%tgam(3)%X)*OP%herm
           OP%tblck(qx)%tgam(8)%X = Transpose(OP%tblck(q)%tgam(2)%X)*OP%herm 
           OP%tblck(qx)%tgam(9)%X = Transpose(OP%tblck(q)%tgam(6)%X)*OP%herm
           OP%tblck(qx)%tgam(3)%X = Transpose(OP%tblck(q)%tgam(7)%X)*OP%herm
           OP%tblck(qx)%tgam(2)%X = Transpose(OP%tblck(q)%tgam(8)%X)*OP%herm 
           OP%tblck(qx)%tgam(6)%X = Transpose(OP%tblck(q)%tgam(9)%X)*OP%herm
        
        end if
      
     end if 
     
  end do

end subroutine
!============================================================
!============================================================
subroutine seed_random_number
  implicit none 
  
  integer :: i,a
  real(8) :: x
  logical :: ext
  
  inquire(file='seed',exist=ext)
  
  if (ext) then 
     open(unit=41,file='seed')
     read(41,*) a
     close(41)
  else 
     a = 2934
  end if 
  
  
  do i = 1, a 
     call random_number(x)
  end do 
  
  a = nint(12349.d0*x)
  
  open(unit=41,file='seed')
  write(41,*)  a
  close(41)

end subroutine
!============================================================
!============================================================
subroutine test_scalar_scalar_commutator(jbas,h1,h2) 
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: AA,BB,OUT,w1,w2
  type(cc_mat) :: AACC,BBCC,WCC
  integer :: a,b,c,d,ja,jb,jc,jd,jmin,jmax,PAR,TZ,Jtot,dick
  integer :: hole,part,iii
  integer,intent(in) :: h1,h2
  real(8) :: val,sm
  real(8) :: vv,xx,yy,zz
    
  call seed_random_number
  
  call allocate_blocks(jbas,AA)
  call duplicate_sq_op(AA,BB)
  call duplicate_sq_op(AA,OUT)
  call duplicate_sq_op(AA,w1) !workspace
  call duplicate_sq_op(AA,w2) !workspace
  call init_ph_mat(AA,AACC,jbas) ! cross coupled ME
  call duplicate_ph_mat(AACC,BBCC) !cross coupled ME
  call init_ph_wkspc(BBCC,WCC) !
  
  call construct_random_rank0(AA,h1,jbas) 
  call construct_random_rank0(BB,h2,jbas) 

  OUT%herm = -1* AA%herm * BB%herm 
  
  print*, 'TESTING SCALAR-SCALAR COMMUTATORS' 
  
  OUT%E0 = commutator_110(AA,BB,jbas) + commutator_220(AA,BB,jbas)
  
  val = scalar_scalar_0body_comm(AA,BB,jbas)
 
  if ( abs( val - OUT%E0) > 1e-10 ) then 
     print*, val, OUT%E0
     STOP 'ZERO BODY FAILURE' 
  end if 
  
  call calculate_cross_coupled(BB,BBCC,jbas)
  call calculate_cross_coupled(AA,AACC,jbas)
  
  call commutator_111(AA,BB,OUT,jbas) 
  call commutator_121(AA,BB,OUT,jbas)
  call commutator_122(AA,BB,OUT,jbas)
 
  call commutator_222_pp_hh(AA,BB,OUT,w1,w2,jbas)
 
  call commutator_221(AA,BB,OUT,w1,w2,jbas)
  call commutator_222_ph(AACC,BBCC,OUT,WCC,jbas)
 
 
!  do a = 1, jbas%total_orbits
 !    do b = 1, jbas%total_orbits
  do iii = 1, 50   
     call random_number(vv)
     call random_number(yy)
   
     a = ceiling(vv*(AA%Nsp))
     b = ceiling(yy*(AA%Nsp))
        
        val = scalar_scalar_1body_comm(AA,BB,a,b,jbas) 
        
        if (abs(val-f_elem(a,b,OUT,jbas)) > 1e-10) then
           print*, 'at: ',a,b
           print*, val, f_elem(a,b,OUT,jbas)
           STOP 'ONE BODY FAILURE'  
        end if 
        
        print*, 'success:', a,b
  !   end do 
  end do 
 
  !do a = 1, jbas%total_orbits
  iii = 0 
  do while (iii < 15)  
     call random_number(vv)
     call random_number(xx)
     call random_number(yy)
     call random_number(zz)
   
     a = ceiling(vv*AA%Nsp)
     b = ceiling(xx*AA%Nsp)
     c = ceiling(yy*AA%Nsp)
     d = ceiling(zz*AA%Nsp)
     
     ja = jbas%jj(a) 
   !  do b = 1, jbas%total_orbits
        jb = jbas%jj(b)
        
        PAR = mod(jbas%ll(a) + jbas%ll(b),2) 
        TZ = jbas%itzp(a) + jbas%itzp(b) 
        
    !    do c = 1, jbas%total_orbits
           jc = jbas%jj(c)
     !      do d = 1, jbas%total_orbits
              jd = jbas%jj(d) 
              
              if (PAR .ne. mod(jbas%ll(c) + jbas%ll(d),2)) cycle 
              if ( TZ .ne.  jbas%itzp(c) + jbas%itzp(d) ) cycle
              iii = iii + 1
              jmin = max( abs(ja-jb) , abs(jc-jd) )
              jmax = min( ja+jb , jc+jd) 
              
              do Jtot = jmin,jmax,2
                 
                 val = scalar_scalar_2body_comm(AA,BB,a,b,c,d,Jtot,jbas)
                 
                 if (abs(val-v_elem(a,b,c,d,Jtot,OUT,jbas)) > 1e-6) then
                    print*, 'at:',a,b,c,d, 'J:', Jtot ,val,v_elem(a,b,c,d,Jtot,OUT,jbas)                
                    print*, val,v_elem(a,b,c,d,Jtot,OUT,jbas)
                    STOP 'TWO BODY FAILURE'  
                 end if
                 
              end do 
              
              print*, 'success:', a,b,c,d
!           end do
 !       end do
  !   end do
  end do
  
  print*, ' COMMUTATOR EXPRESSIONS CONFIRMED '
  
end subroutine
!============================================================
!============================================================
subroutine test_EOM_scalar_scalar_commutator(jbas,h1,h2) 
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: AA,BB,OUT,w1,w2
  type(cc_mat) :: AACC,BBCC,WCC
  integer :: a,b,c,d,q,g,ja,jb,jc,jd,jmin
  integer :: ax,bx,cx,dx,jmax,PAR,TZ,Jtot
  integer,intent(in) :: h1,h2
  real(8) :: val
  
  
  call seed_random_number
  
  call allocate_blocks(jbas,AA)
  call duplicate_sq_op(AA,BB)
  call duplicate_sq_op(AA,OUT)
  call duplicate_sq_op(AA,w1) !workspace
  call duplicate_sq_op(AA,w2) !workspace
  call init_ph_mat(AA,AACC,jbas) ! cross coupled ME
  call duplicate_ph_mat(AACC,BBCC) !cross coupled ME
  call init_ph_wkspc(BBCC,WCC) !
  
  call construct_random_rank0(AA,h1,jbas) 
  call construct_random_rank0(BB,h2,jbas) 
  
  ! make BB look like an excitation operator. 
  BB%fhh = 0.d0 
  BB%fpp = 0.d0 
  
  do q = 1, BB%nblocks
     do g = 1, 6
        if (g==3) cycle
        BB%mat(q)%gam(g)%X = 0.d0 
     end do 
  end do 

  OUT%herm = -1* AA%herm * BB%herm 
  
  print*, 'TESTING EOM SCALAR-SCALAR COMMUTATORS' 
  
  OUT%E0 = EOM_scalar_commutator_110(AA,BB,jbas) + &
       EOM_scalar_commutator_220(AA,BB,jbas)
  
  val = EOM_scalar_scalar_0body_comm(AA,BB,jbas)
 
  if ( abs( val - OUT%E0) > 1e-10 ) then 
     print*, val, OUT%E0
     STOP 'ZERO BODY FAILURE' 
  end if 
  
  call EOM_scalar_cross_coupled(BB,BBCC,jbas)
  call calculate_cross_coupled(AA,AACC,jbas) 
  
  call EOM_scalar_commutator_111(AA,BB,OUT,jbas) 
  call EOM_scalar_commutator_121(AA,BB,OUT,jbas)
  call EOM_scalar_commutator_122(AA,BB,OUT,jbas)
  
  call EOM_scalar_commutator_222_pp_hh(AA,BB,OUT,w1,w2,jbas)
  
  call EOM_scalar_commutator_221(AA,BB,OUT,w1,w2,jbas)
  call EOM_scalar_commutator_222_ph(AACC,BBCC,OUT,WCC,jbas)
  
  do ax = 1, AA%Nsp-AA%belowEF
     a = jbas%parts(ax)
     do bx = 1, AA%belowEF
        b= jbas%holes(bx) 
        
        val = EOM_scalar_scalar_1body_comm(AA,BB,a,b,jbas) 
        
        if (abs(val-f_elem(a,b,OUT,jbas)) > 1e-10) then
           print*, 'at: ',a,b
           print*, val, f_elem(a,b,OUT,jbas)
           STOP 'ONE BODY FAILURE'  
        end if 
        
        print*, 'success:', a,b
     end do 
  end do 
 
  do ax = 1, AA%Nsp-AA%belowEF
     a = jbas%parts(ax)
     ja = jbas%jj(a) 
     do bx = 1, AA%Nsp-AA%belowEF
        b = jbas%parts(bx)
        jb = jbas%jj(b)
        
        PAR = mod(jbas%ll(a) + jbas%ll(b),2) 
        TZ = jbas%itzp(a) + jbas%itzp(b) 
        
        do cx = 1, AA%belowEF
           c = jbas%holes(cx) 
           jc = jbas%jj(c)
           do dx = 1, AA%belowEF
              d = jbas%holes(dx) 
              jd = jbas%jj(d) 
              
              if (PAR .ne. mod(jbas%ll(c) + jbas%ll(d),2)) cycle 
              if ( TZ .ne.  jbas%itzp(c) + jbas%itzp(d) ) cycle
              
              jmin = max( abs(ja-jb) , abs(jc-jd) )
              jmax = min( ja+jb , jc+jd) 
              
              do Jtot = jmin,jmax,2
                 
                 val = EOM_scalar_scalar_2body_comm(AA,BB,a,b,c,d,Jtot,jbas)
                 
                 if (abs(val-v_elem(a,b,c,d,Jtot,OUT,jbas)) > 1e-6) then
                    print*, 'at:',a,b,c,d, 'J:', Jtot ,val,v_elem(a,b,c,d,Jtot,OUT,jbas)                
                    print*, val,v_elem(a,b,c,d,Jtot,OUT,jbas)
                    STOP 'TWO BODY FAILURE'  
                 end if
              end do 
              
              print*, 'success:', a,b,c,d
           end do
        end do
     end do
  end do
  
  print*, ' COMMUTATOR EXPRESSIONS CONFIRMED '
  
end subroutine
!============================================================
!============================================================
subroutine test_scalar_tensor_commutator(jbas,h1,h2,rank,dpar) 
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: AA,BB,OUT,w1,w2
  type(pandya_mat) :: BBCC,WCC
  type(cc_mat) :: AACC
  integer :: a,b,c,d,ja,jb,jc,jd,j1min,j1max
  integer :: j2min,j2max,PAR,TZ,J1,J2,dpar,iii
  integer,intent(in) :: h1,h2,rank
  real(8) :: val,t1,t2,t3,t4,omp_get_wtime
  real(8) :: vv,xx,yy,zz
  
  
  call seed_random_number
  
  BB%rank = rank
  BB%dpar = par
  AA%rank = 0
  BB%pphh_ph = .false.
 
  call allocate_blocks(jbas,AA)
  call allocate_tensor(jbas,BB,AA)
  call construct_random_rank0(AA,h1,jbas) 
  call construct_random_rankX(BB,h2,jbas) 
 
  call duplicate_sq_op(BB,OUT)
  call duplicate_sq_op(BB,w1) !workspace
  call duplicate_sq_op(BB,w2) !workspace
  call init_ph_mat(AA,AACC,jbas) ! cross coupled ME
  call init_ph_mat(BB,BBCC,jbas) !cross coupled ME
  call init_ph_wkspc(BBCC,WCC) !
  

  OUT%herm = -1* AA%herm * BB%herm 
  
  print*, 'TESTING SCALAR-TENSOR COMMUTATORS' 
!  t1 = OMP_get_wtime()
  call calculate_generalized_pandya(BB,BBCC,jbas)
!  t2 = OMP_get_wtime()
  call calculate_cross_coupled(AA,AACC,jbas) 
  
  call TS_commutator_111(AA,BB,OUT,jbas) 
  call TS_commutator_121(AA,BB,OUT,jbas)
  call TS_commutator_211(AACC,BB,OUT,jbas) 
  call TS_commutator_122(AA,BB,OUT,jbas)
  call TS_commutator_212(AA,BB,OUT,jbas)
  
  call TS_commutator_222_pp_hh(AA,BB,OUT,w1,w2,jbas)
  
  call TS_commutator_221(w1,w2,AA%herm*BB%herm,OUT,jbas)
!  t4 = OMP_get_wtime()
  call TS_commutator_222_ph(AACC,BBCC,OUT,WCC,jbas)
!  t3 = OMP_get_wtime()
  
  print*, 'time:', t3-t1,t2-t1,t3-t4
 
!goto 12
!  do a = 1, jbas%total_orbits
 !    do b = 1, jbas%total_orbits
  do iii = 1, 50   
     call random_number(vv)
     call random_number(yy)
   
     a = ceiling(vv*(AA%Nsp))
     b = ceiling(yy*(AA%Nsp))
        
        val = scalar_tensor_1body_comm(AA,BB,a,b,jbas) 
        
        if (abs(val-f_tensor_elem(a,b,OUT,jbas)) > 1e-10) then
           print*, 'at: ',a,b
           print*, val, f_tensor_elem(a,b,OUT,jbas)
           STOP 'ONE BODY FAILURE'  
        end if 
        
        print*, 'success:', a,b
     !end do 
  end do 

 !do a = 12, jbas%total_orbits
     
  iii = 0 
  do while (iii < 15)  
     call random_number(vv)
     call random_number(xx)
     call random_number(yy)
     call random_number(zz)
   
     a = ceiling(vv*AA%Nsp)
     b = ceiling(xx*AA%Nsp)
     c = ceiling(yy*AA%Nsp)
     d = ceiling(zz*AA%Nsp)
     
     ja = jbas%jj(a) 
!     do b = 7, jbas%total_orbits
        jb = jbas%jj(b)
        
        PAR = mod(jbas%ll(a) + jbas%ll(b),2) 
        TZ = jbas%itzp(a) + jbas%itzp(b) 
        
 !       do c = 1, jbas%total_orbits
           jc = jbas%jj(c)
  !         do d = 1, jbas%total_orbits
              jd = jbas%jj(d) 
              
              if (PAR .ne. mod(jbas%ll(c) + jbas%ll(d)+BB%dpar/2,2)) cycle 
              if ( TZ .ne.  jbas%itzp(c) + jbas%itzp(d) ) cycle
              iii = iii+1 
              j1min = abs(ja-jb) 
              j1max = ja+jb 
              j2min = abs(jc-jd) 
              j2max = jc+jd
              
              do J1 = j1min,j1max,2
                 do J2 = j2min,j2max,2
                    
                    if (.not. (triangle(J1,J2,rank))) cycle
                    
                    val = scalar_tensor_2body_comm(AA,BB,a,b,c,d,J1,J2,jbas)
                  
                    if (abs(val-tensor_elem(a,b,c,d,J1,J2,OUT,jbas)) > 1e-8) then
                       print*, 'at:',a,b,c,d, 'J:', J1,J2 ,val,tensor_elem(a,b,c,d,J1,J2,OUT,jbas)                
            
                       STOP 'TWO BODY FAILURE'  
                    end if 
                 end do 
              end do
              
              print*, 'success:', a,b,c,d
   !        end do
    !    end do
     !end do
  end do
  
  print*, ' COMMUTATOR EXPRESSIONS CONFIRMED '
  
end subroutine
!============================================================
!============================================================
subroutine test_EOM_scalar_tensor_commutator(jbas,h1,h2,rank,dpar) 
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: AA,BB,OUT,w1,w2
  type(ex_pandya_mat) :: BBCC,WCC
  type(ex_cc_mat) :: AACC 
  integer :: a,b,c,d,g,q,ja,jb,jc,jd,j1min,j1max,dpar
  integer :: j2min,j2max,PAR,TZ,J1,J2,ax,bx,cx,dx,iii
  integer,intent(in) :: h1,h2,rank
  real(8) :: val,t1,t2,t3,t4,omp_get_wtime
  real(8) :: vv,xx,yy,zz
  
  call seed_random_number
  
  BB%rank = rank
  BB%dpar = dpar
  BB%pphh_ph = .false.
  AA%rank = 0
  call allocate_blocks(jbas,AA)
  call allocate_tensor(jbas,BB,AA)

  BB%herm = 1 
  call construct_random_rank0(AA,h1,jbas) 
  call construct_random_rankX(BB,h2,jbas) 
 
  call duplicate_sq_op(BB,OUT)
  call duplicate_sq_op(BB,w1) !workspace
  call duplicate_sq_op(BB,w2) !workspace
 
  call init_ph_mat(AA,AACC,jbas) ! cross coupled ME
  call init_ph_mat(BB,BBCC,jbas) !cross coupled ME
  call init_ph_wkspc(BBCC,WCC) !
  

  do q = 1, BB%nblocks
     BB%tblck(q)%lam(1) = 1
     OUT%tblck(q)%lam(1) = 1
     do g = 1, 9
        if ( g == 3) cycle
        if ( g == 7) cycle
        BB%tblck(q)%tgam(g)%X = 0.d0 
     end do 
  end do 
  
  OUT%herm = 1
  
  print*, 'TESTING EOM SCALAR-TENSOR COMMUTATORS' 
 ! t1 = OMP_get_wtime()
  call EOM_generalized_pandya(BB,BBCC,jbas)
 ! t2 = OMP_get_wtime()
  call calculate_cross_coupled_pphh(AA,AACC,jbas) 
  
  call EOM_TS_commutator_111(AA,BB,OUT,jbas) 
  call EOM_TS_commutator_121(AA,BB,OUT,jbas)
  call EOM_TS_commutator_211(AACC,BB,OUT,jbas) 
  call EOM_TS_commutator_122(AA,BB,OUT,jbas)
  call EOM_TS_commutator_212(AA,BB,OUT,jbas)
  
  call EOM_TS_commutator_222_pp_hh(AA,BB,OUT,w1,w2,jbas)
  
  call EOM_TS_commutator_221(w1,w2,AA%herm*BB%herm,OUT,jbas)
  ! t4 = OMP_get_wtime()
  call EOM_TS_commutator_222_ph(AACC,BBCC,OUT,WCC,jbas)
  ! t3 = OMP_get_wtime()
  
  print*, 'time:', t3-t1,t2-t1,t3-t4
!goto 12
 ! do ax = 1, AA%Nsp-AA%belowEF
   
  do iii = 1, 50   
     call random_number(vv)
     call random_number(yy)
   
     ax = ceiling(vv*(AA%Nsp-AA%belowEF))
     bx = ceiling(yy*(AA%belowEF))
     
     a = jbas%parts(ax)
  !   do bx = 1, AA%belowEF
        b = jbas%holes(bx)
        
        val = EOM_scalar_tensor_1body_comm(AA,BB,a,b,jbas) 
        

        if (abs(val-f_tensor_elem(a,b,OUT,jbas)) > 1e-10) then
           print*, 'at: ',a,b
           print*, val, f_tensor_elem(a,b,OUT,jbas)
           STOP 'ONE BODY FAILURE'  
        end if 
        
        print*, 'success:', a,b
!     end do 
  end do 

  
  iii = 0 
  do while (iii < 15) 
     call random_number(vv)
     call random_number(xx)
     call random_number(yy)
     call random_number(zz)
   
     ax = ceiling(vv*(AA%Nsp-AA%belowEF))
     bx = ceiling(xx*(AA%Nsp-AA%belowEF))
     cx = ceiling(yy*(AA%belowEF))
     dx = ceiling(zz*(AA%belowEF))
     
     ! do ax = 1, AA%Nsp-AA%belowEF
     a = jbas%parts(ax)
     ja = jbas%jj(a) 
  !   do bx = 1, AA%Nsp-AA%belowEF
        b = jbas%parts(bx)
        jb = jbas%jj(b)
        
        PAR = mod(jbas%ll(a) + jbas%ll(b),2) 
        TZ = jbas%itzp(a) + jbas%itzp(b) 
        
   !     do cx = 1,AA%belowEF
           c = jbas%holes(cx)
           jc = jbas%jj(c)
    !       do dx = 1,AA%belowEF
              d = jbas%holes(dx)
              jd = jbas%jj(d) 
              
              if (PAR .ne. mod(jbas%ll(c) + jbas%ll(d)+BB%dpar/2,2)) cycle 
              if ( TZ .ne.  jbas%itzp(c) + jbas%itzp(d) ) cycle
              iii = iii + 1
              
              j1min = abs(ja-jb) 
              j1max = ja+jb 
              j2min = abs(jc-jd) 
              j2max = jc+jd
              
              do J1 = j1min,j1max,2
                 do J2 = j2min,j2max,2
                    
                    if (.not. (triangle(J1,J2,rank))) cycle
                    
                    val = EOM_scalar_tensor_2body_comm(AA,BB,a,b,c,d,J1,J2,jbas)
!                    print*, a,b,c,d, 'J:', J1,J2 ,val,tensor_elem(a,b,c,d,J1,J2,OUT,jbas)  
                    if (abs(val-tensor_elem(a,b,c,d,J1,J2,OUT,jbas)) > 1e-8) then
                       print*, 'at:',a,b,c,d, 'J:', J1,J2 ,val,tensor_elem(a,b,c,d,J1,J2,OUT,jbas)  
                       print*, tensor_elem(c,d,a,b,J2,J1,OUT,jbas)  
            
                       STOP 'TWO BODY FAILURE'  
                    end if 
                 end do 
              end do
              
              print*, 'success:', a,b,c,d
        !   end do
       ! end do
     !end do
  end do
  
  print*, ' COMMUTATOR EXPRESSIONS CONFIRMED '
  
end subroutine
!============================================================
!============================================================
real(8) function scalar_scalar_1body_comm(AA,BB,a,b,jbas) 
  !returns [AA^0, BB^0]_{ab} 
  ! uses brute force method. 
  implicit none 
  
  integer :: a,b,i,j,k
  integer :: ja,jb,jj,ji,Jtot,JTM,totorb
  type(spd) :: jbas
  type(sq_op) :: AA,BB 
  real(8) :: sm 
  
  sm = 0.d0 
  JTM = jbas%jtotal_max*2
  totorb = jbas%total_orbits
  
  ja = jbas%jj(a) 
  jb = jbas%jj(b) 
  
  do i = 1, totorb
     
     sm = sm + f_elem(a,i,AA,jbas) * f_elem(i,b,BB,jbas) &
          - f_elem(a,i,BB,jbas) * f_elem(i,b,AA,jbas)
  
  end do 
  
  
  If (ja == jb) then 
     ! this if statement is built into the expression.
     do i = 1, totorb
        do j = 1, totorb
        
           do Jtot = 0,JTM,2
              
              sm = sm + (jbas%con(i) -jbas%con(j) ) * (Jtot+1.d0) /(ja +1.d0) * &
                   ( f_elem(i,j,AA,jbas) * v_elem(j,a,i,b,Jtot,BB,jbas) - &
                   f_elem(i,j,BB,jbas) * v_elem(j,a,i,b,Jtot,AA,jbas) )
           end do
        end do
     end do
  
  
     do i =  1, totorb
        do j =  1, totorb
           do k =  1, totorb
              do Jtot = 0,JTM,2 

                 sm = sm + (jbas%con(i)*jbas%con(j)*(1-jbas%con(k)) + &
                      (1-jbas%con(i))*(1-jbas%con(j))*jbas%con(k)) * (Jtot + 1.d0) &
                      / (ja + 1.d0) * ( v_elem(a,k,i,j,Jtot,AA,jbas)*v_elem(i,j,b,k,Jtot,BB,jbas) &
                      - v_elem(a,k,i,j,Jtot,BB,jbas)*v_elem(i,j,b,k,Jtot,AA,jbas))/2.d0
           
              end do
           end do
        end do
     end do
  end if 
  scalar_scalar_1body_comm = sm 
  
end function 
!============================================================
!============================================================
real(8) function EOM_scalar_scalar_1body_comm(AA,BB,a,b,jbas) 
  !returns [AA^0, BB^0]_{ab} 
  ! uses brute force method. 
  implicit none 
  
  integer :: a,b,i,j,k
  integer :: ja,jb,jj,ji,Jtot,JTM,totorb
  type(spd) :: jbas
  type(sq_op) :: AA,BB 
  real(8) :: sm 
  
  sm = 0.d0 
  JTM = jbas%jtotal_max*2
  totorb = jbas%total_orbits
  
  ja = jbas%jj(a) 
  jb = jbas%jj(b) 
  
  do i = 1, totorb
     
     sm = sm + f_elem(a,i,AA,jbas) * ph_elem(i,b,BB,jbas) &
          - ph_elem(a,i,BB,jbas) * f_elem(i,b,AA,jbas)
  
  end do 
  
  
  If (ja == jb) then 
     ! this if statement is built into the expression.
     do i = 1, totorb
        do j = 1, totorb
        
           do Jtot = 0,JTM,2
              
              sm = sm + (jbas%con(i) -jbas%con(j) ) * (Jtot+1.d0) /(ja +1.d0) * &
                   ( f_elem(i,j,AA,jbas) * pphh_elem(j,a,i,b,Jtot,BB,jbas) - &
                   ph_elem(i,j,BB,jbas) * v_elem(j,a,i,b,Jtot,AA,jbas) )
           end do
        end do
     end do
  
  
     do i =  1, totorb
        do j =  1, totorb
           do k =  1, totorb
              do Jtot = 0,JTM,2 

                 sm = sm + (jbas%con(i)*jbas%con(j)*(1-jbas%con(k)) + &
                      (1-jbas%con(i))*(1-jbas%con(j))*jbas%con(k)) * (Jtot + 1.d0) &
                      / (ja + 1.d0) * ( v_elem(a,k,i,j,Jtot,AA,jbas)*pphh_elem(i,j,b,k,Jtot,BB,jbas) &
                      - pphh_elem(a,k,i,j,Jtot,BB,jbas)*v_elem(i,j,b,k,Jtot,AA,jbas))/2.d0
           
              end do
           end do
        end do
     end do
  end if 
  EOM_scalar_scalar_1body_comm = sm 
  
end function 
!============================================================
!============================================================
real(8) function scalar_scalar_0body_comm(AA,BB,jbas) 
  !returns  [AA^0, BB^0]_{0}
  ! uses brute force method. 
  implicit none 
  
  integer :: a,b,i,j,k,l
  integer :: ja,jb,jj,ji,Jtot,JTM,totorb
  type(spd) :: jbas
  type(sq_op) :: AA,BB 
  real(8) :: sm 
  
  sm = 0.d0 
  JTM = jbas%jtotal_max*2
  totorb = jbas%total_orbits
    
  do i = 1, totorb
     ji = jbas%jj(i)
     do j = 1, totorb
        
        sm = sm + (jbas%con(i) - jbas%con(j)) * (ji + 1.d0) * &
             f_elem(i,j,AA,jbas) * f_elem(j,i,BB,jbas) 
     end do 
  end do 

  do i = 1, totorb
     do j = 1, totorb
        do k = 1, totorb
           do l = 1, totorb
              do Jtot = 0,JTM,2
                 
                 sm = sm + (jbas%con(i) * jbas%con(j) * (1-jbas%con(k)) * (1-jbas%con(l)) &
                      - jbas%con(k) * jbas%con(l) * (1-jbas%con(i)) * (1-jbas%con(j)) ) * &
                      (Jtot+1.d0) * v_elem(i,j,k,l,Jtot,AA,jbas)*v_elem(k,l,i,j,Jtot,BB,jbas)/4.d0
              
              end do
           end do 
        end do
     end do
  end do
  scalar_scalar_0body_comm = sm 
  
end function 
!============================================================
!============================================================
real(8) function EOM_scalar_scalar_0body_comm(AA,BB,jbas) 
  !returns  [AA^0, BB^0]_{0}
  ! uses brute force method. 
  implicit none 
  
  integer :: a,b,i,j,k,l
  integer :: ja,jb,jj,ji,Jtot,JTM,totorb
  type(spd) :: jbas
  type(sq_op) :: AA,BB 
  real(8) :: sm 
  
  sm = 0.d0 
  JTM = jbas%jtotal_max*2
  totorb = jbas%total_orbits
    
  do i = 1, totorb
     ji = jbas%jj(i)
     do j = 1, totorb
        
        sm = sm + (jbas%con(i) - jbas%con(j)) * (ji + 1.d0) * &
             f_elem(i,j,AA,jbas) * ph_elem(j,i,BB,jbas) 
     end do 
  end do 

  do i = 1, totorb
     do j = 1, totorb
        do k = 1, totorb
           do l = 1, totorb
              do Jtot = 0,JTM,2
                 
                 sm = sm + (jbas%con(i) * jbas%con(j) * (1-jbas%con(k)) * (1-jbas%con(l)) &
                      - jbas%con(k) * jbas%con(l) * (1-jbas%con(i)) * (1-jbas%con(j)) ) * &
                      (Jtot+1.d0) * v_elem(i,j,k,l,Jtot,AA,jbas)*pphh_elem(k,l,i,j,Jtot,BB,jbas)/4.d0
              
              end do
           end do 
        end do
     end do
  end do
  EOM_scalar_scalar_0body_comm = sm 
  
end function 
!============================================================
!============================================================
real(8) function scalar_scalar_2body_comm(AA,BB,a,b,c,d,Jtot,jbas) 
  !returns  [AA^0, BB^0]_{0}
  ! uses brute force method. 
  implicit none 
  
  integer :: a,b,c,d,i,j,k,l,J1,J2,ji,jj
  integer :: ja,jb,jc,jd,Jtot,JTM,totorb
  type(spd) :: jbas
  type(sq_op) :: AA,BB 
  real(8) :: sm,coef9
  
  sm = 0.d0 
  JTM = jbas%jtotal_max*2
  totorb = jbas%total_orbits
  
  ja = jbas%jj(a) 
  jb = jbas%jj(b) 
  jc = jbas%jj(c)
  jd = jbas%jj(d)

  do i = 1, totorb
     sm = sm + f_elem(a,i,AA,jbas) * v_elem( i,b,c,d,Jtot,BB,jbas) &
          + f_elem(b,i,AA,jbas) * v_elem( a,i,c,d,Jtot,BB,jbas) &
          - f_elem(i,c,AA,jbas) * v_elem( a,b,i,d,Jtot,BB,jbas) &
          - f_elem(i,d,AA,jbas) * v_elem( a,b,c,i,Jtot,BB,jbas) 
     
     sm = sm - f_elem(a,i,BB,jbas) * v_elem( i,b,c,d,Jtot,AA,jbas) &
          - f_elem(b,i,BB,jbas) * v_elem( a,i,c,d,Jtot,AA,jbas) &
          + f_elem(i,c,BB,jbas) * v_elem( a,b,i,d,Jtot,AA,jbas) &
          + f_elem(i,d,BB,jbas) * v_elem( a,b,c,i,Jtot,AA,jbas) 
  
  end do
  

  do i = 1, totorb
     do j = 1, totorb
        
        sm = sm + 0.5*(1- jbas%con(i) - jbas%con(j)) *&
             (v_elem(a,b,i,j,Jtot,AA,jbas)*v_elem(i,j,c,d,Jtot,BB,jbas)   &
             - v_elem(a,b,i,j,Jtot,BB,jbas)*v_elem(i,j,c,d,Jtot,AA,jbas)) 
     end do
  end do

  do i = 1, totorb
     ji =jbas%jj(i)
     do j = 1,totorb
        jj = jbas%jj(j) 
        
        if ((jbas%con(i)-jbas%con(j)) == 0) cycle 
        
        do J1 = 0, JTM,2
           do J2 = 0, JTM,2 
              
              sm = sm + (jbas%con(i)-jbas%con(j)) *  ( &  
                   
                   (-1)** ((J1+J2 + ja-jc)/2) * (J1+1.d0) * (J2+1.d0) &
                   * coef9(jj,J1,ja,J2,ji,jb,jc,jd,Jtot) * v_elem(j,a,i,d,J1,AA,jbas) &
                   * v_elem(i,b,j,c,J2,BB,jbas) &
                   
                   - (-1)** ((J1+J2 + jb-jc)/2) * (J1+1.d0) * (J2+1.d0) &
                   * coef9(jj,J1,jb,J2,ji,ja,jc,jd,Jtot) * v_elem(j,b,i,d,J1,AA,jbas) &
                   * v_elem(i,a,j,c,J2,BB,jbas) *(-1)**((ja+jb-Jtot)/2) &
                   
                   - (-1)** ((J1+J2 + ja-jd)/2) * (J1+1.d0) * (J2+1.d0) &
                   * coef9(jj,J1,ja,J2,ji,jb,jd,jc,Jtot) * v_elem(j,a,i,c,J1,AA,jbas) &
                   * v_elem(i,b,j,d,J2,BB,jbas) * (-1)**((jc+jd-Jtot)/2) &
                   
                    + (-1)** ((J1+J2 + jb-jd)/2) * (J1+1.d0) * (J2+1.d0) &
                    * coef9(jj,J1,jb,J2,ji,ja,jd,jc,Jtot) * v_elem(j,b,i,c,J1,AA,jbas) &
                   * v_elem(i,a,j,d,J2,BB,jbas) *(-1)**((ja+jb+jc+jd)/2)  )
              ! FUCK
           end do
        end do
        
       
        
     end do
  end do
  
  
  ! do J2= 0,4,2
  !    do i = 1, totorb
  !       ji =jbas%jj(i)
  !       do j = 1,totorb
  !          jj = jbas%jj(j) 
           
  !          if ((jbas%con(i)-jbas%con(j)) == 0) cycle 

  !          sm = sm + (jbas%con(i)-jbas%con(j)) * (J2+1.d0)* ( & 
  !               sixj(ja,jb,Jtot,jc,jd,J2) * &
  !                ( Vpandya(a,d,i,j,J2,AA,jbas) * Vpandya(i,j,c,b,J2,BB,jbas) - &
  !               Vpandya(a,d,i,j,J2,BB,jbas) * Vpandya(i,j,c,b,J2,AA,jbas)  ) - &
  !               (-1)**((ja+jb-Jtot)/2) * sixj(jb,ja,Jtot,jc,jd,J2) * &
  !               ( Vpandya(b,d,i,j,J2,AA,jbas) * Vpandya(i,j,c,a,J2,BB,jbas) - &
  !               Vpandya(b,d,i,j,J2,BB,jbas) * Vpandya(i,j,c,a,J2,AA,jbas)  ) )
  !       end do
  !    end do
     
  ! end do
  
  scalar_scalar_2body_comm = sm 
  
end function 
!============================================================
!============================================================
real(8) function EOM_scalar_scalar_2body_comm(AA,BB,a,b,c,d,Jtot,jbas) 
  !returns  [AA^0, BB^0]_{0}
  ! uses brute force method. 
  implicit none 
  
  integer :: a,b,c,d,i,j,k,l,J1,J2,ji,jj
  integer :: ja,jb,jc,jd,Jtot,JTM,totorb
  type(spd) :: jbas
  type(sq_op) :: AA,BB 
  real(8) :: sm,coef9
  
  sm = 0.d0 
  JTM = jbas%jtotal_max*2
  totorb = jbas%total_orbits
  
  ja = jbas%jj(a) 
  jb = jbas%jj(b) 
  jc = jbas%jj(c)
  jd = jbas%jj(d)

  do i = 1, totorb
     sm = sm + f_elem(a,i,AA,jbas) * pphh_elem( i,b,c,d,Jtot,BB,jbas) &
          + f_elem(b,i,AA,jbas) * pphh_elem( a,i,c,d,Jtot,BB,jbas) &
          - f_elem(i,c,AA,jbas) * pphh_elem( a,b,i,d,Jtot,BB,jbas) &
          - f_elem(i,d,AA,jbas) * pphh_elem( a,b,c,i,Jtot,BB,jbas) 
     
     sm = sm - ph_elem(a,i,BB,jbas) * v_elem( i,b,c,d,Jtot,AA,jbas) &
          - ph_elem(b,i,BB,jbas) * v_elem( a,i,c,d,Jtot,AA,jbas) &
          + ph_elem(i,c,BB,jbas) * v_elem( a,b,i,d,Jtot,AA,jbas) &
          + ph_elem(i,d,BB,jbas) * v_elem( a,b,c,i,Jtot,AA,jbas) 
  
  end do
  

  do i = 1, totorb
     do j = 1, totorb
        
        sm = sm + 0.5*(1- jbas%con(i) - jbas%con(j)) *&
             (v_elem(a,b,i,j,Jtot,AA,jbas)*pphh_elem(i,j,c,d,Jtot,BB,jbas)   &
             - pphh_elem(a,b,i,j,Jtot,BB,jbas)*v_elem(i,j,c,d,Jtot,AA,jbas)) 
     end do
  end do

  do i = 1, totorb
     ji =jbas%jj(i)
     do j = 1,totorb
        jj = jbas%jj(j) 
        
        if ((jbas%con(i)-jbas%con(j)) == 0) cycle 
        do J1 = 0, JTM,2
           do J2 = 0, JTM,2 
              
              sm = sm + (jbas%con(i)-jbas%con(j)) *  ( &  
                   
                   (-1)** ((J1+J2 + ja-jc)/2) * (J1+1.d0) * (J2+1.d0) &
                   * coef9(jj,J1,ja,J2,ji,jb,jc,jd,Jtot) * v_elem(j,a,i,d,J1,AA,jbas) &
                   * pphh_elem(i,b,j,c,J2,BB,jbas) &
                   
                   - (-1)** ((J1+J2 + jb-jc)/2) * (J1+1.d0) * (J2+1.d0) &
                   * coef9(jj,J1,jb,J2,ji,ja,jc,jd,Jtot) * v_elem(j,b,i,d,J1,AA,jbas) &
                   * pphh_elem(i,a,j,c,J2,BB,jbas) *(-1)**((ja+jb-Jtot)/2) &
                   
                   - (-1)** ((J1+J2 + ja-jd)/2) * (J1+1.d0) * (J2+1.d0) &
                   * coef9(jj,J1,ja,J2,ji,jb,jd,jc,Jtot) * v_elem(j,a,i,c,J1,AA,jbas) &
                   * pphh_elem(i,b,j,d,J2,BB,jbas) * (-1)**((jc+jd-Jtot)/2) &
                   
                   + (-1)** ((J1+J2 + jb-jd)/2) * (J1+1.d0) * (J2+1.d0) &
                   * coef9(jj,J1,jb,J2,ji,ja,jd,jc,Jtot) * v_elem(j,b,i,c,J1,AA,jbas) &
                   * pphh_elem(i,a,j,d,J2,BB,jbas) *(-1)**((ja+jb+jc+jd)/2)  )
              
           end do
        end do
     end do
  end do
  
  EOM_scalar_scalar_2body_comm = sm 
  
end function 
!============================================================
!============================================================
real(8) function scalar_tensor_1body_comm(AA,BB,a,b,jbas) 
  !returns [AA^0, BB^X]_{ab} 
  ! uses brute force method. 
  implicit none 
  
  integer :: a,b,i,j,k,rank,J1,J2
  integer :: ja,jb,jj,ji,jk,Jtot,JTM,totorb
  type(spd) :: jbas
  type(sq_op) :: AA,BB 
  real(8) :: sm ,d6ji
  
  rank = BB%rank
  sm = 0.d0 
  JTM = jbas%jtotal_max*2
  totorb = jbas%total_orbits
  
  ja = jbas%jj(a) 
  jb = jbas%jj(b) 
  
  do i = 1, totorb
     
     sm = sm + f_elem(a,i,AA,jbas) * f_tensor_elem(i,b,BB,jbas) &
          - f_tensor_elem(a,i,BB,jbas) * f_elem(i,b,AA,jbas)
  
  end do 
  

  do i = 1, totorb
     ji = jbas%jj(i)
     do j = 1, totorb
        jj = jbas%jj(j) 
        
        do J1 = 0,JTM,2
           do J2 = 0, JTM,2 
                 
              sm = sm + (jbas%con(i) -jbas%con(j) )* &
                   f_elem(i,j,AA,jbas) * tensor_elem(j,a,i,b,J1,J2,BB,jbas) &
                   * sqrt( (J1+1.d0)*(J2+1.d0) ) * (-1)**(( J1+ rank +jb +ji)/2) * &
                   d6ji(J1,J2,rank,jb,ja,ji) 
           
           end do
        end do
     end do
  end do
  
  do i = 1, totorb
     ji = jbas%jj(i)
     do j = 1, totorb
        jj = jbas%jj(j) 
        do k = 1, totorb
           jk =jbas%jj(k)
           
           do J1 = 0,JTM,2
              do J2 = 0, JTM,2 
                 
                 sm = sm + 0.5*(jbas%con(k)*jbas%con(j)*(1-jbas%con(i)) + &
                      (1-jbas%con(k))*(1-jbas%con(j))*jbas%con(i) )* &
                      ( v_elem(i,a,j,k,J1,AA,jbas) * tensor_elem(j,k,i,b,J1,J2,BB,jbas) &
                      - tensor_elem(i,a,j,k,J1,J2,BB,jbas) * v_elem(j,k,i,b,J2,AA,jbas) ) &
                      * sqrt( (J1+1.d0)*(J2+1.d0) ) * (-1)**(( J1+ rank +jb +ji)/2) * &
                      d6ji(J1,J2,rank,jb,ja,ji) 
           
              end do
           end do
        end do
     end do
  end do
  

  do i =  1, totorb
     ji =jbas%jj(i)
     do j =  1, totorb
        jj = jbas%jj(j) 
        
        do Jtot = 0,JTM,2 

           sm = sm - (jbas%con(i) - jbas%con(j)) * (Jtot+1.d0) * &
                (-1) ** ((ja-ji+Jtot)/2) *d6ji(ja,jb,rank,jj,ji,Jtot) * &
                f_tensor_elem(j,i,BB,jbas) * v_elem(i,a,j,b,Jtot,AA,jbas)
           

           
        end do
     end do
  end do

  
  scalar_tensor_1body_comm = sm 
  
end function
!============================================================
!============================================================
real(8) function EOM_scalar_tensor_1body_comm(AA,BB,a,b,jbas) 
  !returns [AA^0, BB^X]_{ab} 
  ! uses brute force method. 
  implicit none 
  
  integer :: a,b,i,j,k,rank,J1,J2
  integer :: ja,jb,jj,ji,jk,Jtot,JTM,totorb
  type(spd) :: jbas
  type(sq_op) :: AA,BB 
  real(8) :: sm ,d6ji,sx
  
  rank = BB%rank
  sm = 0.d0 
  JTM = jbas%jtotal_max*2
  totorb = jbas%total_orbits
  
  ja = jbas%jj(a) 
  jb = jbas%jj(b) 
  
  do i = 1, totorb
     
     sm = sm + f_elem(a,i,AA,jbas) * ph_tensor_elem(i,b,BB,jbas) &
          - ph_tensor_elem(a,i,BB,jbas) * f_elem(i,b,AA,jbas)
  
  end do 
  

  do i = 1, totorb
     ji = jbas%jj(i)
     do j = 1, totorb
        jj = jbas%jj(j) 
        
        do J1 = 0,JTM,2
           do J2 = 0, JTM,2 
                 
              sm = sm + (jbas%con(i) -jbas%con(j) )* &
                   f_elem(i,j,AA,jbas) * pphh_tensor_elem(j,a,i,b,J1,J2,BB,jbas) &
                   * sqrt( (J1+1.d0)*(J2+1.d0) ) * (-1)**(( J1+ rank +jb +ji)/2) * &
                   d6ji(J1,J2,rank,jb,ja,ji) 
           
           end do
        end do
     end do
  end do
  
  do i = 1, totorb
     ji = jbas%jj(i)      
     do j = 1, totorb
        jj = jbas%jj(j) 
        do k = 1, totorb
           jk =jbas%jj(k)
           
           do J1 = 0,JTM,2
              do J2 = 0,JTM,2
           
                 sm = sm + 0.5*(jbas%con(k)*jbas%con(j)*(1-jbas%con(i)) + &
                      (1-jbas%con(k))*(1-jbas%con(j))*jbas%con(i) )* &
                      ( v_elem(i,a,j,k,J1,AA,jbas) * pphh_tensor_elem(j,k,i,b,J1,J2,BB,jbas) &
                      - pphh_tensor_elem(i,a,j,k,J1,J2,BB,jbas) * v_elem(j,k,i,b,J2,AA,jbas) ) &
                      * sqrt( (J1+1.d0)*(J2+1.d0) ) * (-1)**(( J1+ rank +jb +ji)/2) * &
                      d6ji(J1,J2,rank,jb,ja,ji) 
           
              end do
           end do
        end do
     end do
  end do
  

  do i =  1, totorb
     ji =jbas%jj(i)
     do j =  1, totorb
        jj = jbas%jj(j) 
        
        do Jtot = 0,JTM,2 

           sm = sm - (jbas%con(i) - jbas%con(j)) * (Jtot+1.d0) * &
                (-1) ** ((ja-ji+Jtot)/2) *d6ji(ja,jb,rank,jj,ji,Jtot) * &
                ph_tensor_elem(j,i,BB,jbas) * v_elem(i,a,j,b,Jtot,AA,jbas)
           

           
        end do
     end do
  end do

  
  EOM_scalar_tensor_1body_comm = sm 
  
end function 
!============================================================
!============================================================
real(8) function scalar_tensor_2body_comm(AA,BB,a,b,c,d,J1,J2,jbas) 
  !returns  [AA^0, BB^0]_{0}
  ! uses brute force method. 
  implicit none 
  
  integer :: a,b,c,d,i,j,k,l,J1,J2,ji,jj,J3,J4,J5,jx
  integer :: ja,jb,jc,jd,Jtot,JTM,totorb,rank
  type(spd) :: jbas
  type(sq_op) :: AA,BB 
  real(8) :: sm,coef9,d6ji

  rank = BB%rank  
  sm = 0.d0 
  JTM = jbas%jtotal_max*2
  totorb = jbas%total_orbits
  
  ja = jbas%jj(a) 
  jb = jbas%jj(b) 
  jc = jbas%jj(c)
  jd = jbas%jj(d)

   do i = 1, totorb
      ji = jbas%jj(i)

     sm = sm + f_elem(a,i,AA,jbas) * tensor_elem( i,b,c,d,J1,J2,BB,jbas) &
          + f_elem(b,i,AA,jbas) * tensor_elem( a,i,c,d,J1,J2,BB,jbas) &
          - f_elem(i,c,AA,jbas) * tensor_elem( a,b,i,d,J1,J2,BB,jbas) &
          - f_elem(i,d,AA,jbas) * tensor_elem( a,b,c,i,J1,J2,BB,jbas) 
     
          
     sm = sm - f_tensor_elem(a,i,BB,jbas) * v_elem( i,b,c,d,J2,AA,jbas) &
          * d6ji(ji,jb,J2,J1,rank,ja) * (-1)**((ja+jb+rank-J2)/2) * &
          sqrt( (J1+1.d0) * (J2+1.d0) ) &
     
          + f_tensor_elem(b,i,BB,jbas) * v_elem( i,a,c,d,J2,AA,jbas) &
          * d6ji(ji,ja,J2,J1,rank,jb) * (-1)**((J1+J2+rank)/2) * &
          sqrt( (J1+1.d0) * (J2+1.d0) ) &

          - f_tensor_elem(i,c,BB,jbas) * v_elem( a,b,d,i,J1,AA,jbas) &
          *d6ji(ji,jc,rank,J2,J1,jd) * (-1)**((J1+J2+rank)/2) *  &
          sqrt( (J1+1.d0) * (J2+1.d0) ) &
         
          
          + f_tensor_elem(i,d,BB,jbas) * v_elem( a,b,c,i,J1,AA,jbas) &
          *d6ji( ji,jd,rank,J2,J1,jc) * (-1)**((jc+jd-J1+rank)/2) * &
          sqrt( (J1+1.d0) * (J2+1.d0) )
     
   end do
  

  do i = 1, totorb
     do j = 1, totorb
        

        sm = sm + 0.5*(1- jbas%con(i) - jbas%con(j)) *&
             (v_elem(a,b,i,j,J1,AA,jbas)*tensor_elem(i,j,c,d,J1,J2,BB,jbas)   &
             - tensor_elem(a,b,i,j,J1,J2,BB,jbas)*v_elem(i,j,c,d,J2,AA,jbas)) 
     end do
  end do
 
!!$OMP PARALLEL DO PRIVATE( ji,jj,i,j,J3,J4,J5,jx) SHARED(AA,BB) REDUCTION(+:sm)
  do i = 1, totorb
     ji =jbas%jj(i)
     do j = 1,totorb
        jj = jbas%jj(j) 
        
        if ((jbas%con(i)-jbas%con(j)) == 0) cycle 
        do J3 = 0, JTM,2
           do J4 = 0, JTM,2 
              do J5 = 0,JTM,2
                 do jx = 1,JTM,2
                    
                    sm = sm + (jbas%con(i)-jbas%con(j)) *  ( &  
                   
                         (-1)** ((J2+J3 + jc - ji )/2) * sqrt( (J1+1.d0) * (J2+1.d0)  &
                         * (J4+1.d0) * (J5+1.d0) ) * (jx+1.d0) * (J3+1.d0)   &
                    * coef9(jj,J3,ja,J4,ji,jb,jx,jd,J1) * d6ji( jj,J4,jx,rank,jc,J5) * &
                         d6ji(J1,jx,jd,jc,J2,rank) * v_elem(a,j,d,i,J3,AA,jbas) *&
                         tensor_elem(i,b,j,c,J4,J5,BB,jbas) &
                   
                         - (-1)** ((J2+J3 + jc - ji )/2) * sqrt( (J1+1.d0) * (J2+1.d0)  &
                         * (J4+1.d0) * (J5+1.d0) ) * (jx+1.d0) * (J3+1.d0)  &
                    * coef9(jj,J3,jb,J4,ji,ja,jx,jd,J1) * d6ji( jj,J4,jx,rank,jc,J5) * &
                         d6ji(J1,jx,jd,jc,J2,rank) * v_elem(b,j,d,i,J3,AA,jbas) *&
                         tensor_elem(i,a,j,c,J4,J5,BB,jbas) *(-1)**((ja+jb-J1)/2) &
                   
                         - (-1)** ((J2+J3 + jd - ji )/2) * sqrt( (J1+1.d0) * (J2+1.d0)  &
                         * (J4+1.d0) * (J5+1.d0) ) * (jx+1.d0) * (J3+1.d0)  &
                    * coef9(jj,J3,ja,J4,ji,jb,jx,jc,J1) * d6ji( jj,J4,jx,rank,jd,J5) * &
                         d6ji(J1,jx,jc,jd,J2,rank) * v_elem(a,j,c,i,J3,AA,jbas) *&
                         tensor_elem(i,b,j,d,J4,J5,BB,jbas) * (-1)**((jc+jd-J2)/2) &
                   
                         + (-1)** ((J2+J3 + jd - ji )/2) * sqrt( (J1+1.d0) * (J2+1.d0)  &
                         * (J4+1.d0) * (J5+1.d0) ) * (jx+1.d0) * (J3+1.d0)  &
                    * coef9(jj,J3,jb,J4,ji,ja,jx,jc,J1) * d6ji( jj,J4,jx,rank,jd,J5) * &
                         d6ji(J1,jx,jc,jd,J2,rank) * v_elem(b,j,c,i,J3,AA,jbas) *&
                         tensor_elem(i,a,j,d,J4,J5,BB,jbas)  *(-1)**((ja+jb+jc+jd+J1+J2)/2)  )
              
                 end do
              end do
           end do
        end do
     end do
  end do
!!$OMP END PARALLEL DO

  scalar_tensor_2body_comm = sm 
  
end function 
!============================================================
!============================================================
real(8) function EOM_scalar_tensor_2body_comm(AA,BB,a,b,c,d,J1,J2,jbas) 
  !returns  [AA^0, BB^0]_{0}
  ! uses brute force method. 
  implicit none 
  
  integer :: a,b,c,d,i,j,k,l,J1,J2,ji,jj,J3,J4,J5,jx
  integer :: ja,jb,jc,jd,Jtot,JTM,totorb,rank
  type(spd) :: jbas
  type(sq_op) :: AA,BB 
  real(8) :: sm,coef9,d6ji

  rank = BB%rank  
  sm = 0.d0 
  JTM = jbas%jtotal_max*2
  totorb = jbas%total_orbits
  
  ja = jbas%jj(a) 
  jb = jbas%jj(b) 
  jc = jbas%jj(c)
  jd = jbas%jj(d)

    do i = 1, totorb
       ji = jbas%jj(i)

      sm = sm + f_elem(a,i,AA,jbas) * pphh_tensor_elem( i,b,c,d,J1,J2,BB,jbas) &
           + f_elem(b,i,AA,jbas) * pphh_tensor_elem( a,i,c,d,J1,J2,BB,jbas) &
           - f_elem(i,c,AA,jbas) * pphh_tensor_elem( a,b,i,d,J1,J2,BB,jbas) &
           - f_elem(i,d,AA,jbas) * pphh_tensor_elem( a,b,c,i,J1,J2,BB,jbas) 
     
          
     sm = sm - ph_tensor_elem(a,i,BB,jbas) * v_elem( i,b,c,d,J2,AA,jbas) &
          * d6ji(ji,jb,J2,J1,rank,ja) * (-1)**((ja+jb+rank-J2)/2) * &
          sqrt( (J1+1.d0) * (J2+1.d0) ) &
     
          + ph_tensor_elem(b,i,BB,jbas) * v_elem( i,a,c,d,J2,AA,jbas) &
          * d6ji(ji,ja,J2,J1,rank,jb) * (-1)**((J1+J2+rank)/2) * &
          sqrt( (J1+1.d0) * (J2+1.d0) ) &

          - ph_tensor_elem(i,c,BB,jbas) * v_elem( a,b,d,i,J1,AA,jbas) &
          *d6ji(ji,jc,rank,J2,J1,jd) * (-1)**((J1+J2+rank)/2) *  &
          sqrt( (J1+1.d0) * (J2+1.d0) ) &
         
          
          + ph_tensor_elem(i,d,BB,jbas) * v_elem( a,b,c,i,J1,AA,jbas) &
          *d6ji( ji,jd,rank,J2,J1,jc) * (-1)**((jc+jd-J1+rank)/2) * &
          sqrt( (J1+1.d0) * (J2+1.d0) )
     
    end do
  

  do i = 1, totorb
     do j = 1, totorb
        

        sm = sm + 0.5*(1- jbas%con(i) - jbas%con(j)) *&
             (v_elem(a,b,i,j,J1,AA,jbas)*pphh_tensor_elem(i,j,c,d,J1,J2,BB,jbas)   &
             - pphh_tensor_elem(a,b,i,j,J1,J2,BB,jbas)*v_elem(i,j,c,d,J2,AA,jbas)) 
     end do
  end do
 
!!$OMP PARALLEL DO PRIVATE( ji,jj,i,j,J3,J4,J5,jx) SHARED(AA,BB) REDUCTION(+:sm)
  do i = 1, totorb
     ji =jbas%jj(i)
     do j = 1,totorb
        jj = jbas%jj(j) 
        
        if ((jbas%con(i)-jbas%con(j)) == 0) cycle 
        do J3 = 0, JTM,2
           do J4 = 0, JTM,2 
              do J5 = 0,JTM,2
                 do jx = 1,JTM,2
                    
                    sm = sm + (jbas%con(i)-jbas%con(j)) *  ( &  
                   
                         (-1)** ((J2+J3 + jc - ji )/2) * sqrt( (J1+1.d0) * (J2+1.d0)  &
                         * (J4+1.d0) * (J5+1.d0) ) * (jx+1.d0) * (J3+1.d0)   &
                    * coef9(jj,J3,ja,J4,ji,jb,jx,jd,J1) * d6ji( jj,J4,jx,rank,jc,J5) * &
                         d6ji(J1,jx,jd,jc,J2,rank) * v_elem(a,j,d,i,J3,AA,jbas) *&
                         pphh_tensor_elem(i,b,j,c,J4,J5,BB,jbas) &
                   
                         - (-1)** ((J2+J3 + jc - ji )/2) * sqrt( (J1+1.d0) * (J2+1.d0)  &
                         * (J4+1.d0) * (J5+1.d0) ) * (jx+1.d0) * (J3+1.d0)  &
                    * coef9(jj,J3,jb,J4,ji,ja,jx,jd,J1) * d6ji( jj,J4,jx,rank,jc,J5) * &
                         d6ji(J1,jx,jd,jc,J2,rank) * v_elem(b,j,d,i,J3,AA,jbas) *&
                         pphh_tensor_elem(i,a,j,c,J4,J5,BB,jbas) *(-1)**((ja+jb-J1)/2) &
                   
                         - (-1)** ((J2+J3 + jd - ji )/2) * sqrt( (J1+1.d0) * (J2+1.d0)  &
                         * (J4+1.d0) * (J5+1.d0) ) * (jx+1.d0) * (J3+1.d0)  &
                    * coef9(jj,J3,ja,J4,ji,jb,jx,jc,J1) * d6ji( jj,J4,jx,rank,jd,J5) * &
                         d6ji(J1,jx,jc,jd,J2,rank) * v_elem(a,j,c,i,J3,AA,jbas) *&
                         pphh_tensor_elem(i,b,j,d,J4,J5,BB,jbas) * (-1)**((jc+jd-J2)/2) &
                   
                         + (-1)** ((J2+J3 + jd - ji )/2) * sqrt( (J1+1.d0) * (J2+1.d0)  &
                         * (J4+1.d0) * (J5+1.d0) ) * (jx+1.d0) * (J3+1.d0)  &
                    * coef9(jj,J3,jb,J4,ji,ja,jx,jc,J1) * d6ji( jj,J4,jx,rank,jd,J5) * &
                         d6ji(J1,jx,jc,jd,J2,rank) * v_elem(b,j,c,i,J3,AA,jbas) *&
                         pphh_tensor_elem(i,a,j,d,J4,J5,BB,jbas)  *(-1)**((ja+jb+jc+jd+J1+J2)/2)  &
                         )
                 end do
              end do
           end do
        end do
     end do
  end do
!!$OMP END PARALLEL DO

  EOM_scalar_tensor_2body_comm = sm 
  
end function 

end module 
