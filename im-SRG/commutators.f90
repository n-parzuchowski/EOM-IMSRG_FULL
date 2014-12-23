module commutators
  use basic_IMSRG
  ! basic commutator functions 
  
  
  
contains
!=========================================================
!=========================================================
real(8) function commutator_110(L,R,jbas) 
  ! zero body part of [L1,R1] 
  implicit none 
  
  type(sq_op) :: L,R 
  type(spd) :: jbas
  integer :: a,i,ji
  real(8) :: sm
  
  sm = 0.d0 
  do i = 1,L%belowEf
     ji = jbas%jj(jbas%holes(i)) 
     
     do a = 1,L%Nsp-L%belowEf 
            
        sm = sm +  ( L%fph(a,i) * R%fph(a,i) * L%herm - &
             R%fph(a,i) * L%fph(a,i) * R%herm ) * (ji + 1) 
    
     end do
  end do 
        
  commutator_110 = sm

end function 
!=========================================================
!=========================================================
subroutine commutator_111(L,R,RES,jbas) 
  ! one body part of [L1,R1] 
  ! VERIFIED CORRECT. 
  implicit none 
  
  type(sq_op) :: L,R,RES
  type(spd) :: jbas
  integer :: p,a,i,ji,hol,par
  real(8) :: sm 
  real(8),dimension(L%belowEF,L%belowEF) :: th1,th2
  real(8),dimension(L%Nsp-L%belowEF,L%belowEF) :: tb1,tb2
  real(8),dimension(L%Nsp-L%belowEF,L%Nsp-L%belowEF) :: tp1,tp2 
 
  hol = L%belowEf
  par = L%Nsp - hol

! dfhh ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ! two matrix mulitplies, one for the sum over holes, 
  ! and one for the sum over particles. Add them together
  call dgemm('N','N',hol,hol,hol,al,L%fhh,hol,R%fhh,hol,bet,th1,hol) 
  call dgemm('T','N',hol,hol,par,al,L%fph,par,R%fph,par,bet,th2,hol)
  
  RES%fhh = th1 + th2*L%herm  - Transpose(th1*L%herm + th2) *R%herm
   
!dfph~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  call dgemm('N','N',par,hol,hol,al,L%fph,par,R%fhh,hol,bet,tb1,par) 
  call dgemm('N','N',par,hol,par,al,L%fpp,par,R%fph,par,bet,tb2,par) 
  
  RES%fph = tb1 + tb2 
  
  call dgemm('N','N',par,hol,hol,al,R%fph,par,L%fhh,hol,bet,tb1,par) 
  call dgemm('N','N',par,hol,par,al,R%fpp,par,L%fph,par,bet,tb2,par) 
  
  RES%fph = RES%fph - tb1 - tb2 
  
!dfpp~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  
  call dgemm('N','T',par,par,hol,al,L%fph,par,R%fph,par,bet,tp1,par) 
  call dgemm('N','N',par,par,par,al,L%fpp,par,R%fpp,par,bet,tp2,par) 
  
  RES%fpp = tp1 * R%herm + tp2 - Transpose(tp1 + tp2 * R%herm) * L%herm
    
end subroutine
!=========================================================
!=========================================================
subroutine commutator_121(L,R,RES,jbas) 
  ! onebody part of [L1,R2] - [R1,L2] 
  ! this one is brute force. 
  ! not sure of a faster way to do this
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: L,R,RES
  integer :: JT,PAR,TZ,ji,ja,jp,jq,a,i,p,q
  integer :: ti,ta,tp,tq,li,la,lp,lq,ak,ik,pk,qk
  real(8) :: sm,smx,smy,smx2,smy2
  
! dfhh ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  do p = 1,L%belowEF
     
     pk = jbas%holes(p)
     jp = jbas%jj(pk) 
     tp = jbas%itzp(pk)
     lp = jbas%ll(pk)
     
     do q = p,L%belowEF
      
        qk = jbas%holes(q) 
        jq = jbas%jj(qk) 
        tq = jbas%itzp(qk)
        lq = jbas%ll(qk)
      
        ! check if this state is allowed
        if (jq .ne. jp) cycle
        if (lq .ne. lp) cycle
        if (tq .ne. tp) cycle 
        
        ! internal sum
        sm = 0.d0
        do a = 1,L%Nsp - L%belowEF
           
           ak = jbas%parts(a)
           ja = jbas%jj(ak) 
           ta = jbas%itzp(ak)
           la = jbas%ll(ak)
           
           do i = 1,L%belowEF
              
              ik = jbas%holes(i)
              ji = jbas%jj(ik) 
              ti = jbas%itzp(ik)
              li = jbas%ll(ik)
              
              ! check if this intermediate exists
              if (ji .ne. ja) cycle
              if (li .ne. la) cycle
              if (ti .ne. ta) cycle 
                 
              smx = 0.d0 
              smy = 0.d0
              smx2 = 0.d0 
              smy2 = 0.d0 
              ! sum over J_total
              do JT = abs(ji - jq),ji+jq,2
                 smx = smx + v_elem(ak,pk,ik,qk,JT,R,jbas)*(JT + 1)
                 smx2 = smx2 + v_same(L)*(JT + 1)
                 smy = smy + v_elem(ik,pk,ak,qk,JT,R,jbas)*(JT + 1) 
                 smy2 = smy2 + v_same(L)*(JT + 1) 
              end do 
              
              sm = sm + L%fph(a,i) * (L%herm*smx - smy) - &
                   R%fph(a,i) * (R%herm*smx2 - smy2)
           end do 
        end do 
        
        RES%fhh(p,q) = RES%fhh(p,q) + sm /(jp + 1.d0) 
        RES%fhh(q,p) = RES%fhh(p,q) * RES%herm
        
     end do 
  end do 
             
!dfpp~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
  do p = 1, L%Nsp - L%belowEF
     pk = jbas%parts(p)
     jp = jbas%jj(pk) 
     tp = jbas%itzp(pk)
     lp = jbas%ll(pk)
     
     do q = p,L%Nsp - L%belowEF
      
        qk = jbas%parts(q) 
        jq = jbas%jj(qk) 
        tq = jbas%itzp(qk)
        lq = jbas%ll(qk)
      
        if (jq .ne. jp) cycle
        if (lq .ne. lp) cycle
        if (tq .ne. tp) cycle 
        
        sm = 0.d0
        do a = 1,L%Nsp - L%belowEF
           ak = jbas%parts(a)
           ja = jbas%jj(ak) 
           ta = jbas%itzp(ak)
           la = jbas%ll(ak)
           
           do i = 1,L%belowEF
              ik = jbas%holes(i)
              ji = jbas%jj(ik) 
              ti = jbas%itzp(ik)
              li = jbas%ll(ik)
              
              if (ji .ne. ja) cycle
              if (li .ne. la) cycle
              if (ti .ne. ta) cycle 
              
              PAR = mod(li + lq,2) 
              TZ = (ti + tq)/2 
              
              smx = 0.d0 
              smy = 0.d0
              smx2 = 0.d0 
              smy2 = 0.d0 
              ! sum over J_total
              do JT = abs(ji - jq),ji+jq,2
                 smx = smx + v_elem(ak,pk,ik,qk,JT,R,jbas)*(JT + 1)
                 smx2 = smx2 + v_same(L)*(JT + 1)
                 smy = smy + v_elem(ik,pk,ak,qk,JT,R,jbas)*(JT + 1) 
                 smy2 = smy2 + v_same(L)*(JT + 1) 
              end do 
              
              sm = sm + L%fph(a,i) * (L%herm*smx - smy)  - &
                   R%fph(a,i) * (R%herm*smx2 - smy2) 
           end do 
        end do 
        
        RES%fpp(p,q) = RES%fpp(p,q) + sm/(jp + 1.d0) 
        RES%fpp(q,p) = RES%fpp(p,q) * RES%herm
        
     end do 
  end do 
         
 !dfph~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  do p = 1, L%Nsp - L%belowEF
     pk = jbas%parts(p)
     jp = jbas%jj(pk) 
     tp = jbas%itzp(pk)
     lp = jbas%ll(pk)
     
     do q = 1, L%belowEF
      
        qk = jbas%holes(q) 
        jq = jbas%jj(qk) 
        tq = jbas%itzp(qk)
        lq = jbas%ll(qk)
      
        if (jq .ne. jp) cycle
        if (lq .ne. lp) cycle
        if (tq .ne. tp) cycle 
        
        sm = 0.d0
        do a = 1,L%Nsp - L%belowEF
           ak = jbas%parts(a)
           ja = jbas%jj(ak) 
           ta = jbas%itzp(ak)
           la = jbas%ll(ak)
           
           do i = 1,L%belowEF
              ik = jbas%holes(i)
              ji = jbas%jj(ik) 
              ti = jbas%itzp(ik)
              li = jbas%ll(ik)
              
              if (ji .ne. ja) cycle
              if (li .ne. la) cycle
              if (ti .ne. ta) cycle 
              
              PAR = mod(li + lq,2) 
              TZ = (ti + tq)/2 
              
              smx = 0.d0 
              smy = 0.d0
              smx2 = 0.d0 
              smy2 = 0.d0 
              ! sum over J_total
              do JT = abs(ji - jq),ji+jq,2
                 smx = smx + v_elem(ak,pk,ik,qk,JT,R,jbas)*(JT + 1)
                 smx2 = smx2 + v_same(L)*(JT + 1)
                 smy = smy + v_elem(ik,pk,ak,qk,JT,R,jbas)*(JT + 1) 
                 smy2 = smy2 + v_same(L)*(JT + 1) 
              end do 
              
              sm = sm + L%fph(a,i) * (L%herm*smx - smy) - &
                   R%fph(a,i) * (R%herm*smx2 - smy2)
           end do 
        end do 
        
        RES%fph(p,q) = RES%fph(p,q) + sm /(jp + 1.d0)
        
     end do 
  end do             
 
end subroutine             
!====================================================================
!====================================================================
subroutine commutator_122(L,R,RES,jbas) 
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: L,R,RES
  integer :: q,IX,JX,nh,np,nb,i,JT
  integer :: a,b,c,d,ja,jb,jc,jd,ji,g_ix,q_sp,i_sp
  integer :: ta,tb,tc,td,ti,la,lb,lc,ld,li,spec
  integer :: jxstart,jxend,ixend,c1,c2,n1,n2
  logical :: square
  real(8) ::  sm
  
  do q = 1, L%nblocks
     
     JT = L%mat(q)%lam(1)
     
     nh = L%mat(q)%nhh
     np = L%mat(q)%npp
     nb = L%mat(q)%nph
  
     do g_ix = 1,6 
   
        ! figure out how big the array is
        n1 = size(L%mat(q)%gam(g_ix)%X(:,1))
        n2 = size(L%mat(q)%gam(g_ix)%X(1,:))
        if ((n1*n2) == 0) cycle 
        
         ! read in information about which 
        ! array we are using from public arrays
        c1 = sea1(g_ix) 
        c2 = sea2(g_ix) 
        square = sqs(g_ix) 
        jxstart = jst(g_ix) 
        
                
     ! main calculation
   
     do IX = 1,n1
        a = L%mat(q)%qn(c1)%Y(IX,1)
        ja = jbas%jj(a)
        la = jbas%ll(a)
        ta = jbas%itzp(a) 
             
        b = L%mat(q)%qn(c1)%Y(IX,2)
        jb = jbas%jj(b)
        lb = jbas%ll(b)
        tb = jbas%itzp(b)
 
        do JX = min(jxstart,IX),n2
           
           c = L%mat(q)%qn(c2)%Y(JX,1)
           jc = jbas%jj(c)
           lc = jbas%ll(c)
           tc = jbas%itzp(c)

           d = L%mat(q)%qn(c2)%Y(JX,2)
           jd = jbas%jj(d)
           ld = jbas%ll(d)
           td = jbas%itzp(d)
                   
           sm = 0.d0 

            ! a is replaced
            q_sp = sp_block_index(ja,la,ta,jbas) 
            do i = 1,size(jbas%states(q_sp)%Z)   
              
               i_sp = jbas%states(q_sp)%Z(i) 
               
               sm = sm + f_elem(a,i_sp,L,jbas)*v_elem(i_sp,b,c,d,JT,R,jbas)&
                    - f_elem(a,i_sp,R,jbas)*v_same(L) 
            end do 
              
            ! b is replaced
            q_sp = sp_block_index(jb,lb,tb,jbas) 
            do i = 1,size(jbas%states(q_sp)%Z)   
              
               i_sp = jbas%states(q_sp)%Z(i) 
               
               sm = sm + f_elem(b,i_sp,L,jbas)*v_elem(a,i_sp,c,d,JT,R,jbas)&
                    - f_elem(b,i_sp,R,jbas)*v_same(L) 
            end do 
            
            ! c is replaced
            q_sp = sp_block_index(jc,lc,tc,jbas) 
            do i = 1,size(jbas%states(q_sp)%Z)   
              
               i_sp = jbas%states(q_sp)%Z(i) 
               
               sm = sm - f_elem(i_sp,c,L,jbas)*v_elem(a,b,i_sp,d,JT,R,jbas)&
                    + f_elem(i_sp,c,R,jbas)*v_same(L) 
            end do 
            
            ! d is replaced
            q_sp = sp_block_index(jd,ld,td,jbas) 
            do i = 1,size(jbas%states(q_sp)%Z)   
              
               i_sp = jbas%states(q_sp)%Z(i) 
               
               sm = sm - f_elem(i_sp,d,L,jbas)*v_elem(a,b,c,i_sp,JT,R,jbas)&
                    + f_elem(i_sp,d,R,jbas)*v_same(L) 
            end do 
          
              sm = sm / sqrt(1.d0 + kron_del(a,b)) /sqrt(1.d0 + kron_del(c,d)) 

           RES%mat(q)%gam(g_ix)%X(IX,JX) = sm 
           if (square) RES%mat(q)%gam(g_ix)%X(JX,IX) = sm * RES%herm

        end do
     end do
   
     end do 
  end do 

end subroutine           
!=============================================
!=============================================
real(8) function commutator_220(L,R,jbas) 
  ! zero body part of [L2,R2] 
  !VERIFIED
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: L,R
  integer :: IX,JX,q,np,nh
  real(8) :: sm,smx
  
  if (R%herm * L%herm == -1) then 
  
  sm = 0.d0 
  do q = 1, L%nblocks
     nh = L%mat(q)%nhh
     np = L%mat(q)%npp
     
     smx = 0 
     do IX = 1, np 
        do JX = 1,nh
           
           smx = smx - L%mat(q)%gam(3)%X(IX,JX) * R%mat(q)%gam(3)%X(IX,JX) 
        end do
     end do 
     
     sm = sm + smx * (L%mat(q)%lam(1) + 1) 
  end do 
  
  commutator_220 = sm * 2.d0 
  else 
  commutator_220 = 0.d0
  end if 

end function              
!===================================================================
!===================================================================
subroutine commutator_221(L,R,RES,w1,w2,jbas) 
  ! verified
  ! THIS NEEDS TO BE RUN AFTER 222_pp_hh 
  ! 222_pp_hh sets up the intermediary matrices (w1,w2) 
  implicit none
  
  type(spd) :: jbas
  type(sq_op) :: L,R,RES,w1,w2
  integer :: i,j,q,Abody,Ntot,nh,np,nb,a,c
  integer :: ik,jk,ck,ji,jj,ti,tj,li,lj,jc,JT,pm
  real(8) :: sm
  
  Abody = L%belowEF
  Ntot = L%Nsp
   
  pm = R%herm*L%herm
  
! fhh
  do i = 1 , Abody
     ik = jbas%holes(i) 
     ji = jbas%jj(ik) 
     li = jbas%ll(ik) 
     ti = jbas%itzp(ik) 
     
     do j = i , Abody
        
        jk = jbas%holes(j) 
        jj = jbas%jj(jk) 
        if (ji .ne. jj)  cycle
        lj = jbas%ll(jk) 
        if (li .ne. lj)  cycle
        tj = jbas%itzp(jk)
        if (tj .ne. ti) cycle 
                
        sm = 0.d0 
        do c = 1, Abody
           ck = jbas%holes(c) 
           jc = jbas%jj(ck)
           do JT = abs(jc - ji),jc+ji,2
              sm = sm + (v_elem(ck,ik,ck,jk,JT,w1,jbas) - pm* &
                   v_elem(ck,jk,ck,ik,JT,w1,jbas))*(JT + 1) 
           end do 
        end do 
        
        do c = 1, Ntot - Abody
           ck = jbas%parts(c) 
           jc = jbas%jj(ck)
           do JT = abs(jc - ji),jc+ji,2
              sm = sm + (v_elem(ck,ik,ck,jk,JT,w2,jbas) - pm* &
                   v_elem(ck,jk,ck,ik,JT,w2,jbas) )*(JT + 1)
           end do 
        end do 
     
        RES%fhh(i,j) = RES%fhh(i,j) + sm / (ji + 1.d0 )
        RES%fhh(j,i) = RES%fhh(i,j) * RES%herm
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
        do c = 1, Abody
           ck = jbas%holes(c) 
           jc = jbas%jj(ck)
           ! w1 matrix results from multiplying the pp channel
           do JT = abs(jc - ji),jc+ji,2
              sm = sm + (v_elem(ck,ik,ck,jk,JT,w1,jbas) - pm* &
                   v_elem(ck,jk,ck,ik,JT,w1,jbas))*(JT + 1)        
           end do 
        end do 
        
        do c = 1, Ntot - Abody
           ck = jbas%parts(c) 
           jc = jbas%jj(ck)
           ! w2 matrix results from multiplying the hh channel
           do JT = abs(jc - ji),jc+ji,2
              sm = sm + (v_elem(ck,ik,ck,jk,JT,w2,jbas) - pm* &
                   v_elem(ck,jk,ck,ik,JT,w2,jbas))*(JT + 1)
           end do 
        end do 
     
        RES%fpp(i,j) = RES%fpp(i,j) + sm / (ji + 1.d0) 
        RES%fpp(j,i) = RES%fpp(i,j) * RES%herm
     end do 
  end do       
  

  ! fph
  do i = 1 , Ntot - Abody
     ik = jbas%parts(i) 
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
        do c = 1, Abody
           ck = jbas%holes(c) 
           jc = jbas%jj(ck)
           do JT = abs(jc - ji),jc+ji,2
              sm = sm + v_elem(ck,ik,ck,jk,JT,w1,jbas)* ( JT + 1) 
           end do 
        end do 
        
        do c = 1, Ntot - Abody
           ck = jbas%parts(c) 
           jc = jbas%jj(ck)
           do JT = abs(jc - ji),jc+ji,2
              sm = sm + v_elem(ck,ik,ck,jk,JT,w2,jbas) * (JT + 1)
           end do 
        end do 
     
        RES%fph(i,j) = RES%fph(i,j) + sm / (ji + 1.d0 )

     end do 
  end do       

end subroutine
!===================================================================
!===================================================================
subroutine commutator_222_pp_hh(L,R,RES,w1,w2,jbas) 
  !VERIFIED
  !NEEDS TO BE RUN BEFORE 221, because it sets up the 
  !intermediary matrices
  implicit none
  
  type(spd) :: jbas
  type(sq_op) ::  L,R,RES,w1,w2
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
     
  
     if (np .ne. 0)  then 
       
     !L_pppp . R_pppp = W1_pppp   
     call dgemm('N','N',np,np,np,al,L%mat(q)%gam(1)%X,np,&
          R%mat(q)%gam(1)%X,np,bet,w1%mat(q)%gam(1)%X,np)
     end if 
     
     if (nh .ne. 0) then 
     !L_hhhh . R_hhhh = W2_hhhh
     call dgemm('N','N',nh,nh,nh,al,L%mat(q)%gam(5)%X,nh,&
          R%mat(q)%gam(5)%X,nh,bet,w2%mat(q)%gam(5)%X,nh)
     end if 
          
     if (np*nh .ne. 0) then 
     !L_hhpp . R_pphh = W1_hhhh
     call dgemm('T','N',nh,nh,np,al,L%mat(q)%gam(3)%X,np,&
          R%mat(q)%gam(3)%X,np,bet,w1%mat(q)%gam(5)%X,nh)
    
     w1%mat(q)%gam(5)%X = w1%mat(q)%gam(5)%X * L%herm 
     
     !L_pphh . R_hhpp = W2_pppp
     call dgemm('N','T',np,np,nh,al,L%mat(q)%gam(3)%X,np,&
          R%mat(q)%gam(3)%X,np,bet,w2%mat(q)%gam(1)%X,np)
     
     w2%mat(q)%gam(1)%X = w2%mat(q)%gam(1)%X * R%herm
  
     bet_off = -1
     !R_pppp . L_pphh = W1_pphh (Transposed)
     call dgemm('N','N',np,nh,np,al,R%mat(q)%gam(1)%X,np,&
          L%mat(q)%gam(3)%X,np,bet,w1%mat(q)%gam(3)%X,np) 
    
     !L_pppp . R_pphh - R_pppp . L_pphh = W1_pphh
     call dgemm('N','N',np,nh,np,al,L%mat(q)%gam(1)%X,np,&
          R%mat(q)%gam(3)%X,np,bet_off,w1%mat(q)%gam(3)%X,np) 
     
     !R_pphh . L_hhhh = W2_pphh (Transposed)
     call dgemm('N','N',np,nh,nh,al,R%mat(q)%gam(3)%X,np,&
          L%mat(q)%gam(5)%X,nh,bet,w2%mat(q)%gam(3)%X,np) 
      
     !L_pphh . R_hhhh - R_pphh . L_hhhh = W2_pphh 
     call dgemm('N','N',np,nh,nh,al,L%mat(q)%gam(3)%X,np,&
          R%mat(q)%gam(5)%X,nh,bet_off,w2%mat(q)%gam(3)%X,np)           
     end if
     
     ! Vpppp
     RES%mat(q)%gam(1)%X = RES%mat(q)%gam(1)%X + &
          w1%mat(q)%gam(1)%X  - w2%mat(q)%gam(1)%X - &
           ( Transpose(w1%mat(q)%gam(1)%X - w2%mat(q)%gam(1)%X ) ) * &
           pm
     
     ! Vhhhh
     RES%mat(q)%gam(5)%X = RES%mat(q)%gam(5)%X + &
          w1%mat(q)%gam(5)%X  - w2%mat(q)%gam(5)%X - &
           ( Transpose(w1%mat(q)%gam(5)%X - w2%mat(q)%gam(5)%X ) ) * &
           pm
     
     ! Vpphh
     RES%mat(q)%gam(3)%X = RES%mat(q)%gam(3)%X + &
          w1%mat(q)%gam(3)%X  - w2%mat(q)%gam(3)%X 
     
     
     
     
     if (nb*nh .ne. 0) then 
     
     !L_phhh . R_hhph = W2_phph
     call dgemm('N','T',nb,nb,nh,al,L%mat(q)%gam(6)%X,nb,&
          R%mat(q)%gam(6)%X,nb,bet,w2%mat(q)%gam(4)%X,nb) 
     
     w2%mat(q)%gam(4)%X = w2%mat(q)%gam(4)%X * R%herm
     
     end if

     if (nb*np .ne. 0) then 
     !L_phpp . R_ppph = W1_phph
     call dgemm('T','N',nb,nb,np,al,L%mat(q)%gam(2)%X,np,&
          R%mat(q)%gam(2)%X,np,bet,w1%mat(q)%gam(4)%X,nb) 
     
     w1%mat(q)%gam(4)%X = w1%mat(q)%gam(4)%X * L%herm 
     end if 
     
     ! Vphph
     RES%mat(q)%gam(4)%X = RES%mat(q)%gam(4)%X + &
          w1%mat(q)%gam(4)%X  - w2%mat(q)%gam(4)%X - &
           ( Transpose(w1%mat(q)%gam(4)%X - w2%mat(q)%gam(4)%X ) ) * &
           pm
     
 
     if (nb*nh .ne. 0)  then 
     !R_phhh . L_hhhh = W2_phhh 
     call dgemm('N','N',nb,nh,nh,al,R%mat(q)%gam(6)%X,nb,&
          L%mat(q)%gam(5)%X,nh,bet,w2%mat(q)%gam(6)%X,nb) 
     
     bet_off = -1
     !L_phhh . R_hhhh - R_phhh . L_hhhh = W2_phhh
     call dgemm('N','N',nb,nh,nh,al,L%mat(q)%gam(6)%X,nb,&
          R%mat(q)%gam(5)%X,nh,bet_off,w2%mat(q)%gam(6)%X,nb) 
     end if 
     
     if (nb*np .ne. 0) then 
     bet_off = -1 
     !R_pppp . L_ppph = W1_ppph
     call dgemm('N','N',np,nb,np,al,L%mat(q)%gam(1)%X,np,&
          R%mat(q)%gam(2)%X,np,bet,w1%mat(q)%gam(2)%X,np) 
          
     !L_pppp . R_ppph + R_pppp . L_ppph = W1_ppph 
     call dgemm('N','N',np,nb,np,al,L%mat(q)%gam(1)%X,np,&
          R%mat(q)%gam(2)%X,np,bet_off,w1%mat(q)%gam(2)%X,np)
     end if 

     !  the following are
     !  slightly messier because they need to be tranposed 
     !  to matrices which i don't have stored
     
     !R_phpp . L_pphh = W1_phhh
     if (nb*np*nh .ne. 0) then 
     call dgemm('T','N',nb,nh,np,al,R%mat(q)%gam(2)%X,np,&
          L%mat(q)%gam(3)%X,np,bet,w1%mat(q)%gam(6)%X,nb)
         
     al_off = L%herm 
     bet_off = -1*R%herm

     !L_phpp . R_pphh - R_phpp . L_pphh = W1_phhh 
     call dgemm('T','N',nb,nh,np,al_off,L%mat(q)%gam(2)%X,np,&
          R%mat(q)%gam(3)%X,np,bet_off,w1%mat(q)%gam(6)%X,nb)

     
     !R_pphh . L_hhph = W2_ppph 
     call dgemm('N','T',np,nb,nh,al,R%mat(q)%gam(3)%X,np,&
          L%mat(q)%gam(6)%X,nb,bet,w2%mat(q)%gam(2)%X,np)
     
     al_off = R%herm 
     bet_off = -1*L%herm 
     !L_pphh . R_hhph + R_pphh . L_hhph  = W2_ppph
     call dgemm('N','T',np,nb,nh,al_off,L%mat(q)%gam(3)%X,np,&
          R%mat(q)%gam(6)%X,nb,bet_off,w2%mat(q)%gam(2)%X,np)
     end if

     ! Vphhh
     RES%mat(q)%gam(6)%X = RES%mat(q)%gam(6)%X + &
          w1%mat(q)%gam(6)%X  - w2%mat(q)%gam(6)%X
     
     ! Vppph
     RES%mat(q)%gam(2)%X = RES%mat(q)%gam(2)%X + &
          w1%mat(q)%gam(2)%X  - w2%mat(q)%gam(2)%X

  end do

end subroutine 
!=================================================================
!=================================================================
 subroutine commutator_222_ph(LCC,RCC,RES,WCC,jbas) 
   ! VERIFIED ph channel 2body commutator. DFWT! 
   implicit none 
  
   type(spd) :: jbas
   type(sq_op) :: RES
   type(cross_coupled_31_mat) :: LCC,RCC,WCC
   integer :: nh,np,nb,q,IX,JX,i,j,k,l,rinx,Tz,PAR,JTM
   integer :: ji,jj,jk,jl,ti,tj,tk,tl,li,lj,lk,ll,n1,n2,c1,c2,jxstart
   integer :: JP, Jtot,Ntot,qx,jmin,jmax,rik,rjl,ril,rjk,g_ix,thread,total_threads
   real(8) :: sm ,pre,pre2,omp_get_wtime ,t1,t2
   logical :: square
   

  Ntot = RES%Nsp
  JTM = jbas%Jtotal_max
  total_threads = size(RES%direct_omp) - 1
   ! construct intermediate matrices
 
   do q = 1,LCC%nblocks
      
      nb = LCC%nph(q)
      
      rinx = LCC%rlen(q)  
      
      if (nb * rinx == 0) cycle
      
      call dgemm('N','N',rinx,rinx,nb,al,LCC%CCX(q)%X,rinx,&
           RCC%CCR(q)%X,nb,bet,WCC%CCX(q)%X,rinx) 
     
      call dgemm('N','N',rinx,rinx,nb,al,RCC%CCX(q)%X,rinx,&
           LCC%CCR(q)%X,nb,bet,WCC%CCR(q)%X,rinx) 
   
   end do 

!$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE), SHARED(RES,WCC)  
   do thread = 1, total_threads
   do q = 1+RES%direct_omp(thread),RES%direct_omp(thread+1) 
     
     Jtot = RES%mat(q)%lam(1)
     
     nh = RES%mat(q)%nhh
     np = RES%mat(q)%npp
     nb = RES%mat(q)%nph
          
     do g_ix = 1,6 
   
        ! figure out how big the array is
        n1 = size(RES%mat(q)%gam(g_ix)%X(:,1))
        n2 = size(RES%mat(q)%gam(g_ix)%X(1,:))
        if ((n1*n2) == 0) cycle 
        
        ! read in information about which 
        ! array we are using from public arrays
        c1 = sea1(g_ix) 
        c2 = sea2(g_ix) 
        square = sqs(g_ix) 
        jxstart = jst(g_ix) 
        
      do  IX =  1, n1 
         pre = 1.d0 

         i = RES%mat(q)%qn(c1)%Y(IX,1)
         j = RES%mat(q)%qn(c1)%Y(IX,2)
 
         if (i == j )  pre  = .70710678118d0
         ji = jbas%jj(i) 
         jj = jbas%jj(j) 
         li = jbas%ll(i) 
         lj = jbas%ll(j)
         ti = jbas%itzp(i) 
         tj = jbas%itzp(j)
         
         do JX = min(jxstart,IX),n2
            pre2 = 1.d0 
            k = RES%mat(q)%qn(c2)%Y(JX,1)
            l = RES%mat(q)%qn(c2)%Y(JX,2)
            
            if (k == l )  pre2 = .70710678118d0
            jk = jbas%jj(k) 
            jl = jbas%jj(l) 
            lk = jbas%ll(k) 
            ll = jbas%ll(l)
            tk = jbas%itzp(k) 
            tl = jbas%itzp(l)
            
            sm = 0.d0 
                       
            jmin = max( abs(jj - jl) , abs(ji - jk )) 
            jmax = min( jj + jl , ji + jk ) 
            
            
            Tz = abs(ti -tk)/2 
            if (abs(tl - tj) .ne. Tz*2)  cycle 
            PAR = mod(li+lk,2) 
            if (mod(ll+lj,2) .ne. PAR) cycle 
            
            
            do JP = jmin,jmax,2
                 
                  qx = JP/2+1 + Tz*(JTM+1) + 2*PAR*(JTM+1)
                  rjl = specific_rval(j,l,Ntot,qx,LCC)
                  rik = specific_rval(i,k,Ntot,qx,LCC)
                  
                  sm = sm - (RES%herm*WCC%CCX(qx)%X(rjl,rik) - &
                       WCC%CCR(qx)%X(rjl,rik) - &
                       RES%herm*WCC%CCR(qx)%X(rik,rjl) + &
                       WCC%CCX(qx)%X(rik,rjl) ) * &
                       sixj(jk,jl,Jtot,jj,ji,JP) * &
                       (-1)**((ji + jl + Jtot)/2) 
            
            end do 

            Tz = abs(ti -tl)/2 
            if (abs(tk - tj) .ne. Tz*2) cycle 
            PAR = mod(li+ll,2) 
            if (mod(lk+lj,2) .ne. PAR) cycle 
            
               jmin = max( abs(ji - jl) , abs(jj - jk )) 
               jmax = min( ji + jl , jj + jk ) 
               
               do JP = jmin,jmax,2
                  
                  !qx = JP/2 + 1
                  qx = JP/2+1 + Tz*(JTM+1) + 2*PAR*(JTM+1)
                
                  ril = specific_rval(i,l,Ntot,qx,LCC)
                  rjk = specific_rval(j,k,Ntot,qx,LCC)
                  
                  sm = sm + ( WCC%CCR(qx)%X(ril,rjk) - &
                       RES%herm*WCC%CCX(qx)%X(ril,rjk) - &
                       WCC%CCX(qx)%X(rjk,ril) + &
                       RES%herm*WCC%CCR(qx)%X(rjk,ril) ) * &
                       sixj(jk,jl,Jtot,ji,jj,JP) * &
                       (-1)**((ji + jl)/2)
            
               end do 

           RES%mat(q)%gam(g_ix)%X(IX,JX) = &
                RES%mat(q)%gam(g_ix)%X(IX,JX) + sm * pre * pre2 
           if (square) RES%mat(q)%gam(g_ix)%X(JX,IX) =  &
                RES%mat(q)%gam(g_ix)%X(IX,JX) * RES%herm
           
         end do 
      end do
      end do 
   end do
   end do 
 
!$OMP END PARALLEL DO 
   
end subroutine 
!=====================================================
!=====================================================      
integer function specific_rval(i,l,Ntot,q,LCC) 
  implicit none 
  
  type(cross_coupled_31_mat) :: LCC
  integer :: i,l,Ntot,x,g,q
  
  x = CCindex(i,l,Ntot)
  g = 1
  do while (LCC%qmap(x)%Z(g) .ne. q )
  
     g = g + 1
  end do
  
  specific_rval = LCC%rmap(x)%Z(g)
end function 
!============================================
!============================================
end module 
  
  
  
