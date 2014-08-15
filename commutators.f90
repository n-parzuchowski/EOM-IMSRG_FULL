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
  
  RES%fpp = tp1 * R%herm + tp2 - Transpose(tp1 - tp2 * R%herm) * L%herm
    
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
                 smy = smy + v_elem(ik,pk,ak,qk,JT,R,jbas)*(JT + 1) 
                 smx2 = smx2 + v_elem(ak,pk,ik,qk,JT,L,jbas)*(JT + 1)
                 smy2 = smy2 + v_elem(ik,pk,ak,qk,JT,L,jbas)*(JT + 1) 
              end do 
              
              sm = sm + L%fph(a,i) * (L%herm*smx - smy) /(ji + 1.d0) - &
                   (R%fph(a,i) * (R%herm*smx2 - smy2) /(ji + 1.d0)) 
           end do 
        end do 
        
        RES%fhh(p,q) = RES%fhh(p,q) + sm 
        RES%fhh(q,p) = sm 
        
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
                 smy = smy + v_elem(ik,pk,ak,qk,JT,R,jbas)*(JT + 1) 
                 smx2 = smx2 + v_elem(ak,pk,ik,qk,JT,L,jbas)*(JT + 1)
                 smy2 = smy2 + v_elem(ik,pk,ak,qk,JT,L,jbas)*(JT + 1) 
              end do 
              
              sm = sm + L%fph(a,i) * (L%herm*smx - smy) /(ji + 1.d0) - &
                   (R%fph(a,i) * (R%herm*smx2 - smy2) /(ji + 1.d0)) 
           end do 
        end do 
        
        RES%fpp(p,q) = RES%fpp(p,q) + sm 
        RES%fpp(q,p) = sm 
        
     end do 
  end do 
         
 !dfph~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  do p = 1, L%Nsp - L%belowEF
     pk = jbas%parts(p)
     jp = jbas%jj(pk) 
     tp = jbas%itzp(pk)
     lp = jbas%ll(pk)
     
     do q = 1, L%belowEF
      
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
                 smy = smy + v_elem(ik,pk,ak,qk,JT,R,jbas)*(JT + 1) 
                 smx2 = smx2 + v_elem(ak,pk,ik,qk,JT,L,jbas)*(JT + 1)
                 smy2 = smy2 + v_elem(ik,pk,ak,qk,JT,L,jbas)*(JT + 1) 
              end do 
              
              sm = sm + L%fph(a,i) * (L%herm*smx - smy) /(ji + 1.d0) - &
                   (R%fph(a,i) * (R%herm*smx2 - smy2) /(ji + 1.d0)) 
           end do 
        end do 
        
        RES%fph(p,q) = RES%fph(p,q) + sm 
        
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
     if (nh*np*nb == 0) cycle
     
     do g_ix = 1,6 
   
        ! figure out how big the array is
        n1 = size(L%mat(q)%gam(g_ix)%X(:,1))
        n2 = size(L%mat(q)%gam(g_ix)%X(1,:))
        if ((n1*n2) == 0) cycle 
        
        ! figure out which type n1 is
        if (n1 == nh) then 
           c1 = 3 
        else if (n1 == np) then 
           c1 = 1
        else 
           c1 = 2
        end if 
        
        
        ! decide if it's a square array or rectangle
        if (n1 == n2) then
           jxstart = 10000
           square = .true.
           c2 = c1 
        else 
           jxstart = 1
           square = .false.
           if (n2 == nh) then 
              c2 = 3 
           else if (n1 == np) then 
              c2 = 1
           else 
              c2 = 2
           end if
        end if 
                
     ! main calculation
     !   call print_matrix( L%mat(q)%gam(1)%X )
   
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
                    - f_elem(a,i_sp,R,jbas)*v_elem(i_sp,b,c,d,JT,L,jbas) 
            end do 
              
            ! b is replaced
            q_sp = sp_block_index(jb,lb,tb,jbas) 
            do i = 1,size(jbas%states(q_sp)%Z)   
              
               i_sp = jbas%states(q_sp)%Z(i) 
               
               sm = sm + f_elem(b,i_sp,L,jbas)*v_elem(a,i_sp,c,d,JT,R,jbas)&
                    - f_elem(b,i_sp,R,jbas)*v_elem(a,i_sp,c,d,JT,L,jbas) 
            end do 
            
            ! c is replaced
            q_sp = sp_block_index(jc,lc,tc,jbas) 
            do i = 1,size(jbas%states(q_sp)%Z)   
              
               i_sp = jbas%states(q_sp)%Z(i) 
               
               sm = sm - f_elem(i_sp,c,L,jbas)*v_elem(a,b,i_sp,d,JT,R,jbas)&
                    + f_elem(i_sp,c,R,jbas)*v_elem(a,b,i_sp,d,JT,L,jbas) 
            end do 
            
            ! d is replaced
            q_sp = sp_block_index(jd,ld,td,jbas) 
            do i = 1,size(jbas%states(q_sp)%Z)   
              
               i_sp = jbas%states(q_sp)%Z(i) 
               
               sm = sm - f_elem(i_sp,d,L,jbas)*v_elem(a,b,c,i_sp,JT,R,jbas)&
                    + f_elem(i_sp,d,R,jbas)*v_elem(a,b,c,i_sp,JT,L,jbas) 
            end do 
          
              
           RES%mat(q)%gam(g_ix)%X(IX,JX) = sm 
           if (square) RES%mat(q)%gam(g_ix)%X(JX,IX) = sm
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
subroutine commutator_221(L,R,RES,jbas) 
  implicit none
  
  type(spd) :: jbas
  type(sq_op) :: L,R,RES
  integer :: i,j,q,Abody,Ntot,a,b,c,d,ix,jx,JT
  real(8) :: sm
  
  Abody = L%belowEF
  Ntot = L%Nsp
  
  do i = 1, Abody
     do j = 1,Abody
        
        ix = jbas%holes(i)
        jx = jbas%holes(j) 
        
        sm = 0.d0 
        
        do a = 1, Ntot
           do b = 1, Ntot
              do c = 1, Ntot
                 
                 do JT = 0,18,2
                 sm = sm + ( v_elem(c,ix,a,b,JT,L,jbas) * v_elem(a,b,c,jx,JT,R,jbas) + &
                 v_elem(c,jx,a,b,JT,L,jbas) * v_elem(a,b,c,ix,JT,R,jbas) ) * &
                 (( 1 - jbas%con(a)) * ( 1 - jbas%con(b)) *jbas%con(c) + &
                 jbas%con(a) * jbas%con(b) *(1 - jbas%con(c)) ) *(JT + 1)
                 end do 
              end do 
           end do
        end do 
        
        RES%fhh(i,j) = RES%fhh(i,j)  + sm /2.d0/(jbas%jj(ix) + 1.d0 ) 
        
     end do 
  end do 
  
  do i = 1, Ntot-Abody
     do j = 1,Ntot-Abody
        
        ix = jbas%parts(i)
        jx = jbas%parts(j) 
        
        sm = 0.d0 
        
        do a = 1, Ntot
           do b = 1, Ntot
              do c = 1, Ntot
                 
                 do JT = 0,18,2
                 sm = sm + ( v_elem(c,ix,a,b,JT,L,jbas) * v_elem(a,b,c,jx,JT,R,jbas) + &
                 v_elem(c,jx,a,b,JT,L,jbas) * v_elem(a,b,c,ix,JT,R,jbas) ) * &
                 (( 1 - jbas%con(a)) * ( 1 - jbas%con(b)) *jbas%con(c) + &
                 jbas%con(a) * jbas%con(b) *(1 - jbas%con(c)) ) *(JT + 1)
                 end do 
              end do 
           end do
        end do 
        
        RES%fpp(i,j) = RES%fpp(i,j)  + sm /2.d0/(jbas%jj(ix) + 1.d0 ) 
        
     end do 
  end do 
  
  do i = 1,Ntot- Abody
     do j = 1,Abody
        
        ix = jbas%parts(i)
        jx = jbas%holes(j) 
        
        
        sm = 0.d0 
        
        do a = 1, Ntot
           do b = 1, Ntot
              do c = 1, Ntot
                 
                 do JT = 0,18,2
                 sm = sm + ( v_elem(c,ix,a,b,JT,L,jbas) * v_elem(a,b,c,jx,JT,R,jbas) + &
                 v_elem(c,jx,a,b,JT,L,jbas) * v_elem(a,b,c,ix,JT,R,jbas) ) * &
                 (( 1 - jbas%con(a)) * ( 1 - jbas%con(b)) *jbas%con(c) + &
                 jbas%con(a) * jbas%con(b) *(1 - jbas%con(c)) ) * ( JT + 1) 
                 end do 
              end do 
           end do
        end do 
        
        RES%fph(i,j) = RES%fph(i,j)  + sm /2.d0/(jbas%jj(ix) + 1.d0 ) 
        
     end do 
  end do 
  
end subroutine
!===================================================================
!===================================================================
subroutine xcommutator_221(L,R,RES,w1,w2,jbas) 
  implicit none
  
  type(spd) :: jbas
  type(sq_op) :: L,R,RES,w1,w2
  integer :: i,j,q,Abody,Ntot,nh,np,nb
  real(8) :: sm
  
  Abody = L%belowEF
  Ntot = L%Nsp
  
  !construct temporary matrices
  do q = 1, L%nblocks
     
     nh = L%mat(q)%nhh
     np = L%mat(q)%npp
     nb = L%mat(q)%nph
     
     !L_hhpp . R_pphh = W1_hhhh
     call dgemm('T','N',nh,nh,np,al,L%mat(q)%gam(3)%X,np,&
          R%mat(q)%gam(3)%X,np,bet,w1%mat(q)%gam(5)%X,nh)
    
     w1%mat(q)%gam(5)%X = w1%mat(q)%gam(5)%X * L%herm 
     
     !L_phhh . R_hhph = W2_phph
     call dgemm('N','T',nb,nb,nh,al,L%mat(q)%gam(6)%X,nb,&
          R%mat(q)%gam(6)%X,nb,bet,w2%mat(q)%gam(4)%X,nb) 
     
     w2%mat(q)%gam(4)%X = w2%mat(q)%gam(4)%X * R%herm
     
     !L_phpp . R_ppph = W1_phph
     call dgemm('T','N',nb,nb,np,al,L%mat(q)%gam(2)%X,np,&
          R%mat(q)%gam(2)%X,np,bet,w1%mat(q)%gam(4)%X,nb) 
     
     w1%mat(q)%gam(2)%X = w1%mat(q)%gam(2)%X * L%herm 
     
     !L_pphh . R_hhpp = W2_pppp
     call dgemm('N','T',np,np,nh,al,L%mat(q)%gam(3)%X,np,&
          R%mat(q)%gam(3)%X,np,bet,w2%mat(q)%gam(1)%X,np)
     
     w2%mat(q)%gam(1)%X = w2%mat(q)%gam(1)%X * R%herm
     
     !L_phpp . R_pphh = W1_phhh 
     call dgemm('T','N',nb,nh,np,al,L%mat(q)%gam(2)%X,np,&
          R%mat(q)%gam(3)%X,np,bet,w1%mat(q)%gam(6)%X,nb)
     
     w1%mat(q)%gam(6)%X = w1%mat(q)%gam(6)%X * L%herm 
     
     !L_pphh . R_hhph = W2_ppph
     call dgemm('N','T',np,nb,nh,al,L%mat(q)%gam(3)%X,np,&
          R%mat(q)%gam(6)%X,nb,bet,w2%mat(q)%gam(2)%X,np)
     
     w2%mat(q)%gam(2)%X = w2%mat(q)%gam(2)%X * R%herm
     
     !R_phpp . L_pphh = W2_phhh (this is transposed) 
     call dgemm('T','N',nb,nh,np,al,R%mat(q)%gam(2)%X,np,&
          L%mat(q)%gam(3)%X,np,bet,w2%mat(q)%gam(6)%X,nb)
     
     ! I think I need to multiply by the opposite herm factor 
     ! in this case because it's all transposed
   
     w2%mat(q)%gam(6)%X = w2%mat(q)%gam(6)%X * L%herm 
     
     !R_pphh . L_hhph = W1_ppph (this is transposed) 
     call dgemm('N','T',np,nb,nh,al,R%mat(q)%gam(3)%X,np,&
          L%mat(q)%gam(6)%X,nb,bet,w1%mat(q)%gam(2)%X,np)
     
     w1%mat(q)%gam(2)%X = w1%mat(q)%gam(2)%X * R%herm
     
  end do 

end subroutine
     
end module 
  
  
  
