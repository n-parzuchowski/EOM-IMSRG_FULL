module TS_commutators
  use basic_IMSRG
  ! tensor-scalar commutator functions 
  
  ! THE TENSOR MUST BE THE SECOND ARGUMENT
  
contains
!=========================================================
!=========================================================
subroutine TS_commutator_111(L,R,RES,jbas) 
  ! one body part of [L1,R1] 
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
  
  RES%fhh = th1 + th2*L%herm  - &
       Transpose(th1*L%herm + th2) * R%herm *phase_hh  
  
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
  
  RES%fpp = tp1 * R%herm*phase_pp + tp2 - &
       Transpose(tp1 + tp2 * R%herm* phase_pp) * L%herm
    
end subroutine
!=========================================================
!=========================================================
subroutine TS_commutator_121(L,R,RES,jbas) 
  ! onebody part of [L1,R2] 
  ! this one is brute force. 
  ! not sure of a faster way to do this
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: L,R,RES
  integer :: J1,J2,PAR,TZ,ji,ja,jp,jq,a,i,p,q,g
  integer :: ti,ta,tp,tq,li,la,lp,lq,ak,ik,pk,qk,rank
  real(8) :: sm,smx,smy,smx2,smy2,d6ji
  
  rank = R%rank 
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
        if ( mod(lq,2) .ne. mod(lp,2)) cycle
        if (tq .ne. tp) cycle 
        if (.not. (triangle(jq,jp,rank))) cycle
        
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
              do J2 = abs(ji - jq),ji+jq,2
                 do J1 = abs(ja - jp),ja+jp,2 
                    
                    smx = smx + tensor_elem(ak,pk,ik,qk,J1,J2,R,jbas)&
                         *sqrt((J1 + 1.d0)*(J2 + 1.d0))*(-1)**(J1/2)*d6ji(J1,J2,rank,jq,jp,ji)             
                 
                 end do
              end do
              
              do J1 = abs(ji - jp),ji+jp,2
                 do J2 = abs(ja - jq),ja+jq,2 
                 
                    smy = smy + tensor_elem(ik,pk,ak,qk,J1,J2,R,jbas)&
                         *sqrt((J1 + 1.d0)*(J2 + 1.d0))*(-1)**(J1/2)*d6ji(J1,J2,rank,jq,jp,ji) 
                                     
                 end do
              end do
                            
              sm = sm + L%fph(a,i) * (L%herm*smx - smy) &
              *(-1)**( (rank + jq + ji )/2 ) 
              
           end do 
        end do 
        
        RES%fhh(p,q) = RES%fhh(p,q) + sm  
        RES%fhh(q,p) = RES%fhh(p,q) * RES%herm * (-1)** ( (jp - jq)/2 ) 
        
     end do 
  end do 
             
!dfpp~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
  do p = 1,L%nsp-L%belowEF
     
     pk = jbas%parts(p)
     jp = jbas%jj(pk) 
     tp = jbas%itzp(pk)
     lp = jbas%ll(pk)
     
     do q = p,L%nsp-L%belowEF
      
        qk = jbas%parts(q) 
        jq = jbas%jj(qk) 
        tq = jbas%itzp(qk)
        lq = jbas%ll(qk)
      
        ! check if this state is allowed
        if ( mod(lq,2) .ne. mod(lp,2)) cycle
        if (tq .ne. tp) cycle 
        if (.not. (triangle(jq,jp,rank))) cycle
        
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
              do J2 = abs(ji - jq),ji+jq,2
                 do J1 = abs(ja - jp),ja+jp,2 
                    
                    smx = smx + tensor_elem(ak,pk,ik,qk,J1,J2,R,jbas)&
                         *sqrt((J1 + 1.d0)*(J2 + 1.d0))*(-1)**(J1/2)&
                         * d6ji(J1,J2,rank,jq,jp,ji)              
                 
                 end do
              end do
              
              do J1 = abs(ji - jp),ji+jp,2
                 do J2 = abs(ja - jq),ja+jq,2 
                 
                    smy = smy + tensor_elem(ik,pk,ak,qk,J1,J2,R,jbas)&
                         *sqrt((J1 + 1.d0)*(J2 + 1.d0))*(-1)**(J1/2) &
                         * d6ji(J1,J2,rank,jq,jp,ji)
                                     
                 end do
              end do
                            
              sm = sm + L%fph(a,i) * (L%herm*smx - smy) &
              *(-1)**( (rank + jq + ji )/2 ) 
               
           end do 
        end do 
        
        RES%fpp(p,q) = RES%fpp(p,q) + sm  
        RES%fpp(q,p) = RES%fpp(p,q) * RES%herm * (-1)** ( (jp - jq)/2 ) 
        
     end do 
  end do 

 !dfph~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  do p = 1,L%nsp-L%belowEF
     
     pk = jbas%parts(p)
     jp = jbas%jj(pk) 
     tp = jbas%itzp(pk)
     lp = jbas%ll(pk)
    
     do q = 1,L%belowEF
        
        qk = jbas%holes(q) 
        jq = jbas%jj(qk) 
        tq = jbas%itzp(qk)
        lq = jbas%ll(qk)
       
        ! check if this state is allowed
        
        if ( mod(lq,2) .ne. mod(lp,2)) cycle
        if (tq .ne. tp) cycle 
        if (.not. (triangle(jq,jp,rank))) cycle
   
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
              ! sum over J_total
              do J2 = abs(ji - jq),ji+jq,2
                 do J1 = abs(ja - jp),ja+jp,2 
                    
                    smx = smx + tensor_elem(ak,pk,ik,qk,J1,J2,R,jbas)&
                         *sqrt((J1 + 1.d0)*(J2 + 1.d0))*(-1)**(J1/2)  &
                          * d6ji(J1,J2,rank,jq,jp,ji)
                    
                 end do
              end do
              
              do J1 = abs(ji - jp),ji+jp,2
                 do J2 = abs(ja - jq),ja+jq,2 
                 
                    smy = smy + tensor_elem(ik,pk,ak,qk,J1,J2,R,jbas)&
                         *sqrt((J1 + 1.d0)*(J2 + 1.d0))*(-1)**(J1/2) & 
                                      * d6ji(J1,J2,rank,jq,jp,ji)
                    
                 end do
              end do
                            
              sm = sm + L%fph(a,i) * (L%herm*smx - smy) &
              *(-1)**( (rank + jq + ji )/2 )
       
           end do 
        end do 
        
        RES%fph(p,q) = RES%fph(p,q) + sm  
        
     end do 
  end do 
 
end subroutine             

!====================================================================
!====================================================================
subroutine TS_commutator_211(LCC,R,RES,jbas) 
  ! onebody part of  - [R1,L2] 
  ! this one is brute force. 
  ! not sure of a faster way to do this
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: R,RES
  type(cross_coupled_31_mat) :: LCC 
  integer :: J1,J2,PAR,TZ,ji,ja,jp,jq,a,i,p,q,g,JTM
  integer :: ti,ta,tp,tq,li,la,lp,lq,ak,ik,pk,qk,rank
  integer :: rai,rpq,rqp,qx,Ntot
  real(8) :: sm,smx,smy,smx2,smy2,d6ji
  
  rank = R%rank 
  JTM = Jbas%Jtotal_max
  Ntot = R%nsp 
! dfhh ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  do p = 1,R%belowEF
     
     pk = jbas%holes(p)
     jp = jbas%jj(pk) 
     tp = jbas%itzp(pk)
     lp = jbas%ll(pk)
     
     do q = p,R%belowEF
      
        qk = jbas%holes(q) 
        jq = jbas%jj(qk) 
        tq = jbas%itzp(qk)
        lq = jbas%ll(qk)
      
        ! check if this state is allowed
        if ( mod(lq,2) .ne. mod(lp,2)) cycle
        if (tq .ne. tp) cycle 
        if (.not. (triangle(jq,jp,rank))) cycle
        
        ! internal sum
        sm = 0.d0
        do a = 1,R%Nsp - R%belowEF
           
           ak = jbas%parts(a)
           ja = jbas%jj(ak) 
           ta = jbas%itzp(ak)
           la = jbas%ll(ak)
           
           do i = 1,R%belowEF
              
              ik = jbas%holes(i)
              ji = jbas%jj(ik) 
              ti = jbas%itzp(ik)
              li = jbas%ll(ik)
              
              ! check if this intermediate exists
              if ( mod(li,2) .ne. mod(la,2)) cycle
              if (ti .ne. ta) cycle 
              if (.not. (triangle(ja,ji,rank))) cycle
               
              smx = 0.d0 
              smy = 0.d0
              smx2 = 0.d0 
              smy2 = 0.d0 

              Tz = abs(ta -ti)/2 
              if (abs(tp - tq) .ne. Tz*2)  cycle 
              PAR = mod(la+li,2) 
              if (mod(lp+lq,2) .ne. PAR) cycle                                 

              qx = rank/2+1 + Tz*(JTM+1) + 2*PAR*(JTM+1)
                
              rai = ph_rval(ak,ik,Ntot,qx,LCC)
              rpq = TS_rval(pk,qk,Ntot,qx,LCC)
              rqp = TS_rval(qk,pk,Ntot,qx,LCC)

              sm = sm + (-1)**(( jp + jq + rank)/2) * R%fph(a,i) & 
                *  ( R%herm *(-1)**((ja-ji)/2)*LCC%CCR(qx)%X(rai,rpq) &
                - LCC%CCX(qx)%X(rpq,rai) ) / sqrt(rank + 1.d0 ) 

              ! the last (rank + 1) is divided out because
              ! the CC matrix elements are scaled by that, 
              ! which is wrong here. 
         
           end do 
        end do 
        
        RES%fhh(p,q) = RES%fhh(p,q) + sm  
        RES%fhh(q,p) = RES%fhh(p,q) * RES%herm * (-1)** ( (jp - jq)/2 ) 
        
     end do 
  end do 

! dfpp ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  do p = 1,R%nsp-R%belowEF
     
     pk = jbas%parts(p)
     jp = jbas%jj(pk) 
     tp = jbas%itzp(pk)
     lp = jbas%ll(pk)
     
     do q = p,R%nsp-R%belowEF
      
        qk = jbas%parts(q) 
        jq = jbas%jj(qk) 
        tq = jbas%itzp(qk)
        lq = jbas%ll(qk)
      
        ! check if this state is allowed
        if ( mod(lq,2) .ne. mod(lp,2)) cycle
        if (tq .ne. tp) cycle 
        if (.not. (triangle(jq,jp,rank))) cycle
        
        ! internal sum
        sm = 0.d0
        do a = 1,R%Nsp - R%belowEF
           
           ak = jbas%parts(a)
           ja = jbas%jj(ak) 
           ta = jbas%itzp(ak)
           la = jbas%ll(ak)
           
           do i = 1,R%belowEF
              
              ik = jbas%holes(i)
              ji = jbas%jj(ik) 
              ti = jbas%itzp(ik)
              li = jbas%ll(ik)
              
              ! check if this intermediate exists
              if ( mod(li,2) .ne. mod(la,2)) cycle
              if (ti .ne. ta) cycle 
              if (.not. (triangle(ja,ji,rank))) cycle
               
              smx = 0.d0 
              smy = 0.d0
              smx2 = 0.d0 
              smy2 = 0.d0 

              Tz = abs(ta -ti)/2 
              if (abs(tp - tq) .ne. Tz*2)  cycle 
              PAR = mod(la+li,2) 
              if (mod(lp+lq,2) .ne. PAR) cycle                                 


              qx = rank/2+1 + Tz*(JTM+1) + 2*PAR*(JTM+1)
                
              rai = ph_rval(ak,ik,Ntot,qx,LCC)
              rpq = TS_rval(pk,qk,Ntot,qx,LCC)
              rqp = TS_rval(qk,pk,Ntot,qx,LCC)
 
              sm = sm + (-1)**(( jp + jq + rank)/2) * R%fph(a,i) & 
                *  ( R%herm *(-1)**((ja-ji)/2)*LCC%CCR(qx)%X(rai,rpq) &
                - LCC%CCX(qx)%X(rpq,rai) ) / sqrt(rank + 1.d0 ) 
              ! the last (rank + 1) is divided out because
              ! the CC matrix elements are scaled by that, 
              ! which is wrong here. 
         
           end do 
        end do 
        
        RES%fpp(p,q) = RES%fpp(p,q) + sm  
        RES%fpp(q,p) = RES%fpp(p,q) * RES%herm * (-1)** ( (jp - jq)/2 ) 
        
     end do 
  end do 

! dfph ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  do p = 1,R%nsp-R%belowEF
     
     pk = jbas%parts(p)
     jp = jbas%jj(pk) 
     tp = jbas%itzp(pk)
     lp = jbas%ll(pk)
     
     do q = 1,R%belowEF
      
        qk = jbas%holes(q) 
        jq = jbas%jj(qk) 
        tq = jbas%itzp(qk)
        lq = jbas%ll(qk)
      
        ! check if this state is allowed
        if ( mod(lq,2) .ne. mod(lp,2)) cycle
        if (tq .ne. tp) cycle 
        if (.not. (triangle(jq,jp,rank))) cycle
        
        ! internal sum
        sm = 0.d0
        do a = 1,R%Nsp - R%belowEF
           
           ak = jbas%parts(a)
           ja = jbas%jj(ak) 
           ta = jbas%itzp(ak)
           la = jbas%ll(ak)
           
           do i = 1,R%belowEF
              
              ik = jbas%holes(i)
              ji = jbas%jj(ik) 
              ti = jbas%itzp(ik)
              li = jbas%ll(ik)
              
              ! check if this intermediate exists
              if ( mod(li,2) .ne. mod(la,2)) cycle
              if (ti .ne. ta) cycle 
              if (.not. (triangle(ja,ji,rank))) cycle
               
              smx = 0.d0 
              smy = 0.d0
              smx2 = 0.d0 
              smy2 = 0.d0 

              Tz = abs(ta -ti)/2 
              if (abs(tp - tq) .ne. Tz*2)  cycle 
              PAR = mod(la+li,2) 
              if (mod(lp+lq,2) .ne. PAR) cycle                                 

              qx = rank/2+1 + Tz*(JTM+1) + 2*PAR*(JTM+1)
                
              rai = ph_rval(ak,ik,Ntot,qx,LCC)
              rpq = TS_rval(pk,qk,Ntot,qx,LCC)
              rqp = TS_rval(qk,pk,Ntot,qx,LCC)

              sm = sm + (-1)**(( jp + jq + rank)/2) * R%fph(a,i) & 
                *  ( R%herm *(-1)**((ja-ji)/2)*LCC%CCR(qx)%X(rai,rpq) &
                - LCC%CCX(qx)%X(rpq,rai) ) / sqrt(rank + 1.d0 ) 

              ! the last (rank + 1) is divided out because
              ! the CC matrix elements are scaled by that, 
              ! which is wrong here. 
         
           end do 
        end do 
        
        RES%fph(p,q) = RES%fph(p,q) + sm  
        
     end do 
  end do 

end subroutine
!==================================================
!==================================================             
subroutine TS_commutator_122(L,R,RES,jbas) 
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: L,R,RES
  integer :: q,IX,JX,nh,np,nb,i,J1,J2
  integer :: a,b,c,d,ja,jb,jc,jd,ji,g_ix,q_sp,i_sp
  integer :: ta,tb,tc,td,ti,la,lb,lc,ld,li,spec
  integer :: jxstart,jxend,ixend,c1,c2,n1,n2
  logical :: square
  real(8) ::  sm
  
  
  do q = 1, R%nblocks
     
     J1 = R%tblck(q)%Jpair(1)
     J2 = R%tblck(q)%Jpair(2) 
    
     do g_ix = 1,9 
   
        ! figure out how big the array is
        n1 = size(R%tblck(q)%tgam(g_ix)%X(:,1))
        n2 = size(R%tblck(q)%tgam(g_ix)%X(1,:))
        if ((n1*n2) == 0) cycle 
        
        ! read in information about which 
        ! array we are using from public arrays
        c1 = sea1(g_ix) 
        c2 = sea2(g_ix) 
        square = sqs(g_ix) 
        jxstart = jst(g_ix) 
        
                
     ! main calculation
   
     do IX = 1,n1
        a = R%tblck(q)%tensor_qn(c1,1)%Y(IX,1)
        ja = jbas%jj(a)
        la = jbas%ll(a)
        ta = jbas%itzp(a) 
             
        b = R%tblck(q)%tensor_qn(c1,1)%Y(IX,2)
        jb = jbas%jj(b)
        lb = jbas%ll(b)
        tb = jbas%itzp(b)
 
        do JX = 1,n2
           
           c = R%tblck(q)%tensor_qn(c2,2)%Y(JX,1)
           jc = jbas%jj(c)
           lc = jbas%ll(c)
           tc = jbas%itzp(c)

           d = R%tblck(q)%tensor_qn(c2,2)%Y(JX,2)
           jd = jbas%jj(d)
           ld = jbas%ll(d)
           td = jbas%itzp(d)
                   
           sm = 0.d0 

            ! a is replaced
            q_sp = sp_block_index(ja,la,ta,jbas) 
            do i = 1,size(jbas%states(q_sp)%Z)   
              
               i_sp = jbas%states(q_sp)%Z(i) 
               
               sm = sm + f_elem(a,i_sp,L,jbas)*tensor_elem(i_sp,b,c,d,J1,J2,R,jbas)

            end do 
              
            ! b is replaced
            q_sp = sp_block_index(jb,lb,tb,jbas) 
            do i = 1,size(jbas%states(q_sp)%Z)   
              
               i_sp = jbas%states(q_sp)%Z(i) 
               
               sm = sm + f_elem(b,i_sp,L,jbas)*tensor_elem(a,i_sp,c,d,J1,J2,R,jbas) 
            end do 
            
            ! c is replaced
            q_sp = sp_block_index(jc,lc,tc,jbas) 
            do i = 1,size(jbas%states(q_sp)%Z)   
              
               i_sp = jbas%states(q_sp)%Z(i) 
               
               sm = sm - f_elem(i_sp,c,L,jbas)*tensor_elem(a,b,i_sp,d,J1,J2,R,jbas)
            end do 
            
            ! d is replaced
            q_sp = sp_block_index(jd,ld,td,jbas) 
            do i = 1,size(jbas%states(q_sp)%Z)   
              
               i_sp = jbas%states(q_sp)%Z(i) 
               
               sm = sm - f_elem(i_sp,d,L,jbas)*tensor_elem(a,b,c,i_sp,J1,J2,R,jbas)
            end do 
          
            sm = sm / sqrt(1.d0 + kron_del(a,b)) /sqrt(1.d0 + kron_del(c,d)) 

           RES%tblck(q)%tgam(g_ix)%X(IX,JX) = sm 

        end do
     end do
   
     end do 
  end do 

end subroutine           
!===================================================================
!===================================================================
subroutine TS_commutator_221(w1,w2,pm,RES,jbas) 
  ! verified
  ! THIS NEEDS TO BE RUN AFTER 222_pp_hh 
  ! 222_pp_hh sets up the intermediary matrices (w1,w2) 
  implicit none
  
  type(spd) :: jbas
  type(sq_op) :: L,R,RES,w1,w2
  integer :: i,j,q,Abody,Ntot,nh,np,nb,a,c
  integer :: ik,jk,ck,ji,jj,ti,tj,li,lj,jc,J1,J2
  integer,intent(in) :: pm
  real(8) :: sm,sm1,sm2,d6ji
  
  Abody = w1%belowEF
  Ntot = w1%Nsp

  ! fpp
  do ik = 1 , Ntot - Abody
     i = jbas%parts(ik) 
     ji = jbas%jj(i) 
     li = jbas%ll(i) 
     ti = jbas%itzp(i) 
     
     do jk = ik , Ntot - Abody
        
        j = jbas%parts(jk) 
        jj = jbas%jj(j) 
        if (.not. (triangle(jj,ji,w1%rank))) cycle
        lj = jbas%ll(j) 
        if (mod(li,2) .ne. mod(lj,2))  cycle
        tj = jbas%itzp(j)
        if (tj .ne. ti) cycle 
                
        sm = 0.d0 
      
        do ck = 1, Abody
           c = jbas%holes(ck) 
           jc = jbas%jj(c)
           ! w1 matrix results from multiplying the pp channel
           sm1 = 0.d0 
           do J1 = abs(jc - ji),jc+ji,2
            
              ! NOTE: 
              ! THESE SUMS HAVE TO BE BROKEN UP SO the J on the left side is 
              ! smaller. I don't have the other matrix multiplication.
              do J2 = abs(jc - jj),min(J1-2,jc+jj),2

                ! use w1, because it sums over the pp indices
                sm1 = sm1 - sqrt((J1+1.d0)*(J2+1.d0))*d6ji(J1,J2,w1%rank,jj,ji,jc) &
                 *tensor_elem(c,j,c,i,J2,J1,w1,jbas)*(-1)**(J2/2) * pm 
                
             end do              
          
             do J2 = max(J1,abs(jc - jj)),jc+jj,2
             
                ! use w1, because it sums over the pp indices
                sm1 = sm1 + sqrt((J1+1.d0)*(J2+1.d0))*d6ji(J1,J2,w1%rank,jj,ji,jc) &
                     *tensor_elem(c,i,c,j,J1,J2,w1,jbas)*(-1)**(J1/2)
                
             end do

          end do
          sm = sm + sm1*(-1)**((jc+jj)/2)
        end do 
        

        do ck = 1, Ntot - Abody
           c = jbas%parts(ck) 
           jc = jbas%jj(c)
           sm2 = 0.d0
           do J1 = abs(jc - ji),jc+ji,2
             do J2 = abs(jc - jj),min(J1-2,jc+jj),2

                ! use w1, because it sums over the pp indices
                sm2 = sm2 - sqrt((J1+1.d0)*(J2+1.d0))*d6ji(J1,J2,w1%rank,jj,ji,jc) &
                 *tensor_elem(c,j,c,i,J2,J1,w2,jbas)*(-1)**(J2/2) * pm 
                
             end do              
          
             do J2 = max(J1,abs(jc - jj)),jc+jj,2
           
                ! use w1, because it sums over the pp indices
                sm2 = sm2 + sqrt((J1+1.d0)*(J2+1.d0))*d6ji(J1,J2,w1%rank,jj,ji,jc) &
                     *tensor_elem(c,i,c,j,J1,J2,w2,jbas)*(-1)**(J1/2)
                
             end do
            
           end do 
           sm = sm + sm2 *  (-1)**((jc+jj)/2)
        end do 
     
        RES%fpp(ik,jk) = RES%fpp(ik,jk) + sm * (-1)**(w1%rank/2)  
        RES%fpp(jk,ik) = RES%fpp(ik,jk) * RES%herm * (-1)**((ji-jj)/2) 
     end do 
  end do       


  ! fhh
  do ik = 1 , Abody
     i = jbas%holes(ik) 
     ji = jbas%jj(i) 
     li = jbas%ll(i) 
     ti = jbas%itzp(i) 
     
     do jk = ik , Abody
        
        j = jbas%holes(jk) 
        jj = jbas%jj(j) 
        if (.not. (triangle(jj,ji,w1%rank))) cycle
        lj = jbas%ll(j) 
        if (mod(li,2) .ne. mod(lj,2))  cycle
        tj = jbas%itzp(j)
        if (tj .ne. ti) cycle 
                
        sm = 0.d0 
      
        do ck = 1, Abody
           c = jbas%holes(ck) 
           jc = jbas%jj(c)
           ! w1 matrix results from multiplying the pp channel
           sm1 = 0.d0 
           do J1 = abs(jc - ji),jc+ji,2
            
              ! NOTE: 
              ! THESE SUMS HAVE TO BE BROKEN UP SO the J on the left side is 
              ! smaller. I don't have the other matrix multiplication.
              do J2 = abs(jc - jj),min(J1-2,jc+jj),2

                ! use w1, because it sums over the pp indices
                sm1 = sm1 - sqrt((J1+1.d0)*(J2+1.d0))*d6ji(J1,J2,w1%rank,jj,ji,jc) &
                 *tensor_elem(c,j,c,i,J2,J1,w1,jbas)*(-1)**(J2/2) * pm 
                
             end do              
          
             do J2 = max(J1,abs(jc - jj)),jc+jj,2
             
                ! use w1, because it sums over the pp indices
                sm1 = sm1 + sqrt((J1+1.d0)*(J2+1.d0))*d6ji(J1,J2,w1%rank,jj,ji,jc) &
                     *tensor_elem(c,i,c,j,J1,J2,w1,jbas)*(-1)**(J1/2)
                
             end do

          end do
          sm = sm + sm1*(-1)**((jc+jj)/2)
        end do 
        

        do ck = 1, Ntot - Abody
           c = jbas%parts(ck) 
           jc = jbas%jj(c)
           sm2 = 0.d0
           do J1 = abs(jc - ji),jc+ji,2
             do J2 = abs(jc - jj),min(J1-2,jc+jj),2

                ! use w1, because it sums over the pp indices
                sm2 = sm2 - sqrt((J1+1.d0)*(J2+1.d0))*d6ji(J1,J2,w1%rank,jj,ji,jc) &
                 *tensor_elem(c,j,c,i,J2,J1,w2,jbas)*(-1)**(J2/2) * pm 
                
             end do              
          
             do J2 = max(J1,abs(jc - jj)),jc+jj,2
           
                ! use w1, because it sums over the pp indices
                sm2 = sm2 + sqrt((J1+1.d0)*(J2+1.d0))*d6ji(J1,J2,w1%rank,jj,ji,jc) &
                     *tensor_elem(c,i,c,j,J1,J2,w2,jbas)*(-1)**(J1/2)
                
             end do
            
           end do 
           sm = sm + sm2 *  (-1)**((jc+jj)/2)
        end do 
     
        RES%fhh(ik,jk) = RES%fhh(ik,jk) + sm * (-1)**(w1%rank/2)  
        RES%fhh(jk,ik) = RES%fhh(ik,jk) * RES%herm * (-1)**((ji-jj)/2) 
     end do 
  end do       
  
  ! fph
  do ik = 1 , Ntot-Abody
     i = jbas%parts(ik) 
     ji = jbas%jj(i) 
     li = jbas%ll(i) 
     ti = jbas%itzp(i) 
    
     do jk = 1 , Abody
        
        j = jbas%holes(jk) 
        jj = jbas%jj(j) 
        if (.not. (triangle(jj,ji,w1%rank))) cycle
        lj = jbas%ll(j) 
        if (mod(li,2) .ne. mod(lj,2))  cycle
        tj = jbas%itzp(j)
        if (tj .ne. ti) cycle 
      
    
        sm = 0.d0 
      
        do ck = 1, Abody
           c = jbas%holes(ck) 
           jc = jbas%jj(c)
           ! w1 matrix results from multiplying the pp channel
           sm1 = 0.d0 
           do J1 = abs(jc - ji),jc+ji,2
            
              ! NOTE: 
              ! THESE SUMS HAVE TO BE BROKEN UP SO the J on the left side is 
              ! smaller. I don't have the other matrix multiplication.
              do J2 = abs(jc - jj),min(J1-2,jc+jj),2

                ! use w1, because it sums over the pp indices
                sm1 = sm1 - sqrt((J1+1.d0)*(J2+1.d0))*d6ji(J1,J2,w1%rank,jj,ji,jc) &
                 *tensor_elem(c,j,c,i,J2,J1,w1,jbas)*(-1)**(J2/2) * pm 
                
             end do              
          
             do J2 = max(J1,abs(jc - jj)),jc+jj,2
             
                ! use w1, because it sums over the pp indices
                sm1 = sm1 + sqrt((J1+1.d0)*(J2+1.d0))*d6ji(J1,J2,w1%rank,jj,ji,jc) &
                     *tensor_elem(c,i,c,j,J1,J2,w1,jbas)*(-1)**(J1/2)
                
             end do

          end do
          sm = sm + sm1*(-1)**((jc+jj)/2)
        end do 
        

        do ck = 1, Ntot - Abody
           c = jbas%parts(ck) 
           jc = jbas%jj(c)
           sm2 = 0.d0
           do J1 = abs(jc - ji),jc+ji,2
             do J2 = abs(jc - jj),min(J1-2,jc+jj),2

                ! use w1, because it sums over the pp indices
                sm2 = sm2 - sqrt((J1+1.d0)*(J2+1.d0))*d6ji(J1,J2,w1%rank,jj,ji,jc) &
                 *tensor_elem(c,j,c,i,J2,J1,w2,jbas)*(-1)**(J2/2) * pm 
                
             end do              
          
             do J2 = max(J1,abs(jc - jj)),jc+jj,2
           
                ! use w1, because it sums over the pp indices
                sm2 = sm2 + sqrt((J1+1.d0)*(J2+1.d0))*d6ji(J1,J2,w1%rank,jj,ji,jc) &
                     *tensor_elem(c,i,c,j,J1,J2,w2,jbas)*(-1)**(J1/2)
                
             end do
            
           end do 
           sm = sm + sm2 *  (-1)**((jc+jj)/2)
        end do 
     
        RES%fph(ik,jk) = RES%fph(ik,jk) + sm * (-1)**(w1%rank/2)  

     end do 
  end do       

end subroutine
!===================================================================
!===================================================================
subroutine TS_commutator_222_pp_hh(L,R,RES,w1,w2,jbas) 
  !VERIFIED
  !NEEDS TO BE RUN BEFORE 221, because it sets up the 
  !intermediary matrices
  implicit none
  
  type(spd) :: jbas
  type(sq_op) ::  L,R,RES,w1,w2
  integer :: q,q1,q2,J1,J2,Tz,Par,phase
  integer :: np1,nb1,nh1,np2,nb2,nh2,pm
  real(8) :: bet_off,al_off
  
  pm = R%herm*L%herm
  
!construct temporary matrices
  do q = 1, R%nblocks
     
     J1 = R%tblck(q)%Jpair(1) 
     J2 = R%tblck(q)%Jpair(2)
     phase = R%tblck(q)%lam(1)
     par = R%tblck(q)%lam(2) 
     Tz = R%tblck(q)%lam(3)
    
     q1 = block_index(J1,Tz,Par) 
     q2 = block_index(J2,Tz,Par) 
     
     nh1 = R%tblck(q)%nhh1
     np1 = R%tblck(q)%npp1
     nb1 = R%tblck(q)%nph1
     nh2 = R%tblck(q)%nhh2
     np2 = R%tblck(q)%npp2
     nb2 = R%tblck(q)%nph2
       

!----------------------------------------------------------------------------
!         Zpppp 
!----------------------------------------------------------------------------
     if (np1*np2 .ne. 0)  then 
     
        !w1pppp = Bpppp.Apppp 
        call dgemm('N','N',np1,np2,np2,al,R%tblck(q)%tgam(1)%X,np1,&
             L%mat(q2)%gam(1)%X,np2,bet,w1%tblck(q)%tgam(1)%X,np1)
       
        !w1pppp = Apppp.Bpppp - Bpppp.Apppp
        bet_off = -1
        call dgemm('N','N',np1,np2,np1,al,L%mat(q1)%gam(1)%X,np1,&
             R%tblck(q)%tgam(1)%X,np1,bet_off,w1%tblck(q)%tgam(1)%X,np1)

        if (nh2 .ne. 0) then 
        
           !w2pppp = -Bpphh.Ahhpp  
           
           al_off = -1*L%herm   ! I don't have the a/h.c. of L so I need
           ! to multiply by L%herm, and transpose in dgemm. 
           call dgemm('N','T',np1,np2,nh2,al_off,R%tblck(q)%tgam(3)%X,np1,&
                L%mat(q2)%gam(3)%X,np2,bet,w2%tblck(q)%tgam(1)%X,np1)
        end if
        
        if (nh1 .ne. 0) then
           ! notice I have only a part of the R matrix being used here
           ! this is because I have both pphh and hhpp parts stored
           ! for the J1,J2 orientation. flipping across the aisle 
           ! gives them for the J2,J1 orientation. It's a pain.  
        
           !w2pppp = Apphh.Bhhpp - Bpphh.Ahhpp  
           bet_off = 1
           call dgemm('N','N',np1,np2,nh1,al,L%mat(q1)%gam(3)%X,np1,&
                R%tblck(q)%tgam(7)%X,nh1,bet_off,w2%tblck(q)%tgam(1)%X,np1)
        end if 
        
     end if
     
     RES%tblck(q)%tgam(1)%X = w1%tblck(q)%tgam(1)%X - w2%tblck(q)%tgam(1)%X


!----------------------------------------------------------------------------
!         Zphph 
!----------------------------------------------------------------------------
         
     if (nb1*nb2 .ne. 0)  then 
        
        if (np2 .ne. 0) then
           !w1phph = -Bphpp.Appph 
           al_off = -1
           call dgemm('N','N',nb1,nb2,np2,al_off,R%tblck(q)%tgam(8)%X,nb1,&
                L%mat(q2)%gam(2)%X,np2,bet,w1%tblck(q)%tgam(4)%X,nb1)
        
        end if
        
        if (np1 .ne. 0 ) then 
           !w1phph = Aphpp.Bppph - Bphpp.Appph
           al_off = L%herm
           bet_off = 1
           call dgemm('T','N',nb1,nb2,np1,al_off,L%mat(q1)%gam(2)%X,np1,&
                R%tblck(q)%tgam(2)%X,np1,bet_off,w1%tblck(q)%tgam(4)%X,nb1)
        end if 
             
         if (nh1 .ne. 0)  then      
            !w2phph = Aphhh.Bhhph
            call dgemm('N','N',nb1,nb2,nh1,al,L%mat(q1)%gam(6)%X,nb1,&
                 R%tblck(q)%tgam(9)%X,nh1,bet,w2%tblck(q)%tgam(4)%X,nb1)        
         end if 
        
        if (nh2 .ne. 0) then
           !w2phph = Aphhh.Bhhph - Bphhh.Ahhph 
           al_off = -1*L%herm
           bet_off = 1
           call dgemm('N','T',nb1,nb2,nh2,al_off,R%tblck(q)%tgam(6)%X,nb1,&
                L%mat(q2)%gam(6)%X,nb2,bet_off,w2%tblck(q)%tgam(4)%X,nb1)           

        end if

        RES%tblck(q)%tgam(4)%X = w1%tblck(q)%tgam(4)%X - w2%tblck(q)%tgam(4)%X
         
    end if 


!----------------------------------------------------------------------------
!         Zhhhh 
!----------------------------------------------------------------------------
     if (nh1*nh2 .ne. 0)  then 
     
        if (np2 .ne. 0) then 
           !w1hhhh = -Bhhpp.Apphh 
           al_off = -1
           call dgemm('N','N',nh1,nh2,np2,al_off,R%tblck(q)%tgam(7)%X,nh1,&
                L%mat(q2)%gam(3)%X,np2,bet,w1%tblck(q)%tgam(5)%X,nh1)
        end if 
        
        if (np1 .ne. 0) then 
           !w1hhhh = Ahhpp.Bpphh - Bhhpp.Apphh
           al_off = L%herm 
           bet_off = 1
           call dgemm('T','N',nh1,nh2,np1,al_off,L%mat(q1)%gam(3)%X,np1,&
                R%tblck(q)%tgam(3)%X,np1,bet_off,w1%tblck(q)%tgam(5)%X,nh1)
        end if 
        
        !w1hhhh = Bhhhh.Ahhhh 
        call dgemm('N','N',nh1,nh2,nh2,al,R%tblck(q)%tgam(5)%X,nh1,&
             L%mat(q2)%gam(5)%X,nh2,bet,w2%tblck(q)%tgam(5)%X,nh1)

        bet_off = -1 
        !w1hhhh = Ahhhh.Bhhhh - Bhhhh.Ahhhh 
        call dgemm('N','N',nh1,nh2,nh1,al,L%mat(q1)%gam(5)%X,nh1,&
             R%tblck(q)%tgam(5)%X,nh1,bet_off,w2%tblck(q)%tgam(5)%X,nh1)
     end if
        
     RES%tblck(q)%tgam(5)%X = w1%tblck(q)%tgam(5)%X - w2%tblck(q)%tgam(5)%X
     

!----------------------------------------------------------------------------
!         Zpphh 
!----------------------------------------------------------------------------
     if (np1*nh2 .ne. 0)  then 
     
        if (np2 .ne. 0) then 
           !w1pphh = Bpppp.Apphh 
           call dgemm('N','N',np1,nh2,np2,al,R%tblck(q)%tgam(1)%X,np1,&
                L%mat(q2)%gam(3)%X,np2,bet,w1%tblck(q)%tgam(3)%X,np1)
        end if 

        
        !w1pphh = Apppp.Bpphh - Bpppp.Apphh
        bet_off = -1
        call dgemm('N','N',np1,nh2,np1,al,L%mat(q1)%gam(1)%X,np1,&
             R%tblck(q)%tgam(3)%X,np1,bet_off,w1%tblck(q)%tgam(3)%X,np1)
        

        
        !w1pphh = -Bpphh.Ahhhh 
        al_off = -1 
        call dgemm('N','N',np1,nh2,nh2,al_off,R%tblck(q)%tgam(3)%X,np1,&
             L%mat(q2)%gam(5)%X,nh2,bet,w2%tblck(q)%tgam(3)%X,np1)
         
        
        if (nh1 .ne. 0 ) then 
           !w1pphh = Apphh.Bhhhh - Bpphh.Ahhhh 
           bet_off = 1
           call dgemm('N','N',np1,nh2,nh1,al,L%mat(q1)%gam(3)%X,np1,&
                R%tblck(q)%tgam(5)%X,nh1,bet_off,w2%tblck(q)%tgam(3)%X,np1)
        end if 
     end if 
       
     RES%tblck(q)%tgam(3)%X = w1%tblck(q)%tgam(3)%X - w2%tblck(q)%tgam(3)%X

!----------------------------------------------------------------------------
!         Zhhpp 
!----------------------------------------------------------------------------
     if (np2*nh1 .ne. 0)  then 
     
        al_off = -1
        !w1hhpp = -Bhhpp.Apppp 
        call dgemm('N','N',nh1,np2,np2,al_off,R%tblck(q)%tgam(7)%X,nh1,&
             L%mat(q2)%gam(1)%X,np2,bet,w1%tblck(q)%tgam(7)%X,nh1)
        

        if (np1 .ne. 0 ) then
           !w1hhpp = Ahhpp.Bpppp - Bhhpp.Apppp
           al_off = L%herm
           bet_off = 1
           call dgemm('T','N',nh1,np2,np1,al_off,L%mat(q1)%gam(3)%X,np1,&
                R%tblck(q)%tgam(1)%X,np1,bet_off,w1%tblck(q)%tgam(7)%X,nh1)
        end if 

        if (nh2 .ne. 0) then 
           !w1hhpp = Bhhhh.Ahhpp
           al_off = L%herm
           call dgemm('N','T',nh1,np2,nh2,al_off,R%tblck(q)%tgam(5)%X,nh1,&
                L%mat(q2)%gam(3)%X,np2,bet,w2%tblck(q)%tgam(7)%X,nh1)
        end if 
        
        bet_off = -1 
        !w1pphh = Ahhhh.Bhhpp - Bhhhh.Ahhpp 
        call dgemm('N','N',nh1,np2,nh1,al,L%mat(q1)%gam(5)%X,nh1,&
             R%tblck(q)%tgam(7)%X,nh1,bet_off,w2%tblck(q)%tgam(7)%X,nh1)

     end if 
       
     RES%tblck(q)%tgam(7)%X = w1%tblck(q)%tgam(7)%X - w2%tblck(q)%tgam(7)%X

!----------------------------------------------------------------------------
!         Zppph 
!----------------------------------------------------------------------------
     if (np1*nb2 .ne. 0)  then 
     
        if (np2 .ne. 0)  then 
           al_off = -1
           !w1ppph = -Bpppp.Appph 
           call dgemm('N','N',np1,nb2,np2,al_off,R%tblck(q)%tgam(1)%X,np1,&
                L%mat(q2)%gam(2)%X,np2,bet,w1%tblck(q)%tgam(2)%X,np1)
        end if


        !w1ppph = Apppp.Bppph - Bpppp.Appph
        bet_off = 1
        call dgemm('N','N',np1,nb2,np1,al,L%mat(q1)%gam(1)%X,np1,&
             R%tblck(q)%tgam(2)%X,np1,bet_off,w1%tblck(q)%tgam(2)%X,np1)


        if (nh2 .ne. 0) then 
           !w2ppph = -Bpphh.Ahhph
           al_off = -1*L%herm
           call dgemm('N','T',np1,nb2,nh2,al_off,R%tblck(q)%tgam(3)%X,np1,&
                L%mat(q2)%gam(6)%X,nb2,bet,w2%tblck(q)%tgam(2)%X,np1)
        end if 
        
        if (nh1 .ne. 0) then 
           bet_off = 1 
           !w2ppph = Apphh.Bhhph - Bpphh.Ahhph 
           call dgemm('N','N',np1,nb2,nh1,al,L%mat(q1)%gam(3)%X,np1,&
                R%tblck(q)%tgam(9)%X,nh1,bet_off,w2%tblck(q)%tgam(2)%X,np1)
        end if 
        
     end if 
       
     RES%tblck(q)%tgam(2)%X = w1%tblck(q)%tgam(2)%X - w2%tblck(q)%tgam(2)%X

!----------------------------------------------------------------------------
!         Zphpp 
!----------------------------------------------------------------------------
     if (nb1*np2 .ne. 0)  then 
     
       
        al_off = -1
        !w1phpp = -Bphpp.Apppp 
        call dgemm('N','N',nb1,np2,np2,al_off,R%tblck(q)%tgam(8)%X,nb1,&
             L%mat(q2)%gam(1)%X,np2,bet,w1%tblck(q)%tgam(8)%X,nb1)
        

        if (np1 .ne. 0) then 
           !w1phpp = Aphpp.Bpppp - Bphpp.Apppp
           bet_off = 1
           al_off = L%herm
           call dgemm('T','N',nb1,np2,np1,al_off,L%mat(q1)%gam(2)%X,np1,&
                R%tblck(q)%tgam(1)%X,np1,bet_off,w1%tblck(q)%tgam(8)%X,nb1)
        end if 

        if (nh2 .ne. 0) then 
           !w2phpp = -Bphhh.Ahhpp
           al_off = -1*L%herm
           call dgemm('N','T',nb1,np2,nh2,al_off,R%tblck(q)%tgam(6)%X,nb1,&
                L%mat(q2)%gam(3)%X,np2,bet,w2%tblck(q)%tgam(8)%X,nb1)
        end if 
        
        if (nh1 .ne. 0) then 
           bet_off = 1 
           !w2phpp = Aphhh.Bhhpp - Bphhh.Ahhpp 
           call dgemm('N','N',nb1,np2,nh1,al,L%mat(q1)%gam(6)%X,nb1,&
                R%tblck(q)%tgam(7)%X,nh1,bet_off,w2%tblck(q)%tgam(8)%X,nb1)
        end if 
        
     end if 
       
     RES%tblck(q)%tgam(8)%X = w1%tblck(q)%tgam(8)%X - w2%tblck(q)%tgam(8)%X


!----------------------------------------------------------------------------
!         Zphhh 
!----------------------------------------------------------------------------
     if (nb1*nh2 .ne. 0)  then 
     
        if ( np2 .ne. 0 ) then 
           al_off = -1
           !w1phhh = -Bphpp.Apphh 
           call dgemm('N','N',nb1,nh2,np2,al_off,R%tblck(q)%tgam(8)%X,nb1,&
                L%mat(q2)%gam(3)%X,np2,bet,w1%tblck(q)%tgam(6)%X,nb1)
        end if 

        if (np1 .ne. 0) then 
           !w1phhh = Aphpp.Bpphh - Bphpp.Apphh
           bet_off = 1
           al_off = L%herm
           call dgemm('T','N',nb1,nh2,np1,al_off,L%mat(q1)%gam(2)%X,np1,&
                R%tblck(q)%tgam(3)%X,np1,bet_off,w1%tblck(q)%tgam(6)%X,nb1)
        end if 

       
        !w2phhh = -Bphhh.Ahhhh
        al_off = -1
        call dgemm('N','N',nb1,nh2,nh2,al_off,R%tblck(q)%tgam(6)%X,nb1,&
             L%mat(q2)%gam(5)%X,nh2,bet,w2%tblck(q)%tgam(6)%X,nb1)
       
        
        if (nh1 .ne. 0) then 
           bet_off = 1 
           !w2phhh = Aphhh.Bhhhh - Bphhh.Ahhhh 
           call dgemm('N','N',nb1,nh2,nh1,al,L%mat(q1)%gam(6)%X,nb1,&
                R%tblck(q)%tgam(5)%X,nh1,bet_off,w2%tblck(q)%tgam(6)%X,nb1)
        end if 
        
     end if 
       
     RES%tblck(q)%tgam(6)%X = w1%tblck(q)%tgam(6)%X - w2%tblck(q)%tgam(6)%X

!----------------------------------------------------------------------------
!         Zhhph 
!----------------------------------------------------------------------------
     if (nh1*nb2 .ne. 0)  then 
     
        if (np2 .ne. 0)  then 
           al_off = -1
           !w1hhph = -Bhhpp.Appph 
           call dgemm('N','N',nh1,nb2,np2,al_off,R%tblck(q)%tgam(7)%X,nh1,&
                L%mat(q2)%gam(2)%X,np2,bet,w1%tblck(q)%tgam(9)%X,nh1)
        end if

        if (np1 .ne. 0) then
           !w1hhph = Ahhpp.Bppph - Bhhpp.Appph
           bet_off = 1
           al_off = L%herm
           call dgemm('T','N',nh1,nb2,np1,al_off,L%mat(q1)%gam(3)%X,np1,&
                R%tblck(q)%tgam(2)%X,np1,bet_off,w1%tblck(q)%tgam(9)%X,nh1)
        end if 

        if (nh2 .ne. 0) then 
           !w2hhph = -Bhhhh.Ahhph
           al_off = -1*L%herm
           call dgemm('N','T',nh1,nb2,nh2,al_off,R%tblck(q)%tgam(5)%X,nh1,&
                L%mat(q2)%gam(6)%X,nb2,bet,w2%tblck(q)%tgam(9)%X,nh1)
        end if 
        
        if (nh1 .ne. 0) then 
           bet_off = 1 
           !w2hhph = Ahhhh.Bhhph - Bhhhh.Ahhph 
           call dgemm('N','N',nh1,nb2,nh1,al,L%mat(q1)%gam(5)%X,nh1,&
                R%tblck(q)%tgam(9)%X,nh1,bet_off,w2%tblck(q)%tgam(9)%X,nh1)
        end if 
        
     end if 
       
     RES%tblck(q)%tgam(9)%X = w1%tblck(q)%tgam(9)%X - w2%tblck(q)%tgam(9)%X
  
  end do
  
end subroutine 
!=================================================================
!=================================================================
 subroutine TS_commutator_222_ph(LCC,RCC,RES,WCC,jbas) 
   ! VERIFIED ph channel 2body TS_commutator. DFWT! 
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
         
         do JX =min(jxstart,IX),n2
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
                  rjl = TS_rval(j,l,Ntot,qx,LCC)
                  rik = TS_rval(i,k,Ntot,qx,LCC)
                  
                  sm = sm - (1.d0*WCC%CCX(qx)%X(rjl,rik) - &
                       WCC%CCR(qx)%X(rjl,rik) - &
                       1.d0*WCC%CCR(qx)%X(rik,rjl) + &
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
                
                  ril = TS_rval(i,l,Ntot,qx,LCC)
                  rjk = TS_rval(j,k,Ntot,qx,LCC)
                  
                  sm = sm + ( WCC%CCR(qx)%X(ril,rjk) - &
                       1.d0*WCC%CCX(qx)%X(ril,rjk) - &
                       WCC%CCX(qx)%X(rjk,ril) + &
                       1.d0*WCC%CCR(qx)%X(rjk,ril) ) * &
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
integer function TS_rval(i,l,Ntot,q,LCC) 
  implicit none 
  
  type(cross_coupled_31_mat) :: LCC
  integer :: i,l,Ntot,x,g,q
  
  x = CCindex(i,l,Ntot)
  g = 1
  do while (LCC%qmap(x)%Z(g) .ne. q )
  
     g = g + 1
  end do
  
  TS_rval = LCC%rmap(x)%Z(g)
end function 
!============================================
!============================================
integer function ph_rval(i,l,Ntot,q,LCC) 
  implicit none 
  
  type(cross_coupled_31_mat) :: LCC
  integer :: i,l,Ntot,x,g,q
  
  x = CCindex(i,l,Ntot)
  g = 1
  do while (LCC%qmap(x)%Z(g) .ne. q )
  
     g = g + 1
  end do
  
  ph_rval = LCC%nbmap(x)%Z(g)
end function 
!=====================================================
!=====================================================      
end module 
  
  
  
