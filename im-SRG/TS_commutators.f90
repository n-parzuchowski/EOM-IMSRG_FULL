module TS_commutators
  use cross_coupled
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
        if ( mod(lq,2) .ne. mod(lp+R%dpar/2,2)) cycle
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
                         *sqrt((J1 + 1.d0)*(J2 + 1.d0))*(-1)**(J1/2)*&
                         xxxsixj(J1,J2,rank,jq,jp,ji)             
                 
                 end do
              end do
              
              do J1 = abs(ji - jp),ji+jp,2
                 do J2 = abs(ja - jq),ja+jq,2 
                 
                    smy = smy + tensor_elem(ik,pk,ak,qk,J1,J2,R,jbas)&
                         *sqrt((J1 + 1.d0)*(J2 + 1.d0))*(-1)**(J1/2)*&
                         xxxsixj(J1,J2,rank,jq,jp,ji) 
                                     
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
        if ( mod(lq,2) .ne. mod(lp+R%dpar/2,2)) cycle
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
                         * xxxsixj(J1,J2,rank,jq,jp,ji)              
                 
                 end do
              end do
              
              do J1 = abs(ji - jp),ji+jp,2
                 do J2 = abs(ja - jq),ja+jq,2 
                 
                    smy = smy + tensor_elem(ik,pk,ak,qk,J1,J2,R,jbas)&
                         *sqrt((J1 + 1.d0)*(J2 + 1.d0))*(-1)**(J1/2) &
                         * xxxsixj(J1,J2,rank,jq,jp,ji)
                                     
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
        
        if ( mod(lq,2) .ne. mod(lp+R%dpar/2,2)) cycle
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
                          * xxxsixj(J1,J2,rank,jq,jp,ji)
                    
                 end do
              end do
              
              do J1 = abs(ji - jp),ji+jp,2
                 do J2 = abs(ja - jq),ja+jq,2 
                 
                    smy = smy + tensor_elem(ik,pk,ak,qk,J1,J2,R,jbas)&
                         *sqrt((J1 + 1.d0)*(J2 + 1.d0))*(-1)**(J1/2) & 
                                      * xxxsixj(J1,J2,rank,jq,jp,ji)
                    
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
subroutine  TS_commutator_211(LCC,R,RES,jbas) 
  ! onebody part of  - [R1,L2] 
  ! this one is brute force. 
  ! not sure of a faster way to do this
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: R,RES
  type(cc_mat) :: LCC 
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
        if ( mod(lq,2) .ne. mod(lp+R%dpar/2,2)) cycle
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
              if ( mod(li,2) .ne. mod(la+R%dpar/2,2)) cycle
              if (ti .ne. ta) cycle 
              if (.not. (triangle(ja,ji,rank))) cycle
               
              smx = 0.d0 
              smy = 0.d0
              smx2 = 0.d0 
              smy2 = 0.d0 
              
            
              Tz = abs(ta -ti)/2 
              if (abs(tp - tq) .ne. Tz*2)  cycle 
              PAR  = mod(la+li,2) 
              if (mod(lp+lq,2) .ne. PAR) cycle                                 
           
              
              qx = rank/2+1 + Tz*(JTM+1) + 2*PAR*(JTM+1)
                
              rai = ph_rval(ak,ik,Ntot,qx,LCC)
              rpq = fetch_rval(pk,qk,Ntot,qx,LCC)
              rqp = fetch_rval(qk,pk,Ntot,qx,LCC)

              sm = sm +  R%fph(a,i) & 
                *  ( (-1)**(( jp - jq + rank)/2) * LCC%CCX(qx)%X(rpq,rai) &
                -R%herm*LCC%herm*(-1)**(rank/2)*LCC%CCX(qx)%X(rqp,rai) ) &
                / sqrt(rank + 1.d0 )

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
        if ( mod(lq,2) .ne. mod(lp+R%dpar/2,2)) cycle
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
              if ( mod(li,2) .ne. mod(la+R%dpar/2,2)) cycle
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
              rpq = fetch_rval(pk,qk,Ntot,qx,LCC)
              rqp = fetch_rval(qk,pk,Ntot,qx,LCC)

              sm = sm +  R%fph(a,i) & 
                *  ( (-1)**(( jp - jq + rank)/2) * LCC%CCX(qx)%X(rpq,rai) &
                -R%herm*LCC%herm*(-1)**(rank/2)*LCC%CCX(qx)%X(rqp,rai) ) &
                / sqrt(rank + 1.d0 )

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
        if ( mod(lq,2) .ne. mod(lp+R%dpar/2,2)) cycle
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
              if ( mod(li,2) .ne. mod(la+R%dpar/2,2)) cycle
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
              rpq = fetch_rval(pk,qk,Ntot,qx,LCC)
              rqp = fetch_rval(qk,pk,Ntot,qx,LCC)

              sm = sm +  R%fph(a,i) & 
                *  ( (-1)**(( jp - jq + rank)/2) * LCC%CCX(qx)%X(rpq,rai) &
                -R%herm*LCC%herm*(-1)**(rank/2)*LCC%CCX(qx)%X(rqp,rai) ) &
                / sqrt(rank + 1.d0 )

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
!==================================================
!==================================================             
subroutine TS_commutator_212(L,R,RES,jbas) 
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: L,R,RES
  integer :: q,IX,JX,nh,np,nb,i,J1,J2
  integer :: a,b,c,d,ja,jb,jc,jd,ji,g_ix,q_sp,i_sp
  integer :: ta,tb,tc,td,ti,modla,modlb,modlc,modld,modli,spec,rank
  integer :: jxstart,jxend,ixend,c1,c2,n1,n2
  logical :: square
  real(8) ::  sm,sm1,sm2,sm3,sm4,d6ji
  
  rank = R%rank
  
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
        modla = mod(jbas%ll(a),2)
        ta = jbas%itzp(a) 
             
        b = R%tblck(q)%tensor_qn(c1,1)%Y(IX,2)
        jb = jbas%jj(b)
        modlb = mod(jbas%ll(b),2)
        tb = jbas%itzp(b)
 
        do JX = 1,n2
           
           c = R%tblck(q)%tensor_qn(c2,2)%Y(JX,1)
           jc = jbas%jj(c)
           modlc = mod(jbas%ll(c),2)
           tc = jbas%itzp(c)

           d = R%tblck(q)%tensor_qn(c2,2)%Y(JX,2)
           jd = jbas%jj(d)
           modld = mod(jbas%ll(d),2)
           td = jbas%itzp(d)
                   
           sm = 0.d0 

           do i = 1,jbas%total_orbits
               
              ji = jbas%jj(i) 
              ti = jbas%itzp(i)
              modli = mod(jbas%ll(i)+R%dpar/2,2) 
              
              sm1 = 0.d0 
              if (ti == ta) then
                 if (modli == modla ) then 
                    if (triangle(ji,ja,rank)) then  
                       
                       sm1 = sm1 - xxxsixj(J1,J2,rank,ji,ja,jb)&
                            *f_tensor_elem(a,i,R,jbas)*v_elem(i,b,c,d,J2,L,jbas)
                     
                    end if
                 end if
              end if
               
              sm1 = sm1*(-1)**((ja + jb-J2)/2) 
               
              
              sm2 = 0.d0 
              if (ti == tb) then
                 if (modli == modlb ) then 
                    if (triangle(ji,jb,rank)) then  
                       
                       sm2 = sm2 + xxxsixj(J1,J2,rank,ji,jb,ja)&
                            *f_tensor_elem(b,i,R,jbas)*v_elem(i,a,c,d,J2,L,jbas)
                       
                    end if
                 end if
              end if
               
              sm2 = sm2*(-1)**((J1+J2)/2) 
              

              sm3 = 0.d0 
              if (ti == td) then
                 if (modli == modld ) then 
                    if (triangle(ji,jd,rank)) then  
                       
                       sm3 = sm3 + xxxsixj(J1,J2,rank,jd,ji,jc)&
                            *f_tensor_elem(i,d,R,jbas)*v_elem(a,b,c,i,J1,L,jbas)
                       
                    end if
                 end if
              end if
               
              sm3 = sm3*(-1)**((jc+jd-J1)/2) 
              

              sm4 = 0.d0 
              if (ti == tc) then
                 if (modli == modlc ) then 
                    if (triangle(ji,jc,rank)) then  
                       
                       sm4 = sm4 -  xxxsixj(J1,J2,rank,jc,ji,jd)&
                            *f_tensor_elem(i,c,R,jbas)*v_elem(a,b,d,i,J1,L,jbas)
                       
                    end if
                 end if
              end if
               
              sm4 = sm4*(-1)**((J2+J1)/2)
              
              sm =  sm + (sm1+sm2+sm3+sm4) 
           end do 
           
           sm = sm * sqrt((J1+1.d0)*(J2+1.d0) / &
              (1.d0 + kron_del(a,b)) /(1.d0 + kron_del(c,d))) * (-1)**(rank/2) 
 
           RES%tblck(q)%tgam(g_ix)%X(IX,JX) = RES%tblck(q)%tgam(g_ix)%X(IX,JX)  +sm 
           
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
  type(sq_op) :: RES,w1,w2
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
        if (mod(li,2) .ne. mod(lj+w1%dpar/2,2))  cycle
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
                sm1 = sm1 - sqrt((J1+1.d0)*(J2+1.d0))*xxxsixj(J1,J2,w1%rank,jj,ji,jc) &
                 *tensor_elem(c,j,c,i,J2,J1,w1,jbas)*(-1)**(J2/2) * pm 
                
             end do              
          
             do J2 = max(J1,abs(jc - jj)),jc+jj,2
             
                ! use w1, because it sums over the pp indices
                sm1 = sm1 + sqrt((J1+1.d0)*(J2+1.d0))*xxxsixj(J1,J2,w1%rank,jj,ji,jc) &
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
                sm2 = sm2 - sqrt((J1+1.d0)*(J2+1.d0))*xxxsixj(J1,J2,w1%rank,jj,ji,jc) &
                 *tensor_elem(c,j,c,i,J2,J1,w2,jbas)*(-1)**(J2/2) * pm 
                
             end do              
          
             do J2 = max(J1,abs(jc - jj)),jc+jj,2
           
                ! use w1, because it sums over the pp indices
                sm2 = sm2 + sqrt((J1+1.d0)*(J2+1.d0))*xxxsixj(J1,J2,w1%rank,jj,ji,jc) &
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
        if (mod(li,2) .ne. mod(lj+w1%dpar/2,2))  cycle
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
                sm1 = sm1 - sqrt((J1+1.d0)*(J2+1.d0))*xxxsixj(J1,J2,w1%rank,jj,ji,jc) &
                 *tensor_elem(c,j,c,i,J2,J1,w1,jbas)*(-1)**(J2/2) * pm 
                
             end do              
          
             do J2 = max(J1,abs(jc - jj)),jc+jj,2
             
                ! use w1, because it sums over the pp indices
                sm1 = sm1 + sqrt((J1+1.d0)*(J2+1.d0))*xxxsixj(J1,J2,w1%rank,jj,ji,jc) &
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
                sm2 = sm2 - sqrt((J1+1.d0)*(J2+1.d0))*xxxsixj(J1,J2,w1%rank,jj,ji,jc) &
                 *tensor_elem(c,j,c,i,J2,J1,w2,jbas)*(-1)**(J2/2) * pm 
                
             end do              
          
             do J2 = max(J1,abs(jc - jj)),jc+jj,2
           
                ! use w1, because it sums over the pp indices
                sm2 = sm2 + sqrt((J1+1.d0)*(J2+1.d0))*xxxsixj(J1,J2,w1%rank,jj,ji,jc) &
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
        if (mod(li,2) .ne. mod(lj+w1%dpar/2,2))  cycle
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
                sm1 = sm1 - sqrt((J1+1.d0)*(J2+1.d0))*xxxsixj(J1,J2,w1%rank,jj,ji,jc) &
                 *tensor_elem(c,j,c,i,J2,J1,w1,jbas)*(-1)**(J2/2) * pm 
                
             end do              
          
             do J2 = max(J1,abs(jc - jj)),jc+jj,2
             
                ! use w1, because it sums over the pp indices
                sm1 = sm1 + sqrt((J1+1.d0)*(J2+1.d0))*xxxsixj(J1,J2,w1%rank,jj,ji,jc) &
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
                sm2 = sm2 - sqrt((J1+1.d0)*(J2+1.d0))*xxxsixj(J1,J2,w1%rank,jj,ji,jc) &
                 *tensor_elem(c,j,c,i,J2,J1,w2,jbas)*(-1)**(J2/2) * pm 
                
             end do              
          
             do J2 = max(J1,abs(jc - jj)),jc+jj,2
           
                ! use w1, because it sums over the pp indices
                sm2 = sm2 + sqrt((J1+1.d0)*(J2+1.d0))*xxxsixj(J1,J2,w1%rank,jj,ji,jc) &
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
  integer :: q,q1,q2,J1,J2,Tz,Par,phase,rank,i
  integer :: np1,nb1,nh1,np2,nb2,nh2,pm
  real(8) :: bet_off,al_off
  
  pm = R%herm*L%herm
  rank = R%rank
!construct temporary matrices
  do q = 1, R%nblocks
     
     J1 = R%tblck(q)%Jpair(1) 
     J2 = R%tblck(q)%Jpair(2)
     phase = R%tblck(q)%lam(1)
     par = R%tblck(q)%lam(2) 
     Tz = R%tblck(q)%lam(3)
    
     q1 = block_index(J1,Tz,Par) 
     q2 = block_index(J2,Tz,mod(Par+R%Dpar/2,2)) 
     
     nh1 = R%tblck(q)%nhh1
     np1 = R%tblck(q)%npp1
     nb1 = R%tblck(q)%nph1
     nh2 = R%tblck(q)%nhh2
     np2 = R%tblck(q)%npp2
     nb2 = R%tblck(q)%nph2
       
     do i = 1, 9 
        w1%tblck(q)%tgam(i)%X=0.d0
        w2%tblck(q)%tgam(i)%X=0.d0
     end do 
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
     
     RES%tblck(q)%tgam(1)%X = RES%tblck(q)%tgam(1)%X + &
          w1%tblck(q)%tgam(1)%X - w2%tblck(q)%tgam(1)%X


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

        RES%tblck(q)%tgam(4)%X = RES%tblck(q)%tgam(4)%X + &
             w1%tblck(q)%tgam(4)%X - w2%tblck(q)%tgam(4)%X
         
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
        
     RES%tblck(q)%tgam(5)%X = RES%tblck(q)%tgam(5)%X + &
          w1%tblck(q)%tgam(5)%X - w2%tblck(q)%tgam(5)%X
     

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
       
     RES%tblck(q)%tgam(3)%X = RES%tblck(q)%tgam(3)%X + &
          w1%tblck(q)%tgam(3)%X - w2%tblck(q)%tgam(3)%X

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
       
     RES%tblck(q)%tgam(7)%X = RES%tblck(q)%tgam(7)%X + &
          w1%tblck(q)%tgam(7)%X - w2%tblck(q)%tgam(7)%X

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
       
     RES%tblck(q)%tgam(2)%X = RES%tblck(q)%tgam(2)%X + &
          w1%tblck(q)%tgam(2)%X - w2%tblck(q)%tgam(2)%X

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
       
     RES%tblck(q)%tgam(8)%X = RES%tblck(q)%tgam(8)%X +&
          w1%tblck(q)%tgam(8)%X - w2%tblck(q)%tgam(8)%X


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
       
     RES%tblck(q)%tgam(6)%X = RES%tblck(q)%tgam(6)%X + &
          w1%tblck(q)%tgam(6)%X - w2%tblck(q)%tgam(6)%X

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
       
     RES%tblck(q)%tgam(9)%X = RES%tblck(q)%tgam(9)%X + &
          w1%tblck(q)%tgam(9)%X - w2%tblck(q)%tgam(9)%X
  
  end do
  
end subroutine 
!=================================================================
!=================================================================
 subroutine TS_commutator_222_ph(LCC,RCC,RES,WCC,jbas) 
   ! VERIFIED ph channel 2body TS_commutator. DFWT! 
   implicit none 
  
   type(spd) :: jbas
   type(sq_op) :: RES
   type(pandya_mat) :: RCC,WCC
   type(cc_mat) :: LCC
   integer :: nh,np,nb1,nb2,q,IX,JX,i,j,k,l,r1,r2,Tz,PAR,JTM,q1,q2,J3,J4,rank
   integer :: ji,jj,jk,jl,ti,tj,tk,tl,li,lj,lk,ll,n1,n2,c1,c2,jxstart,J4min,J4max
   integer :: J1,J2, Jtot,Ntot,qx,J3min,J3max,ril,rjk,rli,rkj,g_ix,thread,total_threads
   integer :: phase1,phase2,phase3,rik,rki,rjl,rlj,PAR2,ass
   real(8) :: sm ,pre,pre2,omp_get_wtime ,t1,t2,coef9,factor,sm_ex
   logical :: square
   
  rank = RES%rank
  Ntot = RES%Nsp
  JTM = jbas%Jtotal_max
  total_threads = size(RES%direct_omp) - 1
   ! construct intermediate matrices

  do q = 1,RCC%nblocks
     if (RCC%jval2(q) > jbas%jtotal_max*2) cycle

     nb2 = size(RCC%CCX(q)%X(1,:))
     nb1 = size(RCC%CCR(q)%X(:,1))
     r1 = size(RCC%CCX(q)%X(:,1))
     r2 = size(RCC%CCR(q)%X(1,:))      

     if (r1 * r2 == 0) cycle
      
     PAR = mod(q-1,2)
     Tz = mod((q-1)/2,2) 
         
     if (nb1 .ne. 0 ) then 
        J1 = RCC%Jval(q) 
        factor = 1.0/sqrt(J1+1.d0)*LCC%herm
        q1 = J1/2+1 + Tz*(JTM+1) + 2*PAR*(JTM+1)
        
         call dgemm('N','N',r1,r2,nb1,factor,LCC%CCX(q1)%X,r1,&
              RCC%CCR(q)%X,nb1,bet,WCC%CCX(q)%X,r1)       
     end if
         
     if (nb2 .ne. 0 ) then 
        
        J2 = RCC%Jval2(q) 
        PAR2 = mod(PAR+RCC%dpar/2,2) 
        q2 = J2/2+1 + Tz*(JTM+1) + 2*PAR2*(JTM+1)
        factor = 1.d0/sqrt(J2+1.d0)
       
        call dgemm('N','T',r1,r2,nb2,factor,RCC%CCX(q)%X,r1,&
             LCC%CCX(q2)%X,r2,bet,WCC%CCR(q)%X,r1) 
     end if
     
  end do


!$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE), SHARED(RES,WCC)  
  do thread = 1, total_threads
     do q = 1+RES%direct_omp(thread),RES%direct_omp(thread+1) 

!   do q = 1,RES%nblocks
      J1 = RES%tblck(q)%jpair(1)
      J2 = RES%tblck(q)%jpair(2)
               
      do g_ix = 1,9 
   
         ! figure out how big the array is
         n1 = size(RES%tblck(q)%tgam(g_ix)%X(:,1))
         n2 = size(RES%tblck(q)%tgam(g_ix)%X(1,:))
         if ((n1*n2) == 0) cycle 
        
         ! read in information about which 
         ! array we are using from public arrays
        c1 = sea1(g_ix) 
        c2 = sea2(g_ix) 
        square = sqs(g_ix) 
        jxstart = jst(g_ix) 
        
      do  IX =  1, n1 
         pre = 1.d0 

         i = RES%tblck(q)%tensor_qn(c1,1)%Y(IX,1)
         j = RES%tblck(q)%tensor_qn(c1,1)%Y(IX,2)
 
         if (i == j )  pre  = .70710678118d0
         ji = jbas%jj(i) 
         jj = jbas%jj(j) 
         li = jbas%ll(i) 
         lj = jbas%ll(j)
         ti = jbas%itzp(i) 
         tj = jbas%itzp(j)
         
         do JX =1,n2
            pre2 = 1.d0 
            k = RES%tblck(q)%tensor_qn(c2,2)%Y(JX,1)
            l = RES%tblck(q)%tensor_qn(c2,2)%Y(JX,2)
            
            if (k == l )  pre2 = .70710678118d0
            jk = jbas%jj(k) 
            jl = jbas%jj(l) 
            lk = jbas%ll(k) 
            ll = jbas%ll(l)
            tk = jbas%itzp(k) 
            tl = jbas%itzp(l)
           
            phase1 = (-1) ** (( ji + jj + jk + jl )/2) 
            
            sm = 0.d0 
            sm_ex = 0.d0 
                       
            J3min = abs(ji - jl) 
            J3max = ji + jl
            
            J4min = abs(jj - jk)
            J4max = jj + jk 
            
            
            Tz = abs(ti -tl)/2                         
            PAR = mod(li+ll,2) 

            if (mod(lk+lj+RCC%dpar/2,2) == PAR) then 
               if (abs(tk - tj) == Tz*2)  then 
             
            do J3 = J3min,J3max,2
            
               q1 = block_index(J3,Tz,PAR)
              
               ril = fetch_rval(i,l,Ntot,q1,RCC)
               rli = fetch_rval(l,i,Ntot,q1,RCC)

                do J4 = max( J3 , J4min ) , J4max,2 
            
                  if (.not. (triangle(J3,J4,rank))) cycle
                  
                  PAR2 = mod(PAR + RCC%dpar/2,2) 
                  q2 = block_index(J4,Tz,PAR2)
               
                  rjk = fetch_rval(j,k,Ntot,q2,RCC)
                  rkj = fetch_rval(k,j,Ntot,q2,RCC)
              
                  qx = CCtensor_block_index(J3,J4,rank,Tz,PAR)
                  sm = sm + sqrt((J3+1.d0)*(J4+1.d0))* &
                       ninej(ji,jl,J3,jj,jk,J4,J1,J2,rank)  * ( &
                       WCC%CCX(qx)%X(rli,rkj)*(-1)**((J3+J4)/2)*(-1)**((jl+ji)/2)  &
                      + WCC%CCR(qx)%X(ril,rkj) * phase1 * LCC%herm &
                      - (-1)**(rank/2)*phase1*RCC%herm*LCC%herm * &
                      WCC%CCX(qx)%X(ril,rjk)*(-1)**((jl+ji)/2) - (-1)**(( J3+J4+rank)/2) &
                      * RCC%herm * WCC%CCR(qx)%X(rli,rjk) )

                end do 
               
                do J4 = J4min , min(J4max,J3-2),2 
                  if (.not. (triangle(J3,J4,rank))) cycle
                  
                  
                  PAR2 = mod(PAR + RCC%dpar/2,2)                  
                  q2 = block_index(J4,Tz,PAR2)
                  
                  rjk = fetch_rval(j,k,Ntot,q2,RCC)     
                  rkj = fetch_rval(k,j,Ntot,q2,RCC)
                  
                  qx = CCtensor_block_index(J4,J3,rank,Tz,PAR2)
                  sm = sm + sqrt((J3+1.d0)*(J4+1.d0))* &
                       ninej(ji,jl,J3,jj,jk,J4,J1,J2,rank)  * ( &
                       WCC%CCR(qx)%X(rjk,rli)*phase1*(-1)**((rank+J3+J4)/2)*LCC%herm &
                       + WCC%CCX(qx)%X(rkj,rli)*(-1)**((jk+jj)/2) * (-1)**(rank/2) &
                        - RCC%herm * WCC%CCR(qx)%X(rkj,ril) &
                     - phase1*(-1)**((J3+J4)/2) * RCC%herm *LCC%herm * &
                          WCC%CCX(qx)%X(rjk,ril)*(-1)**((jk+jj)/2) )
                      
               end do
               
            
            end do 
            
               end if 
            end if 

            ! exchange of 1 set of indeces
            J3min = abs(jj - jl) 
            J3max = jj + jl
            
            J4min = abs(ji - jk)
            J4max = ji + jk 
            
            Tz = abs(tl -tj)/2 
            PAR = mod(ll+lj,2) 

            if (mod(li+lk+RCC%dpar/2,2) == PAR) then 
               if (abs(ti - tk) == Tz*2)  then 

            do J3 = J3min,J3max,2
               q1 = block_index(J3,Tz,PAR)
            
               rjl = fetch_rval(j,l,Ntot,q1,RCC)
               rlj = fetch_rval(l,j,Ntot,q1,RCC)

               do J4 = max( J3 , J4min ) , J4max,2 
                 
                  if (.not. (triangle(J3,J4,rank))) cycle
                 
                  
                  PAR2 = mod(PAR + RCC%dpar/2,2)
                  q2 = block_index(J4,Tz,PAR2)
             
                  rki = fetch_rval(k,i,Ntot,q2,RCC)
                  rik = fetch_rval(i,k,Ntot,q2,RCC)
              
                  qx = CCtensor_block_index(J3,J4,rank,Tz,PAR)
               
                  sm_ex = sm_ex - sqrt((J3+1.d0)*(J4+1.d0))* &
                       ninej(jj,jl,J3,ji,jk,J4,J1,J2,rank) * ( &
                       (-1)**((J3+J4)/2) * WCC%CCX(qx)%X(rlj,rki)*(-1)**((jl+jj)/2)  & 
                       - phase1*RCC%herm*LCC%herm*(-1)**(rank/2)* &
                       WCC%CCX(qx)%X(rjl,rik)*(-1)**((jl+jj)/2) + phase1*LCC%herm * &
                       WCC%CCR(qx)%X(rjl,rki) - (-1)**((J3+J4+rank)/2)* &
                       RCC%herm * WCC%CCR(qx)%X(rlj,rik) )
                  
                end do 
               
                do J4 = J4min , min(J4max,J3-2),2 
                  if (.not. (triangle(J3,J4,rank))) cycle
                  
                  PAR2 = mod(PAR + RCC%dpar/2,2)
                  q2 = block_index(J4,Tz,PAR2)
     
                  rki = fetch_rval(k,i,Ntot,q2,RCC)
                  rik = fetch_rval(i,k,Ntot,q2,RCC)

                  qx = CCtensor_block_index(J4,J3,rank,Tz,PAR2)
                  
                  sm_ex = sm_ex - sqrt((J3+1.d0)*(J4+1.d0))* &
                       ninej(jj,jl,J3,ji,jk,J4,J1,J2,rank) * ( &
                       (-1)**((J3+J4+rank)/2) *phase1*LCC%herm &
                       * WCC%CCR(qx)%X(rik,rlj)  & 
                       - RCC%herm*WCC%CCR(qx)%X(rki,rjl) &
                       + (-1)**(rank/2) * WCC%CCX(qx)%X(rki,rlj)*(-1)**((jk+ji)/2) &
                       - phase1*(-1)**((J3+J4)/2) * LCC%herm * RCC%herm * &
                       WCC%CCX(qx)%X(rik,rjl)*(-1)**((jk+ji)/2) ) 
                  
               end do
               
               
            end do 
               end if 
            end if 

          RES%tblck(q)%tgam(g_ix)%X(IX,JX) = RES%tblck(q)%tgam(g_ix)%X(IX,JX) + ( sm * &
               (-1) ** ((ji+jj+J2)/2) + sm_ex * (-1)**((J1+J2)/2) )&
               *   pre * pre2 *   sqrt((J1+1.d0)*(J2+1.d0))
   
         end do 
      end do
      end do 
   end do

   end do 
!$OMP END PARALLEL DO 
   
end subroutine 
!=====================================================
!=====================================================      
real(8) function TS_commutator_223_single(L,R,ip,iq,ir,is,it,iu,jtot1,jtot2,Jpq,Jst,jbas)
  implicit none 
  
  integer,intent(in) :: ip,iq,ir,is,it,iu,jpq,jst
  integer :: a,b,c,d,Jx,Jy,Jz,J2,J1,J3,phase,rank
  integer :: ja,jb,jc,jd,jp,jq,jr,js,jt,ju,jtot1,jtot2
  integer :: j1min,j1max , j2min,j2max,j3min,j3max
  type(sq_op) :: L,R
  type(spd) :: jbas
  real(8) :: sm,sm_sub,multfact,smtot,d6ji,out,otherfact,dsum
  real(8) :: Vs1,Vs2,sj1,sj2,sj3,sj4
  
  smtot = 0.d0 
  rank = R%rank 
  
  jp = jbas%jj(ip)
  jq = jbas%jj(iq)
  jr = jbas%jj(ir)  
  js = jbas%jj(is)
  jt = jbas%jj(it)
  ju = jbas%jj(iu)  
     
  ! FIRST TERM 
  !changed to q-r instead of q+r
  multfact = (-1)**((jq+jr+jp-jtot2+rank)/2) *sqrt((Jpq+1.d0)*(Jst+1.d0) &
       *(jtot1+1.d0)*(jtot2+1.d0) ) * (rank+1.d0) 
  ! ju isn't in here because I need it to make add with ja
  ! so I get an integer later 

  
  ! Left-Right 
  sm = 0.d0   
  
  j2min = abs(jq - jr)
  j2max = jq+jr 
     
 
  do a = 1,jbas%total_orbits 
     
     ja = jbas%jj(a) 
     
     if (.not. triangle(jp,ja,jst) ) cycle
          
     j3min = abs(ja - ju) 
     j3max = ja+ju 
          
     phase = (-1)**((ja-ju)/2) 
     

     sj1 = v_elem(ip,a,is,it,Jst,L,jbas)*phase

     do J2 = j2min, j2max , 2
        sj2 = sj1*d6ji(jp,jq,Jpq,jr,jtot1,J2)*sqrt(J2+1.d0) 
        do J3 = j3min,j3max,2

           sm = sm - (-1)**(J3/2) * sqrt(J3 + 1.d0) &            
            * sj2 * d6ji(jp,ja,Jst,ju,jtot2,J3) *  &
             d6ji(J2,J3,rank,jtot2,jtot1,jp) * &
             tensor_elem(iq,ir,a,iu,J2,J3,R,jbas)
           
        end do
     end do
  end do
 
  smtot = smtot + sm*multfact
  ! !Right-left  

  do a = 1,jbas%total_orbits 
     
     ja = jbas%jj(a) 
          
     j1min = abs(jp - ja)
     j1max = jp + ja  
     
     j3min = max( abs(jq - jr) , abs(ja - ju) ) 
     j3max = min( jq+jr , ja+ju)
     
     phase = (-1)**((ja - jp)/2)
             
     do J1 = j1min, j1max , 2      

        Vs1 = tensor_elem(ip,a,is,it,J1,Jst,R,jbas)*(-1)**(J1/2) 
        sj1 = sqrt(J1+1.d0)*xxxsixj(J1,Jst,rank,jtot2,jtot1,ju) 

        do J3 = j3min,j3max,2

           sm = sm +  phase *sj1*(J3 + 1.d0) &
                * sixj(jp,jq,Jpq,jr,jtot1,J3) * sixj(jp,ja,J1,ju,jtot1,J3) &
                * Vs1 * v_elem(iq,ir,a,iu,J3,L,jbas)

        end do
     end do
  end do 
  
  
  smtot = smtot + sm*multfact*0 
  
!   !SECOND TERM
!   !Left-Right
!   sm = 0.d0 
  
!   multfact = (-1)**((jp+jq+jr+jt+ju+jtot2+rank)/2) *sqrt((Jpq+1.d0)*(Jst+1.d0)*&
!        (jtot1+1.d0)*(jtot2+1.d0))*(rank+1.d0) 
!   ! added a minus sign
!   do a = 1, jbas%total_orbits

!      ja = jbas%jj(a)

!      j2min = max( abs(jp - ja) , abs(jt - ju) ) 
!      j2max = min( jp+ja , jt+ju) 
     
!      j1min = abs(ja - js) 
!      j1max = ja+js
     
!      j3min = abs(jq - jr) 
!      j3max = jq+jr
     
!      phase = (-1) ** ((ja - js)/2) 
     
!      do J1 = j1min,j1max,2
        
!         sj1 = sqrt(J1+1.d0)*phase  
       
!         do J2 = j2min,j2max,2
           
!            sj2 = sj1 * (J2+1.d0) * (-1)**(J2/2) * sixj(js,jt,Jst,ju,jtot2,J2) * &
!                 sixj(jp,ja,J2,js,jtot2,J1) * v_elem(ip,a,it,iu,J2,L,jbas) 
           
!            do J3 = j3min,j3max,2
           
!               sm = sm - sqrt(J3+1.d0) * sj2 * sixj(jp,jq,Jpq,jr,jtot1,J3) &
!                    * xxxsixj(J1,J3,rank,jtot1,jtot2,jp) * (-1)**(J3/2) * &
!                      tensor_elem(iq,ir,a,is,J3,J1,R,jbas)
!            end do
!         end do         
        
!      end do
 
!   end do


!   do a = 1, jbas%total_orbits

!      ja = jbas%jj(a)
     
!      j1min = abs(jp - ja)
!      j1max = jp+ja 

!      j2min = abs(jt - ju) 
!      j2max = jt+ju 
     
!      j3min = max( abs(jq - jr) , abs(ja - js) ) 
!      j3max = min( jq+jr , ja+js)
     
!      phase = (-1) ** ((ja - jp)/2) 
     
!      do J1 = j1min,j1max,2
        
!         sj1 = sqrt(J1+1.d0) * phase 
        
!         do J2 = j2min,j2max,2 
           
!            sj2 = sj1 * sqrt(J2+1.d0) * sixj(js,jt,Jst,ju,jtot2,J2) * &
!                 xxxsixj(J2,J1,rank,jtot1,jtot2,js) * tensor_elem(ip,a,it,iu,J1,J2,R,jbas) 
           
!            do J3 = j3min,j3max,2
              
!               sm = sm + sj2 * (J3+1.d0) * sixj(jp,jq,Jpq,jr,jtot1,J3) &
!                    * sixj(jp,ja,J1,js,jtot1,J3) * v_elem(iq,ir,a,is,J3,L,jbas)
!            end do 
!         end do 
!      end do
 
!   end do
     
!   smtot = smtot + sm*multfact

!   ! THIRD TERM    
!   !Left-Right
!   sm = 0.d0 
  
!   multfact = (-1)**((jq+jr+jtot2+js+Jst+rank)/2) *sqrt((Jpq+1.d0)*(Jst+1.d0)*&
!        (jtot1+1.d0)*(jtot2+1.d0)) *(rank+1.d0)

!   do a = 1, jbas%total_orbits

!      ja = jbas%jj(a)     
     
!      j2min = max( abs(jp - ja) , abs(js - ju) ) 
!      j2max = min( jp+ja , js+ju) 
     
!      j3min = abs(jq - jr) 
!      j3max = jq+jr 
   
!      j1min = abs(ja - jt) 
!      j1max = ja+jt
     
!      phase = (-1) ** ((ja + jp)/2) ! changed to ja+js rather than ja-js 
     
!      do J3 = j3min,j3max,2
        
!         sj1 = phase*sixj(jp,jq,Jpq,jr,jtot1,J3)*sqrt(J3+1.d0)  
        
!         do J1 = j1min,j1max,2 
           
!            sj2 =  sj1*(-1)**(J1/2)*tensor_elem(iq,ir,a,it,J3,J1,R,jbas)*&
!                 xxxsixj(J1,J3,rank,jtot1,jtot2,jp)*sqrt(J1+1.d0)

!            do J2 = j2min,j2max,2
!               sm = sm - (J2+1.d0) *sj2* sixj(js,jt,Jst,jtot2,ju,J2) &
!                    * sixj(jp,ja,J2,jt,jtot2,J1) * v_elem(ip,a,iu,is,J2,L,jbas)
!            end do
!         end do
!      end do
 
!   end do

!   ! right-left

!   do a = 1, jbas%total_orbits

!      ja = jbas%jj(a)     
    
!      j1min = abs(jp - ja)
!      j1max =  jp+ja 
     
!      j2min =  abs(js - ju) 
!      j2max =  js+ju 
     
!      j3min = max( abs(jq - jr) , abs(ja - jt) ) 
!      j3max = min( jq+jr , ja+jt)
     
!      phase = (-1) ** ((ja + jt)/2) ! changed to ja+js rather than ja-js 
     
!      do J2 = j2min,j2max,2
        
!         sj1 = phase*(-1)**(J2/2)* sixj(js,jt,Jst,jtot2,ju,J2)*sqrt(J2+1.d0) 
        
!         do J1 = j1min,j1max,2 
!            sj2 = sj1* xxxsixj(J1,J2,rank,jtot2,jtot1,jt) *sqrt(J1+1.d0) *&
!                 tensor_elem(ip,a,iu,is,J1,J2,R,jbas)      

!            do J3 = j3min,j3max,2
!               sm = sm + (J3+1.d0) *sj2* sixj(jp,jq,Jpq,jr,jtot1,J3) &
!                    * sixj(jp,ja,J1,jt,jtot1,J3) * v_elem(iq,ir,a,it,J3,L,jbas)
!            end do
        
!         end do
        
!      end do
 
!   end do
     
!   smtot = smtot + sm*multfact



!   ! FOURTH TERM    
!   !Left-Right
!   sm = 0.d0 
  
!   multfact = (-1)**((jp+jtot2+Jpq+rank)/2) *sqrt((Jpq+1.d0)*(Jst+1.d0)*&
!        (jtot1+1.d0)*(jtot2+1.d0)) *(rank+1.d0)

!   do a = 1, jbas%total_orbits

!      ja = jbas%jj(a)     
          
     
!      j1min = abs(jp - jr) 
!      j1max = jp+jr 
   
!      j2min = abs(ja - ju) 
!      j2max = ja+ju
     
!      phase = (-1) ** ((ja - ju)/2) ! changed to ja+js rather than ja-js 
     
!      sj1 = v_elem(iq,a,is,it,Jst,L,jbas)
     
!      do J1 = j1min,j1max,2
        
!         sj2 = phase*sj1*sixj(jp,jq,Jpq,jtot1,jr,J1)*sqrt(J1+1.d0)*(-1)**(J1/2)
        
!         do J2 = j2min,j2max,2 
           
!            sm = sm -  sj2*sqrt(J2+1.d0)*(-1)**(J2/2)*tensor_elem(ir,ip,a,iu,J1,J2,R,jbas)*&
!                 xxxsixj(J1,J2,rank,jtot2,jtot1,jq)*sixj(ja,jq,Jst,jtot2,ju,J2)

!         end do
!      end do
 
!   end do
!   smtot = smtot + sm*multfact
!   ! right-left
     
!   sm = 0.d0 
  
!   multfact = (-1)**((jp+jtot1+Jpq+Jst+rank)/2) *sqrt((Jpq+1.d0)*&
!        (jtot1+1.d0)*(jtot2+1.d0)) *(rank+1.d0)

!   do a = 1, jbas%total_orbits

!      ja = jbas%jj(a)     
          
!      j1min = abs(jq - ja) 
!      j1max = jq+ja 
   
!      j2min = max(abs(ja - ju),abs(jr-jp)) 
!      j2max = min(ja+ju,jr+jp) 
     
!      phase = (-1) ** ((ja - jq)/2) ! changed to ja+js rather than ja-js 
     
! !     sj1 = v_elem(iq,a,is,it,Jst,L,jbas)
     
!      do J1 = j1min,j1max,2
        
!         sj1 = phase*sqrt(J1+1.d0)* &
!              tensor_elem(iq,a,is,it,J1,Jst,R,jbas) 
        
!         do J2 = j2min,j2max,2 
           
!            sm = sm +  sj1*(J2+1.d0)*(-1)**(J2/2)*v_elem(ir,ip,a,iu,J2,L,jbas) *&
!                 xxxsixj(J1,Jst,rank,jtot2,jtot1,ju)*sixj(ja,jq,J1,jtot1,ju,J2) &
!                 *sixj(jp,jq,Jpq,jtot1,jr,J2)
!         end do
!      end do
 
!   end do
  
!   smtot = smtot + sm*multfact

!   ! FIFTH TERM    
!   !Left-Right
!   sm = 0.d0 
  
!   multfact = (-1)**((ju+jt+jtot2+js+Jpq+rank)/2) *sqrt((Jpq+1.d0)*(Jst+1.d0)*&
!        (jtot1+1.d0)*(jtot2+1.d0)) *(rank+1.d0)

!   do a = 1, jbas%total_orbits

!      ja = jbas%jj(a)     
     
!      j1min = max( abs(jq - ja) , abs(jt - ju) ) 
!      j1max = min( jq+ja , jt+ju) 
     
!      j2min = abs(jp - jr) 
!      j2max = jp+jr 
   
!      j3min = abs(ja - js) 
!      j3max = ja+js
     
!      phase = (-1) ** ((ja + jp)/2) ! changed to ja+js rather than ja-js 
     
!      do J1 = j1min,j1max,2
        
!         sj1 = phase*sixj(js,jt,Jst,ju,jtot2,J1)*(J1+1.d0)*(-1)**(J1/2)&
!              *v_elem(iq,a,it,iu,J1,L,jbas) 
        
!         do J3 = j3min,j3max,2 
           
!            sj2 =  sj1*(-1)**(J3/2)*&
!                 sixj(jq,ja,J1,js,jtot2,J3)*sqrt(J3+1.d0)

!            do J2 = j2min,j2max,2
!               sm = sm - sqrt(J2+1.d0)*(-1)**(J2/2)*sj2* sixj(jp,jq,Jpq,jtot1,jr,J2) &
!                    * xxxsixj(J3,J2,rank,jtot1,jtot2,jq) * tensor_elem(ir,ip,a,is,J2,J3,R,jbas)
!            end do
!         end do
!      end do
 
!   end do
!   smtot = smtot + sm*multfact

!   sm = 0.d0 
!   ! right-left
!   multfact = (-1)**((jp+jq+jtot2+ju+Jpq+rank)/2) *sqrt((Jpq+1.d0)*(Jst+1.d0)*&
!        (jtot1+1.d0)*(jtot2+1.d0)) *(rank+1.d0)

!   do a = 1, jbas%total_orbits

!      ja = jbas%jj(a)     
    
!      j1min = abs(jq - ja)
!      j1max =  jq+ja 
     
!      j2min =  abs(jt - ju) 
!      j2max =  jt+ju 
     
!      j3min = max( abs(jp - jr) , abs(ja - js) ) 
!      j3max = min( jp+jr , ja+js)
     
!      phase = (-1) ** ((ja + jt)/2) ! changed to ja+js rather than ja-js 
     
!      do J3 = j3min,j3max,2
        
!         sj1 = phase*(-1)**(J3/2)* sixj(jp,jq,Jpq,jtot1,jr,J3)*(J3+1.d0)*&
!              v_elem(ir,ip,a,is,J3,L,jbas) 
        
!         do J2 = j2min,j2max,2 
!            sj2 = sj1*(-1)**(J2/2)*sixj(js,jt,Jst,ju,jtot2,J2)*sqrt(J2+1.d0)
                    

!            do J1 = j1min,j1max,2
!               sm = sm + sqrt(J1+1.d0) *(-1)**(J1/2) *sj2* sixj(jq,ja,J1,js,jtot1,J3) &
!                    * xxxsixj(J2,J1,rank,jtot1,jtot2,js) * tensor_elem(iq,a,it,iu,J1,J2,R,jbas)
!            end do
        
!         end do
        
!      end do
 
!   end do
     
!   smtot = smtot + sm*multfact 

!   ! SIXTH TERM    
!   !Left-Right
!   sm = 0.d0 
  
!   multfact = (-1)**((jtot2+js+Jpq+Jst+rank)/2) *sqrt((Jpq+1.d0)*(Jst+1.d0)*&
!        (jtot1+1.d0)*(jtot2+1.d0)) *(rank+1.d0)

!   do a = 1, jbas%total_orbits

!      ja = jbas%jj(a)     
     
!      j1min = max( abs(jq - ja) , abs(js - ju) ) 
!      j1max = min( jq+ja , js+ju) 
     
!      j2min = abs(jp - jr) 
!      j2max = jp+jr 
   
!      j3min = abs(ja - jt) 
!      j3max = ja+jt
     
!      phase = (-1) ** ((ja - jp)/2) ! changed to ja+js rather than ja-js 
     
!      do J1 = j1min,j1max,2
        
!         sj1 = phase*sixj(js,jt,Jst,jtot2,ju,J1)*(J1+1.d0) &
!              *v_elem(iq,a,iu,is,J1,L,jbas) 
        
!         do J2 = j2min,j2max,2 
           
!            sj2 =  sj1*(-1)**(J2/2)*&
!                 sixj(jp,jq,Jpq,jtot1,jr,J2)*sqrt(J2+1.d0)

!            do J3 = j3min,j3max,2
!               sm = sm - sqrt(J3+1.d0)*(-1)**(J3/2)*sj2* sixj(jq,ja,J1,jt,jtot2,J3) &
!                    * xxxsixj(J3,J2,rank,jtot1,jtot2,jq) * tensor_elem(ir,ip,a,it,J2,J3,R,jbas)
!            end do
!         end do
!      end do
 
!   end do
!   smtot = smtot + sm*multfact

!   sm = 0.d0 
!   ! right-left
!   multfact = (-1)**((jp+jq+js+jtot1+Jpq+Jst+rank)/2) *sqrt((Jpq+1.d0)*(Jst+1.d0)*&
!        (jtot1+1.d0)*(jtot2+1.d0)) *(rank+1.d0)

!   do a = 1, jbas%total_orbits

!      ja = jbas%jj(a)     
    
!      j1min = abs(jq - ja)
!      j1max =  jq+ja 
     
!      j2min =  abs(js - ju) 
!      j2max =  js+ju 
     
!      j3min = max( abs(jp - jr) , abs(ja - jt) ) 
!      j3max = min( jp+jr , ja+jt)
     
!      phase = (-1) ** ((ja + jt)/2) ! changed to ja+js rather than ja-js 
     
!      do J3 = j3min,j3max,2
        
!         sj1 = phase*(-1)**(J3/2) * sixj(jp,jq,Jpq,jtot1,jr,J3)*(J3+1.d0)*&
!              v_elem(ir,ip,a,it,J3,L,jbas) 
        
!         do J2 = j2min,j2max,2 

!            sj2 = sj1*sixj(js,jt,Jst,jtot2,ju,J2)*sqrt(J2+1.d0)
                    

!            do J1 = j1min,j1max,2
!               sm = sm + sqrt(J1+1.d0) *(-1)**(J1/2) *sj2* sixj(jq,ja,J1,jt,jtot1,J3) &
!                    * xxxsixj(J2,J1,rank,jtot1,jtot2,jt) * tensor_elem(iq,a,iu,is,J1,J2,R,jbas)
!            end do
        
!         end do
        
!      end do
 
!   end do
     
!   smtot = smtot + sm*multfact


!   ! SEVENTH TERM
  
!   ! Left-right
!   sm = 0.d0
!   multfact = (-1)**((Jpq+rank+jr-jtot2)/2) *sqrt((Jst+1.d0)*(jtot1+1.d0)*&
!        (jtot2+1.d0))*(rank+1.d0) 
!   do a = 1, jbas%total_orbits
     
!      ja = jbas%jj(a)

!      if (.not. triangle(jr,ja,jst) ) cycle
     
!      j3min = abs(ja-ju) 
!      j3max = ja + ju 
     
!      sj1 = v_elem(ir,a,is,it,Jst,L,jbas) * (-1)**((ja+ju)/2) 
!      do J3=j3min,j3max,2 
!         sm = sm - sj1* sixj(jr,ja,Jst,ju,jtot2,J3) * (-1)**(J3/2) * sqrt(J3+1.d0) &
!              * xxxsixj(J3,Jpq,rank,jtot1,jtot2,jr) * tensor_elem(ip,iq,a,iu,Jpq,J3,R,jbas)
!      end do 
!   end do 
  
!   smtot = smtot + sm*multfact

!   ! right-left
!   sm = 0.d0
!   multfact = (-1)**((Jpq+Jst+rank)/2) *sqrt((Jpq+1.d0)*(jtot1+1.d0)*&
!        (jtot2+1.d0))*(rank+1.d0) 
!   do a = 1, jbas%total_orbits
     
!      ja = jbas%jj(a)

!      if (.not. triangle(ju,ja,Jpq) ) cycle
     
!      j3min = abs(ja-jr) 
!      j3max = ja + jr 
     
!      sj1 = v_elem(ip,iq,a,iu,Jpq,L,jbas) * (-1)**((ja+jtot2)/2) 
!      do J3=j3min,j3max,2 
!         sm = sm - sj1* sixj(jr,ja,J3,ju,jtot1,Jpq) * sqrt(J3+1.d0) &
!              * xxxsixj(J3,Jst,rank,jtot2,jtot1,ju) * tensor_elem(ir,a,is,it,J3,Jst,R,jbas)
!      end do 
!   end do 
  
!   smtot = smtot + sm*multfact


!   ! EIGHTH TERM 
!   !changed to q-r instead of q+r
!   multfact = (-1)**((jt+jr+js+jtot2+Jpq+rank)/2) *sqrt((Jst+1.d0) &
!        *(jtot1+1.d0)*(jtot2+1.d0) ) * (rank+1.d0) 
!   ! ju isn't in here because I need it to make add with ja
!   ! so I get an integer later 

!   ! Left-Right 
!   sm = 0.d0   
  
!   do a = 1,jbas%total_orbits 
     
!      ja = jbas%jj(a) 
     
!      j2min = abs(ja - js)
!      j2max = ja+js 
          
!      j1min = max(abs(ja - jr),abs(jt-ju)) 
!      j1max = min(ja+jr,jt+ju)  
          
!      phase = (-1)**((ja+ju)/2) 
     
!      do J1 = j1min, j1max , 2
!         sj1 = phase*(-1)**(J1/2)*(J1+1.d0) * sixj(js,jt,Jst,ju,jtot2,J1)&
!              * v_elem(ir,a,it,iu,J1,L,jbas) 
!         do J2 = j2min,j2max,2
        
!            sm = sm - (-1)**(J2/2) * sqrt(J2 + 1.d0) &            
!             * sj1 * sixj(ja,js,J2,jtot2,jr,J1) *  &
!              xxxsixj(Jpq,J2,rank,jtot2,jtot1,jr) * &
!              tensor_elem(ip,iq,a,is,Jpq,J2,R,jbas)
!         end do
!      end do
!   end do
 
!   smtot = smtot + sm*multfact
  
!   ! !Right-left  
!   sm = 0.d0 
!   multfact = (-1)**((jt+jtot2+Jpq+rank)/2) *sqrt((Jst+1.d0)*(Jpq+1.d0) &
!        *(jtot1+1.d0)*(jtot2+1.d0) ) * (rank+1.d0) 
  
!   do a = 1,jbas%total_orbits 
     
!      ja = jbas%jj(a) 
          
!      if (.not. triangle(ja,js,Jpq)) cycle 
     
!      j1min = abs(jr - ja)
!      j1max = jr + ja  
     
!      j2min = abs(jt - ju) 
!      j2max = jt+ju
     
!      phase = (-1)**((ja - ju)/2)
           
!      sj1 = phase * v_elem(ip,iq,a,is,Jpq,L,jbas)
!      do J1 = j1min, j1max , 2      

!         sj2 = sqrt(J1+1.d0)*(-1)**(J1/2)*sixj(jr,ja,J1,js,jtot1,Jpq) * sj1

!         do J2 = j2min,j2max,2

!            sm = sm + sj2* sqrt(J2 + 1.d0) * (-1)**(J2/2) &
!                 * sixj(js,jt,Jst,ju,jtot2,J2) * xxxsixj(J1,J2,rank,jtot2,jtot1,js) &
!                 * tensor_elem(ir,a,it,iu,J1,J2,R,jbas)

!         end do
!      end do
!   end do 
  
  
!   smtot = smtot + sm*multfact
  
! ! NINTH TERM 
!   !changed to q-r instead of q+r
!   multfact = (-1)**((js - jtot2+Jpq+Jst+rank)/2) *sqrt((Jst+1.d0) &
!        *(jtot1+1.d0)*(jtot2+1.d0) ) * (rank+1.d0) 
!   ! ju isn't in here because I need it to make add with ja
!   ! so I get an integer later 

!   ! Left-Right 
!   sm = 0.d0   
  
!   do a = 1,jbas%total_orbits 
     
!      ja = jbas%jj(a) 
     
!      j2min = abs(ja - jt)
!      j2max = ja+jt 
          
!      j1min = max(abs(ja - jr),abs(js-ju)) 
!      j1max = min(ja+jr,js+ju)  
          
!      phase = (-1)**((ja+jr)/2) 
     
!      do J1 = j1min, j1max , 2
!         sj1 = phase*(J1+1.d0) * sixj(js,jt,Jst,jtot2,ju,J1)&
!              * v_elem(ir,a,iu,is,J1,L,jbas) 
        
!         do J2 = j2min,j2max,2
        
!            sm = sm - (-1)**(J2/2) * sqrt(J2 + 1.d0) &            
!             * sj1 * sixj(ja,jr,J1,jtot2,jt,J2) *  &
!              xxxsixj(Jpq,J2,rank,jtot2,jtot1,jr) * &
!              tensor_elem(ip,iq,a,it,Jpq,J2,R,jbas)
!         end do
!      end do
!   end do
 
!   smtot = smtot + sm*multfact
  
!   ! !Right-left  
!   sm = 0.d0 
!   multfact = (-1)**((jt+jtot1+Jpq+Jst+rank)/2) *sqrt((Jst+1.d0)*(Jpq+1.d0) &
!        *(jtot1+1.d0)*(jtot2+1.d0) ) * (rank+1.d0) 
  
!   do a = 1,jbas%total_orbits 
     
!      ja = jbas%jj(a) 
          
!      if (.not. triangle(ja,jt,Jpq)) cycle 
     
!      j1min = abs(jr - ja)
!      j1max = jr + ja  
     
!      j2min = abs(js - ju) 
!      j2max = js+ju
     
!      phase = (-1)**((ja + js)/2)
           
!      sj1 = phase * v_elem(ip,iq,a,it,Jpq,L,jbas)
!      do J1 = j1min, j1max , 2      

!         sj2 = sqrt(J1+1.d0)*(-1)**(J1/2)*sixj(jr,ja,J1,jt,jtot1,Jpq) * sj1

!         do J2 = j2min,j2max,2

!            sm = sm + sj2* sqrt(J2 + 1.d0)  &
!                 * sixj(js,jt,Jst,jtot2,ju,J2) * xxxsixj(J1,J2,rank,jtot2,jtot1,jt) &
!                 * tensor_elem(ir,a,iu,is,J1,J2,R,jbas)

!         end do
!      end do
!   end do 
  
  
!   smtot = smtot + sm*multfact
    

     
  TS_commutator_223_single = smtot

 end function
!=====================================================================================
!=====================================================================================
end module 
  
  
  
