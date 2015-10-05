module EOM_TS_commutators
  use basic_IMSRG
  ! tensor-scalar commutator functions 
  
  ! THE TENSOR MUST BE THE SECOND ARGUMENT
  
contains
!=========================================================
!=========================================================
subroutine EOM_TS_commutator_111(L,R,RES,jbas) 
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

!dfph~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  call dgemm('N','N',par,hol,par,al,L%fpp,par,R%fph,par,bet,tb2,par) 
  call dgemm('N','N',par,hol,hol,al,R%fph,par,L%fhh,hol,bet,tb1,par) 
  
  RES%fph = tb2 - tb1 
    
end subroutine
!=========================================================
!=========================================================
subroutine EOM_TS_commutator_121(L,R,RES,jbas) 
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
        
        if ( mod(lq,2) .ne. mod(lp+rank/2,2)) cycle
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
             
              sm = sm + L%fph(a,i) * L%herm*smx &
              *(-1)**( (rank + jq + ji )/2 )
       
           end do 
        end do 
        
        RES%fph(p,q) = RES%fph(p,q) + sm  
        
     end do 
  end do 
 
end subroutine             

!====================================================================
!====================================================================
subroutine  EOM_TS_commutator_211(LCC,R,RES,jbas) 
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
        if ( mod(lq,2) .ne. mod(lp+rank/2,2)) cycle
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
              if ( mod(li,2) .ne. mod(la+rank/2,2)) cycle
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
                
              rai = EOMph_rval(ak,ik,Ntot,qx,LCC)
              rpq = EOMTS_rval(pk,qk,Ntot,qx,LCC)
              rqp = EOMTS_rval(qk,pk,Ntot,qx,LCC)

              sm = sm - (-1)**(( jp + jq + rank)/2) * R%fph(a,i) & 
                *   LCC%CCX(qx)%X(rpq,rai)  / sqrt(rank + 1.d0 ) 

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
subroutine EOM_TS_commutator_122(L,R,RES,jbas) 
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
subroutine EOM_TS_commutator_212(L,R,RES,jbas) 
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
              modli = mod(jbas%ll(i)+rank/2,2) 
              
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
subroutine EOM_TS_commutator_221(w1,w2,pm,RES,jbas) 
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
        if (mod(li,2) .ne. mod(lj+w1%rank/2,2))  cycle
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
        if (mod(li,2) .ne. mod(lj+w1%rank/2,2))  cycle
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
        if (mod(li,2) .ne. mod(lj+w1%rank/2,2))  cycle
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
subroutine EOM_TS_commutator_222_pp_hh(L,R,RES,w1,w2,jbas) 
  !VERIFIED
  !NEEDS TO BE RUN BEFORE 221, because it sets up the 
  !intermediary matrices
  implicit none
  
  type(spd) :: jbas
  type(sq_op) ::  L,R,RES,w1,w2
  integer :: q,q1,q2,J1,J2,Tz,Par,phase,rank
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
     q2 = block_index(J2,Tz,mod(Par+Rank/2,2)) 
     
     nh1 = R%tblck(q)%nhh1
     np1 = R%tblck(q)%npp1
     nb1 = R%tblck(q)%nph1
     nh2 = R%tblck(q)%nhh2
     np2 = R%tblck(q)%npp2
     nb2 = R%tblck(q)%nph2
       
!----------------------------------------------------------------------------
!         Zpphh 
!----------------------------------------------------------------------------
     if (np1*nh2 .ne. 0)  then 
        
        !w1pphh = Apppp.Bpphh 
 
        call dgemm('N','N',np1,nh2,np1,al,L%mat(q1)%gam(1)%X,np1,&
             R%tblck(q)%tgam(3)%X,np1,bet,w1%tblck(q)%tgam(3)%X,np1)
        
        !w2pphh = -Bpphh.Ahhhh 
        al_off = -1 
        call dgemm('N','N',np1,nh2,nh2,al_off,R%tblck(q)%tgam(3)%X,np1,&
             L%mat(q2)%gam(5)%X,nh2,bet,w2%tblck(q)%tgam(3)%X,np1)
             
     end if 
       
     RES%tblck(q)%tgam(3)%X = RES%tblck(q)%tgam(3)%X + &
          w1%tblck(q)%tgam(3)%X - w2%tblck(q)%tgam(3)%X

!----------------------------------------------------------------------------
!         Zhhpp 
!----------------------------------------------------------------------------
     if (np2*nh1 .ne. 0)  then 
     
        al_off = -1*L%herm
        !w1hhpp = -Bhhpp.Apppp 
        call dgemm('N','N',nh1,np2,np2,al_off,R%tblck(q)%tgam(7)%X,nh1,&
             L%mat(q2)%gam(1)%X,np2,bet,w1%tblck(q)%tgam(7)%X,nh1)
        
        
        al_off = L%herm 
        !w1pphh = Ahhhh.Bhhpp 
        call dgemm('N','N',nh1,np2,nh1,al_off,L%mat(q1)%gam(5)%X,nh1,&
             R%tblck(q)%tgam(7)%X,nh1,bet,w2%tblck(q)%tgam(7)%X,nh1)

     end if 
       
     RES%tblck(q)%tgam(7)%X = RES%tblck(q)%tgam(7)%X - &
          w1%tblck(q)%tgam(7)%X + w2%tblck(q)%tgam(7)%X

!----------------------------------------------------------------------------
!         Zppph 
!----------------------------------------------------------------------------
     if (np1*nb2 .ne. 0)  then 
     
     
        if (nh2 .ne. 0) then 
           !w2ppph = -Bpphh.Ahhph
           al_off = -1*L%herm
           call dgemm('N','T',np1,nb2,nh2,al_off,R%tblck(q)%tgam(3)%X,np1,&
                L%mat(q2)%gam(6)%X,nb2,bet,w2%tblck(q)%tgam(2)%X,np1)
        end if 
        
     end if 
       
!----------------------------------------------------------------------------
!         Zphpp 
!----------------------------------------------------------------------------
     if (nb1*np2 .ne. 0)  then 
        
        if (nh1 .ne. 0) then 
           al_off = -1*L%herm 
           !w2phpp = Aphhh.Bhhpp
           call dgemm('N','N',nb1,np2,nh1,al_off,L%mat(q1)%gam(6)%X,nb1,&
                R%tblck(q)%tgam(7)%X,nh1,bet,w2%tblck(q)%tgam(8)%X,nb1)
        end if 
        
     end if 
       

!----------------------------------------------------------------------------
!         Zphhh 
!----------------------------------------------------------------------------
     if (nb1*nh2 .ne. 0)  then 
     
        if (np1 .ne. 0) then 
           !w1phhh = Aphpp.Bpphh
           al_off = L%herm
           call dgemm('T','N',nb1,nh2,np1,al_off,L%mat(q1)%gam(2)%X,np1,&
                R%tblck(q)%tgam(3)%X,np1,bet,w1%tblck(q)%tgam(6)%X,nb1)
        end if 
 
     end if 

!----------------------------------------------------------------------------
!         Zhhph 
!----------------------------------------------------------------------------
     if (nh1*nb2 .ne. 0)  then 
     
        if (np2 .ne. 0)  then 
           al_off = L%herm
           !w1hhph = -Bhhpp.Appph 
           call dgemm('N','N',nh1,nb2,np2,al_off,R%tblck(q)%tgam(7)%X,nh1,&
                L%mat(q2)%gam(2)%X,np2,bet,w1%tblck(q)%tgam(9)%X,nh1)
        end if
     end if 
  end do
  
end subroutine 
!=================================================================
!=================================================================
 subroutine EOM_TS_commutator_222_ph(LCC,RCC,RES,WCC,jbas) 
   ! VERIFIED ph channel 2body EOM_TS_commutator. DFWT! 
   implicit none 
  
   type(spd) :: jbas
   type(sq_op) :: RES
   type(cross_coupled_31_mat) :: LCC,RCC,WCC
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
        factor = 1.d0/sqrt(J1+1.d0)
        q1 = J1/2+1 + Tz*(JTM+1) + 2*PAR*(JTM+1)
        
        call dgemm('T','N',r1,r2,nb1,factor,LCC%CCR(q1)%X,nb1,&
             RCC%CCR(q)%X,nb1,bet,WCC%CCX(q)%X,r1) 
      
     end if
         
     if (nb2 .ne. 0 ) then 
        
        J2 = RCC%Jval2(q) 
        PAR2 = mod(PAR+rank/2,2) 
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

            if (mod(lk+lj+rank/2,2) == PAR) then 
               if (abs(tk - tj) == Tz*2)  then 
             
            do J3 = J3min,J3max,2
            
               q1 = block_index(J3,Tz,PAR)
              
               ril = EOMTS_rval(i,l,Ntot,q1,RCC)
               rli = EOMTS_rval(l,i,Ntot,q1,RCC)

                do J4 = max( J3 , J4min ) , J4max,2 
            
                  if (.not. (triangle(J3,J4,rank))) cycle
                  
                  PAR2 = mod(PAR + rank/2,2) 
                  q2 = block_index(J4,Tz,PAR2)
               
                  rjk = EOMTS_rval(j,k,Ntot,q2,RCC)
                  rkj = EOMTS_rval(k,j,Ntot,q2,RCC)
              
                  qx = CCtensor_block_index(J3,J4,rank,Tz,PAR)
                  sm = sm + sqrt((J3+1.d0)*(J4+1.d0))* &
                       ninej(ji,jl,J3,jj,jk,J4,J1,J2,rank)  * ( &
                       WCC%CCX(qx)%X(ril,rkj)*(-1)**((J3+J4)/2)  &
                      + WCC%CCR(qx)%X(ril,rkj) * phase1 * LCC%herm &
                      - (-1)**(rank/2)*phase1*RCC%herm*LCC%herm * &
                      WCC%CCX(qx)%X(rli,rjk) - (-1)**(( J3+J4+rank)/2) &
                      * RCC%herm * WCC%CCR(qx)%X(rli,rjk) )

                end do 
               
                do J4 = J4min , min(J4max,J3-2),2 
                  if (.not. (triangle(J3,J4,rank))) cycle
                  
                  
                  PAR2 = mod(PAR + rank/2,2)                  
                  q2 = block_index(J4,Tz,PAR2)
                  
                  rjk = EOMTS_rval(j,k,Ntot,q2,RCC)     
                  rkj = EOMTS_rval(k,j,Ntot,q2,RCC)
                  
                  qx = CCtensor_block_index(J4,J3,rank,Tz,PAR2)
                  sm = sm + sqrt((J3+1.d0)*(J4+1.d0))* &
                       ninej(ji,jl,J3,jj,jk,J4,J1,J2,rank)  * ( &
                       WCC%CCR(qx)%X(rjk,rli)*phase1*(-1)**((rank+J3+J4)/2)*LCC%herm &
                       + WCC%CCX(qx)%X(rjk,rli) * (-1)**(rank/2) &
                        - RCC%herm * WCC%CCR(qx)%X(rkj,ril) &
                     - phase1*(-1)**((J3+J4)/2) * RCC%herm *LCC%herm * &
                          WCC%CCX(qx)%X(rkj,ril) )
                      
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

            if (mod(li+lk+rank/2,2) == PAR) then 
               if (abs(ti - tk) == Tz*2)  then 

            do J3 = J3min,J3max,2
               q1 = block_index(J3,Tz,PAR)
            
               rjl = EOMTS_rval(j,l,Ntot,q1,RCC)
               rlj = EOMTS_rval(l,j,Ntot,q1,RCC)

               do J4 = max( J3 , J4min ) , J4max,2 
                 
                  if (.not. (triangle(J3,J4,rank))) cycle
                 
                  
                  PAR2 = mod(PAR + rank/2,2)
                  q2 = block_index(J4,Tz,PAR2)
             
                  rki = EOMTS_rval(k,i,Ntot,q2,RCC)
                  rik = EOMTS_rval(i,k,Ntot,q2,RCC)
              
                  qx = CCtensor_block_index(J3,J4,rank,Tz,PAR)
               
                  sm_ex = sm_ex - sqrt((J3+1.d0)*(J4+1.d0))* &
                       ninej(jj,jl,J3,ji,jk,J4,J1,J2,rank) * ( &
                       (-1)**((J3+J4)/2) * WCC%CCX(qx)%X(rjl,rki)  & 
                       - phase1*RCC%herm*LCC%herm*(-1)**(rank/2)* &
                       WCC%CCX(qx)%X(rlj,rik) + phase1*LCC%herm * &
                       WCC%CCR(qx)%X(rjl,rki) - (-1)**((J3+J4+rank)/2)* &
                       RCC%herm * WCC%CCR(qx)%X(rlj,rik) )
                  
                end do 
               
                do J4 = J4min , min(J4max,J3-2),2 
                  if (.not. (triangle(J3,J4,rank))) cycle
                  
                  PAR2 = mod(PAR + rank/2,2)
                  q2 = block_index(J4,Tz,PAR2)
     
                  rki = EOMTS_rval(k,i,Ntot,q2,RCC)
                  rik = EOMTS_rval(i,k,Ntot,q2,RCC)

                  qx = CCtensor_block_index(J4,J3,rank,Tz,PAR2)
                  
                  sm_ex = sm_ex - sqrt((J3+1.d0)*(J4+1.d0))* &
                       ninej(jj,jl,J3,ji,jk,J4,J1,J2,rank) * ( &
                       (-1)**((J3+J4+rank)/2) *phase1*LCC%herm &
                       * WCC%CCR(qx)%X(rik,rlj)  & 
                       - RCC%herm*WCC%CCR(qx)%X(rki,rjl) &
                       + (-1)**(rank/2) * WCC%CCX(qx)%X(rik,rlj) &
                       - phase1*(-1)**((J3+J4)/2) * LCC%herm * RCC%herm * &
                       WCC%CCX(qx)%X(rki,rjl) ) 
                  
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
integer function EOMTS_rval(i,l,Ntot,q,LCC) 
  implicit none 
  
  type(cross_coupled_31_mat) :: LCC
  integer :: i,l,Ntot,x,g,q
  
  x = CCindex(i,l,Ntot)
  g = 1
  do while (LCC%qmap(x)%Z(g) .ne. q )
  
     g = g + 1
  end do
  
  EOMTS_rval = LCC%rmap(x)%Z(g)
end function 
!============================================
!============================================
integer function EOMph_rval(i,l,Ntot,q,LCC) 
  implicit none 
  
  type(cross_coupled_31_mat) :: LCC
  integer :: i,l,Ntot,x,g,q
  
  x = CCindex(i,l,Ntot)
  g = 1

  do while (LCC%qmap(x)%Z(g) .ne. q )
     g = g + 1
  end do
  
  EOMph_rval = LCC%nbmap(x)%Z(g)
end function 
!=====================================================
!=====================================================      
end module 
  
  
  
