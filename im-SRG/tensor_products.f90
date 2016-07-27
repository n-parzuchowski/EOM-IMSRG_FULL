module tensor_products
  use cross_coupled
  ! tensor-scalar commutator functions 
  
  ! THE TENSOR MUST BE THE SECOND ARGUMENT
  
contains

!===================================================================
!===================================================================
subroutine tensor_product_222_pp_hh(AA,BB,CC,jbas) 
  !VERIFIED
  !NEEDS TO BE RUN BEFORE 221, because it sets up the 
  !intermediary matrices
  implicit none
  
  type(spd) :: jbas
  type(sq_op) ::  AA,BB,CC
  integer :: q,q1,q2,J1,J2,Tz,Par,phase,rank_A,rank_B,rank_c,i,Cpar
  integer :: np1,nb1,nh1,np2,nb2,nh2,pm,P_a,P_b,J3min,J3max,J3
  real(8) :: bet_off,al_off,d6ji
  real(8),allocatable,dimension(:,:) :: W 
  
!  pm = AA%herm*BB%herm
  rank_A = AA%rank
  rank_B = BB%rank
  rank_C = CC%rank 
!construct temporary matrices
  do q = 1, CC%nblocks
     
     J1 = CC%tblck(q)%Jpair(1) 
     J2 = CC%tblck(q)%Jpair(2)

     phase = CC%tblck(q)%lam(1)
     par = CC%tblck(q)%lam(2) 
     Tz = CC%tblck(q)%lam(3)
             
     J3Min = max(abs(rank_a-J1),abs(rank_b-J2))
     J3Max = min(rank_a + J1,rank_b + J2)
     np1 = CC%tblck(q)%npp1
     nh2 = CC%tblck(q)%nhh2

     do J3 = J3Min, J3Max, 2
!        if ( J3 > 2*jbas%jtotal_max) cycle
        If ( J3 .ge. J1 )  then

           IF (J2 .ge. J3 ) then 
           
              P_a = CC%tblck(q)%lam(2)
              Tz = CC%tblck(q)%lam(3) 
              P_b = mod(P_a + AA%dpar/2,2)
              
              q1 = tensor_block_index(J1,J3,rank_a,Tz,P_a) 
              q2 = tensor_block_index(J3,J2,rank_b,Tz,P_b) 
     
!----------------------------------------------------------------------------
!         Zpphh 
!----------------------------------------------------------------------------

              !DICK 
              np2 = AA%tblck(q1)%npp2
                            
              if (np1*nh2 .ne. 0)  then                       
                 if (np2 .ne. 0) then 
                    !w1pphh = Bpppp.Apphh
                    al_off = d6ji(rank_a,rank_b,rank_c,J2,J1,J3)*&
                         (-1)**((J1+J2+rank_c)/2)*sqrt(rank_c+1.d0)
                    bet_off = 1.d0

                   call dgemm('N','N',np1,nh2,np2,al_off,AA%tblck(q1)%tgam(1)%X,np1,&
                        BB%tblck(q2)%tgam(3)%X,np2,bet_off,CC%tblck(q)%tgam(3)%X,np1)
                 end if
              end if
                   
                              
           else ! ( J2< J3) 
                                      
              P_a = CC%tblck(q)%lam(2)
              Tz = CC%tblck(q)%lam(3) 
              P_b = mod(P_a + CC%dpar/2,2) 
              
              q1 = tensor_block_index(J1,J3,rank_a,Tz,P_a) 
              q2 = tensor_block_index(J2,J3,rank_b,Tz,P_b) 
     
!----------------------------------------------------------------------------
!         Zpphh 
!----------------------------------------------------------------------------
              np2 = AA%tblck(q1)%npp2
              nh1 = BB%tblck(q2)%nhh1
              
              if (np1*nh2 .ne. 0)  then                       
                 if (np2*nh1 .ne. 0) then 
                    !w1pphh = Bpppp.Apphh
                    al_off = d6ji(rank_a,rank_b,rank_c,J2,J1,J3)*&
                         (-1)**((J1+J2+rank_c)/2)*sqrt(rank_c+1.d0)
                    bet_off = 1.d0
                    call dgemm('N','T',np1,nh2,np2,al_off,AA%tblck(q1)%tgam(1)%X,np1,&
                         BB%tblck(q2)%tgam(7)%X,nh1,bet_off,CC%tblck(q)%tgam(3)%X,np1)
                 end if
              end if
              
              
           end if

        else ! ( J3 < J1 .le. J2 )
           if ( J2 .ge. J1) then

              P_a = mod(CC%tblck(q)%lam(2)+AA%dpar/2,2) 
              Tz = CC%tblck(q)%lam(3) 
              P_b = P_a
              
              q1 = tensor_block_index(J3,J1,rank_a,Tz,P_a) 
              q2 = tensor_block_index(J3,J2,rank_b,Tz,P_b) 
     
!----------------------------------------------------------------------------
!         Zpphh 
!----------------------------------------------------------------------------

              !DICK 
              np2 = AA%tblck(q1)%npp1
                            
              if (np1*nh2 .ne. 0)  then                       
                 if (np2 .ne. 0) then 
                    !w1pphh = Bpppp.Apphh
                    al_off = d6ji(rank_a,rank_b,rank_c,J2,J1,J3)*AA%tblck(q1)%lam(1)*&
                         (-1)**((J1+J2+rank_c)/2)*sqrt(rank_c+1.d0)
                    bet_off = 1.d0

                   call dgemm('T','N',np1,nh2,np2,al_off,AA%tblck(q1)%tgam(1)%X,np2,&
                        BB%tblck(q2)%tgam(3)%X,np2,bet_off,CC%tblck(q)%tgam(3)%X,np1)
                 end if
              end if
                   
                   
              
            end if

              
           
        end if

     end do

     J3Min = max(abs(rank_b-J1),abs(rank_a-J2))
     J3Max = min(rank_b + J1,rank_a + J2)

     
     do J3 = J3Min, J3Max, 2
!        if ( J3 > 2*jbas%jtotal_max) cycle        
        If ( J3 .ge. J1 )  then

           IF (J2 .ge. J3 ) then 


              Tz = CC%tblck(q)%lam(3) 

              !----------------------------------------------------------------------------
              !         Zpphh 
              !----------------------------------------------------------------------------

              !DICK 

              P_b = CC%tblck(q)%lam(2)
              P_a = mod(P_b + BB%dpar/2,2)

              q1 = tensor_block_index(J1,J3,rank_b,Tz,P_b) 
              q2 = tensor_block_index(J3,J2,rank_a,Tz,P_a) 

              nh1 = BB%tblck(q1)%nhh2     
              !w1pphh = -Bpphh.Ahhhh 

              if (np1*nh2 .ne. 0)  then                       
                 if (nh1 .ne. 0) then 
                    al_off = d6ji(rank_a,rank_b,rank_c,J1,J2,J3)*&
                         (-1)**((J1+J2+rank_a+rank_b)/2)*sqrt(rank_c+1.d0)
                    bet_off = 1.d0
                    call dgemm('N','N',np1,nh2,nh1,al_off,BB%tblck(q1)%tgam(3)%X,np1,&
                         AA%tblck(q2)%tgam(5)%X,nh1,bet_off,CC%tblck(q)%tgam(3)%X,np1)
                 end if
              end if

           else ! ( J2< J3) 

              Tz = CC%tblck(q)%lam(3) 

              !----------------------------------------------------------------------------
              !         Zpphh 
              !----------------------------------------------------------------------------

              P_b = CC%tblck(q)%lam(2)
              P_a = mod(P_b + CC%dpar/2,2)

              q1 = tensor_block_index(J1,J3,rank_b,Tz,P_b) 
              q2 = tensor_block_index(J2,J3,rank_a,Tz,P_a) 

              nh1 = BB%tblck(q1)%nhh2     
              !w1pphh = -Bpphh.Ahhhh 

              if (np1*nh2 .ne. 0)  then                       
                 if (nh1 .ne. 0) then 
                    al_off = d6ji(rank_a,rank_b,rank_c,J1,J2,J3)*AA%tblck(q2)%lam(1)*&
                         (-1)**((J1+J2+rank_a+rank_b)/2)*sqrt(rank_c+1.d0)
                    bet_off = 1.d0
                    call dgemm('N','T',np1,nh2,nh1,al_off,BB%tblck(q1)%tgam(3)%X,np1,&
                         AA%tblck(q2)%tgam(5)%X,nh2,bet_off,CC%tblck(q)%tgam(3)%X,np1)
                 end if
              end if


           end if

        else ! ( J3 < J1 .le. J2 )
           if ( J2 .ge. J1) then
              

              Tz = CC%tblck(q)%lam(3) 
              
           !----------------------------------------------------------------------------
           !         Zpphh 
           !----------------------------------------------------------------------------

              !DICK 
              
              P_b = mod(CC%tblck(q)%lam(2)+BB%dpar/2,2)
              P_a = P_b

              q1 = tensor_block_index(J3,J1,rank_b,Tz,P_b) 
              q2 = tensor_block_index(J3,J2,rank_a,Tz,P_a) 

              nh1 = BB%tblck(q1)%nhh1     
              !w1pphh = -Bpphh.Ahhhh 

              if (np1*nh2 .ne. 0)  then                       
                 if (nh1 .ne. 0) then 
                    al_off = d6ji(rank_a,rank_b,rank_c,J1,J2,J3)*&
                         (-1)**((J1+J2+rank_a+rank_b)/2)*sqrt(rank_c+1.d0)
                    bet_off = 1.d0
                    call dgemm('T','N',np1,nh2,nh1,al_off,BB%tblck(q1)%tgam(7)%X,nh1,&
                         AA%tblck(q2)%tgam(5)%X,nh1,bet_off,CC%tblck(q)%tgam(3)%X,np1)
                 end if
              end if

           end if



        end if

     end do

     !cocks

     ! NOW DO THE OPPOSITE SIDE OF THE MATRIX (J1 > J2)

     ! just switching the names so that I don't have to change
     ! variables in the expression below. 
     J3 = J1
     J1 = J2
     J2 = J3 

     if (J1 > jbas%jtotal_max*2) cycle
     
     J3Min = max(abs(rank_a-J1),abs(rank_b-J2))
     J3Max = min(rank_a + J1,rank_b + J2)
     Cpar= mod(CC%tblck(q)%lam(2)+ CC%dpar/2,2)

     np1 = CC%tblck(q)%npp2
     nh2 = CC%tblck(q)%nhh1

     allocate(W(np1,nh2)) 
     W = 0.d0
     do J3 = J3Min, J3Max, 2
        if ( J3 > 2*jbas%jtotal_max) cycle
        If ( J3 .ge. J1 )  then
           IF (J1 .ge. J2 ) then 
           
              
              Tz = CC%tblck(q)%lam(3) 
              
              P_a = CC%tblck(q)%lam(2)
              P_b = mod(P_b + CC%dpar/2,2)

              q1 = tensor_block_index(J1,J3,rank_a,Tz,P_a) 
              q2 = tensor_block_index(J2,J3,rank_b,Tz,P_b) 
!----------------------------------------------------------------------------
!         Zpphh 
!----------------------------------------------------------------------------

              !DICK 
              np2 = AA%tblck(q1)%npp2
              nh1 = BB%tblck(q2)%nhh1
                            
              if (np1*nh2 .ne. 0)  then                       
                 if (np2*nh1 .ne. 0) then 
                    !w1pphh = Bpppp.Apphh
                    al_off = d6ji(rank_a,rank_b,rank_c,J2,J1,J3)*&
                         (-1)**((J1+J2+rank_c)/2)*sqrt(rank_c+1.d0)
                    bet_off = 1.d0
                    
                    call dgemm('N','T',np1,nh2,np2,al_off,AA%tblck(q1)%tgam(1)%X,np1,&
                         BB%tblck(q2)%tgam(7)%X,nh1,bet_off,W,np1)
                    
                 end if
              end if
              
           end if
        else ! J3 < J1 
           if (J2 .ge. J3)  then 
                                                    
!               Tz = CC%tblck(q)%lam(3) 
              
!               P_a = mod(CC%tblck(q)%lam(2)+AA%dpar/2,2)
!               P_b = P_a

!               q1 = tensor_block_index(J3,J1,rank_a,Tz,P_a) 
!               q2 = tensor_block_index(J3,J2,rank_b,Tz,P_b) 
     
! !----------------------------------------------------------------------------
! !         Zpphh 
! !----------------------------------------------------------------------------
!               np2 = AA%tblck(q1)%npp1
              
!               if (np1*nh2 .ne. 0)  then                       
!                  if (np2 .ne. 0) then 
!                     !w1pphh = Bpppp.Apphh
!                     al_off = d6ji(rank_a,rank_b,rank_c,J2,J1,J3)*AA%tblck(q1)%lam(1)*&
!                          (-1)**((J1+J2+rank_c)/2)*sqrt(rank_c+1.d0)
!                     bet_off = 1.d0
!                     call dgemm('T','N',np1,nh2,np2,al_off,AA%tblck(q1)%tgam(1)%X,np2,&
!                          BB%tblck(q2)%tgam(3)%X,np2,bet_off,W,np1)
!                  end if
!               end if
              
              
              
           else ! ( J3 < J1 ,  J3 >= J2  )
                 
!               P_a = mod(Cpar+AA%dpar/2,2) 
!               P_b = mod(Cpar+CC%dpar/2,2) 
!               Tz = CC%tblck(q)%lam(3) 
              
!               q1 = tensor_block_index(J3,J1,rank_a,Tz,P_a) 
!               q2 = tensor_block_index(J2,J3,rank_b,Tz,P_b) 
     
! !----------------------------------------------------------------------------
! !         Zpphh 
! !----------------------------------------------------------------------------

!               !DICK 
!               np2 = AA%tblck(q1)%npp1
!               nh1 = BB%tblck(q2)%nhh1
                            
!               if (np1*nh2 .ne. 0)  then                       
!                  if (np2*nh1 .ne. 0) then 
!                     !w1pphh = Bpppp.Apphh
!                     al_off = d6ji(rank_a,rank_b,rank_c,J2,J1,J3)*AA%tblck(q1)%lam(1)*&
!                          (-1)**((J1+J2+rank_c)/2)*sqrt(rank_c+1.d0)
!                     bet_off = 1.d0

!                    call dgemm('T','T',np1,nh2,np2,al_off,AA%tblck(q1)%tgam(1)%X,np2,&
!                         BB%tblck(q2)%tgam(7)%X,nh1,bet_off,W,np1)
!                  end if
!               end if
                   
                   
              
            end if

              
           
         end if
         
      end do

      ! J3Min = max(abs(rank_b-J1),abs(rank_a-J2))
     ! J3Max = min(rank_b + J1,rank_a + J2)

     
     ! do J3 = J3Min, J3Max, 2
        
     !    If ( J3 .ge. J1 )  then

     !       IF (J2 .ge. J3 ) then 


     !          Tz = CC%tblck(q)%lam(3) 

     !          !----------------------------------------------------------------------------
     !          !         Zpphh 
     !          !----------------------------------------------------------------------------

     !          !DICK 

     !          P_b = Cpar
     !          P_a = mod(P_b + BB%dpar/2,2)

     !          q1 = tensor_block_index(J1,J3,rank_b,Tz,P_b) 
     !          q2 = tensor_block_index(J3,J2,rank_a,Tz,P_a) 

     !          nh1 = BB%tblck(q1)%nhh2     
     !          !w1pphh = -Bpphh.Ahhhh 

     !          if (np1*nh2 .ne. 0)  then                       
     !             if (nh1 .ne. 0) then 
     !                al_off = d6ji(rank_a,rank_b,rank_c,J1,J2,J3)*&
     !                     (-1)**((J1+J2+rank_a+rank_b)/2)*sqrt(rank_c+1.d0)
     !                bet_off = 1.d0
     !                call dgemm('N','N',np1,nh2,nh1,al_off,BB%tblck(q1)%tgam(3)%X,np1,&
     !                     AA%tblck(q2)%tgam(5)%X,nh1,bet_off,W,np1)
     !             end if
     !          end if

     !       else ! ( J2< J3) 

     !          Tz = CC%tblck(q)%lam(3) 

     !          q1 = tensor_block_index(J1,J3,rank_a,Tz,P_a) 
     !          q2 = tensor_block_index(J2,J3,rank_b,Tz,P_b) 

     !          !----------------------------------------------------------------------------
     !          !         Zpphh 
     !          !----------------------------------------------------------------------------

     !          P_b = Cpar
     !          P_a = mod(P_b + CC%dpar/2,2)

     !          q1 = tensor_block_index(J1,J3,rank_b,Tz,P_b) 
     !          q2 = tensor_block_index(J2,J3,rank_a,Tz,P_a) 

     !          nh1 = BB%tblck(q1)%nhh2     
     !          !w1pphh = -Bpphh.Ahhhh 

     !          if (np1*nh2 .ne. 0)  then                       
     !             if (nh1 .ne. 0) then 
     !                al_off = d6ji(rank_a,rank_b,rank_c,J1,J2,J3)*AA%tblck(q2)%lam(1)*&
     !                     (-1)**((J1+J2+rank_a+rank_b)/2)*sqrt(rank_c+1.d0)
     !                bet_off = 1.d0
     !                call dgemm('N','T',np1,nh2,nh1,al_off,BB%tblck(q1)%tgam(3)%X,np1,&
     !                     AA%tblck(q2)%tgam(5)%X,nh2,bet_off,W,np1)
     !             end if
     !          end if


     !       end if

     !    else ! ( J3 < J1 .le. J2 )
     !       if ( J2 .ge. J1) then
              

     !          Tz = CC%tblck(q)%lam(3) 
              
     !       !----------------------------------------------------------------------------
     !       !         Zpphh 
     !       !----------------------------------------------------------------------------

     !          !DICK 
              
     !          P_b = mod(Cpar+BB%dpar/2,2)
     !          P_a = P_b

     !          q1 = tensor_block_index(J3,J1,rank_b,Tz,P_b) 
     !          q2 = tensor_block_index(J3,J2,rank_a,Tz,P_a) 

     !          nh1 = BB%tblck(q1)%nhh1     
     !          !w1pphh = -Bpphh.Ahhhh 

     !          if (np1*nh2 .ne. 0)  then                       
     !             if (nh1 .ne. 0) then 
     !                al_off = d6ji(rank_a,rank_b,rank_c,J1,J2,J3)*&
     !                     (-1)**((J1+J2+rank_a+rank_b)/2)*sqrt(rank_c+1.d0)
     !                bet_off = 1.d0
     !                call dgemm('T','N',np1,nh2,nh1,al_off,BB%tblck(q1)%tgam(7)%X,nh1,&
     !                     AA%tblck(q2)%tgam(5)%X,nh1,bet_off,W,np1)
     !             end if
     !          end if

     !       end if



     !    end if
  
     ! end do

     CC%tblck(q)%tgam(7)%X = CC%tblck(q)%tgam(7)%X + Transpose(W) 
     deallocate(W) 
  end do
  
end subroutine tensor_product_222_pp_hh
!=================================================================
!=================================================================
 subroutine tensor_product_222_ph(LCC,RCC,R,RES,jbas) 
   ! VERIFIED ph channel 2body tensor_product. DFWT! 
   implicit none 
  
   type(spd) :: jbas
   type(sq_op) :: RES,R
   type(pandya_mat) :: RCC
   real(8),allocatable,dimension(:,:) :: Wx,Wy 
   type(cc_mat) :: LCC
   integer :: nh,np,nb1,nb2,q,IX,JX,i,j,k,l,r1,r2,Tz,PAR,JTM,q1,q2,J3,J4,rank,a,b,c,d
   integer :: ji,jj,jk,jl,ti,tj,tk,tl,li,lj,lk,ll,n1,n2,c1,c2,jxstart,J4min,J4max,ja,jb,jc,jd
   integer :: J1,J2, Jtot,Ntot,qx,J3min,J3max,ril,rjk,rli,rkj,g_ix,thread,total_threads
   integer :: phase1,phase2,phase3,rik,rki,rjl,rlj,PAR2,J1min,J2min,J1max,J2max
   integer :: phase_34,phase_abcd,phase_ac,phase_bc
   integer :: phase_bd,phase_ad,nj_perm,full_int_phase  
   real(8) :: sm ,pre,pre2,omp_get_wtime ,t1,t2,coef9,factor,sm_ex, nj1,nj2  
   real(8) :: prefac_34,prefac_134,prefac_1234,Xelem,Yelem,V
   logical :: square
   
   rank = RES%rank
   Ntot = RES%Nsp
   JTM = jbas%Jtotal_max
   total_threads = size(RES%direct_omp) - 1
   ! construct intermediate matrices

   do qx = 1,RCC%nblocks
      if (RCC%jval2(qx) > jbas%jtotal_max*2) cycle
      
      nb2 = RCC%nb2(qx)
      nb1 = RCC%nb1(qx)
      r1 = size(RCC%qn1(qx)%Y(:,1))
      r2 = size(RCC%qn2(qx)%Y(:,1))      
      
      if (r1 * r2 == 0) cycle
      
      allocate(Wx(r1,r2),Wy(r1,r2)) 
      Wx = 0.d0 
      Wy = 0.d0
      PAR = mod(qx-1,2)
      Tz = mod((qx-1)/2,2) 
      J3 = RCC%Jval(qx) 
      J4 = RCC%Jval2(qx)          
      
      if (nb1 .ne. 0 ) then 

         allocate(RCC%CCR(qx)%X(nb1,r2))
         call calculate_single_pandya(R,RCC,jbas,qx,2)  
         factor = 1.0/sqrt(J3+1.d0)*LCC%herm
         q1 = J3/2+1 + Tz*(JTM+1) + 2*PAR*(JTM+1) 
         call dgemm('N','N',r1,r2,nb1,factor,LCC%CCX(q1)%X,r1,&
              RCC%CCR(qx)%X,nb1,bet,Wx,r1)       
         deallocate(RCC%CCR(qx)%X)
      end if
         
      if ((nb2 .ne. 0) .and. (J3 .ne. J4)) then 
         allocate(RCC%CCX(qx)%X(r1,nb2))
         call calculate_single_pandya(R,RCC,jbas,qx,1)
         PAR2 = mod(PAR+RCC%dpar/2,2) 
         q2 = J4/2+1 + Tz*(JTM+1) + 2*PAR2*(JTM+1)
         factor = 1.d0/sqrt(J4+1.d0)
         
        call dgemm('N','T',r1,r2,nb2,factor,RCC%CCX(qx)%X,r1,&
             LCC%CCX(q2)%X,r2,bet,Wy,r1) 
        deallocate(RCC%CCX(qx)%X)
      end if

      
      prefac_34 = sqrt((J3+1.d0)*(J4+1.d0))
      phase_34 = (-1)**((J3+J4)/2) 

!!!      !$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) SHARED(RCC,LCC,RES,Wx,Wy)
      do IX = 1,r1
         
         ! GET BRA
         d = RCC%qn1(qx)%Y(IX,1)
         a = RCC%qn1(qx)%Y(IX,2)

         ja = jbas%jj(a)
         jd = jbas%jj(d)

         do JX = 1, r2 
                       
            ! GET KET 
            c = RCC%qn2(qx)%Y(JX,1) 
            b = RCC%qn2(qx)%Y(JX,2)            


            jc = jbas%jj(c)
            jb = jbas%jj(b)
            
            ! CALCULATE X CONTRIBUTIONS
            
            J1min = abs(ja-jb)
            J1max = ja+jb
            
            J2min = abs(jc-jd) 
            J2max = jc+jd 
            
            ! these are the results of the Matmuls 
            Xelem = Wx(IX,JX)
            Yelem = Wy(IX,JX)

            phase_abcd= (-1)**((ja+jb+jc+jd)/2)
            
            if (abs(Xelem) > 1e-6) then 
               if (b .ge. a) then 
                  if (d .ge. c) then 
                     
                     ! CALCULATE V^{J1 J2}_{abcd} and V^{J1 J2}_{cdab}
                     
                     do J1 = J1min,J1max,2
                        if ((a==b).and.(mod(J1/2,2)==1)) cycle
                        prefac_134 = prefac_34 * sqrt((J1+1.d0))
                        do J2 = max(J1,J2min,abs(rank-J1)),min(J2max,rank+J1),2 
                           if ((c==d).and.(mod(J2/2,2)==1)) cycle
                           prefac_1234 = prefac_134 * sqrt((J2+1.d0))

                           ! V^{J1 J2}_{abcd} 
                           V = prefac_1234* ninej(ja,jd,J3,jb,jc,J4,J1,J2,rank) &
                                * (-1)**((jb+jd+J2)/2) * phase_34 * LCC%herm &
                                * Xelem

                           call add_elem_to_tensor(V,a,b,c,d,J1,J2,RES,jbas) 

                        end do
                     end do

                     do J1 = J2min,J2max,2
                        if ((c==d).and.(mod(J1/2,2)==1)) cycle
                        prefac_134 = prefac_34 * sqrt((J1+1.d0))
                        do J2 = max(J1,J1min,abs(J1-rank)),min(J1max,J1+rank),2 
                           if ((a==b).and.(mod(J2/2,2)==1)) cycle
                           prefac_1234 = prefac_134 * sqrt((J2+1.d0))

                           ! V^{J1 J2}_{cdab}
                           V = -1*prefac_1234* ninej(jd,ja,J3,jc,jb,J4,J1,J2,rank) &
                                * (-1)**((jc+ja+J1+rank)/2) * RCC%herm &
                                * Xelem
                           call add_elem_to_tensor(V,c,d,a,b,J1,J2,RES,jbas) 

                        end do
                     end do

                  else 

                     ! CALCULATE V^{J1 J2}_{abdc} and V^{J1 J2}_{dcab}

                     do J1 = J1min,J1max,2
                        if ((a==b).and.(mod(J1/2,2)==1)) cycle
                        prefac_134 = prefac_34 * sqrt((J1+1.d0))
                        do J2 = max(J1,J2min,abs(J1-rank)),min(J2max,J1+rank),2 
                           prefac_1234 = prefac_134 * sqrt((J2+1.d0))

                           !V^{J1 J2}_{abdc}
                           V = prefac_1234 * ninej(ja,jd,J3,jb,jc,j4,J1,J2,rank)&
                                *(-1)**((jb+jc)/2)*phase_34*LCC%herm*Xelem

                           call add_elem_to_tensor(V,a,b,d,c,J1,J2,RES,jbas)                        
                        end do
                     end do

                     do J1 = J2min,J2max,2
                        prefac_134 = prefac_34 * sqrt((J1+1.d0))
                        do J2 = max(J1,J1min,abs(J1-rank)),min(J1max,J1+rank),2 
                           if ((a==b).and.(mod(J2/2,2)==1)) cycle
                           prefac_1234 = prefac_134 * sqrt((J2+1.d0))

                           !V^{J1 J2}_{dcab}
                           V = prefac_1234 * ninej(jd,ja,J3,jc,jb,J4,J1,J2,rank) &
                                * (-1)**((ja-jd+rank)/2) * RCC%herm *Xelem
                           ! the mapping here inverts the indeces to {cdba} so you need 
                           ! an additional factor of phase(a+b+c+d+J1+J2) 
                           call add_elem_to_tensor(V,d,c,a,b,J1,J2,RES,jbas)

                        end do
                     end do
                  end if

               else 
                  if (d .ge. c) then 

                     ! CALCULATE V^{J1 J2}_{bacd} and V^{J1 J2}_{cdba}
                     do J1 = J1min,J1max,2
                        prefac_134 = prefac_34 * sqrt((J1+1.d0))
                        do J2 = max(J1,J2min,abs(J1-rank)),min(J2max,J1+rank),2 
                           if ((c==d).and.(mod(J2/2,2)==1)) cycle
                           prefac_1234 = prefac_134 * sqrt((J2+1.d0))

                           !V^{J1 J2}_{bacd}

                           V = prefac_1234 * ninej(ja,jd,J3,jb,jc,J4,J1,J2,rank) &
                                *(-1)** ((ja+jd+J1+J2)/2) * phase_34 * LCC%herm * Xelem

                           call add_elem_to_tensor(V,b,a,c,d,J1,J2,RES,jbas)                        

                        end do
                     end do

                     do J1 = J2min,J2max,2
                        if ((c==d).and.(mod(J1/2,2)==1)) cycle
                        prefac_134 = prefac_34 * sqrt((J1+1.d0))
                        do J2 = max(J1,J1min,abs(J1-rank)),min(J1max,J1+rank),2 
                           prefac_1234 = prefac_134 * sqrt((J2+1.d0))

                           !V^{J1 J2}_{cdba}

                           V = prefac_1234 * ninej(jd,ja,J3,jc,jb,J4,J1,J2,rank) &
                                *(-1)**((jc-jb+J1+J2+rank)/2) * RCC%herm * Xelem
                           ! the mapping here inverts the indeces to {dcab} so you need 
                           ! an additional factor of phase(a+b+c+d+J1+J2)
                           call add_elem_to_tensor(V,c,d,b,a,J1,J2,RES,jbas)               

                        end do
                     end do

                  else

                     ! CALCULATE V^{J1 J2}_{badc} and V^{J1 J2}_{dcba}
                     do J1 = J1min,J1max,2
                        prefac_134 = prefac_34 * sqrt((J1+1.d0))
                        do J2 = max(J1,J2min,abs(J1-rank)),min(J2max,J1+rank),2 
                           prefac_1234 = prefac_134 * sqrt((J2+1.d0))

                           !V^{J1 J2}_{badc}
                           V = prefac_1234 * ninej(ja,jd,J3,jb,jc,J4,J1,J2,Rank) &
                                * (-1) ** ((ja+jc+J1)/2) *phase_34 * LCC%herm * Xelem 

                           call add_elem_to_tensor(V,b,a,d,c,J1,J2,RES,jbas)               

                        end do
                     end do

                     do J1 = J2min,J2max,2
                        prefac_134 = prefac_34 * sqrt((J1+1.d0))
                        do J2 = max(J1,J1min,abs(J1-rank)),min(J1max,J1+rank),2 
                           prefac_1234 = prefac_134 * sqrt((J2+1.d0))

                           !V^{J1 J2}_{dcba}
                           V = prefac_1234 * ninej(jd,ja,J3,jc,jb,J4,J1,J2,Rank) &
                                * (-1)**((jb-jd+J2 + rank )/2) *RCC%herm * Xelem 

                           call add_elem_to_tensor(V,d,c,b,a,J1,J2,RES,jbas)               

                        end do
                     end do
                  end if
               end if
            end if
            
            
            J1min = abs(jd-jb)
            J1max = jd+jb
            
            J2min = abs(jc-ja) 
            J2max = jc+ja 
            
            if ( abs(Yelem) > 1e-6) then 
               if (b .ge. d ) then                
                  if ( a .ge. c ) then 

                     do J1 = J1min,J1max,2
                        if ((d==b).and.(mod(J1/2,2)==1)) cycle
                        prefac_134 = prefac_34 * sqrt((J1+1.d0))
                        do J2 = max(J1,J2min,abs(rank-J1)),min(J2max,rank+J1),2 
                           if ((c==a).and.(mod(J2/2,2)==1)) cycle
                           prefac_1234 = prefac_134 * sqrt((J2+1.d0))

                           !V^{J1 J2}_{dbca}
                           V = prefac_1234* ninej(jd,ja,J3,jb,jc,J4,J1,J2,rank)&
                                * (-1) ** ((jc+ja+J2)/2) *LCC%herm * Yelem 
                           call add_elem_to_tensor(V,d,b,c,a,J1,J2,RES,jbas)               

                        end do
                     end do

                     !V^{J1 J2}_{cadb}
                     do J1 = J2min,J2max,2
                        if ((c==a).and.(mod(J1/2,2)==1)) cycle
                        prefac_134 = prefac_34 * sqrt((J1+1.d0))
                        do J2 = max(J1,J1min,abs(J1-rank)),min(J1max,J1+rank),2 
                           if ((b==d).and.(mod(J2/2,2)==1)) cycle
                           prefac_1234 = prefac_134 * sqrt((J2+1.d0))

                           V = prefac_1234* ninej(ja,jd,J3,jc,jb,J4,J1,J2,rank)&
                                * (-1)** ((jd-jb+J1+rank)/2) * phase_34 *RCC%herm * Yelem
                           call add_elem_to_tensor(V,c,a,d,b,J1,J2,RES,jbas)
                        end do
                     end do

                  else
                     do J1 = J1min,J1max,2
                        if ((d==b).and.(mod(J1/2,2)==1)) cycle
                        prefac_134 = prefac_34 * sqrt((J1+1.d0))
                        do J2 = max(J1,J2min,abs(rank-J1)),min(J2max,rank+J1),2 
                           if ((c==a).and.(mod(J2/2,2)==1)) cycle
                           prefac_1234 = prefac_134 * sqrt((J2+1.d0))

                           !V^{J1 J2}_{dbca}
                           V = -1*prefac_1234* ninej(jd,ja,J3,jb,jc,J4,J1,J2,rank)&
                                * LCC%herm * Yelem 
                           call add_elem_to_tensor(V,d,b,a,c,J1,J2,RES,jbas)               

                        end do
                     end do

                     do J1 = J2min,J2max,2
                        if ((a==c).and.(mod(J1/2,2)==1)) cycle
                        prefac_134 = prefac_34 * sqrt((J1+1.d0))
                        do J2 = max(J1,J1min,abs(J1-rank)),min(J1max,J1+rank),2 
                           if ((b==d).and.(mod(J2/2,2)==1)) cycle
                           prefac_1234 = prefac_134 * sqrt((J2+1.d0))

                           !V^{J1 J2}_{dbca}
                           V = prefac_1234* ninej(ja,jd,J3,jc,jb,J4,J1,J2,rank)&
                                *phase_abcd*phase_34*(-1)**(rank/2)* RCC%herm * Yelem
                           call add_elem_to_tensor(V,a,c,d,b,J1,J2,RES,jbas)               

                        end do
                     end do
                     !do nothing
                  end if
               else 
                  if ( a .ge. c ) then 

                     do J1 = J1min,J1max,2
                        if ((b==d).and.(mod(J1/2,2)==1)) cycle
                        prefac_134 = prefac_34 * sqrt((J1+1.d0))
                        do J2 = max(J1,J2min,abs(rank-J1)),min(J2max,rank+J1),2 
                           if ((a==c).and.(mod(J2/2,2)==1)) cycle
                           prefac_1234 = prefac_134 * sqrt((J2+1.d0))

                           V = -1*prefac_1234  * ninej(jd,ja,J3,jb,jc,J4,J1,J2,rank) &
                                * phase_abcd*(-1)**((J1+J2)/2) * LCC%herm * Yelem 

                           call add_elem_to_tensor(V,b,d,c,a,J1,J2,RES,jbas)               
                        end do
                     end do

                     do J1 = J2min,J2max,2
                        if ((c==a).and.(mod(J1/2,2)==1)) cycle
                        prefac_134 = prefac_34 * sqrt((J1+1.d0))
                        do J2 = max(J1,J1min,abs(J1-rank)),min(J1max,J1+rank),2 
                           if ((b==d).and.(mod(J2/2,2)==1)) cycle
                           prefac_1234 = prefac_134 * sqrt((J2+1.d0))

                           V = prefac_1234  * ninej(ja,jd,J3,jc,jb,J4,J1,J2,rank) &
                                * (-1)**((J1+J2+rank)/2) * phase_34 * RCC%herm * Yelem

                           call add_elem_to_tensor(V,c,a,b,d,J1,J2,RES,jbas)               
                        end do
                     end do

                  else
                     do J1 = J1min,J1max,2
                        if ((b==d).and.(mod(J1/2,2)==1)) cycle
                        prefac_134 = prefac_34 * sqrt((J1+1.d0))
                        do J2 = max(J1,J2min,abs(rank-J1)),min(J2max,rank+J1),2 
                           if ((a==c).and.(mod(J2/2,2)==1)) cycle
                           prefac_1234 = prefac_134 * sqrt((J2+1.d0))

                           !V^{J1 J2}_{bdac}
                           V = prefac_1234 * ninej(jd,ja,J3,jb,jc,J4,J1,J2,rank) &
                                * (-1)**((jb+jd+J1)/2) * LCC%herm *Yelem
                           call add_elem_to_tensor(V,b,d,a,c,J1,J2,RES,jbas)               
                        end do
                     end do

                     do J1 = J2min,J2max,2
                        if ((a==c).and.(mod(J1/2,2)==1)) cycle
                        prefac_134 = prefac_34 * sqrt((J1+1.d0))
                        do J2 = max(J1,J1min,abs(J1-rank)),min(J1max,J1+rank),2 
                           if ((d==b).and.(mod(J2/2,2)==1)) cycle
                           prefac_1234 = prefac_134 * sqrt((J2+1.d0))

                           !V^{J1 J2}_{acbd}
                           V = prefac_1234 * ninej(ja,jd,J3,jc,jb,J4,J1,J2,rank) &
                                * (-1)**((ja-jc+J2+rank)/2) * phase_34 * RCC%herm *Yelem

                           call add_elem_to_tensor(V,a,c,b,d,J1,J2,RES,jbas)               
                        end do
                     end do


                  end if
               end if
            end if
         end do
      end do
!!!!      !$OMP END PARALLEL DO
      deallocate(Wx,Wy)
   end do
   
   !if ( (d == 11).and.(c==18).and.(b==12).and.(a==17).and.(J1==0).and.(J2==4)) then 
   !   !if (abs(V)>1e-6) print*, J3,J4,V, '8'  
   !end if
                    
 end subroutine tensor_product_222_ph
!=================================================================
!=================================================================
real(8) function tensor_product_223_single(L,R,ip,iq,ir,is,it,iu,jtot1,jtot2,Jpq,Jst,jbas)
  implicit none 
  
  integer,intent(in) :: ip,iq,ir,is,it,iu,jpq,jst
  integer :: a,b,c,d,Jx,Jy,Jz,J2,J1,J3,phase,rank,ia
  integer :: tp,tq,tr,ts,tt,tu,lp,lq,lr,ls,lt,lu,astart
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

  tp = jbas%itzp(ip)
  tq = jbas%itzp(iq)
  tr = jbas%itzp(ir)  
  ts = jbas%itzp(is)
  tt = jbas%itzp(it)
  tu = jbas%itzp(iu)       

  lp = jbas%ll(ip)
  lq = jbas%ll(iq)
  lr = jbas%ll(ir)  
  ls = jbas%ll(is)
  lt = jbas%ll(it)
  lu = jbas%ll(iu)  

  ! FIRST TERM 
  !changed to q-r instead of q+r
  multfact = (-1)**((jq+jr+jp-jtot2+rank)/2) *sqrt((Jpq+1.d0)*(Jst+1.d0) &
       *(jtot1+1.d0)*(jtot2+1.d0) )  
  ! ju isn't in here because I need it to make add with ja
  ! so I get an integer later 

  ! Left-Right 
  sm = 0.d0   
  
  j2min = max(abs(jq - jr),abs(jp-jtot1))
  j2max = min(jq+jr,jp+jtot1)
      
  do a=1,jbas%total_orbits 
     

     if ( (tp+jbas%itzp(a)) .ne. (ts+tt) ) cycle
     if ( mod(lp+jbas%ll(a),2) .ne. mod(ls+lt,2) ) cycle 
     
     ja = jbas%jj(a) 
     
     if (.not. triangle(jp,ja,jst) ) cycle
                    
     phase = (-1)**((ja-ju)/2) 
     

     sj1 = v_elem(ip,a,is,it,Jst,L,jbas)*phase

     do J2 = j2min, j2max , 2
        
        sj2 = sj1*sixj(jp,jq,Jpq,jr,jtot1,J2)*sqrt(J2+1.d0) 
        
        j3min = min(abs(ja - ju),abs(jp-jtot2),abs(J2-rank))
        j3max = max(ja+ju,jp+jtot2,J2+rank) 
     
        do J3 = j3min,j3max,2

           sm = sm - (-1)**(J3/2) * sqrt(J3 + 1.d0) &            
            * sj2 * sixj(jp,ja,Jst,ju,jtot2,J3) *  &
             xxxsixj(J2,J3,rank,jtot2,jtot1,jp) * &
             tensor_elem(iq,ir,a,iu,J2,J3,R,jbas)
           
        end do
     end do
  end do
 
  smtot = smtot + sm*multfact
  
  !Right-left  

  multfact = (-1)**((jq+jr+jp-jtot2+rank)/2) *sqrt((Jpq+1.d0) &
       *(jtot1+1.d0)*(jtot2+1.d0) ) 
  sm = 0.d0 
  do a=1,jbas%total_orbits 
     
     
     if ( (tu+jbas%itzp(a)) .ne. (tq+tr) ) cycle
     if ( mod(lu+jbas%ll(a),2) .ne. mod(lq+lr,2) ) cycle 
     
     ja = jbas%jj(a) 
          
     j1min = min(abs(jp - ja),abs(Jst-rank))
     j1max = max(jp + ja,Jst+rank) 
     
     j3min = max( abs(jq - jr) , abs(ja - ju) , abs(jp-jtot1)) 
     j3max = min( jq+jr , ja+ju , jp+jtot1) 
     
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
  
  
  smtot = smtot + sm*multfact 
  
  !SECOND TERM
  !Left-Right
  sm = 0.d0 
  
  multfact = (-1)**((jp+jq+jr+jt+ju+jtot2+rank)/2) *sqrt((Jpq+1.d0)*(Jst+1.d0)*&
       (jtot1+1.d0)*(jtot2+1.d0)) 
  ! added a minus sign
  do a=1,jbas%total_orbits
     
     
     if ( (tp+jbas%itzp(a)) .ne. (tt+tu) ) cycle
     if ( mod(lp+jbas%ll(a),2) .ne. mod(lt+lu,2) ) cycle 
     
     
     ja = jbas%jj(a)

     j2min = max(abs(jp - ja) , abs(jt - ju) ,abs(js-jtot2)) 
     j2max = min(jp+ja,jt+ju,js+jtot2) 
     
     j1min = max(abs(ja - js),abs(jp-jtot2)) 
     j1max = min(ja+js,jp+jtot2) 
     
     phase = (-1) ** ((ja - js)/2) 
     
     do J1 = j1min,j1max,2
        
        sj1 = sqrt(J1+1.d0)*phase  
       
        do J2 = j2min,j2max,2
           
           sj2 = sj1 * (J2+1.d0) * (-1)**(J2/2) * sixj(js,jt,Jst,ju,jtot2,J2) * &
                sixj(jp,ja,J2,js,jtot2,J1) * v_elem(ip,a,it,iu,J2,L,jbas) 
     
           j3min = min(abs(jq - jr),abs(J1-rank)) 
           j3max = max(jq+jr,J1+rank) 
     
           do J3 = j3min,j3max,2
           
              sm = sm - sqrt(J3+1.d0) * sj2 * sixj(jp,jq,Jpq,jr,jtot1,J3) &
                   * xxxsixj(J1,J3,rank,jtot1,jtot2,jp) * (-1)**(J1/2) * &
                     tensor_elem(iq,ir,a,is,J3,J1,R,jbas)
           end do
        end do         
        
     end do
 
  end do


  do a=1,jbas%total_orbits
     
     
     if ( (ts+jbas%itzp(a)) .ne. (tq+tr) ) cycle
     if ( mod(ls+jbas%ll(a),2) .ne. mod(lq+lr,2) ) cycle 
     
     ja = jbas%jj(a)
     
     j1min = min(abs(jp - ja),abs(js-jtot1))
     j1max = max(jp+ja ,js+jtot1) 
     
     j3min = max( abs(jq - jr) , abs(ja - js) ,abs(jp-jtot1)) 
     j3max = min( jq+jr,ja+js,jp+jtot1)
     
     phase = (-1) ** ((ja - jp)/2) 
     
     do J1 = j1min,j1max,2
        
        sj1 = sqrt(J1+1.d0) * phase *(-1)**(J1/2)
        
        j2min = min(abs(jt - ju),abs(J1-rank)) 
        j2max = max(jt+ju ,J1+rank) 
        
        do J2 = j2min,j2max,2 
           
           sj2 = sj1 * sqrt(J2+1.d0) * sixj(js,jt,Jst,ju,jtot2,J2) * (-1)**(J2/2) *&
                xxxsixj(J2,J1,rank,jtot1,jtot2,js) * tensor_elem(ip,a,it,iu,J1,J2,R,jbas) 
           
           do J3 = j3min,j3max,2
              
              sm = sm + sj2 * (J3+1.d0) * sixj(jp,jq,Jpq,jr,jtot1,J3) &
                   * sixj(jp,ja,J1,js,jtot1,J3) * v_elem(iq,ir,a,is,J3,L,jbas)
           end do 
        end do 
     end do
 
  end do
     
  smtot = smtot + sm*multfact

  ! THIRD TERM    
  !Left-Right
  sm = 0.d0 
  
  multfact = (-1)**((jq+jr+jtot2+js+Jst+rank)/2) *sqrt((Jpq+1.d0)*(Jst+1.d0)*&
       (jtot1+1.d0)*(jtot2+1.d0)) 

  do a=1,jbas%total_orbits
     
          
     if ( (tp+jbas%itzp(a)) .ne. (ts+tu) ) cycle
     if ( mod(lp+jbas%ll(a),2) .ne. mod(ls+lu,2) ) cycle 
     
     ja = jbas%jj(a)     
     
     j2min = max( abs(jp - ja) , abs(js - ju),abs(jt-jtot2) ) 
     j2max = min( jp+ja , js+ju,jt+jtot2) 
     
     j3min = max(abs(jq - jr) ,abs(jp-jtot1))
     j3max = min(jq+jr,jp+jtot1) 
   
     
     phase = (-1) ** ((ja + jp)/2) ! changed to ja+js rather than ja-js 
     
     do J3 = j3min,j3max,2
        
        sj1 = phase*sixj(jp,jq,Jpq,jr,jtot1,J3)*sqrt(J3+1.d0)  
        
        j1min = max(abs(ja - jt),abs(rank-J3))
        j1max = min(ja+jt,rank+J3) 
             
        do J1 = j1min,j1max,2 
           
           sj2 =  sj1*(-1)**(J1/2)*tensor_elem(iq,ir,a,it,J3,J1,R,jbas)*&
                xxxsixj(J1,J3,rank,jtot1,jtot2,jp)*sqrt(J1+1.d0)

           do J2 = j2min,j2max,2
              sm = sm - (J2+1.d0) *sj2* sixj(js,jt,Jst,jtot2,ju,J2) &
                   * sixj(jp,ja,J2,jt,jtot2,J1) * v_elem(ip,a,iu,is,J2,L,jbas)
           end do
        end do
     end do
 
  end do

  ! right-left

  do a=1,jbas%total_orbits
     
     
     if ( (tt+jbas%itzp(a)) .ne. (tq+tr) ) cycle
     if ( mod(lt+jbas%ll(a),2) .ne. mod(lq+lr,2) ) cycle 

     ja = jbas%jj(a)     
     
     j2min =  max(abs(js - ju),abs(jt-jtot2))
     j2max =  min(js+ju,jt+jtot2)
     
     j3min = max( abs(jq - jr) , abs(ja - jt) ,abs(jp-jtot1)) 
     j3max = min( jq+jr , ja+jt,jp+jtot1)
     
     phase = (-1) ** ((ja + jt)/2) ! changed to ja+js rather than ja-js 
     
     do J2 = j2min,j2max,2
        
        sj1 = phase* sixj(js,jt,Jst,jtot2,ju,J2)*sqrt(J2+1.d0) 
    
        j1min = max(abs(jp - ja),abs(rank-J2))
        j1max = min(jp+ja ,rank+J2)
    
        do J1 = j1min,j1max,2 
           sj2 = sj1*(-1)**(J1/2)*xxxsixj(J1,J2,rank,jtot2,jtot1,jt) *sqrt(J1+1.d0) *&
                tensor_elem(ip,a,iu,is,J1,J2,R,jbas)      

           do J3 = j3min,j3max,2
              sm = sm + (J3+1.d0) *sj2* sixj(jp,jq,Jpq,jr,jtot1,J3) &
                   * sixj(jp,ja,J1,jt,jtot1,J3) * v_elem(iq,ir,a,it,J3,L,jbas)
           end do
        
        end do
        
     end do
 
  end do
     
  smtot = smtot + sm*multfact



!  FOURTH TERM    
!  Left-Right
  sm = 0.d0 
  
  multfact = (-1)**((jp+jtot2+Jpq+rank)/2) *sqrt((Jpq+1.d0)*(Jst+1.d0)*&
       (jtot1+1.d0)*(jtot2+1.d0)) 

  do a=1,jbas%total_orbits
     
     
     if ( (tq+jbas%itzp(a)) .ne. (ts+tt) ) cycle
     if ( mod(lq+jbas%ll(a),2) .ne. mod(ls+lt,2) ) cycle 

     ja = jbas%jj(a)     
     
     
     j1min = max(abs(jp - jr),abs(jtot1-jq))
     j1max = min(jp+jr ,jtot1+jq) 
        
     phase = (-1) ** ((ja - ju)/2) ! changed to ja+js rather than ja-js 
     
     sj1 = v_elem(iq,a,is,it,Jst,L,jbas)*phase
     
     do J1 = j1min,j1max,2
        
        sj2 = sj1*sixj(jp,jq,Jpq,jtot1,jr,J1)*sqrt(J1+1.d0)*(-1)**(J1/2)
        
        j2min = max(abs(ja - ju),abs(rank-J1),abs(jq-jtot2))
        j2max = min(ja+ju,rank+J1,jq+jtot2)
        
        do J2 = j2min,j2max,2 
           
           sm = sm -  sj2*sqrt(J2+1.d0)*(-1)**(J2/2)*tensor_elem(ir,ip,a,iu,J1,J2,R,jbas)*&
                xxxsixj(J1,J2,rank,jtot2,jtot1,jq)*sixj(ja,jq,Jst,jtot2,ju,J2)

        end do
     end do

  end do

  smtot = smtot + sm*multfact

  ! right-left
      
  sm = 0.d0 
  
  multfact = (-1)**((jp-jtot2+Jpq+rank)/2) *sqrt((Jpq+1.d0)*&
       (jtot1+1.d0)*(jtot2+1.d0)) 

  do a=1,jbas%total_orbits
     

     if ( (tu+jbas%itzp(a)) .ne. (tr+tp) ) cycle
     if ( mod(lu+jbas%ll(a),2) .ne. mod(lr+lp,2) ) cycle 

     ja = jbas%jj(a)     
          
     j1min = max(abs(jq - ja),abs(rank-Jst))
     j1max = min(jq+ja,rank+Jst) 
   
     j2min = max(abs(ja - ju),abs(jr-jp),abs(jq-jtot1)) 
     j2max = min(ja+ju,jr+jp,jq+jtot1) 
     
     phase = (-1) ** ((ja + jq)/2) ! changed to ja+js rather than ja-js 
          
     do J1 = j1min,j1max,2
        
        sj1 = phase*sqrt(J1+1.d0)* &
             tensor_elem(iq,a,is,it,J1,Jst,R,jbas)*(-1)**(J1/2)
        
        do J2 = j2min,j2max,2 
           
           sm = sm +  sj1*(J2+1.d0)*(-1)**(J2/2)*v_elem(ir,ip,a,iu,J2,L,jbas) *&
                xxxsixj(J1,Jst,rank,jtot2,jtot1,ju)*sixj(ja,jq,J1,jtot1,ju,J2) &
                *sixj(jp,jq,Jpq,jtot1,jr,J2)
        end do
     end do

  end do

  smtot = smtot + sm*multfact
 
  ! FIFTH TERM    
  !Left-Right
  sm = 0.d0 

  multfact = (-1)**((ju+jt+jtot2+js+Jpq+rank)/2) *sqrt((Jpq+1.d0)*(Jst+1.d0)*&
       (jtot1+1.d0)*(jtot2+1.d0)) 

  do a=1,jbas%total_orbits
     
     
     
     if ( (tq+jbas%itzp(a)) .ne. (tt+tu) ) cycle
     if ( mod(lq+jbas%ll(a),2) .ne. mod(lt+lu,2) ) cycle 

     ja = jbas%jj(a)     
     
     j1min = max( abs(jq - ja) , abs(jt - ju) ) 
     j1max = min( jq+ja , jt+ju) 
     
     j3min = max(abs(ja - js),abs(jq-jtot2))
     j3max = min(ja+js,jq+jtot2) 
     
     phase = (-1) ** ((ja + jp)/2) ! changed to ja+js rather than ja-js 
     
     do J1 = j1min,j1max,2
        
        sj1 = phase*sixj(js,jt,Jst,ju,jtot2,J1)*(J1+1.d0)*(-1)**(J1/2)&
             *v_elem(iq,a,it,iu,J1,L,jbas) 
        
        do J3 = j3min,j3max,2 
           
           sj2 =  sj1*(-1)**(J3/2)*&
                sixj(jq,ja,J1,js,jtot2,J3)*sqrt(J3+1.d0)

           j2min = max(abs(jp - jr),abs(rank-J3),abs(jq-jtot1))
           j2max = min(jp+jr,rank+J3,jq+jtot1) 
   
           do J2 = j2min,j2max,2
              sm = sm - sqrt(J2+1.d0)*(-1)**(J2/2)*sj2* sixj(jp,jq,Jpq,jtot1,jr,J2) &
                   * xxxsixj(J3,J2,rank,jtot1,jtot2,jq) * tensor_elem(ir,ip,a,is,J2,J3,R,jbas)
           end do
        end do
     end do
 
  end do
  smtot = smtot + sm*multfact

  sm = 0.d0 
  ! right-left
  multfact = (-1)**((jp+jq+jtot2+ju+Jpq+rank)/2) *sqrt((Jpq+1.d0)*(Jst+1.d0)*&
       (jtot1+1.d0)*(jtot2+1.d0)) 

  do a=1,jbas%total_orbits
     
     if ( (ts+jbas%itzp(a)) .ne. (tr+tp) ) cycle
     if ( mod(ls+jbas%ll(a),2) .ne. mod(lr+lp,2) ) cycle 

     ja = jbas%jj(a)     
    
     
     j2min =  max(abs(jt - ju),abs(js-jtot2))
     j2max =  max(jt+ju,js+jtot2)
     
     j3min = max( abs(jp - jr) , abs(ja - js) ) 
     j3max = min( jp+jr , ja+js)
     
     phase = (-1) ** ((ja + jt)/2) ! changed to ja+js rather than ja-js 
     
     do J3 = j3min,j3max,2
        
        sj1 = phase*(-1)**(J3/2)* sixj(jp,jq,Jpq,jtot1,jr,J3)*(J3+1.d0)*&
             v_elem(ir,ip,a,is,J3,L,jbas) 
        
        do J2 = j2min,j2max,2 
           sj2 = sj1*(-1)**(J2/2)*sixj(js,jt,Jst,ju,jtot2,J2)*sqrt(J2+1.d0)

           j1min = max(abs(jq - ja),abs(rank-J2),abs(js-jtot1)) 
           j1max = min(jq+ja,rank+J2,js+jtot1)
                    
           do J1 = j1min,j1max,2
              sm = sm + sqrt(J1+1.d0) *(-1)**(J1/2) *sj2* sixj(jq,ja,J1,js,jtot1,J3) &
                   * xxxsixj(J2,J1,rank,jtot1,jtot2,js) * tensor_elem(iq,a,it,iu,J1,J2,R,jbas)
           end do
        
        end do
        
     end do
 
  end do
     
  smtot = smtot + sm*multfact 

  ! SIXTH TERM    
  !Left-Right
  sm = 0.d0 
  
  multfact = (-1)**((jtot2+js+Jpq+Jst+rank)/2) *sqrt((Jpq+1.d0)*(Jst+1.d0)*&
       (jtot1+1.d0)*(jtot2+1.d0)) 

  do a=1,jbas%total_orbits
        
     
          
     if ( (tq+jbas%itzp(a)) .ne. (tu+ts) ) cycle
     if ( mod(lq+jbas%ll(a),2) .ne. mod(lu+ls,2) ) cycle 

     ja = jbas%jj(a)     
     
     j1min = max( abs(jq - ja) , abs(js - ju) ,abs(jt-jtot2) ) 
     j1max = min( jq+ja , js+ju , jt+jtot2 ) 
     
     j2min = max(abs(jp - jr),abs(jq-jtot1))
     j2max = min(jp+jr ,jq+jtot1) 
       
     phase = (-1) ** ((ja - jp)/2) ! changed to ja+js rather than ja-js 
     
     do J1 = j1min,j1max,2
        
        sj1 = phase*sixj(js,jt,Jst,jtot2,ju,J1)*(J1+1.d0) &
             *v_elem(iq,a,iu,is,J1,L,jbas) 
        
        do J2 = j2min,j2max,2 
           
           sj2 =  sj1*(-1)**(J2/2)*&
                sixj(jp,jq,Jpq,jtot1,jr,J2)*sqrt(J2+1.d0)
           
           j3min = max(abs(ja - jt),abs(rank-J2)) 
           j3max = min(ja+jt,rank+J2)
     
           do J3 = j3min,j3max,2
              sm = sm - sqrt(J3+1.d0)*(-1)**(J3/2)*sj2* sixj(jq,ja,J1,jt,jtot2,J3) &
                   * xxxsixj(J3,J2,rank,jtot1,jtot2,jq) * tensor_elem(ir,ip,a,it,J2,J3,R,jbas)
           end do
        end do
     end do
 
  end do
  
  smtot = smtot + sm*multfact

  sm = 0.d0 
  ! right-left
  multfact = (-1)**((jp+jq+js+jtot2+Jpq+Jst+rank)/2) *sqrt((Jpq+1.d0)*(Jst+1.d0)*&
       (jtot1+1.d0)*(jtot2+1.d0)) 

  do a=1,jbas%total_orbits
     
     if ( (tt+jbas%itzp(a)) .ne. (tr+tp) ) cycle
     if ( mod(lt+jbas%ll(a),2) .ne. mod(lr+lp,2) ) cycle 

     ja = jbas%jj(a)     
    
     j2min = max(abs(js - ju),abs(jt-jtot2)) 
     j2max = min(js+ju ,jt+jtot2) 
     
     j3min = max( abs(jp - jr) , abs(ja - jt) ) 
     j3max = min( jp+jr , ja+jt)
     
     phase = (-1) ** ((ja + jt)/2) ! changed to ja+js rather than ja-js 
     
     do J3 = j3min,j3max,2
        
        sj1 = phase*(-1)**(J3/2) * sixj(jp,jq,Jpq,jtot1,jr,J3)*(J3+1.d0)*&
             v_elem(ir,ip,a,it,J3,L,jbas) 
        
        do J2 = j2min,j2max,2 

           sj2 = sj1*sixj(js,jt,Jst,jtot2,ju,J2)*sqrt(J2+1.d0)
                    
           j1min = max(abs(jq - ja),abs(rank-J2),abs(jt-jtot1))
           j1max = min(jq+ja,rank+J2,jt+jtot1)
           
           do J1 = j1min,j1max,2
              sm = sm + sqrt(J1+1.d0) *(-1)**(J1/2) *sj2* sixj(jq,ja,J1,jt,jtot1,J3) &
                   * xxxsixj(J2,J1,rank,jtot1,jtot2,jt) * tensor_elem(iq,a,iu,is,J1,J2,R,jbas)
           end do
        
        end do
        
     end do
 
  end do
  smtot = smtot + sm*multfact


  ! SEVENTH TERM
  
  ! Left-right
  sm = 0.d0
  multfact = (-1)**((Jpq+rank+jr-jtot2)/2) *sqrt((Jst+1.d0)*(jtot1+1.d0)*&
       (jtot2+1.d0)) 

  do a=1,jbas%total_orbits

     
     
     if ( (tr+jbas%itzp(a)) .ne. (ts+tt) ) cycle
     if ( mod(lr+jbas%ll(a),2) .ne. mod(ls+lt,2) ) cycle 
     
     ja = jbas%jj(a)

     if (.not. triangle(jr,ja,jst) ) cycle
     
     j3min = max(abs(ja-ju),abs(Jpq-rank),abs(jr-jtot2))
     j3max = min(ja+ju,Jpq+rank,jr+jtot2)  
     
     sj1 = v_elem(ir,a,is,it,Jst,L,jbas) * (-1)**((ja+ju)/2) 

     do J3=j3min,j3max,2 
        sm = sm - sj1* sixj(jr,ja,Jst,ju,jtot2,J3) * (-1)**(J3/2) * sqrt(J3+1.d0) &
             * xxxsixj(J3,Jpq,rank,jtot1,jtot2,jr) * tensor_elem(ip,iq,a,iu,Jpq,J3,R,jbas)
     end do 

  end do 
  
  smtot = smtot + sm*multfact

  ! right-left
  sm = 0.d0
  multfact = (-1)**((Jpq+rank)/2) *sqrt((Jpq+1.d0)*(jtot1+1.d0)*&
       (jtot2+1.d0)) 

  do a=1,jbas%total_orbits
     
     if ( (tu+jbas%itzp(a)) .ne. (tp+tq) ) cycle
     if ( mod(lu+jbas%ll(a),2) .ne. mod(lp+lq,2) ) cycle 

     ja = jbas%jj(a)

     if (.not. triangle(ju,ja,Jpq) ) cycle
     
     j3min = max(abs(ja-jr),abs(rank-Jst),abs(ju-jtot1))
     j3max = min(ja+jr,rank+Jst,ju+jtot1)  
     
     sj1 = v_elem(ip,iq,a,iu,Jpq,L,jbas) * (-1)**((ja+jtot2)/2) 

     do J3=j3min,j3max,2 

        sm = sm + sj1* sixj(jr,ja,J3,ju,jtot1,Jpq) * sqrt(J3+1.d0) *(-1)**(J3/2) &
             * xxxsixj(J3,Jst,rank,jtot2,jtot1,ju) * tensor_elem(ir,a,is,it,J3,Jst,R,jbas)

     end do 

  end do 
  
  smtot = smtot + sm*multfact


  ! EIGHTH TERM 
  !changed to q-r instead of q+r
  multfact = (-1)**((jt+jr+js+jtot2+Jpq+rank)/2) *sqrt((Jst+1.d0) &
       *(jtot1+1.d0)*(jtot2+1.d0) )  
  ! ju isn't in here because I need it to make add with ja
  ! so I get an integer later 

  ! Left-Right 
  sm = 0.d0   
  
  do a=1,jbas%total_orbits 

     
          
     if ( (tr+jbas%itzp(a)) .ne. (tt+tu) ) cycle
     if ( mod(lr+jbas%ll(a),2) .ne. mod(lt+lu,2) ) cycle 

     ja = jbas%jj(a) 
     
     j2min = max(abs(ja - js),abs(rank-Jpq),abs(jtot2-jr))
     j2max = min(ja+js ,rank+Jpq,jtot2+jr)
          
     j1min = max(abs(ja - jr),abs(jt-ju)) 
     j1max = min(ja+jr,jt+ju)  
          
     phase = (-1)**((ja+ju)/2) 
     
     do J1 = j1min, j1max , 2
        sj1 = phase*(-1)**(J1/2)*(J1+1.d0) * sixj(js,jt,Jst,ju,jtot2,J1)&
             * v_elem(ir,a,it,iu,J1,L,jbas) 
     
        do J2 = j2min,j2max,2
        
           sm = sm - (-1)**(J2/2) * sqrt(J2 + 1.d0) &            
            * sj1 * sixj(ja,js,J2,jtot2,jr,J1) *  &
             xxxsixj(Jpq,J2,rank,jtot2,jtot1,jr) * &
             tensor_elem(ip,iq,a,is,Jpq,J2,R,jbas)
        end do
     end do
  end do
 
  smtot = smtot + sm*multfact
  
  ! !Right-left  
  sm = 0.d0 
  multfact = (-1)**((jt+jtot2+Jpq+rank)/2) *sqrt((Jst+1.d0)*(Jpq+1.d0) &
       *(jtot1+1.d0)*(jtot2+1.d0) )  
  
  do a=1,jbas%total_orbits 
     
     if ( (ts+jbas%itzp(a)) .ne. (tp+tq) ) cycle
     if ( mod(ls+jbas%ll(a),2) .ne. mod(lp+lq,2) ) cycle 

     ja = jbas%jj(a) 
          
     if (.not. triangle(ja,js,Jpq)) cycle 
     
     j1min = max(abs(jr - ja),abs(js-jtot1))
     j1max = min(jr + ja  ,js+jtot1) 
     
     
     phase = (-1)**((ja - ju)/2)
           
     sj1 = phase * v_elem(ip,iq,a,is,Jpq,L,jbas)
     do J1 = j1min, j1max , 2      

        sj2 = sqrt(J1+1.d0)*(-1)**(J1/2)*sixj(jr,ja,J1,js,jtot1,Jpq) * sj1

        j2min = max(abs(jt - ju),abs(J1-rank),abs(js-jtot2))
        j2max = min(jt+ju,J1+rank,js+jtot2)
     
        do J2 = j2min,j2max,2

           sm = sm + sj2* sqrt(J2 + 1.d0) * (-1)**(J2/2) &
                * sixj(js,jt,Jst,ju,jtot2,J2) * xxxsixj(J1,J2,rank,jtot2,jtot1,js) &
                * tensor_elem(ir,a,it,iu,J1,J2,R,jbas)

        end do
     end do
  end do 
  
  
  smtot = smtot + sm*multfact
  
! NINTH TERM 
  !changed to q-r instead of q+r
  multfact = (-1)**((jr+jtot2+Jpq+Jst+rank)/2) *sqrt((Jst+1.d0) &
       *(jtot1+1.d0)*(jtot2+1.d0) )  
  ! ju isn't in here because I need it to make add with ja
  ! so I get an integer later 

  ! Left-Right 
  sm = 0.d0   
  
  do a=1,jbas%total_orbits 
     
     
     if ( (tr+jbas%itzp(a)) .ne. (tu+ts) ) cycle
     if ( mod(lr+jbas%ll(a),2) .ne. mod(lu+ls,2) ) cycle 

     ja = jbas%jj(a) 
     
     j2min = max(abs(ja - jt),abs(rank-Jpq))
     j2max = min(ja+jt,rank+Jpq) 
          
     j1min = max(abs(ja - jr),abs(js-ju),abs(jt-jtot2)) 
     j1max = min(ja+jr,js+ju,jt+jtot2)  
          
     phase = (-1)**((ja-js)/2) 
     
     do J1 = j1min, j1max , 2
        sj1 = phase*(J1+1.d0) * sixj(js,jt,Jst,jtot2,ju,J1)&
             * v_elem(ir,a,iu,is,J1,L,jbas) 
        
        do J2 = j2min,j2max,2
        
           sm = sm - (-1)**(J2/2) * sqrt(J2 + 1.d0) &            
            * sj1 * sixj(ja,jr,J1,jtot2,jt,J2) *  &
             xxxsixj(Jpq,J2,rank,jtot2,jtot1,jr) * &
             tensor_elem(ip,iq,a,it,Jpq,J2,R,jbas)
        end do
        
     end do
  end do
 
  smtot = smtot + sm*multfact
  
  ! !Right-left  
  sm = 0.d0 
  multfact = (-1)**((jt-jtot2+Jpq+Jst+rank)/2) *sqrt((Jst+1.d0)*(Jpq+1.d0) &
       *(jtot1+1.d0)*(jtot2+1.d0) )  
  
  do a=1,jbas%total_orbits 
     
     if ( (tt+jbas%itzp(a)) .ne. (tp+tq) ) cycle
     if ( mod(lt+jbas%ll(a),2) .ne. mod(lp+lq,2) ) cycle 

     ja = jbas%jj(a) 
          
     if (.not. triangle(ja,jt,Jpq)) cycle 
     
     j1min = max(abs(jr - ja),abs(jt-jtot1))
     j1max = max(jr + ja,jt+jtot1) 
     
     
     phase = (-1)**((ja + js)/2)
           
     sj1 = phase * v_elem(ip,iq,a,it,Jpq,L,jbas)
     do J1 = j1min, j1max , 2      

        sj2 = sqrt(J1+1.d0)*(-1)**(J1/2)*sixj(jr,ja,J1,jt,jtot1,Jpq) * sj1

        j2min = max(abs(js - ju),abs(rank-J1),abs(jt-jtot2))
        j2max = min(js+ju,rank+J1,jt+jtot2)
     
        do J2 = j2min,j2max,2

           sm = sm + sj2* sqrt(J2 + 1.d0) &
                * sixj(js,jt,Jst,jtot2,ju,J2) * xxxsixj(J1,J2,rank,jtot2,jtot1,jt) &
                * tensor_elem(ir,a,iu,is,J1,J2,R,jbas)

        end do
     end do
  end do 
  
  
  smtot = smtot + sm*multfact
    

     
  tensor_product_223_single = smtot

end function tensor_product_223_single
!=====================================================================================
!=====================================================================================
end module 
  
  
  
