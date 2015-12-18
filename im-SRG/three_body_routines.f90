module three_body_routines 
  use basic_IMSRG
  
  
  contains


real(8) function overlap_3b(a1,b1,c1,Jab1,Tab1,a2,b2,c2,Jab2,Tab2,jtot,ttot,jbas) 
  ! INCLUDES ANTI-SYMMETRY FACTOR
  implicit none 
  
  type(spd) :: jbas
  integer :: a1,b1,c1,a2,b2,c2
  integer :: Jab1,Jab2,Tab1,Tab2,jtot,ttot
  integer :: ja,jb,jc
  logical :: zero
 
  ! check right away that this thing can exist.
  zero = .true. 
  if ( (a1 == a2) .or. (a1 == b2) .or. (a1==c2) ) then 
     if ( (b1 == a2) .or. (b1 == b2) .or. (b1==c2) ) then  
        if ( (c1 == a2) .or. (c1 == b2) .or. (c1==c2) ) then 
           zero = .false.
        end if
     end if
  end if
  
  if ( zero ) then 
     overlap_3b = 0.d0 
     return
  end if
          
  if ( (a1 == b1).neqv.(a1==c1) ) then 
        
     ! we've got two like indeces. 
     
     if ( ( a1 == b1) .and. (a2 == b2) ) then 

        !  < [(aac) Jaa Taa] jt |  [(aac) Jaa Taa] jt >

        If (Jab1 .ne. Jab2) then 
           overlap_3b = 0.d0 
           return
        end if

        If (Tab1 .ne. Tab2) then 
           overlap_3b = 0.d0 
           return
        end if

        overlap_3b = (1.d0 - (-1)**((Jab1 +Tab1)/2) )/2.d0   !!! normalization? 
        return
     end if

     if ( ( a1 == b1 ) .and. ( a2 == c2 ) ) then 

        !  < [(aac) Jaa Taa] jt |  [(aca) Jaa Taa] jt >           

        ja = jbas%jj(a1)
        jc = jbas%jj(c1) 

        overlap_3b = ((-1) **((Jab1+Tab1)/2)-1) &
             * sqrt( (Jab1+1.d0)*(Jab2+1.d0)*(Tab1+1.d0)*(Tab2+1.d0) ) &
             * sixj(ja,ja,Jab1,jc,jtot,Jab2) * sixj(1,1,Tab1,1,ttot,Tab2)&
             * (-1) ** ((ja+jc - Jab2 - Tab2)/2) / 2.d0 ! normalization?  
        return
     end if

     if ( ( a1 == c1 ) .and. ( a2 == b2 ) ) then 
        !  < [(aca) Jac Tac] jt |  [(aac) Jaa Taa] jt >           

        ja = jbas%jj(a1)
        jc = jbas%jj(b1) 

        overlap_3b = ((-1) **((Jab1+Tab1)/2)-1) &
             * sqrt( (Jab1+1.d0)*(Jab2+1.d0)*(Tab1+1.d0)*(Tab2+1.d0) ) &
             * sixj(ja,ja,Jab1,jc,jtot,Jab2) * sixj(1,1,Tab1,1,ttot,Tab2)&
             * (-1) ** ((ja+jc - Jab2 - Tab2)/2)/2.d0!normalization 
        return
     end if

     if ( ( a1 == b1 ) .and. ( b2 == c2 ) ) then 

        !  < [(aac) Jaa Taa] jt |  [(caa) Jaa Taa] jt >                            
        ja = jbas%jj(a1)
        jc = jbas%jj(c1) 

        overlap_3b = ((-1) **((Jab1+Tab1)/2)-1) &
             * sqrt( (Jab1+1.d0)*(Jab2+1.d0)*(Tab1+1.d0)*(Tab2+1.d0) ) &
             * sixj(ja,ja,Jab1,jc,jtot,Jab2) * sixj(1,1,Tab1,1,ttot,Tab2)/2.d0!normalization

        return
     end if


     if ( ( b1 == c1 ) .and. ( a2 == b2 ) ) then 

        !  < [(caa) Jca Tca] jt |  [(aac) Jaa Taa] jt >                            
        ja = jbas%jj(b1)
        jc = jbas%jj(a1) 

        overlap_3b = ((-1) **((Jab1+Tab1)/2)-1) &
             * sqrt( (Jab1+1.d0)*(Jab2+1.d0)*(Tab1+1.d0)*(Tab2+1.d0) ) &
             * sixj(ja,ja,Jab1,jc,jtot,Jab2) * sixj(1,1,Tab1,1,ttot,Tab2)/2.d0!normalization

        return
     end if

     
     if ( ( a1 == c1 ) .and. ( a2 == c2 ) ) then 

        !  < [(aca) Jac Tac] jt |  [(aca) Jac Tac] jt >                            
        ja = jbas%jj(a1)
        jc = jbas%jj(b1) 

        
        overlap_3b = kron_del(Jab1,Jab2)*kron_del(Tab1,Tab2) - &
              sqrt( (Jab1+1.d0)*(Jab2+1.d0)*(Tab1+1.d0)*(Tab2+1.d0) ) &
             * sixj(ja,jc,Jab1,ja,jtot,Jab2) * sixj(1,1,Tab1,1,ttot,Tab2)! not sure if this needs it

        return
     end if
     
     if ( ( a1 == c1 ) .and. ( b2 == c2 ) ) then 

        !  < [(aca) Jac Tac] jt |  [(caa) Jac Tac] jt >                            
        ja = jbas%jj(a1)
        jc = jbas%jj(b1) 

        overlap_3b = (kron_del(Jab1,Jab2)*kron_del(Tab1,Tab2) - &
              sqrt( (Jab1+1.d0)*(Jab2+1.d0)*(Tab1+1.d0)*(Tab2+1.d0) ) &
             * sixj(ja,jc,Jab1,ja,jtot,Jab2) * sixj(1,1,Tab1,1,ttot,Tab2)) &
             * ( -1) ** ((ja + jc -Jab2 - Tab2)/2 )
        
        return
     end if
                 
     if ( ( b1 == c1 ) .and. ( a2 == c2 ) ) then 

        !  < [(caa) Jac Tac] jt |  [(aca) Jac Tac] jt >                            
        ja = jbas%jj(b1)
        jc = jbas%jj(a1) 

        
        overlap_3b = (kron_del(Jab1,Jab2)*kron_del(Tab1,Tab2) - &
              sqrt( (Jab1+1.d0)*(Jab2+1.d0)*(Tab1+1.d0)*(Tab2+1.d0) ) &
             * sixj(ja,jc,Jab1,ja,jtot,Jab2) * sixj(1,1,Tab1,1,ttot,Tab2)) &
             * ( -1) ** ((ja + jc -Jab1 - Tab1)/2 )
        
        return
     end if   
     
     if ( ( b1 == c1 ) .and. ( b2 == c2 ) ) then 

        !  < [(caa) Jac Tac] jt |  [(caa) Jac Tac] jt >                            
        ja = jbas%jj(b1)
        jc = jbas%jj(a1) 

        
        overlap_3b = (kron_del(Jab1,Jab2)*kron_del(Tab1,Tab2) - &
              sqrt( (Jab1+1.d0)*(Jab2+1.d0)*(Tab1+1.d0)*(Tab2+1.d0) ) &
             * sixj(ja,jc,Jab1,ja,jtot,Jab2) * sixj(1,1,Tab1,1,ttot,Tab2)) &
             * ( -1) ** ((Jab1+Tab1+Jab2+Tab2)/2 )
        
        return
     end if
           
  end if


  if ( (a1 == b1).and.(a1==c1) ) then 
     ! < (aaa) Jaa Taa jt |   (aaa) Jaa Taa jt > 


     ja= jbas%jj(a1) 
     
     overlap_3b = (kron_del(Jab1,Jab2)*kron_del(Tab1,Tab2) - &
              (sqrt( (Jab1+1.d0)*(Jab2+1.d0)*(Tab1+1.d0)*(Tab2+1.d0) ) &
             * sixj(ja,jc,Jab1,ja,jtot,Jab2) * sixj(1,1,Tab1,1,ttot,Tab2) &
             * (1-( -1) ** ((Jab2+Tab2)/2 )))) * (1- (-1)**((Jab1+Tab1)/2)) /6.d0 ! normalization? 
     return
  
  end if

! otherwise there are no duplicate indeces. 
  if ( a1 == a2 ) then 
     
     if ( b1 == b2) then 
               
        if (c1 == c2) then 
           
           if ( Jab1 == Jab2 ) then 
              if (Tab1 == Tab2) then 
                 overlap_3b = 1.d0 
              else 
                 overlap_3b = 0.d0 
              end if
           else
              overlap_3b=0.d0 
           end if 
           
        else
           
           overlap_3b = 0.d0
        
        end if 
  
     else if (b1 == c2 )  then 
        
        if (c1 == b2 ) then 
           
           ja=jbas%jj(a1) 
           jb=jbas%jj(b1)
           jc=jbas%jj(c1) 
           
           overlap_3b = (-1) **(( jc + jb + Jab1 + Jab2 +Tab1 + Tab2)/2 ) &
                * sqrt( (Jab1+1.d0)*(Jab2+1.d0)*(Tab1+1.d0)*(Tab2+1.d0) ) &
                * sixj(jb,ja,Jab1,jc,jtot,Jab2) * sixj(1,1,Tab1,1,ttot,Tab2) 
           
        else 
           
           overlap_3b = 0.d0 
        end if 
     else 
        overlap_3b = 0.d0 
     end if 
  else if (a1 == b2) then 
     
     if (b1 == a2 ) then 
        
        if ( c1 == c2)  then 
           
           if (Jab1 == Jab2) then
              if (Tab1 == Tab2) then 
                 
                 ja = jbas%jj(a1) 
                 jb = jbas%jj(b1) 
           
                 overlap_3b =  (-1)**((ja+jb-Jab1-Tab1)/2) 
              else 
                 overlap_3b = 0.d0 
              end if 
           else 
              overlap_3b = 0.d0 
           end if 
        else
           overlap_3b = 0.d0 
        end if 
        
     else if (b1 == c2) then 
        if (c1 == a2 ) then 
           
           ja=jbas%jj(a1) 
           jb=jbas%jj(b1)
           jc=jbas%jj(c1) 
           
           overlap_3b = (-1) **(( jc + jb + Jab1 + Jab2 +Tab1 + Tab2)/2 ) &
                * sqrt( (Jab1+1.d0)*(Jab2+1.d0)*(Tab1+1.d0)*(Tab2+1.d0) ) &
                * sixj(jb,ja,Jab1,jc,jtot,Jab2) * sixj(1,1,Tab1,1,ttot,Tab2) &
                * (-1) ** (( ja + jc  - Jab2 -Tab2 )/2) 
           
        else 
           overlap_3b = 0.d0 
        end if 
     else 
        overlap_3b = 0.d0 
     end if 
     
  else if (a1 == c2)  then 
     
     if (b1 == b2) then 
        
        if ( c1 == a2) then
           
           ja=jbas%jj(a1) 
           jb=jbas%jj(b1)
           jc=jbas%jj(c1) 
           
           overlap_3b = -1*sqrt( (Jab1+1.d0)*(Jab2+1.d0)*(Tab1+1.d0)*(Tab2+1.d0) ) &
                * sixj(ja,jb,Jab1,jc,jtot,Jab2) * sixj(1,1,Tab1,1,ttot,Tab2)
        else
           overlap_3b = 0.d0
        end if 
        
     else if (b1 == a2) then 
        if ( c1 == b2) then 
           
           ja=jbas%jj(a1) 
           jb=jbas%jj(b1)
           jc=jbas%jj(c1) 
           
           overlap_3b = -1*sqrt( (Jab1+1.d0)*(Jab2+1.d0)*(Tab1+1.d0)*(Tab2+1.d0) ) &
                * sixj(ja,jb,Jab1,jc,jtot,Jab2) * sixj(1,1,Tab1,1,ttot,Tab2) &
                * (-1)**((jb +jc-Jab2-Tab2)/2)
           
        else 
           overlap_3b = 0.d0 
        end if 
     end if 
  else 
     overlap_3b = 0.d0
  end if 

  
end function
!====================================================================
!====================================================================
subroutine allocate_three_body_storage(jbas,store_3b)
  implicit none 
  
  type(spd) :: jbas
  type(three_body_force) :: store_3b
  integer :: jtot, Tz,PAR,Jij,Jlm,i,j,k,l,m,n,nsp_iso,Tij
  integer :: ji,jj,jl,jk,jm,jn,a,b,c,d,e,f,qx,Tab_indx,TTab_indx
  integer :: jtot_max,Nsp,num_3b,num_blocks,x1,x2
  integer :: j_half,j_half_min,j_half_max,NN
  integer :: Jij_min,Jij_max,q,q1Tij,ttot
  integer :: ti,tj,tk,li,lj,lk,num_included 
  integer :: nlj1,nlj2,nlj3,nnlj1,nnlj2,nnlj3,aux,aux1,aux2,aux3,aux4
  integer :: nnlj2_end,nnlj3_end,twoTMin,twoTmax,twoJCMin,twoJCMax
  integer :: twoJCMindown,twoJCMaxup,twoJCMindownket,twoJCMindownbra
  integer :: twoJCMaxupket,twoJCMaxupbra,la,lb,lc,ld,le,lf
  integer :: ja,jb,jc,jd,je,jf,iblock,Jab,JJab,Tab,TTab,elems
  integer :: ea,eb,ec,ed,ef,ee,e1max,E3max,JabMax,JabMin,JJabMax,JJabMin
  real(8) :: mem ,ass

  jtot_max = 3*(jbas%jtotal_max)
  Nsp = jbas%total_orbits
  nsp_iso = Nsp/2  ! isospin coupled
  num_3b = (Nsp_iso+Nsp_iso**2)/2 + (Nsp_iso**3 - Nsp_iso)/6
  num_blocks = ((jtot_max-1)/2+1)*4
 
  allocate(store_3b%mat(num_blocks),store_3b%hashmap(num_3b))
  allocate(store_3b%lam(num_blocks,3),store_3b%kets(num_blocks))
  allocate(store_3b%Nsize(num_blocks), store_3b%bras(num_blocks)) 
  
  store_3b%kets = 0 
  store_3b%bras = 0 
  
  l = 1 
  do i = 1,Nsp_iso
     ji = jbas%jj(2*i) ! isospin coupled so only looking at Tz=1 states in jbas
     do j= 1,i
        jj = jbas%jj(2*j)
        do k= 1,j 
           jk = jbas%jj(2*k)
           
           if ( l .ne. threebody_index(i,j,k) ) print*, 'fuck' 
           
           Jij_min = abs(ji-jj)
           Jij_max = ji+jj
           aux = (Jij_max-Jij_min)/2+1
           
           store_3b%hashmap(l)%Jij_start= Jij_min
           
           allocate(store_3b%hashmap(l)%position(aux,2))  ! the 2 are the Tab projections 
           allocate(store_3b%hashmap(l)%jhalf_start(aux)) 
           allocate(store_3b%hashmap(l)%halfsize(aux)) 
           
           do Jij = Jij_min , Jij_max, 2
              
              aux = (Jij-Jij_min)/2+1
              
              j_half_min = abs(Jij-jk)
              j_half_max = Jij+jk 
              
              store_3b%hashmap(l)%jhalf_start(aux) = j_half_min
              
              aux2 = (j_half_max-j_half_min)/2+1 
              
              store_3b%hashmap(l)%halfsize(aux) = aux2 
              do Tij = 1,2! not actual Tab value
                 if (Tij ==2) aux2 = aux2 * 2 ! two possible thalf
                 allocate( store_3b%hashmap(l)%position(aux,Tij)%Y(aux2,3) ) 
                 store_3b%hashmap(l)%position(aux,Tij)%Y = 0
              end do
              
           end do
       
           l=l+1
        end do
     end do
  end do
  
  E3max = 12
  
  elems = 0
  mem = 0.d0 

  do nlj1 = 1, nsp_iso 
     la = jbas%ll(2*nlj1)
     ja = jbas%jj(2*nlj1)
     ea = 2*jbas%nn(2*nlj1)+la 
     if (ea > E3Max) exit
     
     do nlj2 = 1,nlj1
        lb = jbas%ll(2*nlj2)
        jb = jbas%jj(2*nlj2)
        eb = 2*jbas%nn(2*nlj2)+lb 
        if (ea + eb > E3Max) exit
        
        do nlj3 = 1,nlj2
           lc = jbas%ll(2*nlj3)
           jc = jbas%jj(2*nlj3)
           ec = 2*jbas%nn(2*nlj3)+lc 
           if (ea + eb + ec > E3Max) exit
           
     do nnlj1 = 1, nlj1 
        ld = jbas%ll(2*nnlj1)
        jd = jbas%jj(2*nnlj1)
        ed = 2*jbas%nn(2*nnlj1)+ld 
        if (ed > E3Max) exit
        
        if ( nlj1 == nnlj1 ) then 
           nnlj2_end = nlj2
        else
           nnlj2_end = nnlj1
        end if
                  
        do nnlj2 = 1, nnlj2_end            
           le = jbas%ll(2*nnlj2)
           je = jbas%jj(2*nnlj2)
           ee = 2*jbas%nn(2*nnlj2)+le 
           if (ed + ee > E3Max) exit
        
           if ( (nlj1 == nnlj1).and.(nlj2==nnlj2)) then 
              nnlj3_end = nlj3
           else
              nnlj3_end = nnlj2
           end if 
           
           do nnlj3 = 1,nnlj3_end
              lf = jbas%ll(2*nnlj3)
              jf = jbas%jj(2*nnlj3)
              ef = 2*jbas%nn(2*nnlj3)+lf 
              if (ed + ee + ef > E3Max) exit

              
              ! check parity 
              PAR = mod(la+lb+lc,2)
              if ( mod(ld+le+lf,2) .ne. PAR ) cycle
              
              JabMax =  ja+jb 
              JabMin = abs(ja-jb) 
              
              JJabMax = jd+je
              JJabMin = abs(jd-je) 
              
              !determine roughly the two*J bounds
              
              if (abs(ja-jb) > jc) then 
                 twoJCMindownbra = abs(ja-jb)-jc
              else if ( jc < ja+jb) then 
                 twoJCMindownbra = 1
              else
                 twoJCMindownbra = jc - ja -jb
              end if 
              
              if (abs(jd-je) > jf) then 
                 twoJCMindownket = abs(jd-je)-jf
              else if ( jf < jd+je) then 
                 twoJCMindownket = 1
              else
                 twoJCMindownket = jf - jd -je
              end if 
              
              twoJCMaxupbra = ja+jb+jc
              twoJCMaxupket = jd+je+jf
              
              twoJCMindown = max(twoJCMindownket,twoJCMindownbra) 
              twoJCMaxup = min(twoJCMaxupket,twoJCMaxupbra)
              
              if (twoJCMindown > twoJCMaxup) cycle
              
              do Jab = JabMin,JabMax,2
                 do JJab=JJabMin,JJabMax,2
                    
                    twoJCMin = max(abs(Jab-jc),abs(JJab-jf))
                    twoJCMax = min(Jab+jc,JJab+jf) 
                    
                    do jtot = twoJCMin,twoJCMax,2
                       
                       do Tab = 0,1
                          do TTab = 0,1
                             
                             twoTMin = max(abs(2*Tab-1),abs(2*TTab-1))
                             twoTMax = min(2*Tab+1,2*TTab+1) 
                             
                             do ttot = twoTMin,twoTmax,2
                               
                                ! what block are we in?   
                                q = block_index_3b(jtot,ttot,PAR) 
                                    
                                x1=threebody_index(nlj1,nlj2,nlj3)
                                x2=threebody_index(nnlj1,nnlj2,nnlj3)
                              
                                aux1 = (Jab-store_3b%hashmap(x1)%Jij_start)/2 + 1  
                                aux2 = (JJab-store_3b%hashmap(x2)%Jij_start)/2 + 1  
                                aux3 = (jtot - store_3b%hashmap(x1)%jhalf_start(aux1))/2+1 & 
                                     + store_3b%hashmap(x1)%halfsize(aux1)*(ttot-1)/2   
                                aux4 = (jtot - store_3b%hashmap(x2)%jhalf_start(aux2))/2+1 & 
                                     + store_3b%hashmap(x2)%halfsize(aux2)*(ttot-1)/2
                                
                                Tab_indx = (2*Tab+2)/2
                                TTab_indx = (2*TTab+2)/2

                                ! ONLY ADD INFORMATION IF IT HASN'T BEEN ADDED YET.
                                if (store_3b%hashmap(x1)%position(aux1,Tab_indx)%Y(aux3,2) == 0 ) then  
                                   store_3b%bras(q) = store_3b%bras(q) + 1 
                                   store_3b%hashmap(x1)%position(aux1,Tab_indx)%Y(aux3,1)= q
                                   store_3b%hashmap(x1)%position(aux1,Tab_indx)%Y(aux3,2)=store_3b%bras(q)     
                                end if    

                                if (store_3b%hashmap(x2)%position(aux2,TTab_indx)%Y(aux4,3) == 0 ) then
                                   store_3b%kets(q) = store_3b%kets(q) + 1  
                                   store_3b%hashmap(x2)%position(aux2,TTab_indx)%Y(aux4,1) = q
                                   store_3b%hashmap(x2)%position(aux2,TTab_indx)%Y(aux4,3) = store_3b%kets(q) 
                                end if 
                                                     
                                
                                elems = elems + 1                               
                                
                                
                             end do !ttot
                          end do !TTab
                       end do ! Tab
                    end do !jtot
                 end do !JJab
              end do !Jab
           end do !nnlj3
        end do !nnlj2
     end do !nnlj1 
        end do !nlj3
     end do !nlj2
  end do !nlj1
   
  
  print*, 'number of matrix elements: ',elems
  store_3b%num_elems = elems
  do q = 1, num_blocks     
     NN = store_3b%kets(q) 
     allocate(store_3b%mat(q)%XX((NN*NN+NN)/2))
     store_3b%Nsize(q) = (NN*NN+NN)/2 
     mem = mem +sizeof(store_3b%mat(q)%XX)
  end do
  
  mem = mem/1024.d0/1024.d0/1024.d0
  print*, 'MEMORY OF 3 BODY STORAGE IS: ',mem,'GB' 
end subroutine allocate_three_body_storage

end module
        
           
           
     
     
        
        
  
  
  
  
  



