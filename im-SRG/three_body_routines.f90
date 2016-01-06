module three_body_routines 
  use basic_IMSRG
  
type mono_3b
   type(real_vec),dimension(8) :: mat
   integer,dimension(8,2) :: lam
   integer,allocatable,dimension(:,:) :: hash
   integer,dimension(8) :: dm
end type mono_3b  
 
contains

subroutine get_mono_3b_elems(STOR,monoSTOR,jbas) 
   implicit none 
  
   type(three_body_force) :: STOR
   type(mono_3b) :: monoSTOR
   type(spd) :: jbas
   real(8) :: mem,sm
   integer :: ta,tb,tc,la,lb,lc,N,x1,x2
   integer :: ja,jb,jc,jaa,jbb,jcc,II,JJ,jtot,Jab
   integer :: Jab_min,Jab_max,j_min,j_max,aux
   integer :: q,Tz,PAR,a,b,c,aa,bb,cc,items
   ! i'm too lazy right now. 
   ! monopoles --- 
   
   N = size(jbas%con) 
   allocate(monoSTOR%hash(N*N*N,2))
   mem = 0.d0 
   q = 1
   do Tz=-3,3,2
      do PAR=0,1
         
         monoSTOR%lam(q,1) = Tz 
         monoSTOR%lam(q,2) = PAR 

         items = 0 
         do a=1,N
            ta = jbas%itzp(a) 
            la = jbas%ll(a)
            do b=1,N
               tb = jbas%itzp(a) 
               lb = jbas%ll(a)
               do c=1,N
                  tc = jbas%itzp(a) 
                  lc = jbas%ll(a)
   
                  if (ta+tb+tc .ne. Tz) cycle
                  if (mod(la+lb+lc,2) .ne.PAR) cycle
                  
                  !!! member of this block 
                  items=items + 1
                  monoSTOR%hash(N*N*(a-1)+N*(b-1)+c,1) = q
                  monoSTOR%hash(N*N*(a-1)+N*(b-1)+c,2) = items
               end do 
            end do 
         end do 
      
         allocate(monoSTOR%mat(q)%XX(items*(items+1)/2)) 
         monoSTOR%dm(q) = items
         mem = mem +sizeof(monoSTOR%mat(q)%XX)
         q=q+1
      end do 
   end do 
   mem = mem+sizeof(monoSTOR%hash)
   mem = mem/1024./1024./1024.
   print*, 'monopole storage: ',mem,'GB'

 !$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) SHARED(STOR,monoSTOR,jbas)
   do x1 = 1, N**3 
      c = mod(x1-1,N)+1
      a = (x1-1)/(N**2)+1
      b = (x1-c-N*N*(a-1))/N+1 

      q   = monoSTOR%hash(x1,1)
      II  = monoSTOR%hash(x1,2)
      
      Tz  = monoSTOR%lam(q,1)      
      PAR = monoSTOR%lam(q,2)
      
      ja = jbas%jj(a) 
      jb = jbas%jj(b)
      jc = jbas%jj(c) 
      print*, x1,N**3
      do x2 = x1,N**3 
         
         if (monoSTOR%hash(x2,1) .ne. q) cycle
         JJ  = monoSTOR%hash(x2,2)
         
         cc = mod(x2-1,N)+1          
         aa = (x2-1)/(N**2)+1
         bb = (x2-cc-N*N*(aa-1))/N+1
         
         jaa = jbas%jj(aa) 
         jbb = jbas%jj(bb)
         jcc = jbas%jj(cc) 
      
         aux = bosonic_tp_index(II,JJ,monoSTOR%dm(q)) 
         
         sm = 0.d0 
         Jab_min = max(abs(ja-jb),abs(jaa-jbb))
         Jab_max = min(ja+jb,jaa+jbb) 
         
         do Jab = Jab_min,Jab_max,2
            
            j_min = max(abs(Jab-jc),abs(Jab-jcc))
            j_max = min(Jab+jc,Jab+jcc) 
            
            do jtot = j_min,j_max,2
               
               sm = sm + (jtot+1.d0)&
                    *GetME_pn(Jab,Jab,jtot,a,b,c,aa,bb,cc,STOR,jbas)
            end do
         end do 
         
         monoSTOR%mat(q)%XX(aux) = sm 
      end do
   end do
!$OMP END PARALLEL DO
end subroutine           
                  
   
  

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
  integer :: Jij_min,Jij_max,q,q1Tij,ttot,mores
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
  store_3b%E3max = 12
  
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
  elems = 0
  do q = 1, num_blocks     
     if (store_3b%kets(q) .ne. store_3b%bras(q)) print*, 'shit'
     NN = store_3b%kets(q) 
     elems = elems + (NN*NN+NN)/2
     allocate(store_3b%mat(q)%XX((NN*NN+NN)/2))
     store_3b%Nsize(q) = NN
     mem = mem +sizeof(store_3b%mat(q)%XX)
  end do
  
  mem = mem/1024.d0/1024.d0/1024.d0
  print*, 'MEMORY OF 3 BODY STORAGE IS: ',mem,'GB' 
  print*, 'slots:', elems
end subroutine allocate_three_body_storage

 
real(8) function GetME_pn(Jab,Jde,jtot,a,b,c,d,e,f,STOR,jbas) 
  implicit none 
  
  type(spd) :: jbas
  type(three_body_force) :: STOR
  integer :: Jab,Jde,jtot,a,b,c,d,e,f
  integer :: Tab,Tde,ttot,Tmin,Tmax
  integer :: ta,tb,tc,td,te,tf
  real(8) :: sm,CG1,CG2,CG3,CG4,d6ji
  
  ta = jbas%itzp(a)
  tb = jbas%itzp(b)  
  tc = jbas%itzp(c)
  td = jbas%itzp(d)  
  te = jbas%itzp(e)
  tf = jbas%itzp(f)
  
  if ((ta+tb+tc-td-te-tf).ne.0) then 
     GetME_pn = 0.d0 
     return
  end if 
  
  sm = 0.d0  
  Tmin = min(abs(ta+tb+tc),abs(td+te+tf)) 
  do Tab = abs(ta+tb),2,2 
     CG1 = iso_clebsch_halfhalf(1,1,Tab,ta,tb,ta+tb)        
     do Tde = abs(td+te),2,2        
        CG2 = iso_clebsch_halfhalf(1,1,Tde,td,te,td+te)        
        do ttot = Tmin,3,2
           
           CG3 = iso_clebsch_wholehalf(Tab,1,ttot,ta+tb,tc,ta+tb+tc)
           CG4 = iso_clebsch_wholehalf(Tde,1,ttot,td+te,tf,td+te+tf)
           if (CG3*CG4 == 0) cycle
           sm = sm + CG1*CG2*CG3*CG4 *&
                GetME(Jab,Jde,jtot,Tab,Tde,ttot,a,b,c,d,e,f,STOR,jbas) 
        end do
     end do
  end do
  
  GetME_pn = sm 
end function GetME_pn
  
subroutine SetME(Jab_in,Jde_in,jtot,Tab_in,Tde_in,ttot,a_in,b_in,c_in,d_in,e_in,f_in,V_in,STOR,jbas)
  ! here we accept pn basis indeces, but treat them as iso-spin coupled. 
  ! e.g. 1==2, 3==4, 5==6 ...   
  implicit none 
  
  integer :: Jab_in,Jde_in,jtot,Tab_in,Tde_in,ttot,a_in,b_in,c_in,d_in,e_in,f_in
  real(8) :: V_in,x
  type(three_body_force) :: STOR
  type(spd) :: jbas
  
  x = AddToME(Jab_in,Jde_in,jtot,Tab_in,Tde_in,ttot,a_in,b_in,c_in,d_in,e_in,f_in,V_in,STOR,jbas)
  !x is not important here
end subroutine SetME

real(8) function GetME(Jab_in,Jde_in,jtot,Tab_in,Tde_in,ttot,a_in,b_in,c_in,d_in,e_in,f_in,STOR,jbas) 
  ! here we accept pn basis indeces, but treat them as iso-spin coupled. 
  implicit none 
  
  integer,intent(in) :: Jab_in,Jde_in,jtot,Tab_in,Tde_in,ttot,a_in,b_in,c_in,d_in,e_in,f_in
  integer :: Jab_i,Jde_i,jt,Tab_i,Tde_i,tt,a_i,b_i,c_i,d_i,e_i,f_i
  real(8) :: x
  type(three_body_force) :: STOR
  type(spd) :: jbas
  
  Jab_i = Jab_in; Jde_i = Jde_in; jt = jtot;Tab_i = Tab_in; Tde_i = Tde_in; tt = ttot
  a_i=a_in;  b_i=b_in;  c_i=c_in;  d_i=d_in;  e_i=e_in;  f_i=f_in

  GetME = AddToME(Jab_i,Jde_i,jt,Tab_i,Tde_i,tt,a_i,b_i,c_i,d_i,e_i,f_i,0.d0,STOR,jbas) 
end function GetME
  
!!! These functions are ported from Ragnar Stroberg's code: 
real(8) function &
AddToME(Jab_in,Jde_in,jtot,Tab_in,Tde_in,ttot,a_in,b_in,c_in,d_in,e_in,f_in,V_in,STOR,jbas) 
  implicit none 

  type(spd) :: jbas 
  type(three_body_force) :: STOR
  integer :: a,b,c,d,e,f,abc_recoup,def_recoup,jmin
  integer :: Tde_min,ttot,jmax
  integer :: ja,jb,jc,jd,je,jf,x1,x2,q,II,JJ,Jindx
  integer :: aux1,aux2,aux3,aux4,Tab_indx,TTab_indx
  integer :: Jab,Jde,Tab,Tde,aux,Jab_max,Jab_min
  integer :: Jde_max,Jde_min,Tab_max,Tab_min,Tde_max
  real(8) :: V_out,Cj_abc,Cj_def,Ct_abc,Ct_def
 
  integer,intent(in) :: a_in,b_in,c_in,d_in,e_in,f_in
  integer,intent(in) :: Jab_in,Jde_in,jtot,Tab_in,Tde_in  
  real(8),intent(in) :: V_in
  
  abc_recoup = SortOrbits(a_in,b_in,c_in,a,b,c)
  def_recoup = SortOrbits(d_in,e_in,f_in,d,e,f)
    
  if ((d>a).or.((d==a).and.(e>b)).or.((d==a).and.(e==b).and.(f>c))) then 
     call swap(a,d)
     call swap(b,e)
     call swap(c,f)
     call swap(Jab_in,Jde_in)
     call swap(Tab_in,Tde_in)
     call swap(abc_recoup,def_recoup)
  end if 
  
  ! make sure the element isn't above the E3max truncation
  if ( (2* ( jbas%nn(a)+jbas%nn(b) + jbas%nn(c)) + &
       jbas%ll(a) + jbas%ll(b) +jbas%ll(c) ) > STOR%E3max ) then 
     AddToME = 0.d0 
     return
  end if 

  if ( (2* ( jbas%nn(d)+jbas%nn(e) + jbas%nn(f)) + &
       jbas%ll(d) + jbas%ll(e) +jbas%ll(f) ) > STOR%E3max ) then 
     AddToME = 0.d0 
     return
  end if      
  
  ja = jbas%jj(a) 
  jb = jbas%jj(b) 
  jc = jbas%jj(c) 
  jd = jbas%jj(d) 
  je = jbas%jj(e) 
  jf = jbas%jj(f)  
  
  Jab_min = abs(ja-jb)
  Jab_max = ja+jb
  Jde_min = abs(jd-je)
  Jde_max = jd+je 
  
  if (ttot == 3) then 
     Tab_min = 2
     Tde_min = 2
  else
     Tab_min = 0
     Tde_min = 0
  end if 
  Tab_max = 2
  Tde_max = 2
  
  !!! Find indeces:::  
  x1=threebody_index((a-1)/2+1,(b-1)/2+1,(c-1)/2+1)
  x2=threebody_index((d-1)/2+1,(e-1)/2+1,(f-1)/2+1) 
    

  V_out = 0.d0 
  
   if (abc_recoup == 0) then 
      Jab_min = Jab_in;Jab_max = Jab_in
      Tab_min = Tab_in;Tab_max = Tab_in
   end if 
  
  if (def_recoup == 0) then 
     Jde_min = Jde_in;Jde_max = Jde_in
     Tde_min = Tde_in;Tde_max = Tde_in
  end if 
  
  do Jab = Jab_min,Jab_max,2
     Cj_abc = Recoupling_Coef(abc_recoup,ja,jb,jc,Jab_in,Jab,jtot) 
     if (abc_recoup>2) Cj_abc = Cj_abc*(-1) 

     do Jde = Jde_min,Jde_max,2
        Cj_def = Recoupling_Coef(def_recoup,jd,je,jf,Jde_in,Jde,jtot) 
        if (def_recoup>2) Cj_def = Cj_def*(-1) 
        
        jmin = max( abs( Jab-jc) , abs(Jde-jf) ) 
        jmax = min( Jab+jc , Jde+jf) 

        if ( jmin > jmax) cycle
        
        if ( (jtot .ge. jmin).and.(jtot.le.jmax) ) then

           do Tab = Tab_min,Tab_max,2
              Ct_abc = Recoupling_Coef(abc_recoup,1,1,1,Tab_in,Tab,ttot) 
              do Tde= Tde_min,Tde_max,2
                 Ct_def = Recoupling_Coef(def_recoup,1,1,1,Tde_in,Tde,ttot)          
                
                 aux1 = (Jab-STOR%hashmap(x1)%Jij_start)/2 + 1  
                 aux2 = (Jde-STOR%hashmap(x2)%Jij_start)/2 + 1  
                 aux3 = (jtot - STOR%hashmap(x1)%jhalf_start(aux1))/2+1 & 
                      + STOR%hashmap(x1)%halfsize(aux1)*(ttot-1)/2   
                 aux4 = (jtot - STOR%hashmap(x2)%jhalf_start(aux2))/2+1 & 
                      + STOR%hashmap(x2)%halfsize(aux2)*(ttot-1)/2

                 Tab_indx = (Tab+2)/2
                 TTab_indx = (Tde+2)/2

                 q=STOR%hashmap(x1)%position(aux1,Tab_indx)%Y(aux3,1)
                 II=STOR%hashmap(x1)%position(aux1,Tab_indx)%Y(aux3,2)
                 JJ=STOR%hashmap(x2)%position(aux2,TTab_indx)%Y(aux4,3)
                 
                 if ( II >= JJ) then 
                    aux = bosonic_tp_index(JJ,II,STOR%Nsize(q))
                    STOR%mat(q)%XX(aux) = STOR%mat(q)%XX(aux) + &
                         Cj_abc * Cj_def * Ct_abc * Ct_def * V_in

                 else 
                    
                    aux = bosonic_tp_index(II,JJ,STOR%Nsize(q))                    
                 end if
                  
                 V_out = V_out + Cj_abc * Cj_def * Ct_abc * Ct_def &
                      *STOR%mat(q)%XX(aux)
         
              end do
           end do
        end if
     end do
  end do
 ! end if 
  AddToME = V_out 

end function 
!=====================================================================================
!=====================================================================================
integer function SortOrbits(a_in,b_in,c_in,a,b,c) 
  ! order so that a>=b>=c
  implicit none 
  
  integer :: a_in,b_in,c_in,a,b,c
  integer :: recoupling_case
  integer :: ax,bx,cx
  
  ! this puts us in isospin basis
  a_in = a_in - mod(a_in-1,2)
  b_in = b_in - mod(b_in-1,2) 
  c_in = c_in - mod(c_in-1,2)
  
  a=a_in
  b=b_in
  c=c_in 

  if (a<b) call swap(a,b)
  if (b<c) call swap(b,c)
  if (a<b) call swap(a,b) 

  ! any even perm is less than 3
  if (a_in == a) then 
     if ( b_in==b) then 
        recoupling_case = 0
     else
        recoupling_case = 3
     end if 
  else if (a_in==b) then 
     if (b_in == a) then 
        recoupling_case = 4
     else 
        recoupling_case = 1
     end if 
  else 
     if (b_in==a) then 
        recoupling_case = 2
     else
        recoupling_case = 5
     end if 
  end if 

  SortOrbits = recoupling_case 
  
end function
!=====================================================================================
!=====================================================================================
real(8) function Recoupling_Coef(recoup_case,ja,jb,jc,Jab_in,Jab,jtot)
  ! these are just overlaps. 
  ! the return value of SortOrbits determines the case.
  implicit none
  
  integer :: recoup_case,ja,jb,jc,Jab_in,Jab,jtot
  real(8) :: coeff,d6ji
 
  
  select case (recoup_case)
   
    case (0)
       if( Jab==Jab_in) then
          coeff = 1.d0 
       else
          coeff = 0.d0 
       end if 
    
    case (1) 
       coeff= (-1)**((jb-jc+Jab_in)/2) * sqrt((Jab_in+1.d0)*(Jab+1.d0)) * sixj(ja, jb, Jab, jc, jtot, Jab_in)
    case (2) 
       coeff= (-1)**((ja-jb-Jab)/2) * sqrt((Jab_in+1.d0)*(Jab+1.d0)) * sixj(jb, ja, Jab, jc, jtot, Jab_in)
    case (3) 
       coeff = (-1)**((jb+jc+Jab_in-Jab)/2) * sqrt((Jab_in+1.d0)*(Jab+1.d0)) * sixj(jb, ja, Jab, jc, jtot, Jab_in)
    case (4) 
       if (Jab==Jab_in) then 
         coeff = (-1)**((ja+jb-Jab)/2)  
       else
         coeff = 0.d0 
       end if 
    case (5) 
       coeff = -sqrt((Jab_in+1.d0)*(Jab+1.d0)) * sixj(ja, jb, Jab, jc,jtot, Jab_in)
    case default
       coeff = 0.d0 
    end select
    
    Recoupling_Coef = coeff
end function
!==================================================================
!==================================================================
real(8) function iso_clebsch_halfhalf(T1,T2,T3,m1,m2,m3) 
  !dangerous, this does not check very many things. 
  ! don't actually use these for anything important. They aren't nice. 
  ! they are only used in the iso to pn conversion routines. 
  implicit none 
  
  integer :: T1,T2,T3,m1,m2,m3 
  
  if ((m1 + m2).ne.m3) then 
     iso_clebsch_halfhalf =0.d0
     return
  end if 
  
  if (T3 == 0) then 
     iso_clebsch_halfhalf = sign(1,m1)*1.d0/sqrt(2.d0)    
  else 
     select case (abs(m3)) 
     case (0) 
        iso_clebsch_halfhalf = 1.d0 /sqrt(2.d0) 
     case (2) 
        iso_clebsch_halfhalf = 1.d0 
     case default
        iso_clebsch_halfhalf = 0.d0 
     end select
  end if
  
end function 
!==================================================================
!==================================================================
real(8) function iso_clebsch_wholehalf(T1,T2,T3,m1,m2,m3) 
  ! don't actually use these for anything important. They aren't nice. 
  ! they are only used in the iso to pn conversion routines. 
  implicit none 
  
  integer :: T1,T2,T3,m1,m2,m3 
  
  if ((m1 + m2).ne.m3) then 
     iso_clebsch_wholehalf =0.d0
     return
  end if 

     !Tab 1/2 t
  if (T1 == 0) then 
     if (T3 == 1) then 
        iso_clebsch_wholehalf = 1.d0 
     else 
        iso_clebsch_wholehalf = 0.d0 
     end if
  else if (T1 == 2) then 
  
     if (T3 == 1) then  
        select case (abs(m1))   
        case (0) 
           iso_clebsch_wholehalf = -1.d0*sign(1,m2)/sqrt(3.d0) 
        case (2) 
           iso_clebsch_wholehalf = -sqrt(2.d0/3.d0)*sign(1,m2)
        case default 
           iso_clebsch_wholehalf = 0.d0
        end select
     else if (T3 == 3) then 

        select case (abs(m3) ) 
        case (1) 
           if ( abs(m1) == 2) then 
              iso_clebsch_wholehalf = sqrt(1.d0/3.d0)
           else
              iso_clebsch_wholehalf = sqrt(2.d0/3.d0)
           end if 
        case (3) 
           iso_clebsch_wholehalf = 1.d0
        case default
           iso_clebsch_wholehalf = 0.d0 
        end select
     else
        iso_clebsch_wholehalf = 0.d0
     end if
  end if
  
end function 
           
  
end module
        
           
           
     
     
        
        
  
   
  
  



