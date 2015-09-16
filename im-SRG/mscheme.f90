module mscheme
  use basic_IMSRG
  implicit none 


contains
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
real(8) function V_mscheme(a,ma,b,mb,c,mc,d,md,Op,jbas) 
  implicit none 
  
  integer :: a,b,c,d,ma,mb,mc,md
  integer :: Jmin,Jmax,ja,jb,jc,jd,JT,MT
  type(sq_op) :: Op
  type(spd) :: jbas
  real(8) :: dcgi , sm,dcgi00,x
 
  if ( (jbas%itzp(a) + jbas%itzp(b) ).ne. (jbas%itzp(c) + jbas%itzp(d)) ) then
     V_mscheme = 0.d0 
     return
  end if 
  
  if ( mod(jbas%ll(a) + jbas%ll(b),2) .ne. mod(jbas%ll(c) + jbas%ll(d),2) ) then
     V_mscheme = 0.d0 
     return
  end if 
  MT = ma + mb 
  if ((mc + md) .ne. MT ) then 
     V_mscheme = 0.d0 
     return
  end if 
  
  ja = jbas%jj(a)
  jb = jbas%jj(b)  
  jc = jbas%jj(c)  
  jd = jbas%jj(d)  
  
  Jmin = max(abs(ja-jb),abs(jc-jd))
  Jmax = min(ja+jb,jc+jd)
  
  sm = 0.d0
 
  do JT = Jmin,Jmax,2
      
        sm = sm + dcgi( ja, ma , jb, mb, JT, MT) * &
             dcgi( jc, mc , jd, md, JT, MT) * &
             v_elem(a,b,c,d,JT,Op,jbas) ! this is independent of M

  end do 
  
  v_mscheme = sm 
  
end function
!==========================================================================
subroutine allocate_3body_storage( XR , jbas) 
  ! stores three-body matrix elements in M-scheme  
  implicit none 
  
  type(mscheme_3body) :: XR 
  type(spd) :: jbas 
  integer :: occ,unocc, dimo,dimu,dimfull,q,MTOT,Tz,PAR  
  integer :: i,j,k,a,b,c,ix,jx,kx,ax,bx,cx
  real(8) :: memry
  character(1) :: choice
  
 
  XR%nblocks = 8*3*maxval(jbas%jj)+8
  
  allocate(XR%Wmat(XR%nblocks),XR%lam(XR%nblocks))
  allocate(XR%qn_p(XR%nblocks),XR%qn_h(XR%nblocks)) 
 
  memry = 0.d0
  occ = 0
  unocc = 0
! fill XR%qnm arrays
  do i = 1,jbas%total_orbits 
     occ = occ +jbas%con(i)*(jbas%jj(i)+1)
     unocc = unocc +(1-jbas%con(i))*(jbas%jj(i)+1)
  end do 
  
  allocate( XR%qnm_holes( occ, 2) ) ! stores ( jscheme nums, m ) 
  allocate( XR%qnm_parts( unocc, 2)) 
  
  XR%oc_m = occ
  XR%un_m = unocc 
  occ = 0
  unocc = 0
! fill XR%qnm arrays
  do i = 1,jbas%total_orbits 
     
     if (jbas%con(i) == 1) then 
        
        do k = -jbas%jj(i) ,jbas%jj(i) , 2 
           occ = occ + 1
           XR%qnm_holes(occ,1) = i
           XR%qnm_holes(occ,2) = k 
        end do 
        
     else 
     
        do k = -jbas%jj(i) ,jbas%jj(i) , 2 
           unocc = unocc + 1
           XR%qnm_parts(unocc,1) = i
           XR%qnm_parts(unocc,2) = k 
        end do
        
     end if 
     
  end do 

! now iterate over the blocks of conserved quantities
  q = 1
  do MTOT = -3*maxval(jbas%jj) , 3*maxval(jbas%jj) , 2 
     do Tz = - 3, 3 , 2 
        do PAR = 0,1 
           
           allocate(XR%lam(q)%Z(3)) 
           XR%lam(q)%Z(1) = MTOT
           XR%lam(q)%Z(2) = PAR
           XR%lam(q)%Z(3) = Tz 
           
           unocc = 0 
           occ = 0
  
           do i = 1, XR%oc_m
              ix = XR%qnm_holes(i,1)
              do j = i+1, XR%oc_m
                 jx = XR%qnm_holes(j,1)
                 do k = j+1, XR%oc_m
                    kx = XR%qnm_holes(k,1)
                    
                    if ((XR%qnm_holes(i,2)+XR%qnm_holes(j,2) &
                         +XR%qnm_holes(k,2)) .ne. MTOT ) cycle
                    
                    if ((jbas%itzp(ix)+jbas%itzp(jx) + &
                          jbas%itzp(kx)) .ne. Tz )  cycle
                     
                    if (mod(jbas%ll(ix)+jbas%ll(jx)&
                         +jbas%ll(kx),2) .ne. PAR) cycle
                     
                     
                    occ = occ + 1 
                 end do 
              end do 
           end do
            
           do i = 1, XR%un_m
              ix = XR%qnm_parts(i,1)
              do j = i+1, XR%un_m
                 jx = XR%qnm_parts(j,1)
                 do k = j+1, XR%un_m
                    kx = XR%qnm_parts(k,1)
                     
                    if ((XR%qnm_parts(i,2)+XR%qnm_parts(j,2) &
                          +XR%qnm_parts(k,2)) .ne. MTOT ) cycle
                     
                    if ((jbas%itzp(ix)+jbas%itzp(jx) + &
                          jbas%itzp(kx)) .ne. Tz )  cycle
                     
                    if (mod(jbas%ll(ix)+jbas%ll(jx)&
                          +jbas%ll(kx),2) .ne. PAR) cycle
                     
                     
                    unocc = unocc + 1 
                 end do 
              end do 
           end do
            
           allocate(XR%qn_h(q)%Y(occ,3),XR%qn_p(q)%Y(unocc,3))
        
           memry = memry + unocc*occ*8.d0
           unocc = 0 
           occ = 0
           
           do i = 1, XR%oc_m
              ix = XR%qnm_holes(i,1)
              do j = i+1, XR%oc_m
                 jx = XR%qnm_holes(j,1)
                 do k = j+1, XR%oc_m
                    kx = XR%qnm_holes(k,1)
                    
                    if ((XR%qnm_holes(i,2)+XR%qnm_holes(j,2) &
                         +XR%qnm_holes(k,2)) .ne. MTOT ) cycle
                    
                    if ((jbas%itzp(ix)+jbas%itzp(jx) + &
                          jbas%itzp(kx)) .ne. Tz )  cycle
                     
                    if (mod(jbas%ll(ix)+jbas%ll(jx)&
                         +jbas%ll(kx),2) .ne. PAR) cycle
                     
                     
                    occ = occ + 1 
                    XR%qn_h(q)%Y(occ,1)=i
                    XR%qn_h(q)%Y(occ,2)=j
                    XR%qn_h(q)%Y(occ,3)=k
                 end do 
              end do 
           end do
            
           do i = 1, XR%un_m
              ix = XR%qnm_parts(i,1)
              do j = i+1, XR%un_m
                 jx = XR%qnm_parts(j,1)
                 do k = j+1, XR%un_m
                    kx = XR%qnm_parts(k,1)
                     
                    if ((XR%qnm_parts(i,2)+XR%qnm_parts(j,2) &
                          +XR%qnm_parts(k,2)) .ne. MTOT ) cycle
                     
                    if ((jbas%itzp(ix)+jbas%itzp(jx) + &
                          jbas%itzp(kx)) .ne. Tz )  cycle
                     
                    if (mod(jbas%ll(ix)+jbas%ll(jx)&
                          +jbas%ll(kx),2) .ne. PAR) cycle
                     
                     
                    unocc = unocc + 1
                    XR%qn_p(q)%Y(unocc,1)=i
                    XR%qn_p(q)%Y(unocc,2)=j
                    XR%qn_p(q)%Y(unocc,3)=k
                 end do 
              end do 
           end do
 
           q = q + 1
        end do 
     end do 
  end do 
  
  call divide_work_threebody(XR) 
!  memry = memry/ 1024. /1024. 
  
!  IF ( memry > 100.d0 ) then  ! ask user if the memory is getting
                              ! out of control
     
  ! print*, '3 - body ALLOCATION will require in MB: ', memry
  ! print*, 'CONTINUE?? (Y/N)'
  ! read*, choice 
  
  ! if (choice .ne. 'Y') STOP 'aborted...'
  ! else 
  
  ! print*, 'allocating',memry,'MB'
  ! end if 
  
  ! do q = 1, XR%nblocks
  !   allocate(XR%Wmat(q)%X(size(XR%qn_p(q)%Y(:,1)) &
  !        ,size(XR%qn_h(q)%Y(:,1))))
  ! end do 
    
end subroutine     
!==================================================================
subroutine divide_work_threebody(XR) 
  implicit none 
  
  type(mscheme_3body) :: XR
  integer :: A,N,threads,omp_get_num_threads
  integer :: i ,g,q,k,b,j
  
!$omp parallel
  threads=omp_get_num_threads() 
!$omp end parallel
!threads = 1
  b = 0.d0
  do q = 1, XR%nblocks
     b = b + size(XR%qn_p(q)%Y(:,1))*size(XR%qn_h(q)%Y(:,1))
  end do
  
  allocate( XR%direct_omp(threads+1) )
  
  g = 0
  k = 0
  q = 1
  do i = 1,threads
     
     g = g+k 
     j = (b-g)/(threads-i+1) 
     
     k = 0
     
     do while ( q .le. XR%nblocks) 
        k = k + size(XR%qn_p(q)%Y(:,1))*size(XR%qn_h(q)%Y(:,1))
        q = q + 1
        
        if (k .ge. j) then 
           XR%direct_omp(i+1) = q-1
           exit
        end if
     end do 
  end do 
  XR%direct_omp(threads+1) = XR%nblocks
  XR%direct_omp(1) = 0 

end subroutine     
!==================================================================           
!===============================================================
real(8) function fourth_order_restore(OMEGA,H,XR,jbas) 
  ! OKAY THIS ISN'T GENERAL LIKE THE OTHER COMMUTATORS, 
  ! IT ASSUMES THAT YOU ARE USING THE ANTI-HERMITIAN OMEGA
  ! OPERATOR IN MAGNUS AND THE HAMILTONIAN 
  
  ! only W_abcijk type elements are calculated, and they are stored
  ! in XR
  implicit none 
  
  type(sq_op) :: omega,H 
  type(mscheme_3body) :: XR
  type(spd) :: jbas 
  integer :: p,a,b,c,i,j,k , II,AA,q,ompspot
  integer :: ax,bx,cx,ix,jx,kx,tot_threads,thread
  real(8) :: sm,correction,denom
  
  
  print*, 'calculating approximate 4th order triples' 
  print*, 'this will probably take a while...' 
  
  tot_threads = size(XR%direct_omp)-1
  correction = 0.d0 

!$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE), SHARED(OMEGA,H) &
!$OMP& REDUCTION(+:correction)
do thread = 1,tot_threads
  do q = 1+XR%direct_omp(thread), XR%direct_omp(thread+1)
     
    ! print*, q,XR%nblocks
     do II = 1, size(XR%qn_h(q)%Y(:,1))
     
        i = XR%qn_h(q)%Y(II,1)
        j = XR%qn_h(q)%Y(II,2)
        k = XR%qn_h(q)%Y(II,3)
     
        ix = XR%qnm_holes(i,1)
        jx = XR%qnm_holes(j,1)
        kx = XR%qnm_holes(k,1)
        
        do AA = 1, size(XR%qn_p(q)%Y(:,1))
           
           a = XR%qn_p(q)%Y(AA,1)
           b = XR%qn_p(q)%Y(AA,2)
           c = XR%qn_p(q)%Y(AA,3)
        
           ax = XR%qnm_parts(a,1)
           bx = XR%qnm_parts(b,1)
           cx = XR%qnm_parts(c,1)
        
 ! derivative expression ( LONG AS FORK )                     
                    sm = 0.d0 
                    do p = 1, XR%oc_m  ! unforuntately I have to do this FORKING sum for
                       !both particles and holes seperately 
                      
                       sm = sm + & ! each of these are a single matrix element 
( v_mscheme( ax, XR%qnm_parts(a,2) , bx, XR%qnm_parts(b,2), &
  ix, XR%qnm_holes(i,2) , XR%qnm_holes(p,1), XR%qnm_holes(p,2),omega,jbas) &
 
* v_mscheme( XR%qnm_holes(p,1), XR%qnm_holes(p,2) , cx, XR%qnm_parts(c,2), &
  jx, XR%qnm_holes(j,2) , kx, XR%qnm_holes(k,2),H,jbas) - &

( v_mscheme( ax, XR%qnm_parts(a,2) , bx, XR%qnm_parts(b,2), &
  ix, XR%qnm_holes(i,2) , XR%qnm_holes(p,1), XR%qnm_holes(p,2),H,jbas) &
 
* v_mscheme( XR%qnm_holes(p,1), XR%qnm_holes(p,2) , cx, XR%qnm_parts(c,2), &
  jx, XR%qnm_holes(j,2) , kx, XR%qnm_holes(k,2),omega,jbas) ) )
         
!!!  - P_ac 
                       sm = sm - & ! each of these are a single matrix element 
( v_mscheme( cx, XR%qnm_parts(c,2) , bx, XR%qnm_parts(b,2), &
  ix, XR%qnm_holes(i,2) , XR%qnm_holes(p,1), XR%qnm_holes(p,2),omega,jbas) &
 
* v_mscheme( XR%qnm_holes(p,1), XR%qnm_holes(p,2) , ax, XR%qnm_parts(a,2), &
  jx, XR%qnm_holes(j,2) , kx, XR%qnm_holes(k,2),H,jbas) - &

( v_mscheme( ax, XR%qnm_parts(a,2) , bx, XR%qnm_parts(b,2), &
  ix, XR%qnm_holes(i,2) , XR%qnm_holes(p,1), XR%qnm_holes(p,2),H,jbas) &
 
* v_mscheme( XR%qnm_holes(p,1), XR%qnm_holes(p,2) , cx, XR%qnm_parts(c,2), &
  jx, XR%qnm_holes(j,2) , kx, XR%qnm_holes(k,2),omega,jbas) ) )

!!!  - P_bc
                       sm = sm - & ! each of these are a single matrix element 
( v_mscheme( ax, XR%qnm_parts(a,2) , cx, XR%qnm_parts(c,2), &
  ix, XR%qnm_holes(i,2) , XR%qnm_holes(p,1), XR%qnm_holes(p,2),omega,jbas) &
 
* v_mscheme( XR%qnm_holes(p,1), XR%qnm_holes(p,2) , bx, XR%qnm_parts(b,2), &
  jx, XR%qnm_holes(j,2) , kx, XR%qnm_holes(k,2),H,jbas) - &

( v_mscheme( ax, XR%qnm_parts(a,2) , cx, XR%qnm_parts(c,2), &
  ix, XR%qnm_holes(i,2) , XR%qnm_holes(p,1), XR%qnm_holes(p,2),H,jbas) &
 
* v_mscheme( XR%qnm_holes(p,1), XR%qnm_holes(p,2) , bx, XR%qnm_parts(b,2), &
  jx, XR%qnm_holes(j,2) , kx, XR%qnm_holes(k,2),omega,jbas) ) )

!!! - P_ij 
                       sm = sm - & ! each of these are a single matrix element 
( v_mscheme( ax, XR%qnm_parts(a,2) , bx, XR%qnm_parts(b,2), &
  jx, XR%qnm_holes(j,2) , XR%qnm_holes(p,1), XR%qnm_holes(p,2),omega,jbas) &
 
* v_mscheme( XR%qnm_holes(p,1), XR%qnm_holes(p,2) , cx, XR%qnm_parts(c,2), &
  ix, XR%qnm_holes(i,2) , kx, XR%qnm_holes(k,2),H,jbas) - &

( v_mscheme( ax, XR%qnm_parts(a,2) , bx, XR%qnm_parts(b,2), &
  jx, XR%qnm_holes(j,2) , XR%qnm_holes(p,1), XR%qnm_holes(p,2),H,jbas) &
 
* v_mscheme( XR%qnm_holes(p,1), XR%qnm_holes(p,2) , cx, XR%qnm_parts(c,2), &
  ix, XR%qnm_holes(i,2) , kx, XR%qnm_holes(k,2),omega,jbas) ) )

!!!  - P_ik
                       sm = sm - & ! each of these are a single matrix element 
( v_mscheme( ax, XR%qnm_parts(a,2) , bx, XR%qnm_parts(b,2), &
  kx, XR%qnm_holes(k,2) , XR%qnm_holes(p,1), XR%qnm_holes(p,2),omega,jbas) &
 
* v_mscheme( XR%qnm_holes(p,1), XR%qnm_holes(p,2) , cx, XR%qnm_parts(c,2), &
  jx, XR%qnm_holes(j,2) , ix, XR%qnm_holes(i,2),H,jbas) - &

( v_mscheme( ax, XR%qnm_parts(a,2) , bx, XR%qnm_parts(b,2), &
  kx, XR%qnm_holes(k,2) , XR%qnm_holes(p,1), XR%qnm_holes(p,2),H,jbas) &
 
* v_mscheme( XR%qnm_holes(p,1), XR%qnm_holes(p,2) , cx, XR%qnm_parts(c,2), &
  jx, XR%qnm_holes(j,2) , ix, XR%qnm_holes(i,2),omega,jbas) ) )

!!!  + P_ac P_ij
                       sm = sm + & ! each of these are a single matrix element 
( v_mscheme( cx, XR%qnm_parts(c,2) , bx, XR%qnm_parts(b,2), &
  jx, XR%qnm_holes(j,2) , XR%qnm_holes(p,1), XR%qnm_holes(p,2),omega,jbas) &
 
* v_mscheme( XR%qnm_holes(p,1), XR%qnm_holes(p,2) , ax, XR%qnm_parts(a,2), &
  ix, XR%qnm_holes(i,2) , kx, XR%qnm_holes(k,2),H,jbas) - &

( v_mscheme( cx, XR%qnm_parts(c,2) , bx, XR%qnm_parts(b,2), &
  jx, XR%qnm_holes(j,2) , XR%qnm_holes(p,1), XR%qnm_holes(p,2),H,jbas) &
 
* v_mscheme( XR%qnm_holes(p,1), XR%qnm_holes(p,2) , ax, XR%qnm_parts(a,2), &
  ix, XR%qnm_holes(i,2) , kx, XR%qnm_holes(k,2),omega,jbas) ) )

!!!  + P_ac P_ik
                       sm = sm + & ! each of these are a single matrix element 
( v_mscheme( cx, XR%qnm_parts(c,2) , bx, XR%qnm_parts(b,2), &
  kx, XR%qnm_holes(k,2) , XR%qnm_holes(p,1), XR%qnm_holes(p,2),omega,jbas) &
 
* v_mscheme( XR%qnm_holes(p,1), XR%qnm_holes(p,2) , ax, XR%qnm_parts(a,2), &
  jx, XR%qnm_holes(j,2) , ix, XR%qnm_holes(i,2),H,jbas) - &

( v_mscheme( cx, XR%qnm_parts(c,2) , bx, XR%qnm_parts(b,2), &
  kx, XR%qnm_holes(k,2) , XR%qnm_holes(p,1), XR%qnm_holes(p,2),H,jbas) &
 
* v_mscheme( XR%qnm_holes(p,1), XR%qnm_holes(p,2) , ax, XR%qnm_parts(a,2), &
  jx, XR%qnm_holes(j,2) , ix, XR%qnm_holes(i,2),omega,jbas) ) )

!!!  + P_bc P_ij 
                       sm = sm + & ! each of these are a single matrix element 
( v_mscheme( ax, XR%qnm_parts(a,2) , cx, XR%qnm_parts(c,2), &
  jx, XR%qnm_holes(j,2) , XR%qnm_holes(p,1), XR%qnm_holes(p,2),omega,jbas) &
 
* v_mscheme( XR%qnm_holes(p,1), XR%qnm_holes(p,2) , bx, XR%qnm_parts(b,2), &
  ix, XR%qnm_holes(i,2) , kx, XR%qnm_holes(k,2),H,jbas) - &

( v_mscheme( ax, XR%qnm_parts(a,2) , cx, XR%qnm_parts(c,2), &
  jx, XR%qnm_holes(j,2) , XR%qnm_holes(p,1), XR%qnm_holes(p,2),H,jbas) &
 
* v_mscheme( XR%qnm_holes(p,1), XR%qnm_holes(p,2) , bx, XR%qnm_parts(b,2), &
  ix, XR%qnm_holes(i,2) , kx, XR%qnm_holes(k,2),omega,jbas) ) )

!!! + P_bc P_ik 
                       sm = sm + & ! each of these are a single matrix element 
( v_mscheme( ax, XR%qnm_parts(a,2) , cx, XR%qnm_parts(c,2), &
  kx, XR%qnm_holes(k,2) , XR%qnm_holes(p,1), XR%qnm_holes(p,2),omega,jbas) &
 
* v_mscheme( XR%qnm_holes(p,1), XR%qnm_holes(p,2) , bx, XR%qnm_parts(b,2), &
  jx, XR%qnm_holes(j,2) , ix, XR%qnm_holes(i,2),H,jbas) - &

( v_mscheme( ax, XR%qnm_parts(a,2) , cx, XR%qnm_parts(c,2), &
  kx, XR%qnm_holes(k,2) , XR%qnm_holes(p,1), XR%qnm_holes(p,2),H,jbas) &
 
* v_mscheme( XR%qnm_holes(p,1), XR%qnm_holes(p,2) , bx, XR%qnm_parts(b,2), &
  jx, XR%qnm_holes(j,2) , ix, XR%qnm_holes(i,2),omega,jbas) ) )

                       end do 






                    do p = 1,XR%un_m  ! unforuntately I have to do this FORKING sum for
                       !both particles and holes seperately 
                      
                       sm = sm + & ! each of these are a single matrix element 
( v_mscheme( ax, XR%qnm_parts(a,2) , bx, XR%qnm_parts(b,2), &
  ix, XR%qnm_holes(i,2) , XR%qnm_parts(p,1), XR%qnm_parts(p,2),omega,jbas) &
 
* v_mscheme( XR%qnm_parts(p,1), XR%qnm_parts(p,2) , cx, XR%qnm_parts(c,2), &
  jx, XR%qnm_holes(j,2) , kx, XR%qnm_holes(k,2),H,jbas) - &

( v_mscheme( ax, XR%qnm_parts(a,2) , bx, XR%qnm_parts(b,2), &
  ix, XR%qnm_holes(i,2) , XR%qnm_parts(p,1), XR%qnm_parts(p,2),H,jbas) &
 
* v_mscheme( XR%qnm_parts(p,1), XR%qnm_parts(p,2) , cx, XR%qnm_parts(c,2), &
  jx, XR%qnm_holes(j,2) , kx, XR%qnm_holes(k,2),omega,jbas) ) )
         
!!!  - P_ac 
                       sm = sm - & ! each of these are a single matrix element 
( v_mscheme( cx, XR%qnm_parts(c,2) , bx, XR%qnm_parts(b,2), &
  ix, XR%qnm_holes(i,2) , XR%qnm_parts(p,1), XR%qnm_parts(p,2),omega,jbas) &
 
* v_mscheme( XR%qnm_parts(p,1), XR%qnm_parts(p,2) , ax, XR%qnm_parts(a,2), &
  jx, XR%qnm_holes(j,2) , kx, XR%qnm_holes(k,2),H,jbas) - &

( v_mscheme( ax, XR%qnm_parts(a,2) , bx, XR%qnm_parts(b,2), &
  ix, XR%qnm_holes(i,2) , XR%qnm_parts(p,1), XR%qnm_parts(p,2),H,jbas) &
 
* v_mscheme( XR%qnm_parts(p,1), XR%qnm_parts(p,2) , cx, XR%qnm_parts(c,2), &
  jx, XR%qnm_holes(j,2) , kx, XR%qnm_holes(k,2),omega,jbas) ) )

!!!  - P_bc
                       sm = sm - & ! each of these are a single matrix element 
( v_mscheme( ax, XR%qnm_parts(a,2) , cx, XR%qnm_parts(c,2), &
  ix, XR%qnm_holes(i,2) , XR%qnm_parts(p,1), XR%qnm_parts(p,2),omega,jbas) &
 
* v_mscheme( XR%qnm_parts(p,1), XR%qnm_parts(p,2) , bx, XR%qnm_parts(b,2), &
  jx, XR%qnm_holes(j,2) , kx, XR%qnm_holes(k,2),H,jbas) - &

( v_mscheme( ax, XR%qnm_parts(a,2) , cx, XR%qnm_parts(c,2), &
  ix, XR%qnm_holes(i,2) , XR%qnm_parts(p,1), XR%qnm_parts(p,2),H,jbas) &
 
* v_mscheme( XR%qnm_parts(p,1), XR%qnm_parts(p,2) , bx, XR%qnm_parts(b,2), &
  jx, XR%qnm_holes(j,2) , kx, XR%qnm_holes(k,2),omega,jbas) ) )

!!! - P_ij 
                       sm = sm - & ! each of these are a single matrix element 
( v_mscheme( ax, XR%qnm_parts(a,2) , bx, XR%qnm_parts(b,2), &
  jx, XR%qnm_holes(j,2) , XR%qnm_parts(p,1), XR%qnm_parts(p,2),omega,jbas) &
 
* v_mscheme( XR%qnm_parts(p,1), XR%qnm_parts(p,2) , cx, XR%qnm_parts(c,2), &
  ix, XR%qnm_holes(i,2) , kx, XR%qnm_holes(k,2),H,jbas) - &

( v_mscheme( ax, XR%qnm_parts(a,2) , bx, XR%qnm_parts(b,2), &
  jx, XR%qnm_holes(j,2) , XR%qnm_parts(p,1), XR%qnm_parts(p,2),H,jbas) &
 
* v_mscheme( XR%qnm_parts(p,1), XR%qnm_parts(p,2) , cx, XR%qnm_parts(c,2), &
  ix, XR%qnm_holes(i,2) , kx, XR%qnm_holes(k,2),omega,jbas) ) )

!!!  - P_ik
                       sm = sm - & ! each of these are a single matrix element 
( v_mscheme( ax, XR%qnm_parts(a,2) , bx, XR%qnm_parts(b,2), &
  kx, XR%qnm_holes(k,2) , XR%qnm_parts(p,1), XR%qnm_parts(p,2),omega,jbas) &
 
* v_mscheme( XR%qnm_parts(p,1), XR%qnm_parts(p,2) , cx, XR%qnm_parts(c,2), &
  jx, XR%qnm_holes(j,2) , ix, XR%qnm_holes(i,2),H,jbas) - &

( v_mscheme( ax, XR%qnm_parts(a,2) , bx, XR%qnm_parts(b,2), &
  kx, XR%qnm_holes(k,2) , XR%qnm_parts(p,1), XR%qnm_parts(p,2),H,jbas) &
 
* v_mscheme( XR%qnm_parts(p,1), XR%qnm_parts(p,2) , cx, XR%qnm_parts(c,2), &
  jx, XR%qnm_holes(j,2) , ix, XR%qnm_holes(i,2),omega,jbas) ) )

!!!  + P_ac P_ij
                       sm = sm + & ! each of these are a single matrix element 
( v_mscheme( cx, XR%qnm_parts(c,2) , bx, XR%qnm_parts(b,2), &
  jx, XR%qnm_holes(j,2) , XR%qnm_parts(p,1), XR%qnm_parts(p,2),omega,jbas) &
 
* v_mscheme( XR%qnm_parts(p,1), XR%qnm_parts(p,2) , ax, XR%qnm_parts(a,2), &
  ix, XR%qnm_holes(i,2) , kx, XR%qnm_holes(k,2),H,jbas) - &

( v_mscheme( cx, XR%qnm_parts(c,2) , bx, XR%qnm_parts(b,2), &
  jx, XR%qnm_holes(j,2) , XR%qnm_parts(p,1), XR%qnm_parts(p,2),H,jbas) &
 
* v_mscheme( XR%qnm_parts(p,1), XR%qnm_parts(p,2) , ax, XR%qnm_parts(a,2), &
  ix, XR%qnm_holes(i,2) , kx, XR%qnm_holes(k,2),omega,jbas) ) )

!!!  + P_ac P_ik
                       sm = sm + & ! each of these are a single matrix element 
( v_mscheme( cx, XR%qnm_parts(c,2) , bx, XR%qnm_parts(b,2), &
  kx, XR%qnm_holes(k,2) , XR%qnm_parts(p,1), XR%qnm_parts(p,2),omega,jbas) &
 
* v_mscheme( XR%qnm_parts(p,1), XR%qnm_parts(p,2) , ax, XR%qnm_parts(a,2), &
  jx, XR%qnm_holes(j,2) , ix, XR%qnm_holes(i,2),H,jbas) - &

( v_mscheme( cx, XR%qnm_parts(c,2) , bx, XR%qnm_parts(b,2), &
  kx, XR%qnm_holes(k,2) , XR%qnm_parts(p,1), XR%qnm_parts(p,2),H,jbas) &
 
* v_mscheme( XR%qnm_parts(p,1), XR%qnm_parts(p,2) , ax, XR%qnm_parts(a,2), &
  jx, XR%qnm_holes(j,2) , ix, XR%qnm_holes(i,2),omega,jbas) ) )

!!!  + P_bc P_ij 
                       sm = sm + & ! each of these are a single matrix element 
( v_mscheme( ax, XR%qnm_parts(a,2) , cx, XR%qnm_parts(c,2), &
  jx, XR%qnm_holes(j,2) , XR%qnm_parts(p,1), XR%qnm_parts(p,2),omega,jbas) &
 
* v_mscheme( XR%qnm_parts(p,1), XR%qnm_parts(p,2) , bx, XR%qnm_parts(b,2), &
  ix, XR%qnm_holes(i,2) , kx, XR%qnm_holes(k,2),H,jbas) - &

( v_mscheme( ax, XR%qnm_parts(a,2) , cx, XR%qnm_parts(c,2), &
  jx, XR%qnm_holes(j,2) , XR%qnm_parts(p,1), XR%qnm_parts(p,2),H,jbas) &
 
* v_mscheme( XR%qnm_parts(p,1), XR%qnm_parts(p,2) , bx, XR%qnm_parts(b,2), &
  ix, XR%qnm_holes(i,2) , kx, XR%qnm_holes(k,2),omega,jbas) ) )

!!! + P_bc P_ik 
                       sm = sm + & ! each of these are a single matrix element 
( v_mscheme( ax, XR%qnm_parts(a,2) , cx, XR%qnm_parts(c,2), &
  kx, XR%qnm_holes(k,2) , XR%qnm_parts(p,1), XR%qnm_parts(p,2),omega,jbas) &
 
* v_mscheme( XR%qnm_parts(p,1), XR%qnm_parts(p,2) , bx, XR%qnm_parts(b,2), &
  jx, XR%qnm_holes(j,2) , ix, XR%qnm_holes(i,2),H,jbas) - &

( v_mscheme( ax, XR%qnm_parts(a,2) , cx, XR%qnm_parts(c,2), &
  kx, XR%qnm_holes(k,2) , XR%qnm_parts(p,1), XR%qnm_parts(p,2),H,jbas) &
 
* v_mscheme( XR%qnm_parts(p,1), XR%qnm_parts(p,2) , bx, XR%qnm_parts(b,2), &
  jx, XR%qnm_holes(j,2) , ix, XR%qnm_holes(i,2),omega,jbas) ) )

                       end do 
                       

                       ix = jbas%holesb4(ix) + 1
                       jx = jbas%holesb4(jx) + 1
                       kx = jbas%holesb4(kx) + 1
                       ax = jbas%partsb4(ax) + 1
                       bx = jbas%partsb4(bx) + 1                       
                       cx = jbas%partsb4(cx) + 1
                       
                       denom = H%fhh(ix,ix) + H%fhh(jx,jx) + H%fhh(kx,kx) - &
                            H%fpp(ax,ax) - H%fpp(bx,bx) - H%fpp(cx,cx)
                      
                       correction = correction + sm**2/denom
                       
           end do
        end do
     end do 
  end do 
!$OMP END PARALLEL DO

     fourth_order_restore = correction
end function

  
real(8) function slow_triples_restore(OMEGA,H,XR,jbas) 
  ! OKAY THIS ISN'T GENERAL LIKE THE OTHER COMMUTATORS, 
  ! IT ASSUMES THAT YOU ARE USING THE ANTI-HERMITIAN OMEGA
  ! OPERATOR IN MAGNUS AND THE HAMILTONIAN 
  
  ! only W_abcijk type elements are calculated, and they are stored
  ! in XR
  implicit none 
  
  type(sq_op) :: omega,H 
  type(mscheme_3body) :: XR
  type(spd) :: jbas 
  integer :: p,a,b,c,i,j,k , II,AA,q,ompspot
  integer :: ma,mb,mc,mi,mj,mk,mp
  integer :: ja,jb,jc,ji,jj,jk,jp
  integer :: ax,bx,cx,ix,jx,kx,tot_threads,thread
  real(8) :: sm,correction,denom,smtot
  

  
  smtot = 0.d0 
  do ax = 1, H%nsp - H%belowEF 
     a = jbas%parts(ax) 
     ja = jbas%jj(a) 
     
     do ma  = -1*ja, ja, 2 
     
        do bx = 1, H%nsp - H%belowEF
           b = jbas%parts(bx) 
           jb = jbas%jj(b)
        
           do mb  = -1*jb, jb, 2 
              
              do cx = 1, H%nsp - H%belowEF
                 c = jbas%parts(cx) 
                 jc = jbas%jj(c)
              
                 do mc  = -1*jc, jc, 2 

  do ix = 1, H%belowEF 
     i = jbas%holes(ix) 
     ji = jbas%jj(i) 
     
     do mi  = -1*ji, ji, 2 
     
        do jx = 1, H%belowEF
           j = jbas%holes(jx) 
           jj = jbas%jj(j)
        
           do mj  = -1*jj, jj, 2 
              
              do kx = 1, H%belowEF
                 k = jbas%holes(kx) 
                 jk = jbas%jj(k)
              
                 do mk  = -1*jk, jk, 2 

                 
                    sm = 0.d0 
                    !obtain matrix element
                    do p = 1, jbas%total_orbits
                       jp = jbas%jj(p) 
                       
                       do mp = -1*jp, jp, 2 
                          
                    
                    
                    sm = sm + ( v_mscheme(a,ma,p,mp,i,mi,j,mj,H,jbas)*v_mscheme(b,mb,c,mc,p,mp,k,mk,omega,jbas)  &
                    - v_mscheme(a,ma,p,mp,i,mi,j,mj,omega,jbas)*v_mscheme(b,mb,c,mc,p,mp,k,mk,H,jbas) ) 

                    sm = sm + ( v_mscheme(a,ma,p,mp,j,mj,k,mk,H,jbas)*v_mscheme(b,mb,c,mc,p,mp,i,mi,omega,jbas)  &
                    - v_mscheme(a,ma,p,mp,j,mj,k,mk,omega,jbas)*v_mscheme(b,mb,c,mc,p,mp,i,mi,H,jbas) ) 
                    
                    sm = sm + ( v_mscheme(a,ma,p,mp,k,mk,i,mi,H,jbas)*v_mscheme(b,mb,c,mc,p,mp,j,mj,omega,jbas)  &
                    - v_mscheme(a,ma,p,mp,k,mk,i,mi,omega,jbas)*v_mscheme(b,mb,c,mc,p,mp,j,mj,H,jbas) ) 
                    
                    sm = sm + ( v_mscheme(b,mb,p,mp,i,mi,j,mj,H,jbas)*v_mscheme(c,mc,a,ma,p,mp,k,mk,omega,jbas)  &
                    - v_mscheme(b,mb,p,mp,i,mi,j,mj,omega,jbas)*v_mscheme(c,mc,a,ma,p,mp,k,mk,H,jbas) ) 
                    
                    sm = sm + ( v_mscheme(b,mb,p,mp,j,mj,k,mk,H,jbas)*v_mscheme(c,mc,a,ma,p,mp,i,mi,omega,jbas)  &
                    - v_mscheme(b,mb,p,mp,j,mj,k,mk,omega,jbas)*v_mscheme(c,mc,a,ma,p,mp,i,mi,H,jbas) ) 
                    
                    sm = sm + ( v_mscheme(b,mb,p,mp,k,mk,i,mi,H,jbas)*v_mscheme(c,mc,a,ma,p,mp,j,mj,omega,jbas)  &
                    - v_mscheme(b,mb,p,mp,k,mk,i,mi,omega,jbas)*v_mscheme(c,mc,a,ma,p,mp,j,mj,H,jbas) )
                    
                    sm = sm + ( v_mscheme(c,mc,p,mp,i,mi,j,mj,H,jbas)*v_mscheme(a,ma,b,mb,p,mp,k,mk,omega,jbas)  &
                    - v_mscheme(c,mc,p,mp,i,mi,j,mj,omega,jbas)*v_mscheme(a,ma,b,mb,p,mp,k,mk,H,jbas) )
                    
                    sm = sm + ( v_mscheme(c,mc,p,mp,j,mj,k,mk,H,jbas)*v_mscheme(a,ma,b,mb,p,mp,i,mi,omega,jbas)  &
                    - v_mscheme(c,mc,p,mp,j,mj,k,mk,omega,jbas)*v_mscheme(a,ma,b,mb,p,mp,i,mi,H,jbas) )
                    
                    sm = sm + ( v_mscheme(c,mc,p,mp,k,mk,i,mi,H,jbas)*v_mscheme(a,ma,b,mb,p,mp,j,mj,omega,jbas)  &
                    - v_mscheme(c,mc,p,mp,k,mk,i,mi,omega,jbas)*v_mscheme(a,ma,b,mb,p,mp,j,mj,H,jbas) )
                    
                    
                    end do 
                 end do 
                 
                 denom = f_elem(i,i,H,jbas)+f_elem(j,j,H,jbas)+f_elem(k,k,H,jbas)
                 denom = denom - (f_elem(a,a,H,jbas)+f_elem(b,b,H,jbas)+f_elem(c,c,H,jbas))
                 
                 smtot = smtot + sm**2/denom 
                 
                 end do 
                 end do 
                 end do 
                 end do 
                 end do 
                 end do  

                 end do 
                 end do 
                 end do 
                 end do 
                 end do 
                 end do  


                 slow_triples_restore = smtot /36.d0 
end function

end module

! TYPE :: TPD 
!    integer,allocatable,dimension(:,:) ppp,hhh 
! END TYPE TPD


! subroutine enumerate_three_body(threebas,jbas) 
!   implicit none 
  
!   type(spd) :: jbas
!   type(tpd),allocatable,dimension(:) :: threebas
!   integer :: a,b,c,i,j,k,JTot, PAR, Tz,q
!   integer :: Jmaxx,blocks,q,holes,parts
!   integer :: ix,jx,kx,ax,bx,cx
!   integer :: ja,jb,jc,ji,jj,jk 
!   integer :: ta,tb,tc,ti,tj,tk 
!   integer :: la,lb,lc,li,lj,lk 
  
!   holes = sum(jbas%con) 
!   parts = size(jbas%con) - holes
  
!   Jmaxx = maxval(jbas%jj)*3
  
!   blocks= (Jmaxx + 1)/2 * 8 
!   allocate(threebas(blocks))
  
!   q = 1
!   do Jtot = 1,Jmaxx,2 ! Jtot is odd 
     
!      do Tz = -3,3,2 
        
!         do PAR = 0, 1
       
           
           
! ! sum over all possible un-symmetrized hhh combinations 
!    do ix = 1,holes
!       i =jbas%holes(ix)
!       ji = jbas%jj(i)
!       li = jbas%ll(i)      
!       ti = jbas%itzp(i)
      
!       do jx = 1,holes
!          j =jbas%holes(jx)
!          jj = jbas%jj(j)
!          lj = jbas%ll(j)      
!          tj = jbas%itzp(j)
                  
!          Jij_min = abs(ji - jj) 
!          Jij_max = ji + jj 
         
!          do kx= 1,holes
            
!             k =jbas%holes(kx)
!             jk = jbas%jj(k)
!             lk = jbas%ll(k)      
!             tk = jbas%itzp(k)
            
!             if ( .not. triangle(J_
           
              
         
        

  
  
  
  
  
  
  
