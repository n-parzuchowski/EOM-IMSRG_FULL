module me2j_format
  use basic_IMSRG
  implicit none
  
  contains

subroutine get_me2j_spfile(eMaxchr)
  ! This constructs the sps file for
  ! me2j's format, after it's been converted to a pn-basis 
  implicit none 
  
  integer :: e,eMax,l,jj,tz,q,n
  character(2) :: eMaxchr
  character(13) :: fmt
  

  read(eMaxchr,'(I2)') eMax
  eMaxchr = adjustl(eMaxchr) 
  open(unit=24,file=trim(SP_DIRECTORY_LIST(1))//'hk'//trim(eMaxchr)//'Lmax10.sps')
  
  q = 1
  
  fmt = '(5(I5),e17.7)'
  do  e = 0,eMax
     
     do l = mod(e,2),min(e,10),2
        
        n = (e-l)/2
  
        do jj = abs(2*l-1),2*l+1,2
           
           do tz = -1,1,2
              
              write(24,fmt) q, n , l , jj , tz, float(e) 
              q = q+1
           end do 
        end do 
     end do 
  end do 
  
  close(24) 
  print*, trim(SP_DIRECTORY_LIST(1))//'hk'//trim(eMaxchr)//'Lmax10.sps'
  STOP
end subroutine


subroutine read_me2j_interaction(H,jbas,htype,hw,rr,pp) 
  use gzipmod
  implicit none 
  
  integer :: nlj1,nlj2,nnlj1,nnlj2,j,T,Mt,nljMax,endpoint,j_min,j_max,htype
  integer :: l1,l2,ll1,ll2,j1,j2,jj1,jj2,Ntot,i,q,bospairs,qx,ta,tb,tc,td
  integer :: eMax,iMax,jmax,jmin,JT,a,b,c,d,C1,C2,i1,i2,pre,COM,x,PAR,endsz
  integer,allocatable,dimension(:) :: indx 
  real(8),allocatable,dimension(:) :: ME,MEpp,MErr,me_fromfile,ppff,rrff
  real(8) :: V,g1,g2,g3,hw,pre2
  type(spd) :: jbas 
  type(sq_op) :: H
  type(sq_op),optional :: pp,rr
  logical :: pp_calc,rr_calc
  character(1) :: rem
  character(2) :: eMaxchr
  type(c_ptr) :: buf,buf2,buf3
  integer(c_int) :: hndle,hndle2,hndle3,sz,sz2,sz3,rx
  character(kind=C_CHAR,len=200) :: buffer,buffer2,buffer3
  
  pp_calc=.false.
  rr_calc=.false.
  COM = 0
  
  if (present(pp)) pp_calc=.true.
  if (present(rr)) rr_calc=.true.
  if (htype == 1) COM = 1

  Ntot = jbas%total_orbits/2
  
  i = 0
  q = 0
  ! allocate array to store positions of matrix elements
  bospairs = bosonic_tp_index(Ntot,Ntot,Ntot) 
  allocate(indx(bospairs**2)) 
  indx = 0

  ! move in increments of two, because I use a pn basis,
  !  heiko's is isospin coupled (half the states) 

  eMax = maxval(jbas%e)   
  
  ! counting the states, and labeling them
  do nlj1=1, 2*Ntot,2 
     l1= jbas%ll(nlj1)
     j1= jbas%jj(nlj1)
     
     do nlj2 = 1, nlj1,2
        l2= jbas%ll(nlj2)
        j2= jbas%jj(nlj2)
        
       ! if (nint(jbas%e(nlj1) + jbas%e(nlj2)) > eMax )  exit
      
        do nnlj1 = 1 ,nlj1 , 2
           ll1= jbas%ll(nnlj1)
           jj1= jbas%jj(nnlj1)
           
           endpoint = nnlj1 
           if ( nlj1==nnlj1 ) endpoint = nlj2
         
           do nnlj2 = 1, endpoint , 2 
              ll2= jbas%ll(nnlj2)
              jj2= jbas%jj(nnlj2)
              
             ! if (nint(jbas%e(nlj1) + jbas%e(nlj2)) > eMax )  exit 
             
              if (mod(l1+l2,2) .ne. mod(ll1+ll2,2)) cycle
              jmin = max( abs(j1-j2) , abs(jj1-jj2) ) 
              jmax = min( j1+j2  , jj1+jj2) 
              
              if (jmin > jmax) cycle 
      
              indx(bosonic_tp_index((nlj2+1)/2,(nlj1+1)/2,Ntot)&
                   +bospairs*(bosonic_tp_index((nnlj2+1)/2,(nnlj1+1)/2,Ntot)-1)) = i+1 
        
              do JT = jmin,jmax,2
                 i = i + 4
              end do 
         
           
           end do
        end do 
     end do 
  end do 
  iMax = i 
  
  allocate(me(iMax)) 
  allocate(mepp(iMax))
  allocate(me_fromfile(10)) 
  allocate(ppff(10)) 
  if (rr_calc) then 
     allocate(rrff(10)) 
     allocate(merr(iMax)) 
  end if 
  

  write(eMaxchr,'(I2)') eMax 
  eMaxchr = adjustl(eMaxchr)  

  ! using zlib c library, which is bound with fortran in file "gzipmod.f90" 
  
  ! I don't know why you have to tack on those //achars(0) but it seems nessecary 
  hndle=gzOpen(trim(TBME_DIR)//trim(adjustl(intfile))//achar(0),"r"//achar(0)) 
  
  ! opening the pipj and rirj files 
  if (len(trim(eMaxchr)) == 1) then 
     hndle2=gzOpen(trim(TBME_DIR)//"tpp_eMax0"//trim(eMaxchr)//".me2j.gz"//achar(0),"r"//achar(0)) 
     if (rr_calc) then 
        hndle3=gzOpen(trim(TBME_DIR)//"r1r2_eMax0"//trim(eMaxchr)//".me2j.gz"//achar(0),"r"//achar(0)) 
     end if
  else
      hndle2=gzOpen(trim(TBME_DIR)//"tpp_eMax"//trim(eMaxchr)//".me2j.gz"//achar(0),"r"//achar(0)) 
     if (rr_calc) then 
        hndle3=gzOpen(trim(TBME_DIR)//"r1r2_eMax"//trim(eMaxchr)//".me2j.gz"//achar(0),"r"//achar(0)) 
     end if
  end if 
  
  
  sz=200;sz2=200;sz3=200 !c_ints, don't reuse them 
  
  buf=gzGets(hndle,buffer,sz) 
  buf2=gzGets(hndle2,buffer2,sz2)
  if (rr_calc) buf3=gzGets(hndle3,buffer3,sz3)
  
  endpoint = 10 
  write(rem,'(I1)') endpoint-1
  endsz = 130 
  
  do i = 1,iMax,10
  
     if (i+10 > iMax) then 
        deallocate(me_fromfile)
        deallocate(ppff) 
        allocate(me_fromfile( iMax - i + 1) ) 
        allocate(ppff(iMax-i + 1)) 
        if (rr_calc) then 
           deallocate(rrff)
           allocate(rrff(iMax-i + 1)) 
        end if 
        endpoint = iMax-i + 1
        endsz = 13+(endpoint-1)*13 
        write(rem,'(I1)') endpoint-1
     end if
  
     buf = gzGets(hndle,buffer(1:sz),sz)
     buf2 = gzGets(hndle2,buffer2(1:sz2),sz2)
     
  
     read(buffer(1:endsz),'(f12.7,'//rem//'(f13.7))') me_fromfile 
     read(buffer2(1:endsz),'(f12.7,'//rem//'(f13.7))') ppff 
   
     if (rr_calc) then 
        
        buf3 = gzGets(hndle3,buffer3(1:sz3),sz3)
        read(buffer3(1:endsz),'(f12.7,'//rem//'(f13.7))') rrff 
        
        do j = 1,endpoint 
           ME(i+j-1) = me_fromfile(j)
           MEpp(i+j-1) = ppff(j) 
           MErr(i+j-1) = rrff(j)
        end do
        
     else
        
        do j = 1,endpoint 
           ME(i+j-1) = me_fromfile(j)
           MEpp(i+j-1) = ppff(j) 
        end do
     end if 
  end do

  rx = gzClose(hndle) 
  rx = gzClose(hndle2)
  if(rr_calc) then 
     rx = gzClose(hndle3)
     deallocate(rrff)
     allocate(rrff(4))
  end if 
  
  deallocate(me_fromfile)
  deallocate(ppff)
  allocate(me_fromfile(4))
  allocate(ppff(4)) 
  
  ! redo this loop to put everything in pn basis
  
  i=0
  do nlj1=1, 2*Ntot,2 
     l1= jbas%ll(nlj1)
     j1= jbas%jj(nlj1)
     
     do nlj2 = 1, nlj1,2
        l2= jbas%ll(nlj2)
        j2= jbas%jj(nlj2)
         
        do nnlj1 = 1 ,nlj1 , 2
           ll1= jbas%ll(nnlj1)
           jj1= jbas%jj(nnlj1)
           
           endpoint = nnlj1 
           if ( nlj1==nnlj1 ) endpoint = nlj2
         
           do nnlj2 = 1, endpoint , 2 
              ll2= jbas%ll(nnlj2)
              jj2= jbas%jj(nnlj2)
                  
              if (mod(l1+l2,2) .ne. mod(ll1+ll2,2)) cycle
              jmin = max( abs(j1-j2) , abs(jj1-jj2) ) 
              jmax = min( j1+j2  , jj1+jj2) 
              PAR = mod(l1+l2,2) 
              if (jmin > jmax) cycle 
            
              do JT = jmin,jmax,2
                 me_fromfile=ME(i+1:i+4)
                 ppff = MEpp(i+1:i+4)
                 if (rr_calc) rrff = MErr(i+1:i+4)
                 i = i + 4 ! four different TMt qnums
 

 
!sum over all isospin configs
do a = nlj1,nlj1+1
   do b = nlj2,nlj2+1
      do c = nnlj1,nnlj1+1 
         do d = nnlj2,nnlj2+1  
           
            ! conversion factor to mT scheme 
            pre2 = 1.d0 
            if ( a == b ) pre2 = pre2*sqrt(0.5d0) 
            if ( c == d ) pre2 = pre2*sqrt(0.5d0) 
            
            ! heikos labeling is backwards
            ta = -jbas%itzp(a)
            tb = -jbas%itzp(b)
            tc = -jbas%itzp(c)
            td = -jbas%itzp(d)
            
            T = ta+tb
            if (tc+td .ne. T) cycle
            T = -T/2
            q = block_index(JT,T,Par)

     ! convert to pn matrix element       
     V =  0.125d0*(ta-tb)*(tc-td)*me_fromfile(1)+&   ! 00 clebsch
          kron_del(ta+tb,-2)*kron_del(tc+td,-2)*me_fromfile(2)+& ! 1-1 
          kron_del(ta+tb,2)*kron_del(tc+td,2)*me_fromfile(4)+& !11 
          0.125d0*abs((ta-tb)*(tc-td))*me_fromfile(3) !10 

     ! pipj 
     g3 = 0.125d0*(ta-tb)*(tc-td)*ppff(1)+&   ! 00 clebsch
          kron_del(ta+tb,-2)*kron_del(tc+td,-2)*ppff(2)+& ! 1-1 
          kron_del(ta+tb,2)*kron_del(tc+td,2)*ppff(4)+& !11 
          0.125d0*abs((ta-tb)*(tc-td))*ppff(3) !10 

     if (rr_calc) then 
        g2 = 0.125d0*(ta-tb)*(tc-td)*rrff(1)+&   ! 00 clebsch
          kron_del(ta+tb,-2)*kron_del(tc+td,-2)*rrff(2)+& ! 1-1 
          kron_del(ta+tb,2)*kron_del(tc+td,2)*rrff(4)+& !11 
          0.125d0*abs((ta-tb)*(tc-td))*rrff(3) !10 
     end if 
     
     ! getting rid of weird mass scaling 
     g3 = -2.d0*g3/hbarc2_over_mc2 

     
     ! center of mass subtraction
     V = (V - g3*COM*hw/(H%Aneut+H%Aprot)) *pre2 
  

     
     
     g3 = g3*pre2
     g2 = g2*pre2
     
     C1 = jbas%con(a)+jbas%con(b) + 1 !ph nature
     C2 = jbas%con(c)+jbas%con(d) + 1
    
     qx = C1*C2
     qx = qx + adjust_index(qx)   !Vpppp nature  

     ! get the indeces in the correct order
     pre = 1
     if ( a > b )  then 
        
        x = bosonic_tp_index(b,a,Ntot*2) 
        j_min = H%xmap(x)%Z(1)  
        i1 = H%xmap(x)%Z( (JT-j_min)/2 + 2) 
        pre = (-1)**( 1 + (jbas%jj(a) + jbas%jj(b) -JT)/2 ) 
     else
       ! if (a == b) pre = pre * sqrt( 2.d0 )
       
        x = bosonic_tp_index(a,b,Ntot*2) 
        j_min = H%xmap(x)%Z(1)  
        i1 = H%xmap(x)%Z( (JT-j_min)/2 + 2) 
     end if
  
     if (c > d)  then     
        
        x = bosonic_tp_index(d,c,Ntot*2) 
        j_min = H%xmap(x)%Z(1)  
        i2 = H%xmap(x)%Z( (JT-j_min)/2 + 2) 
        
        pre = pre * (-1)**( 1 + (jbas%jj(c) + jbas%jj(d) -JT)/2 ) 
     else 
       ! if (c == d) pre = pre * sqrt( 2.d0 )
      
        x = bosonic_tp_index(c,d,Ntot*2) 
        j_min = H%xmap(x)%Z(1)  
        i2 = H%xmap(x)%Z( (JT-j_min)/2 + 2) 
     end if

     ! kets/bras are pre-scaled by sqrt(2) if they 
     ! have two particles in the same sp-shell
        
     ! get the units right. I hope 
     
     if ((qx == 1) .or. (qx == 5) .or. (qx == 4)) then 
        H%mat(q)%gam(qx)%X(i2,i1)  = V *pre
        H%mat(q)%gam(qx)%X(i1,i2)  = V *pre
        
        if (rr_calc) then 
           rr%mat(q)%gam(qx)%X(i2,i1)  = hw*g2*pre/(H%Aneut + H%Aprot)
           rr%mat(q)%gam(qx)%X(i1,i2)  = hw*g2*pre/(H%Aneut + H%Aprot)
        end if 

        if (pp_calc) then 
           pp%mat(q)%gam(qx)%X(i2,i1)  = hw*g3*pre/(H%Aneut + H%Aprot)
           pp%mat(q)%gam(qx)%X(i1,i2)  = hw*g3*pre/(H%Aneut + H%Aprot)
        end if 

        
     else if (C1>C2) then
        H%mat(q)%gam(qx)%X(i2,i1)  = V *pre
        
        if (rr_calc) then 
           rr%mat(q)%gam(qx)%X(i2,i1)  = hw*g2*pre/(H%Aneut + H%Aprot) 
        end if
        
        if (pp_calc) then 
           pp%mat(q)%gam(qx)%X(i2,i1)  = hw*g3*pre/(H%Aneut + H%Aprot) 
        end if

     else
        H%mat(q)%gam(qx)%X(i1,i2) = V * pre
        
        if (rr_calc) then 
           rr%mat(q)%gam(qx)%X(i1,i2)  = hw*g2*pre/(H%Aneut + H%Aprot) 
        end if

        if (pp_calc) then 
           pp%mat(q)%gam(qx)%X(i1,i2)  = hw*g3*pre/(H%Aneut + H%Aprot) 
        end if

     end if 
     ! I shouldn't have to worry about hermiticity here, input is assumed to be hermitian
    
 end do;end do; end do; end do !end sums over isospin  
     
             end do ! end sum over j 
         
           
           end do  !end sums over Tcoupled lables
        end do 
     end do 
  end do   
 
end subroutine
!==========================================================
subroutine export_to_nushellX(H,jbas) 
  implicit none 
  
  type(sq_op) :: H
  type(spd) :: jbas
  integer :: i,j,k,l,J_tot,T_tot,Jx,Tx
  integer :: tot_pp_elems,tot_nn_elems,tot_pn_elems
  real(8) :: pre,mat_elem
  real(8),allocatable,dimension(:) :: sp_ens

  

  open(unit=34,file='../../hamiltonians/'// &
       trim(adjustl(prefix))//'_nushell_1bd.sps') 

  write(34,'(A3)') 'iso' ! SHERPA only works for isospin coupled right now
  write(34,'(I2)') jbas%total_orbits/2 ! isospin coupled so devide by 2
  
  ! write .sps file
  do i = 1,jbas%total_orbits,2
     
     write(34,'(f10.1,2(f5.1),I3)') float(jbas%nn(i)),float(jbas%ll(i)), &
          float(jbas%jj(i))/2.0, 2 
     
  end do 
  close(34)
  
  
  open(unit=33,file='t1') ! pn 
  open(unit=34,file='t3') ! pp  
  open(unit=35,file='t4') ! nn
  
  tot_pn_elems = 0 
  tot_nn_elems = 0 
  tot_pp_elems = 0 
  
  do T_tot = 0,1
     do J_tot = 0,jbas%Jtotal_max
        
        
        do i = 1, jbas%total_orbits/2
           do j = i, jbas%total_orbits/2
              
              do k = 1, jbas%total_orbits/2
                 do l = k,jbas%total_orbits/2
                                      
                    pre = 1.d0 
                    if ( i == j)  pre = pre * .70710681186  
                    if ( k == l)  pre = pre * .70710681186 
                    
                    Jx = J_tot*2
                    Tx = T_tot*2
          
                    mat_elem = 0.5 * ( v_elem( 2*i-1 , 2*j , 2*k-1, 2*l , Jx,H,jbas ) + &
                         v_elem( 2*i , 2*j-1 , 2*k , 2*l-1 , Jx,H,jbas ) - (-1)**(T_tot) * &
                         ( v_elem( 2*i-1 , 2*j , 2*k , 2*l-1 , Jx,H,jbas ) +  &
                         v_elem( 2*i , 2*j-1 , 2*k-1 , 2*l , Jx,H,jbas ) ) ) !* pre**4
                     
                         mat_elem = mat_elem +0.5 * ( T_twobody( 2*i-1 , 2*j , 2*k-1,2*l , Jx,Tx,H,jbas ) + &
                              T_twobody( 2*i , 2*j-1 , 2*k,2*l-1 , Jx,Tx,H,jbas ) - (-1)**(T_tot) * &
                              ( T_twobody( 2*i-1 , 2*j , 2*k,2*l-1 , Jx,Tx,H,jbas ) +  &
                              T_twobody( 2*i , 2*j-1 , 2*k-1,2*l , Jx,Tx,H,jbas ) ) ) !* pre**4
                    
                    if (abs(mat_elem) > 1e-10) then 
                       write(33,'(4(I3),I5,I3,f15.7)') i,j,k,l,J_tot,T_tot,mat_elem
                       tot_pn_elems = tot_pn_elems + 1 
                    end if 
  
                    if (T_tot == 0)  cycle ! pp and nn terms are not included in T=0
                    
                    ! pp terms  Mt = -1
                    mat_elem = v_elem( 2*i-1 , 2*j-1 , 2*k-1, 2*l-1 , Jx,H,jbas )*pre
                     
                    if (abs(mat_elem) > 1e-10) then 
                       write(34,'(4(I3),I5,I3,f15.7)') i,j,k,l,J_tot,T_tot,mat_elem
                       tot_pp_elems = tot_pp_elems + 1 
                    end if 
                   
                    ! nn terms Mt = +1
                    mat_elem = v_elem( 2*i , 2*j , 2*k, 2*l , Jx,H,jbas )*pre
                   
                    if (abs(mat_elem) > 1e-10) then 
                       write(35,'(4(I3),I5,I3,f15.7)') i,j,k,l,J_tot,T_tot,mat_elem
                       tot_nn_elems = tot_nn_elems + 1 
                    end if 
                   
                    
                  end do 
               end do 
            end do 
         end do 
     end do 
  end do 
  close(33)
  close(34)
  close(35)
 
  allocate(sp_ens(jbas%total_orbits/2)) 
  do i = 1, jbas%total_orbits/2
     sp_ens(i) = T_elem(2*i,2*i,H,jbas)  
  end do
!=====write pn to file=======
  open(unit=33,file='t2') 
  write(33,*) tot_pn_elems, sp_ens,1.d0,1.d0,0.d0
  close(33) 
  
  call system('cat t2 t1 > '//'../../hamiltonians/'// &
       trim(adjustl(prefix))//'_nushell_TBME_Tz0.int && rm t1 && rm t2') 
!=====write pp to file=======
  open(unit=33,file='t2') 
  write(33,*) tot_pp_elems, sp_ens,1.d0,1.d0,0.d0
  close(33) 
  
  call system('cat t2 t3 > '//'../../hamiltonians/'// &
       trim(adjustl(prefix))//'_nushell_TBME_TzM.int && rm t3 && rm t2') 
!=====write nn to file=======
  open(unit=33,file='t2') 
  write(33,*) tot_nn_elems, sp_ens,1.d0,1.d0,0.d0
  close(33) 
  
  call system('cat t2 t4 > '//'../../hamiltonians/'// &
       trim(adjustl(prefix))//'_nushell_TBME_Tz1.int && rm t4 && rm t2') 
  
end subroutine

subroutine read_me2b_interaction(H,jbas,htype,hw,rr,pp,Lawson) 
  use gzipmod
  implicit none 
  
  integer :: nlj1,nlj2,nnlj1,nnlj2,j,T,Mt,nljMax,endpoint,j_min,j_max,htype,Lmax
  integer :: l1,l2,ll1,ll2,j1,j2,jj1,jj2,Ntot,i,q,bospairs,qx,ta,tb,tc,td,bMax
  integer :: eMax,iMax,jmax,jmin,JT,a,b,c,d,C1,C2,i1,i2,COM,x,PAR,endsz,aMax
  integer :: t1,t2,lj1,lj2,n1,n2,Pi,Tz,AA,BB,qq,iq,jq,a_hh,a_ph,a_pp
  integer,allocatable,dimension(:) :: indx , nMax_lj
  real(8),allocatable,dimension(:) :: ME,MEpp,MErr,me_fromfile,ppff,rrff
  real(8) :: V,g1,g2,g3,hw,pre2,pre,sm
  type(spd) :: jbas 
  type(sq_op) :: H,stors
  type(sq_op),optional :: pp,rr
  logical :: pp_calc,rr_calc,file_there,noteffedup
  character(1),optional :: Lawson
  character(2) :: eMaxchr
  character(200) :: me1bfile
  integer :: lj,twol,twoj,ljMax,idx,idxx
  integer,allocatable,dimension(:,:) :: SPBljs 
  type(c_ptr) :: buf,buf2,buf3
  integer(c_int) :: hndle,hndle2,hndle3,sz,sz2,sz3,rx
  character(kind=C_CHAR,len=200) :: buffer,buffer2,buffer3

  
  rr_calc = .false.
  pp_calc = .false. 
  Ntot = jbas%total_orbits
  Lmax = maxval(jbas%ll) 
  eMax = maxval(jbas%e)
! populate lj array
  lj = 0
  do twol = 0, 2 * Lmax , 2
     do  twoj = abs(twol - 1) , twol+1 , 2
        lj=lj+1
     end do 
  end do 
  ljMax = lj 
  allocate(SPBljs(lj,2)) 
  allocate(nMax_lj(lj))
  
  lj = 0
  do twol = 0, 2 * Lmax , 2
     do  twoj = abs(twol - 1) , twol+1 , 2
        lj=lj+1
        SPBljs(lj,1) = twol
        sPBljs(lj,2) = twoj
        nMax_lj(lj) = (eMax - twol/2)/2
     end do
  end do
  
  allocate(stors%mat(H%nblocks))

  
 
  ! using zlib c library, which is bound with fortran in file "gzipmod.f90" 
  
  ! I don't know why you have to tack on those //achars(0) but it seems nessecary 
  if ( present(Lawson) ) then 
     hndle=gzOpen(trim(TBME_DIR)//"O16_Hcm_eMax10_hwHO020.ham0.me2b.gz"//achar(0),"r"//achar(0)) 
  else 
     hndle=gzOpen(trim(TBME_DIR)//trim(adjustl(intfile))//achar(0),"r"//achar(0)) 
  end if 
  
! here is where we start dealing with the two body piece
  
  sz=200
  
  buf=gzGets(hndle,buffer,sz) 
  buf=gzGets(hndle,buffer,sz) 
  buf=gzGets(hndle,buffer,sz) 
  
  read(buffer(6:9),'(I4)',iostat=qq) bMax 
  if (qq .ne. 0 ) then 
     read(buffer(6:8),'(I3)',iostat=qq) bMax 
     if (qq .ne. 0 ) then 
        read(buffer(6:7),'(I2)',iostat=qq) bMax
     end if 
  end if

 q = 0
! heiko's code calls protons 1 and neutrons 0

do Tz = 1 , -1, -1  
  do Pi = 0,1
     do JT = 0, 2*jbas%Jtotal_max,2 
        if ((Lmax == eMax) .and. (JT == 2*jbas%Jtotal_max)&
             .and. (Abs(Tz)==1)) cycle
        if ((JT == 2*jbas%Jtotal_max) .and. (Pi==1)) cycle
        q = q+1
     
        stors%mat(q)%lam(1) = JT
        stors%mat(q)%lam(2) = Pi
        stors%mat(q)%lam(3) = Tz
        
        select case ( Tz)
           case ( -1 ) 
              t1 = -1
              t2 = -1
           case ( 0 ) 
              t1 = 1 
              t2 = -1
           case ( 1 ) 
              t1 = 1
              t2 = 1 
        end select
                 
        a = 0
        a_hh = 0
        a_ph = 0 
        a_pp = 0
        do lj1 = 1, ljMax
           do lj2 = 1, ljMax

              j1 = SPBljs(lj1,2) 
              j2 = SPBljs(lj2,2)
              l1 = SPBljs(lj1,1)/2
              l2 = SPBljs(lj2,1)/2
           
              if ( ( JT < abs(j1-j2) ) .or. (JT > j1 + j2) ) cycle
              if ( mod(l1 + l2 ,2 ) .ne.Pi ) cycle 

              
              do n1 = 0,nMax_lj(lj1)
                 idx = (lj1-1) * (nMax_lj(1) +1 ) +n1 
                 do n2 = 0,nMax_lj(lj2) 
                    idxx = (lj2-1) * (nMax_lj(1) +1 ) +n2                 
                 
                    if ( (Tz .ne. 0) .and. (idx > idxx) ) cycle
                    if ( (mod(JT/2,2) == 1) .and. (lj1==lj2) .and. &
                         (n1==n2) .and. (Tz .ne. 0) ) cycle
                  
                    ! now search for sp labels
                    do i = 1, jbas%total_orbits 
                       if ( jbas%jj(i) .ne. j1 ) cycle
                       if ( jbas%nn(i) .ne. n1 ) cycle
                       if ( jbas%ll(i) .ne. l1 ) cycle
                       if ( jbas%itzp(i) .ne. t1 ) cycle                     
                       exit
                    end do
                  
                    do j = 1, jbas%total_orbits 
                       if ( jbas%jj(j) .ne. j2 ) cycle
                       if ( jbas%nn(j) .ne. n2 ) cycle
                       if ( jbas%ll(j) .ne. l2 ) cycle
                       if ( jbas%itzp(j) .ne. t2 ) cycle                     
                       exit
                    end do
                  
                    
                    a = a + 1
                    select case(jbas%con(i) + jbas%con(j))
                       case(0)
                          a_pp = a_pp + 1
                       case(1)
                          a_ph = a_ph + 1
                       case(2)
                          a_hh = a_hh + 1
                    end select

                 end do
              end do
           end do
        end do
        
    
        stors%mat(q)%npp = a_pp 
        stors%mat(q)%nph = a_ph
        stors%mat(q)%nhh = a_hh
        stors%mat(q)%ntot = a 

    
     end do
  end do
end do
  
sz = 20
 q = 0
! heiko's code calls protons 1 and neutrons 0

do Tz = 1 , -1, -1  
  do Pi = 0,1
     do JT = 0, 2*jbas%Jtotal_max,2 
        if ((Lmax == eMax) .and. (JT == 2*jbas%Jtotal_max)&
             .and. (Abs(Tz)==1)) cycle
        if ((JT == 2*jbas%Jtotal_max) .and. (Pi==1)) cycle
        q = q+1
     
        buf=gzGets(hndle,buffer,sz) 
  
        read(buffer(10:16),'(I6)') aMax 
        
 
        allocate(stors%mat(q)%qn(1)%Y( aMax+1, 2) ) 
  
        select case ( Tz)
           case ( -1 ) 
              t1 = -1
              t2 = -1
           case ( 0 ) 
              t1 = 1 
              t2 = -1
           case ( 1 ) 
              t1 = 1
              t2 = 1 
        end select
                 
        a = 0
        a_hh = 0
        a_ph = 0 
        a_pp = 0
        do lj1 = 1, ljMax
           do lj2 = 1, ljMax

              j1 = SPBljs(lj1,2) 
              j2 = SPBljs(lj2,2)
              l1 = SPBljs(lj1,1)/2
              l2 = SPBljs(lj2,1)/2
           
              if ( ( JT < abs(j1-j2) ) .or. (JT > j1 + j2) ) cycle
              if ( mod(l1 + l2 ,2 ) .ne.Pi ) cycle 

              
              do n1 = 0,nMax_lj(lj1)
                 idx = (lj1-1) * (nMax_lj(1) +1 ) +n1 
                 do n2 = 0,nMax_lj(lj2) 
                    idxx = (lj2-1) * (nMax_lj(1) +1 ) +n2                 
                 
                    if ( (Tz .ne. 0) .and. (idx > idxx) ) cycle
                    if ( (mod(JT/2,2) == 1) .and. (lj1==lj2) .and. &
                         (n1==n2) .and. (Tz .ne. 0) ) cycle
                  
                    ! now search for sp labels
                    do i = 1, jbas%total_orbits 
                       if ( jbas%jj(i) .ne. j1 ) cycle
                       if ( jbas%nn(i) .ne. n1 ) cycle
                       if ( jbas%ll(i) .ne. l1 ) cycle
                       if ( jbas%itzp(i) .ne. t1 ) cycle                     
                       exit
                    end do
                  
                    do j = 1, jbas%total_orbits 
                       if ( jbas%jj(j) .ne. j2 ) cycle
                       if ( jbas%nn(j) .ne. n2 ) cycle
                       if ( jbas%ll(j) .ne. l2 ) cycle
                       if ( jbas%itzp(j) .ne. t2 ) cycle                     
                       exit
                    end do
                  
                    
                    select case(jbas%con(i) + jbas%con(j))
                       case(0)
                          a_pp = a_pp + 1   
                          a = stors%mat(q)%nhh + &
                               stors%mat(q)%nph + a_pp
                       case(1)
                          a_ph = a_ph + 1
                          a = stors%mat(q)%nhh + a_ph
                       case(2)
                          a_hh = a_hh + 1
                          a = a_hh 
                    end select
                      
                    stors%mat(q)%qn(1)%Y(a,1) = i
                    stors%mat(q)%qn(1)%Y(a,2) = j
                    
                      
                 end do
              end do
           end do
        end do
            
     end do
  end do
end do



! okay for now on there is a space, then a line that specifies the block
! and aMax for the block. We already know that stuff so we will just ignore
! it and read in the matrix elements
sz = 200
qq = 0
do Tz = 1 , -1, -1  
  do Pi = 0,1
     do JT = 0, 2*jbas%Jtotal_max,2 
        if ((Lmax == eMax) .and. (JT == 2*jbas%Jtotal_max)&
             .and. (Abs(Tz)==1)) cycle
        if ((JT == 2*jbas%Jtotal_max) .and. (Pi==1)) cycle
        qq = qq+1 ! heikos block index
      !  print*, qq
        q = block_index(JT,Tz,Pi) ! my block index
        
        ! space then label
        buf=gzGets(hndle,buffer,sz) 
        
        if (noteffedup) then  
           buf=gzGets(hndle,buffer,sz)
        end if
        noteffedup=.true. 
        ! ignore
        sz = 200
    
      do 
            
      buf=gzGets(hndle,buffer,sz)
      
      read(buffer(1:1),'(I1)',iostat=iq) AA
      if (iq .ne. 0 ) then 
         ! because it's damned near impossible to 
         ! figure out where the block ends
         noteffedup=.false. 
         exit
      end if
         ! figure out where the spaces are that separate things 
      i = 1
      ! first space
       do 
         if ( buffer(i:i) == ' ' ) then
            read(buffer(1:i-1),'(I5)')  AA 
            i = i + 2
            j = i 
            exit
         end if
         i = i + 1
      end do
      ! second space
      do 
         if ( buffer(i:i) == ' ' ) then 
            read(buffer(j:i-1),'(I5)') BB 
            i = i + 2
            exit
         end if
         i = i + 1
      end do
      
      
      AA = AA + 1
      BB = BB + 1

      
      if ( buffer(i:i) == '-' ) then 
         ! negative number
         read(buffer(i:i+10), '( f11.8 )' )  V 
      else 
         ! positive
         read(buffer(i:i+9), '( f10.8 )' )  V 
      end if 
      

      ! oTay should have the matrix element now. 
     
      ! indeces     
      a = stors%mat(qq)%qn(1)%Y(AA,1)
      b = stors%mat(qq)%qn(1)%Y(AA,2)      
      c = stors%mat(qq)%qn(1)%Y(BB,1)
      d = stors%mat(qq)%qn(1)%Y(BB,2)      
      
      ! i think the scaling and COM subtraction have already been done
      ! I HOpe. 

!=========================================================================
      ! start the classical method of sorting these into my arrays now
!=========================================================================     
     C1 = jbas%con(a)+jbas%con(b) + 1 !ph nature
     C2 = jbas%con(c)+jbas%con(d) + 1
    
     qx = C1*C2
     qx = qx + adjust_index(qx)   !Vpppp nature  

     ! get the indeces in the correct order
     pre = 1.d0
     if ( a > b )  then 
        
        x = bosonic_tp_index(b,a,Ntot) 
        j_min = H%xmap(x)%Z(1)  
        i1 = H%xmap(x)%Z( (JT-j_min)/2 + 2) 
        pre = pre * (-1.)**( 1 + (jbas%jj(a) + jbas%jj(b) -JT)/2 ) 
     
     else
        if (a == b) pre = pre / sqrt( 2.d0 )
        
        x = bosonic_tp_index(a,b,Ntot) 
        
        j_min = H%xmap(x)%Z(1)  
        i1 = H%xmap(x)%Z( (JT-j_min)/2 + 2) 
     end if
  
     if (c > d)  then     
        
        x = bosonic_tp_index(d,c,Ntot) 
        j_min = H%xmap(x)%Z(1)  
        i2 = H%xmap(x)%Z( (JT-j_min)/2 + 2) 
        
        pre = pre * (-1.)**( 1 + (jbas%jj(c) + jbas%jj(d) -JT)/2 ) 
    
     else 
        if (c == d) pre = pre / sqrt( 2.d0 )
       
        x = bosonic_tp_index(c,d,Ntot) 
        j_min = H%xmap(x)%Z(1)  
        i2 = H%xmap(x)%Z( (JT-j_min)/2 + 2) 
     end if

     ! kets/bras are pre-scaled by sqrt(2) if they 
     ! have two particles in the same sp-shell
        
     ! get the units right. I hope 
 

     if ((qx == 1) .or. (qx == 5) .or. (qx == 4)) then 
        H%mat(q)%gam(qx)%X(i2,i1)  = V *pre
        H%mat(q)%gam(qx)%X(i1,i2)  = V *pre
        
        if (rr_calc) then 
           STOP 'fuck, this is not implemented yet' 
           rr%mat(q)%gam(qx)%X(i2,i1)  = hw*g2*pre/(H%Aneut + H%Aprot)
           rr%mat(q)%gam(qx)%X(i1,i2)  = hw*g2*pre/(H%Aneut + H%Aprot)
        end if 

        if (pp_calc) then 
           STOP 'fuck, this is not implemented yet' 
           pp%mat(q)%gam(qx)%X(i2,i1)  = hw*g3*pre/(H%Aneut + H%Aprot)
           pp%mat(q)%gam(qx)%X(i1,i2)  = hw*g3*pre/(H%Aneut + H%Aprot)
        end if 

        
     else if (C1>C2) then
        H%mat(q)%gam(qx)%X(i2,i1)  = V *pre
        
        if (rr_calc) then 
           STOP 'fuck, this is not implemented yet' 
           rr%mat(q)%gam(qx)%X(i2,i1)  = hw*g2*pre/(H%Aneut + H%Aprot) 
        end if
        
        if (pp_calc) then 
           STOP 'fuck, this is not implemented yet' 
           pp%mat(q)%gam(qx)%X(i2,i1)  = hw*g3*pre/(H%Aneut + H%Aprot) 
        end if

     else
        H%mat(q)%gam(qx)%X(i1,i2) = V * pre
        
        if (rr_calc) then 
           STOP 'fuck, this is not implemented yet' 
           rr%mat(q)%gam(qx)%X(i1,i2)  = hw*g2*pre/(H%Aneut + H%Aprot) 
        end if

        if (pp_calc) then 
           STOP 'fuck, this is not implemented yet' 
           pp%mat(q)%gam(qx)%X(i1,i2)  = hw*g3*pre/(H%Aneut + H%Aprot) 
        end if

     end if 
     ! I shouldn't have to worry about hermiticity here, input is assumed to be hermitian

            if (AA == stors%mat(qq)%ntot) then 
               if (BB == stors%mat(qq)%ntot) then 
                  exit
               end if 
            end if

        end do 
      end do   ! end loops over conserved quantities
   end do 
end do 
      

      
! i guess we are done with the two body piece

14 me1bfile = intfile(1:len(trim(intfile))-5)//'1b.gz'
  
  if ( present(Lawson) ) then 
     hndle=gzOpen(trim(TBME_DIR)//"O16_Hcm_eMax10_hwHO020.ham0.me1b.gz"//achar(0),"r"//achar(0)) 
  else 
     hndle=gzOpen(trim(TBME_DIR)//trim(adjustl(me1bfile))//achar(0),"r"//achar(0)) 
  end if 

!    print*, trim(TBME_DIR)//trim(adjustl(me1bfile))//achar(0)
  sz=200

  ! read verion line, and then some integer
  buf=gzGets(hndle,buffer,sz) 
  ! the integer probably has to do with the file size
  buf=gzGets(hndle,buffer,sz) 
 
  
  read(buffer(1:4),'(I4)',iostat=qq) aMax 
  if (qq .ne. 0 ) then 
     read(buffer(1:3),'(I3)',iostat=qq) aMax 
     if (qq .ne. 0 ) then 
        read(buffer(1:2),'(I2)',iostat=qq) aMax
     end if 
  end if
  sz = 20
  ! I assume this is the zero body piece right here
  buf=gzGets(hndle,buffer,sz) 
 
  ! lets get it right since it can have up to 4 digits before the decimal
  if (buffer(4:4) == '.') then 
     read(buffer(1:10),'(f10.6)') H%E0 
  else if (buffer(5:5) == '.') then 
     read(buffer(1:11),'(f11.6)') H%E0 
  else if (buffer(6:6) == '.') then 
     read(buffer(1:12),'(f12.6)') H%E0
  else 
     print*, 'what the fuck is going on with the me1b file? ' 
  end if 
  
! now lets read the 1 body piece
  

sz=200

do a= 1, aMax
  

    buf=gzGets(hndle,buffer,sz) 
    
    read(buffer(2:2),'(I1)') t1
    read(buffer(4:5),'(I2)') lj
    read(buffer(8:9),'(I2)') n1
    read(buffer(11:12),'(I2)') n2
      
    t1 = (-2*t1+ 1)  
    
    lj = lj + 1
    l1 = SPBljs(lj,1)/2
    j1 = SPBljs(lj,2)
      
          
    read(buffer(15:28),'(f14.10)') V
  

    ! V is now the one body matrix element
             ! now search for sp labels
                  do i = 1, jbas%total_orbits 
                     if ( jbas%jj(i) .ne. j1 ) cycle
                     if ( jbas%nn(i) .ne. n1 ) cycle
                     if ( jbas%ll(i) .ne. l1 ) cycle
                     if ( jbas%itzp(i) .ne. t1 ) cycle                     
                     exit
                  end do 
                  
                  do j = 1, jbas%total_orbits 
                     if ( jbas%jj(j) .ne. j1 ) cycle
                     if ( jbas%nn(j) .ne. n2 ) cycle
                     if ( jbas%ll(j) .ne. l1 ) cycle
                     if ( jbas%itzp(j) .ne. t1 ) cycle                     
                     exit
                  end do
                  
                  ! okay now I have my indeces 
           
                 
                  if (jbas%con(i) + jbas%con(j) == 2 ) then 
                     
                     !fhh 
                     
                     H%fhh( jbas%holesb4(i)+1 , jbas%holesb4(j)+1 ) = V 
                     
                  else if (jbas%con(i) + jbas%con(j) == 0 ) then 
                     
                     !fpp 
                     
                     H%fpp( jbas%partsb4(i)+1 , jbas%partsb4(j)+1 ) = V 
                  
                  else if ((jbas%con(i)==0) .and. (jbas%con(j) == 1) ) then 
                     !fph
                     H%fph( jbas%partsb4(i)+1 , jbas%holesb4(j)+1 ) = V
                  end if 
                  
           
   end do 

   rx = gzClose(hndle)

end subroutine
!==========================================================
subroutine read_me3j(store_3b,jbas) 
  use three_body_routines
  implicit none 
  
  type(spd) :: jbas
  type(three_body_force) :: store_3b
  integer :: nlj1,nlj2,nlj3,nnlj1,nnlj2,nnlj3,aux,aux1,aux2,aux3,aux4
  integer :: nnlj2_end,nnlj3_end,twoTMin,twoTmax,twoJCMin,twoJCMax
  integer :: twoJCMindown,twoJCMaxup,twoJCMindownket,twoJCMindownbra
  integer :: twoJCMaxupket,twoJCMaxupbra,la,lb,lc,ld,le,lf
  integer :: nsp,ja,jb,jc,jd,je,jf,iblock,Jab,JJab,Tab,TTab,ttot,jtot
  integer :: ea,eb,ec,ed,ef,ee,e1max,E3max,JabMax,JabMin,JJabMax,JJabMin
  integer :: Tij_indx,Tlm_indx,blocksize,endpoint,endsz,endsz2
  integer :: lmax3,jtot_max,jtot_max_1,jtot_max_2,jtot_min,jtot_min_1
  integer :: jtot_min_2,i,II,JJ,Jab_max,Jab_min,jc_max,jc_min
  integer :: Jde_min,Jde_max,x1,x2,q,NN,nsp_iso,tc_min,tc_max
  integer :: spot_in_gz_file,iMax,PAR,Tab_indx,TTab_indx,r,w
  real(8) :: szofblock,V
  character(1)::rem
  real(8),allocatable,dimension(:) :: xblock,me_fromfile
  type(c_ptr) :: buf,buf2,buf3
  integer(c_int) :: hndle,hndle2,hndle3,sz,sz2,sz3,rx
  character(kind=C_CHAR,len=200) :: buffer,buffer2,buffer3
  logical :: autozero ,thing
  
  E3max = 12 
  iMax = store_3b%num_elems 
  allocate(me_fromfile(10))
  nsp = jbas%total_orbits
  nsp_iso = nsp/2
  !open file 
  hndle=gzOpen(trim(TBME_DIR)//trim(adjustl(threebody_file))//achar(0),"r"//achar(0)) 
  sz = 200.
  ! read first line, it's not important to me.
  buf = gzGets(hndle,buffer(1:sz),sz)
  rem = '9'
  endsz=130
  
  i = 1
  r = 0 
  w = 0
  spot_in_gz_file = 0 
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
                                
                                q = block_index_3b(jtot,ttot,PAR)                                 
                                
                                if (spot_in_gz_file == 0) then 
                                   ! we need to read a line
                                   if (i+10 > iMax) then 
                                      deallocate(me_fromfile)
                                      allocate(me_fromfile( iMax - i + 1) ) 
                                      endpoint = iMax-i + 1
                                      endsz = 13+(endpoint-1)*13 
                                      write(rem,'(I1)') endpoint-1
                                   end if
                                   buf = gzGets(hndle,buffer(1:sz),sz)
                                   !print*, buffer(1:endsz)
                                   read(buffer(1:endsz),'(f12.7,'//rem//'(f13.7))') me_fromfile 
                                end if 
                                
                                call SetME(Jab,JJab,jtot,2*Tab,2*TTab,ttot,2*nlj1,2*nlj2,2*nlj3,2*nnlj1&
                                     ,2*nnlj2,2*nnlj3,me_fromfile(spot_in_gz_file+1),store_3b,jbas)
                                
                                i = i + 1
                                spot_in_gz_file = mod(spot_in_gz_file + 1, 10) 
                                
                                 
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
               
  rx = gzClose(hndle)  
end subroutine



















end module



