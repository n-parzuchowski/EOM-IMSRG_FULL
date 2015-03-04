module heiko_format
  use basic_IMSRG
  implicit none
  
  contains

subroutine get_heiko_spfile(eMaxchr)
  ! This constructs the sps file for
  ! heiko's format, after it's been converted to a pn-basis 
  implicit none 
  
  integer :: e,eMax,l,jj,tz,q,n
  character(2) :: eMaxchr
  character(13) :: fmt
  
  read(eMaxchr,'(I2)') eMax
  eMaxchr = adjustl(eMaxchr) 
  open(unit=24,file='hk'//trim(eMaxchr)//'.sps')
  
  q = 1
  fmt = '(5(I5),e17.7)'
  do  e = 0,eMax
     
     do l = mod(e,2),e,2
        
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
end subroutine


subroutine read_heiko_interaction(H,jbas,htype,hw) 

  integer :: nlj1,nlj2,nnlj1,nnlj2,j,T,Mt,nljMax,endpoint,j_min,j_max,htype
  integer :: l1,l2,ll1,ll2,j1,j2,jj1,jj2,Ntot,i,q,bospairs,qx,ta,tb,tc,td
  integer :: eMax,iMax,jmax,jmin,JT,a,b,c,d,C1,C2,i1,i2,pre,COM,x,PAR
  integer,allocatable,dimension(:) :: indx 
  real(8),allocatable,dimension(:) :: ME,MEpp,me_fromfile,ppff,rrff
  real(8) :: V,g1,g2,g3,hw,pre2
  type(spd) :: jbas 
  type(sq_op) :: H,pp,rr
  logical :: pp_calc,rr_calc
  
  pp_calc=.false.
  rr_calc=.false.
  COM = 1


  Ntot = jbas%total_orbits/2
  
  i = 0
  q = 0
  ! allocate array to store positions of matrix elements
  bospairs = bosonic_tp_index(Ntot,Ntot,Ntot) 
  allocate(indx(bospairs**2)) 
  indx = 0

  ! move in increments of two, because I use a pn basis,
  !  heikos is isospin coupled (half the states) 
  
  eMax = 4
  
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
  
  if (pp_calc) allocate(rrff(10)) 
  
  open(unit=25,file='../../TBME_input/chi2b_srg0200_eMax04_hwHO024.me2j') 
  open(unit=24,file='../../TBME_input/tpp_eMax04.me2j') 
  read(25,*) ! first line is garbage
  read(24,*) 
  
  endpoint = 10 
 
  do i = 1,iMax,10
     
     if (i+10 > iMax) then 
        deallocate(me_fromfile)
        deallocate(ppff) 
        allocate(me_fromfile( iMax - i - 1) ) 
        allocate(ppff(iMax-i-1)) 
        endpoint = iMax-i-1
     end if 
     
     read(25,*) me_fromfile
     read(24,*) ppff 
    
     do j = 1,endpoint 
        ME(i+j-1) = me_fromfile(j)
        MEpp(i+j-1) = ppff(j) 
     end do 
     
  end do
  
  close(25);close(24)
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
            
            ta = jbas%itzp(a)
            tb = jbas%itzp(b)
            tc = jbas%itzp(c)
            td = jbas%itzp(d)
            
            T = ta+tb
            if (tc+td .ne. T) cycle
            T = T/2
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


     ! center of mass subtraction
     V = (V + g3*2.d0/hbarc2_over_mc2*COM*hw/(H%Aneut+H%Aprot)) *pre2 
    
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

end module
