module generators 
  use basic_IMSRG
  implicit none 
  
  
contains
!==========================================================
!==========================================================
subroutine build_gs_white(H,ETA,jbas) 
  ! calculates the traditional white generator for
  ! ground state decoupling
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: H,ETA 
  integer :: a,b,i,j,ji,ja,ti,ta,li,la,JT,TZ,PAR
  integer :: q,IX,JX,jj,jb,lb,lj,tb,tj,ik,ak
  real(8) :: Eden,sm,Javerage
  
  ETA%herm = -1 ! anti-hermitian operator
  
  ! one body part
  do a = 1,H%Nsp - H%belowEF
     ak = jbas%parts(a) 
     ja = jbas%jj(ak)
     la = jbas%ll(ak)
     ta = jbas%itzp(ak)
     
     do i = 1,H%belowEF
       
        ik = jbas%holes(i)
        ji = jbas%jj(ik)       
        li = jbas%ll(ik)        
        ti = jbas%itzp(ik)
             
        ! the generator is zero if these are true: 
        if ( ji .ne. ja) cycle
        if ( li .ne. la) cycle
        if ( ti .ne. ta) cycle 
        
        PAR = 0 ! both states have same l
        TZ = ti ! both states have same tz
        
        ! energy denominator has a sum over J  to factor out m dep. 
        Eden = 0.0 
        
        do JT = 0, 2*ji , 2
           Eden = Eden + (JT + 1) * v_elem(ak,ik,ak,ik,JT,TZ,PAR,H,jbas) 
        end do 
        
        ! sum is averaged over ji ** 2  
        Eden = Eden / (ji + 1.d0)/(ji + 1.d0) 
        
        Eden = Eden + H%fpp(a,a) - H%fhh(i,i) 
        
        ETA%fph(a,i) = H%fph(a,i) / Eden 
        
     end do 
  end do 
  
  ! two body part 
  
  do  q = 1, H%nblocks
         
     do IX = 1,H%mat(q)%npp 

        ! figure out which sp states compose IX
        a = H%mat(q)%qn(1)%Y(IX,1)  ! pp descriptor qn(1)%Y 
        b = H%mat(q)%qn(1)%Y(IX,2)

        ja = jbas%jj(a)
        jb = jbas%jj(b)
        la = jbas%ll(a)
        lb = jbas%ll(b)
        ta = jbas%itzp(a)
        tb = jbas%itzp(b)
        
        do JX = 1,H%mat(q)%nhh 

           i = H%mat(q)%qn(3)%Y(JX,1) !hh descriptor qn(3)%Y
           j = H%mat(q)%qn(3)%Y(JX,2)

           ji = jbas%jj(i)
           jj = jbas%jj(j)
           li = jbas%ll(i)
           lj = jbas%ll(j)
           ti = jbas%itzp(i)
           tj = jbas%itzp(j)
           
           Eden = 0.d0 
           
           ! constructing the App'hh' term is rather codey... 
          
           !pp'pp' 

           Eden = Eden + Javerage(a,b,ja,jb,la,lb,ta,tb,H,jbas) 
           Eden = Eden + Javerage(i,j,ji,jj,li,lj,ti,tj,H,jbas) 
           Eden = Eden - Javerage(a,i,ja,ji,la,li,ta,ti,H,jbas) 
           Eden = Eden - Javerage(a,j,ja,jj,la,lj,ta,tj,H,jbas) 
           Eden = Eden - Javerage(i,b,ji,jb,li,lb,ti,tb,H,jbas) 
           Eden = Eden - Javerage(j,b,jj,jb,lj,lb,tj,tb,H,jbas) 
           
           Eden = Eden + f_elem(a,a,H,jbas) + f_elem(b,b,H,jbas)  - &
                f_elem(i,i,H,jbas) - f_elem(j,j,H,jbas) 
           
           ETA%mat(q)%gam(3)%X(IX,JX) = H%mat(q)%gam(3)%X(IX,JX)/Eden 
           
        end do 
     end do 
  end do 

end subroutine 
!==========================================================
!==========================================================
end module
!==========================================================
!==========================================================
real(8) function Javerage(a,b,ja,jb,la,lb,ta,tb,H,jbas) 
  ! average over J used a lot in white generator
  use basic_IMSRG
  implicit none 
  
  integer :: ja,jb,la,lb,ta,tb,PAR,TZ,JT,a,b
  type(sq_op) :: H 
  type(spd) :: jbas 
  real(8) :: sm
  
  PAR = mod(la+lb,2) 
  TZ = (ta + tb)/2 
           
  sm = 0.d0 
  do JT = abs(ja-jb),ja+jb,2
     sm = sm + (JT + 1) * v_elem(a,b,a,b,JT,TZ,PAR,H,jbas) 
  end do 
  
  Javerage = sm /(ja + 1.d0) / (jb + 1.d0) 
end function 
           
