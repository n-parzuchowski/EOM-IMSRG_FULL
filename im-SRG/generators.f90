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
     
        ! energy denominator has a sum over J  to factor out m dep. 
        Eden = 0.0 
        
        do JT = 0, 2*ji , 2
           Eden = Eden - (JT + 1) * v_elem(ak,ik,ak,ik,JT,H,jbas) 
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
        
        do JX = 1,H%mat(q)%nhh 

           i = H%mat(q)%qn(3)%Y(JX,1) !hh descriptor qn(3)%Y
           j = H%mat(q)%qn(3)%Y(JX,2)

           ji = jbas%jj(i)
           jj = jbas%jj(j)
         
           
           Eden = 0.d0 
           
           ! constructing the App'hh' term is rather codey... 
          
           !pp'pp' 

           Eden = Eden + Javerage(a,b,ja,jb,H,jbas) 
           Eden = Eden + Javerage(i,j,ji,jj,H,jbas) 
           Eden = Eden - Javerage(a,i,ja,ji,H,jbas) 
           Eden = Eden - Javerage(a,j,ja,jj,H,jbas) 
           Eden = Eden - Javerage(i,b,ji,jb,H,jbas) 
           Eden = Eden - Javerage(j,b,jj,jb,H,jbas) 
           
           Eden = Eden + f_elem(a,a,H,jbas) + f_elem(b,b,H,jbas)  - &
                f_elem(i,i,H,jbas) - f_elem(j,j,H,jbas) 
           
           ETA%mat(q)%gam(3)%X(IX,JX) = H%mat(q)%gam(3)%X(IX,JX)/Eden 
           
        end do 
     end do 
  end do 

end subroutine 
!==========================================================
!==========================================================
subroutine build_ex_white(H,ETA,jbas) 
  ! calculates the traditional white generator for
  ! ground state decoupling
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: H,ETA 
  integer :: a,b,i,j,ji,ja,ti,ta,li,la,JT,TZ,PAR
  integer :: q,IX,JX,jj,jb,lb,lj,tb,tj,ik,ak,p,hl,jp,jh
  real(8) :: Eden,sm,Javerage
  real(8),parameter :: dcut = .2
  
  ETA%herm = -1 ! anti-hermitian operator
  ETA%fph = 0.d0
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
     
        ! energy denominator has a sum over J  to factor out m dep. 
        Eden = 0.0 
        
        do JT = 0, 2*ji , 2
           Eden = Eden - (JT + 1) * v_elem(ak,ik,ak,ik,JT,H,jbas) 
        end do 
        
        ! sum is averaged over ji ** 2  
        Eden = Eden / (ji + 1.d0)/(ji + 1.d0) 
        
        Eden = Eden + H%fpp(a,a) - H%fhh(i,i) 
        
         
        if (abs(Eden) > dcut) ETA%fph(a,i) = H%fph(a,i) / Eden
        
     end do 
  end do 
  
  ! two body part 
  
  do  q = 1, H%nblocks
         
     ETA%mat(q)%gam(2)%X = 0.d0
     ETA%mat(q)%gam(6)%X = 0.d0
     
     do IX = 1,H%mat(q)%npp 

        ! figure out which sp states compose IX
        a = H%mat(q)%qn(1)%Y(IX,1)  ! pp descriptor qn(1)%Y 
        b = H%mat(q)%qn(1)%Y(IX,2)

        ja = jbas%jj(a)
        jb = jbas%jj(b)   
        
        do JX = 1,H%mat(q)%nph 

           i = H%mat(q)%qn(2)%Y(JX,1) !hh descriptor qn(3)%Y
           j = H%mat(q)%qn(2)%Y(JX,2)

           ji = jbas%jj(i)
           jj = jbas%jj(j)
         
           hl = i*jbas%con(i) + j * jbas%con(j) 
           jh =  ji*jbas%con(i) + jj * jbas%con(j)
                  
           Eden = 0.d0 
           
           ! constructing the App'hh' term is rather codey... 
          
           !pp'pp' 

           Eden = Eden + Javerage(a,b,ja,jb,H,jbas) 
           Eden = Eden - Javerage(a,hl,ja,jh,H,jbas) 
           Eden = Eden - Javerage(b,hl,jb,jh,H,jbas) 
                   
           Eden = Eden + f_elem(a,a,H,jbas) + f_elem(b,b,H,jbas)  - &
                f_elem(i,i,H,jbas) - f_elem(j,j,H,jbas) 
           
           if (abs(Eden) > dcut)  then 
              ETA%mat(q)%gam(2)%X(IX,JX) = H%mat(q)%gam(2)%X(IX,JX)/Eden 
           end if 
           
        end do 
     end do
     
     do IX = 1,H%mat(q)%nhh 

        ! figure out which sp states compose IX
        a = H%mat(q)%qn(3)%Y(IX,1)  ! pp descriptor qn(1)%Y 
        b = H%mat(q)%qn(3)%Y(IX,2)

        ja = jbas%jj(a)
        jb = jbas%jj(b)   
        
        do JX = 1,H%mat(q)%nph 

           i = H%mat(q)%qn(2)%Y(JX,1) !hh descriptor qn(3)%Y
           j = H%mat(q)%qn(2)%Y(JX,2)

           ji = jbas%jj(i)
           jj = jbas%jj(j)
         
           p = i*(1-jbas%con(i)) + j * (1-jbas%con(j)) 
           jp =  ji*(1-jbas%con(i)) + jj *(1- jbas%con(j))
                  
           Eden = 0.d0 
           
           ! constructing the App'hh' term is rather codey... 
          
           !pp'pp' 

           Eden = Eden + Javerage(a,b,ja,jb,H,jbas) 
           Eden = Eden - Javerage(a,p,ja,jp,H,jbas) 
           Eden = Eden - Javerage(b,p,jb,jp,H,jbas) 
                   
           Eden = Eden - f_elem(a,a,H,jbas) - f_elem(b,b,H,jbas)  + &
                f_elem(i,i,H,jbas) + f_elem(j,j,H,jbas) 
           
           if (abs(Eden) > dcut) then  
              ETA%mat(q)%gam(6)%X(JX,IX) = H%mat(q)%gam(6)%X(JX,IX)/Eden 
           end if 

        end do 
     end do 
     
  end do 

end subroutine
!==========================================================
!==========================================================
subroutine build_ex_imtime(H,ETA,jbas) 
  ! calculates the traditional white generator for
  ! ground state decoupling
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: H,ETA 
  integer :: a,b,i,j,ji,ja,ti,ta,li,la,JT,TZ,PAR
  integer :: q,IX,JX,jj,jb,lb,lj,tb,tj,ik,ak,p,hl,jp,jh
  real(8) :: Eden,sm,Javerage
  real(8),parameter :: dcut = .2
  
  ETA%herm = -1 ! anti-hermitian operator
  ETA%fph = 0.d0
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
     
        ! energy denominator has a sum over J  to factor out m dep. 
        Eden = 0.0 
        
        do JT = 0, 2*ji , 2
           Eden = Eden - (JT + 1) * v_elem(ak,ik,ak,ik,JT,H,jbas) 
        end do 
        
        ! sum is averaged over ji ** 2  
        Eden = Eden / (ji + 1.d0)/(ji + 1.d0) 
        
        Eden = Eden + H%fpp(a,a) - H%fhh(i,i) 
        
         
        ETA%fph(a,i) = H%fph(a,i)*sign(1.d0,Eden)*abs(Eden)**.0001 
        
     end do 
  end do 
  
  ! two body part 
  
  do  q = 1, H%nblocks
         
     ETA%mat(q)%gam(2)%X = 0.d0
     ETA%mat(q)%gam(6)%X = 0.d0
     
     do IX = 1,H%mat(q)%npp 

        ! figure out which sp states compose IX
        a = H%mat(q)%qn(1)%Y(IX,1)  ! pp descriptor qn(1)%Y 
        b = H%mat(q)%qn(1)%Y(IX,2)

        ja = jbas%jj(a)
        jb = jbas%jj(b)   
        
        do JX = 1,H%mat(q)%nph 

           i = H%mat(q)%qn(2)%Y(JX,1) !hh descriptor qn(3)%Y
           j = H%mat(q)%qn(2)%Y(JX,2)

           if (i > 12) cycle
           if (j > 12) cycle
           if (i < 3) cycle 
           if (j < 3) cycle 
           ji = jbas%jj(i)
           jj = jbas%jj(j)
         
           hl = i*jbas%con(i) + j * jbas%con(j) 
           !if (hl < 3) cycle
           jh =  ji*jbas%con(i) + jj * jbas%con(j)
                  
           Eden = 0.d0 
           
           ! constructing the App'hh' term is rather codey... 
          
           !pp'pp' 

           Eden = Eden + Javerage(a,b,ja,jb,H,jbas) 
           Eden = Eden - Javerage(a,hl,ja,jh,H,jbas) 
           Eden = Eden - Javerage(b,hl,jb,jh,H,jbas) 
                   
           Eden = Eden + f_elem(a,a,H,jbas) + f_elem(b,b,H,jbas)  - &
                f_elem(i,i,H,jbas) - f_elem(j,j,H,jbas) 
           
         
           ETA%mat(q)%gam(2)%X(IX,JX) = &
                H%mat(q)%gam(2)%X(IX,JX)*sign(1.d0,Eden)*abs(Eden)**.0001  
         
           
        end do 
     end do
     
     do IX = 1,H%mat(q)%nhh 

        ! figure out which sp states compose IX
        a = H%mat(q)%qn(3)%Y(IX,1)  ! pp descriptor qn(1)%Y 
        b = H%mat(q)%qn(3)%Y(IX,2)

        ja = jbas%jj(a)
        jb = jbas%jj(b)   
        
        do JX = 1,H%mat(q)%nph 

           i = H%mat(q)%qn(2)%Y(JX,1) !hh descriptor qn(3)%Y
           j = H%mat(q)%qn(2)%Y(JX,2)

           ji = jbas%jj(i)
           jj = jbas%jj(j)
         
           p = i*(1-jbas%con(i)) + j * (1-jbas%con(j)) 
           jp =  ji*(1-jbas%con(i)) + jj *(1- jbas%con(j))
                  
           Eden = 0.d0 
           
           ! constructing the App'hh' term is rather codey... 
          
           !pp'pp' 

           Eden = Eden + Javerage(a,b,ja,jb,H,jbas) 
           Eden = Eden - Javerage(a,p,ja,jp,H,jbas) 
           Eden = Eden - Javerage(b,p,jb,jp,H,jbas) 
                   
           Eden = Eden - f_elem(a,a,H,jbas) - f_elem(b,b,H,jbas)  + &
                f_elem(i,i,H,jbas) + f_elem(j,j,H,jbas) 
           
         
           ETA%mat(q)%gam(6)%X(JX,IX) = &
                H%mat(q)%gam(6)%X(JX,IX) * sign(1.d0,Eden)*abs(Eden)**.0001  
         

        end do 
     end do 
     
  end do 

end subroutine 
!==========================================================
!==========================================================
subroutine build_wegner(H,ETA,jbas) 
  ! calculates the traditional white generator for
  ! ground state decoupling
  use commutators
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: H,ETA,HD,w1,w2
  type(cross_coupled_31_mat) :: WCC,HDCC,HCC
  integer :: a,b,i,j,ji,ja,ti,ta,li,la,JT,TZ,PAR
  integer :: q,IX,JX,jj,jb,lb,lj,tb,tj,ik,ak
  real(8) :: Eden,sm,Javerage
  
  ETA%herm = -1 ! anti-hermitian operator
  
  call duplicate_sq_op(H,HD) 
  call duplicate_sq_op(H,w1) !workspace
  call duplicate_sq_op(H,w2) !workspace
  call allocate_CCMAT(H,HCC,jbas) ! cross coupled ME
  call duplicate_CCMAT(HCC,HDCC) !cross coupled ME
  call allocate_CC_wkspc(HCC,WCC) ! workspace for CCME
  
  HD%fhh = H%fhh
  HD%fpp = H%fpp
  
  do q = 1, H%nblocks
     HD%mat(q)%gam(1)%X = H%mat(q)%gam(1)%X
     !HD%mat(q)%gam(2)%X = H%mat(q)%gam(2)%X
     HD%mat(q)%gam(4)%X = H%mat(q)%gam(4)%X
     HD%mat(q)%gam(5)%X = H%mat(q)%gam(5)%X
     !HD%mat(q)%gam(6)%X = H%mat(q)%gam(6)%X
  end do 

  call calculate_cross_coupled(H,HCC,jbas,.true.)
  call calculate_cross_coupled(HD,HDCC,jbas,.false.) 
  
  call commutator_111(HD,H,ETA,jbas) 
  call commutator_121(HD,H,ETA,jbas)
  call commutator_122(HD,H,ETA,jbas)
  
  call commutator_222_pp_hh(HD,H,ETA,w1,w2,jbas)
  
  call commutator_221(HD,H,ETA,w1,w2,jbas)
  call commutator_222_ph(HDCC,HCC,ETA,WCC,jbas)

end subroutine 
!==========================================================
!==========================================================
end module
!==========================================================
!==========================================================
real(8) function Javerage(a,b,ja,jb,H,jbas) 
  ! average over J used a lot in white generator
  use basic_IMSRG
  implicit none 
  
  integer :: ja,jb,JT,a,b
  type(sq_op) :: H 
  type(spd) :: jbas 
  real(8) :: sm
            
  sm = 0.d0 
  do JT = abs(ja-jb),ja+jb,2
     sm = sm + (JT + 1) * v_elem(a,b,a,b,JT,H,jbas) 
  end do 
  
  Javerage = sm /(ja + 1.d0) / (jb + 1.d0) 
end function 
           
