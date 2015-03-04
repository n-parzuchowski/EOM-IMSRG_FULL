module generators 
  use basic_IMSRG
  implicit none 
  
  
contains

!==========================================================
!==========================================================
subroutine build_gs_wegner(H,ETA,jbas,HCC,HODCC,WCC,w1,w2) 
  use commutators
  ! calculates the traditional white generator for
  ! ground state decoupling
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: H,ETA,HOD,w1,w2
  type(cross_coupled_31_mat) :: HCC,HODCC,WCC
  integer :: a,b,i,j,ji,ja,ti,ta,li,la,JT,TZ,PAR
  integer :: q,IX,JX,jj,jb,lb,lj,tb,tj,ik,ak
  real(8) :: Eden,sm,Javerage
  
  ETA%herm = -1

  call duplicate_sq_op(H,HOD) 
  
  HOD%fph = H%fph

  do q = 1,HOD%nblocks
     HOD%mat(q)%gam(3)%X = H%mat(q)%gam(3)%X
  end do 

  call calculate_cross_coupled(HOD,HODCC,jbas,.true.)
  call calculate_cross_coupled(H,HCC,jbas,.false.) 
 
  call commutator_111(H,HOD,ETA,jbas) 
  call commutator_121(H,HOD,ETA,jbas)
  call commutator_122(H,HOD,ETA,jbas)    

  call commutator_222_pp_hh(H,HOD,ETA,w1,w2,jbas)
  
  call commutator_221(H,HOD,ETA,w1,w2,jbas)
  call commutator_222_ph(HCC,HODCC,ETA,WCC,jbas)

end subroutine
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
subroutine build_gs_imtime(H,ETA,jbas) 
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
        
        ETA%fph(a,i) = H%fph(a,i)*sign(1.d0,Eden)*abs(Eden)**.0001
        
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
           
           ETA%mat(q)%gam(3)%X(IX,JX) = H%mat(q)%gam(3)%X(IX,JX)* &
                sign(1.d0,Eden)*abs(Eden)**.0001
           
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
 
  ETA%fpp = 0.d0 
   do a = 1,H%Nsp - H%belowEF
     ak = jbas%parts(a)
     if (ak < H%valcut+1) cycle
     ja = jbas%jj(ak)
     la = jbas%ll(ak)
     ta = jbas%itzp(ak)
     
     do i = a+1, H%Nsp - H%belowEF
       
        ik = jbas%parts(i)
        if (ik > H%valcut) cycle
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
        
         
        ETA%fpp(a,i) = H%fpp(a,i)*sign(1.d0,Eden)*abs(Eden)**.0001 
        ETA%fpp(i,a) = -1*ETA%fpp(a,i)
     end do 
  end do 
  
  ! two body part 
  
  do  q = 1, H%nblocks
         
     ETA%mat(q)%gam(2)%X = 0.d0
     ETA%mat(q)%gam(6)%X = 0.d0
     ETA%mat(q)%gam(4)%X = 0.d0 
     
     ! Vppph 
     do IX = 1,H%mat(q)%npp 

        ! figure out which sp states compose IX
        a = H%mat(q)%qn(1)%Y(IX,1)  ! pp descriptor qn(1)%Y 
        b = H%mat(q)%qn(1)%Y(IX,2)

        ja = jbas%jj(a)
        jb = jbas%jj(b)   
        
        do JX = 1,H%mat(q)%nph 

           i = H%mat(q)%qn(2)%Y(JX,1) !hh descriptor qn(3)%Y
           j = H%mat(q)%qn(2)%Y(JX,2)

           if (i > H%valcut) cycle
           if (j > H%valcut) cycle
      
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
           
         
           ETA%mat(q)%gam(2)%X(IX,JX) = &
                H%mat(q)%gam(2)%X(IX,JX)*sign(1.d0,Eden)*abs(Eden)**.0001  
         
           
        end do 
     end do
     
     ! Vphph 
     do IX = 1,H%mat(q)%nph 

        ! figure out which sp states compose IX
        a = H%mat(q)%qn(2)%Y(IX,1)  ! pp descriptor qn(1)%Y 
        b = H%mat(q)%qn(2)%Y(IX,2)

        if ( a > H%valcut )  cycle 
        if ( b > H%valcut )  cycle
             
        ja = jbas%jj(a)
        jb = jbas%jj(b)   
        
        do JX = IX + 1,H%mat(q)%nph 

           i = H%mat(q)%qn(2)%Y(JX,1) !hh descriptor qn(3)%Y
           j = H%mat(q)%qn(2)%Y(JX,2)
 
           ji = jbas%jj(i)
           jj = jbas%jj(j)
         
           !hl = i*jbas%con(i) + j * jbas%con(j) 
           !jh =  ji*jbas%con(i) + jj * jbas%con(j)
            
           
           p = i*(1-jbas%con(i)) + j * (1-jbas%con(j)) 
           jp =  ji*(1-jbas%con(i)) + jj * (1-jbas%con(j))
           
           if (p < H%valcut+1) cycle
           Eden = 0.d0 
           
           ! constructing the App'hh' term is rather codey... 
          
           !pp'pp' 

           !Eden = Eden + Javerage(a,b,ja,jb,H,jbas) 
           !Eden = Eden - Javerage(a,hl,ja,jh,H,jbas) 
           !Eden = Eden - Javerage(b,hl,jb,jh,H,jbas) 
                   
           Eden = Eden + f_elem(a,a,H,jbas) + f_elem(b,b,H,jbas)  - &
                f_elem(i,i,H,jbas) - f_elem(j,j,H,jbas) 
           
         
           ETA%mat(q)%gam(4)%X(IX,JX) = &
                H%mat(q)%gam(4)%X(IX,JX)*sign(1.d0,Eden)*abs(Eden)**.0001  
           ETA%mat(q)%gam(4)%X(JX,IX) = -1*ETA%mat(q)%gam(4)%X(IX,JX)
           
        end do 
     end do
     
     !Vphhh is decoupled in full...  
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
subroutine build_specific_space(H,ETA,jbas) 
  ! calculates imaginary time generator for a specific space 
  ! defined by the TDA operator for specified qnums
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: H,ETA 
  integer :: a,b,i,j,ji,ja,ti,ta,li,la,JT,TZ,PAR,hspf,pspf,sz 
  integer :: q,IX,JX,jj,jb,lb,lj,tb,tj,ik,ak,p,hl,jp,jh,pos
  real(8) :: Eden,sm,Javerage
  logical :: in, inSD
 
  
  ETA%herm = -1 ! anti-hermitian operator
  ETA%fph = 0.d0
  sz = size(H%exlabels(:,1))
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
 
  ETA%fpp = 0.d0 
   do a = 1,H%Nsp - H%belowEF
     ak = jbas%parts(a)
     if (ak < H%valcut+1) cycle
     ja = jbas%jj(ak)
     la = jbas%ll(ak)
     ta = jbas%itzp(ak)
     
     do i = a+1, H%Nsp - H%belowEF
       
        ik = jbas%parts(i)        
        if (.not. in(ik,H%exlabels(:,2),pos,sz)) cycle
        
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
        
         
        ETA%fpp(a,i) = H%fpp(a,i)*sign(1.d0,Eden)*abs(Eden)**.0001 
        ETA%fpp(i,a) = -1*ETA%fpp(a,i)
     end do 
  end do 
  
  ! two body part 
  
  do  q = 1, H%nblocks
         
     ETA%mat(q)%gam(2)%X = 0.d0
     ETA%mat(q)%gam(6)%X = 0.d0
     ETA%mat(q)%gam(4)%X = 0.d0 
     
     ! Vppph 
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

           p= i*(1-jbas%con(i)) + j * (1-jbas%con(j)) 
           jp =  ji*(1-jbas%con(i)) + jj * (1-jbas%con(j))
           
           if (.not. in(p,H%exlabels(:,2),pos,sz) ) cycle 

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
     ! Vphph 
     do IX = 1,H%mat(q)%nph 

        ! figure out which sp states compose IX
        a = H%mat(q)%qn(2)%Y(IX,1)  ! pp descriptor qn(1)%Y 
        b = H%mat(q)%qn(2)%Y(IX,2)

        !if ( a > H%valcut )  cycle 
        !if ( b > H%valcut )  cycle
             
        ja = jbas%jj(a)
        jb = jbas%jj(b)   
                
        p= a*(1-jbas%con(i)) + b * (1-jbas%con(j)) 
        jp =  ja*(1-jbas%con(i)) + jb * (1-jbas%con(j))
        
        do JX = IX + 1,H%mat(q)%nph 

           i = H%mat(q)%qn(2)%Y(JX,1) !hh descriptor qn(3)%Y
           j = H%mat(q)%qn(2)%Y(JX,2)
 
           ji = jbas%jj(i)
           jj = jbas%jj(j)
         
           hl = i*jbas%con(i) + j * jbas%con(j) 
           jh =  ji*jbas%con(i) + jj * jbas%con(j)
           
           if (.not. inSD(hl,p,H%exlabels,pos,sz)) cycle
           
           p = i*(1-jbas%con(i)) + j * (1-jbas%con(j)) 
           jp =  ji*(1-jbas%con(i)) + jj * (1-jbas%con(j))
           
           if (p < H%valcut+1) cycle
           Eden = 0.d0 
           
           ! constructing the App'hh' term is rather codey... 
          
           !pp'pp' 

           !Eden = Eden + Javerage(a,b,ja,jb,H,jbas) 
           !Eden = Eden - Javerage(a,hl,ja,jh,H,jbas) 
           !Eden = Eden - Javerage(b,hl,jb,jh,H,jbas) 
                   
           Eden = Eden + f_elem(a,a,H,jbas) + f_elem(b,b,H,jbas)  - &
                f_elem(i,i,H,jbas) - f_elem(j,j,H,jbas) 
           
         
           ETA%mat(q)%gam(4)%X(IX,JX) = &
                H%mat(q)%gam(4)%X(IX,JX)*sign(1.d0,Eden)*abs(Eden)**.0001  
           ETA%mat(q)%gam(4)%X(JX,IX) = -1*ETA%mat(q)%gam(4)%X(IX,JX)
           
        end do 
     end do
  
     !Vphhh  
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
           
           hl = i*(jbas%con(i)) + j * (jbas%con(j)) 
           jh =  ji*(jbas%con(i)) + jj *(jbas%con(j))
                  
           
           if (.not. in(hl,H%exlabels(:,1),pos,sz)) cycle
           
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
     23 continue
  end do 

end subroutine 
!==========================================================
!==========================================================
subroutine build_valence_decouple(H,ETA,jbas) 
  ! calculates the traditional white generator for
  ! ground state decoupling
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: H,ETA 
  integer :: a,b,i,j,ji,ja,ti,ta,li,la,JT,TZ,PAR
  integer :: q,IX,JX,jj,jb,lb,lj,tb,tj,ik,ak,p,hl,jp,jh
  real(8) :: Eden,sm,Javerage
  
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
  

ETA%fpp = 0.d0
! F_QV

  do a = 1,H%Nsp - H%belowEF
     ak = jbas%parts(a)
     if (ak < H%valcut+1 ) cycle
     ja = jbas%jj(ak)
     la = jbas%ll(ak)
     ta = jbas%itzp(ak)
     
     do i = a,H%Nsp - H%belowEF
       
        ik = jbas%parts(i)
        if (ik > H%valcut )  cycle 
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
        
         
        ETA%fpp(a,i) = H%fpp(a,i)*sign(1.d0,Eden)*abs(Eden)**.0001 
        ETA%fpp(i,a) = H%fpp(a,i)*(-1)
     end do 
  end do



  ! two body part 
  
  do  q = 1, H%nblocks
         
     ETA%mat(q)%gam(2)%X = 0.d0
     ETA%mat(q)%gam(1)%X = 0.d0
     
     do IX = 1,H%mat(q)%npp 

        ! figure out which sp states compose IX
        a = H%mat(q)%qn(1)%Y(IX,1)  ! pp descriptor qn(1)%Y 
        b = H%mat(q)%qn(1)%Y(IX,2)

        ja = jbas%jj(a)
        jb = jbas%jj(b)   
        
        do JX = 1,H%mat(q)%nph 

           i = H%mat(q)%qn(2)%Y(JX,1) !hh descriptor qn(3)%Y
           j = H%mat(q)%qn(2)%Y(JX,2)

           if (i > H%valcut) cycle
           if (j > H%valcut) cycle
            
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
           
         
           ETA%mat(q)%gam(2)%X(IX,JX) = &
                H%mat(q)%gam(2)%X(IX,JX)*sign(1.d0,Eden)*abs(Eden)**.0001  
         
           
        end do 
     end do
     
     do IX = 1,H%mat(q)%npp

        ! figure out which sp states compose IX
        a = H%mat(q)%qn(1)%Y(IX,1)  ! pp descriptor qn(1)%Y 
        b = H%mat(q)%qn(1)%Y(IX,2)

        if ( (a < H%valcut+1) .and. (b < H%valcut+1) )  cycle 
        
        ja = jbas%jj(a)
        jb = jbas%jj(b)   
        
        
        do JX = IX,H%mat(q)%npp 

           i = H%mat(q)%qn(1)%Y(JX,1) !hh descriptor qn(3)%Y
           j = H%mat(q)%qn(1)%Y(JX,2)
           
           if (i > H%valcut) cycle
           if (j > H%valcut) cycle 
           
           ji = jbas%jj(i)
           jj = jbas%jj(j)
         
           p = i*(1-jbas%con(i)) + j * (1-jbas%con(j)) 
           jp =  ji*(1-jbas%con(i)) + jj *(1- jbas%con(j))
                  
           Eden = 0.d0 
           
           ! constructing the App'hh' term is rather codey... 
          
           !pp'pp' 
           Eden = Eden - f_elem(a,a,H,jbas) - f_elem(b,b,H,jbas)  + &
                f_elem(i,i,H,jbas) + f_elem(j,j,H,jbas) 
           
         
           ETA%mat(q)%gam(1)%X(IX,JX) = &
                H%mat(q)%gam(1)%X(IX,JX) * sign(1.d0,Eden)*abs(Eden)**.0001  
           ETA%mat(Q)%gam(1)%X(JX,IX) = -1.d0 * ETA%mat(q)%gam(1)%X(IX,JX) 

        end do 
     end do 
     
  end do 

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
!=====================================================
!=====================================================
logical function in(element,list,position,sz) 
  implicit none 
  
  integer :: sz
  integer,dimension(sz) :: list
  integer :: element,position,i
  logical :: dum
 
  dum = .false.
  position = -1 
  do i = 1,size(list)
     if ( list(i) == element ) then 
        dum = .true. 
        exit
        position = i 
     end if 
  end do 

  in = dum 
end function
!=====================================================
!=====================================================
logical function inSD(i1,i2,list,position,sz) 
  implicit none 
  
  integer :: sz
  integer,dimension(sz,2) :: list
  integer :: position,i,i1,i2
  logical :: dum
  
  dum = .false. 
  position = -1 
  do i = 1,size(list(:,1)) 
     if ( list(i,1) == i1) then 
        if ( list(i,2) == i2 ) then 
           position = i 
           dum = .true. 
           exit
        end if 
     end if 
  end do 
  
  inSD = dum 
end function
