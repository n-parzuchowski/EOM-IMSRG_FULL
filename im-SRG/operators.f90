module operators
  use basic_IMSRG
  implicit none 
  

contains
!===================================================================
!===================================================================
subroutine initialize_TDA(TDA,jbas,Jtarget,PARtarget,cut) 
  ! figure out how big the TDA matrix has to be
  ! allocate it
  ! uses full_sp_block_mat, just because it's got the right shape
  implicit none 
  
  type(spd) :: jbas
  type(full_sp_block_mat) :: TDA
  integer :: JT,PI,Jmax,q,i,a,r,nh,np,ix,ax
  integer :: Jtarget, PARtarget,cut
 
  Jmax = jbas%Jtotal_max
  nh = sum(jbas%con) 
  np = jbas%total_orbits - nh 
  TDA%blocks = 1    ! (Jmax+1)*2
  !allocate(TDA%blkM((Jmax+1)*2)) 
  !allocate(TDA%map((Jmax+1)*2))
  allocate(TDA%blkM(1))
  allocate(TDA%map(1)) 
  
  q = 1
  
  !do PI = 0,1
   !  do JT = 0,2*Jmax,2
    
  PI = PARtarget 
  JT = Jtarget 
      
        r = 0
        ! count the states
        do ix = 1,nh
           do ax = 1,np 
              i =jbas%holes(ix) 
              a =jbas%parts(ax)
              if (a > cut) cycle 
              if (jbas%itzp(i) .ne. jbas%itzp(a) ) cycle ! cannot change from p to n
            
              if (.not. (triangle(jbas%jj(i),jbas%jj(a),JT))) cycle 
              if (mod(jbas%ll(i)+jbas%ll(a),2) .ne. PI )  cycle ! technically l_a - l_i  
              
              r = r + 1 
           end do 
        end do 
              
        ! allocate
        TDA%map(q) = r
        allocate(TDA%blkM(q)%matrix(r,r)) 
        allocate(TDA%blkM(q)%extra(10*r))
        allocate(TDA%blkM(q)%eigval(r)) 
        allocate(TDA%blkM(q)%labels(r,2)) 
        allocate(TDA%blkM(q)%states(r))
        TDA%blkM(q)%states = 0.d0
        
        ! fill the states
        r = 0
        do ix = 1,nh
           do ax = 1,np 
              i =jbas%holes(ix) 
              a =jbas%parts(ax)              
              if (a > cut) cycle 
              if (jbas%itzp(i) .ne. jbas%itzp(a) ) cycle ! cannot change from p to n
            
              if (.not. (triangle(jbas%jj(i),jbas%jj(a),JT))) cycle 
              if (mod(jbas%ll(i)+jbas%ll(a),2) .ne. PI )  cycle ! technically l_a - l_i 
              r = r + 1 
              TDA%blkM(q)%labels(r,1) = i 
              TDA%blkM(q)%labels(r,2) = a 
           end do 
        end do 
        TDA%blkM(q)%lmda(1) = JT
        TDA%blkM(q)%lmda(2) = PI 
        TDA%blkM(q)%lmda(3) = 0 
 !       q = q + 1
 !    end do 
 ! end do 

      ! print*, TDA%blkM(1)%labels(:,1)
      ! print*, TDA%blkM(1)%labels(:,2)
end subroutine 
!==========================================
!==========================================
subroutine calc_TDA(TDA,HS,HSCC,jbas) 
  implicit none 
  
  type(spd) :: jbas
  type(full_sp_block_mat) :: TDA
  type(sq_op) :: HS 
  type(cross_coupled_31_mat) :: HSCC
  integer :: q,JT,r1,r2,a,b,i,j,x,g,NBindx,Rindx,q1,Tz,PAR,JTM
  
  JTM = jbas%Jtotal_max
  do q = 1, TDA%blocks
     JT = TDA%blkM(q)%lmda(1) 
     Tz = 0
     PAR = TDA%blkM(q)%lmda(2)
  
     q1 = JT/2+1 + Tz*(JTM+1) + 2*PAR*(JTM+1)
     
     do r1 = 1, TDA%map(q) 
        i = TDA%blkM(q)%labels(r1,1)
        a = TDA%blkM(q)%labels(r1,2)
        
        ! get CC index for this pair
        x = CCindex(a,i,HS%Nsp)
        g = 1
        do while (HSCC%qmap(x)%Z(g) .ne. q1 )
           g = g + 1
        end do

        NBindx = HSCC%nbmap(x)%Z(g) 
           
        do r2 = r1, TDA%map(q)
           j = TDA%blkM(q)%labels(r2,1)
           b = TDA%blkM(q)%labels(r2,2)
           
           ! get CCindex for this pair
           x = CCindex(j,b,HS%Nsp) 
           g = 1
           do while (HSCC%qmap(x)%Z(g) .ne. q1 )
              g = g + 1
           end do
              
           Rindx = HSCC%rmap(x)%Z(g)  
           
           ! one body piece
           TDA%blkM(q)%matrix(r1,r2) = f_elem(a,b,HS,jbas) * &
                kron_del(i,j) - f_elem(j,i,HS,jbas) * kron_del(a,b)
           
           ! two body piece
           TDA%blkM(q)%matrix(r1,r2) = TDA%blkM(q)%matrix(r1,r2) - &
                HSCC%CCR(q1)%X(NBindx,Rindx) * &
                (-1) ** ((JT)/2) /sqrt(JT+1.d0) 
           ! I have to account for the fact that the CCME are scaled by 
           ! Sqrt[2J+1] * (-1)^(j_i + j_a) 
        
           ! hermiticity 
           TDA%blkM(q)%matrix(r2,r1) = TDA%blkM(q)%matrix(r1,r2) 
        end do 
     end do 
  end do 
        
end subroutine 
!===================================================
!===================================================
subroutine TDA_expectation_value(TDA_HAM,TDA_OP) 
  ! takes expectation values of TDA_OP for the eigenvectors defined in TDA_HAM
  implicit none 
  
  type(full_sp_block_mat) :: TDA_HAM, TDA_OP,AUX  
  integer :: q,dm,i 
  
  call duplicate_sp_mat(TDA_HAM,AUX) 
 
  do q = 1, TDA_OP%blocks
        
     dm = TDA_OP%map(q)
     if (dm == 0)  cycle
     
     ! transform operator into basis defined by TDA vectors
     call dgemm('N','N',dm,dm,dm,al,TDA_OP%blkM(q)%matrix&
          ,dm,TDA_HAM%blkM(q)%matrix,dm,bet,AUX%blkM(q)%matrix,dm) 
     call dgemm('T','N',dm,dm,dm,al,TDA_HAM%blkM(q)%matrix&
          ,dm,AUX%blkM(q)%matrix,dm,bet,TDA_OP%blkM(q)%matrix,dm)
     
     ! diagonal elements are the expectation values
     do i = 1, dm 
        TDA_OP%blkM(q)%eigval(i) = TDA_OP%blkM(q)%matrix(i,i) 
     end do    
     
  end do

end subroutine 
!=========================================================
subroutine initialize_rms_radius(rms,rr,jbas)
  implicit none 
  
  type(sq_op) :: rr, rms 
  type(spd) :: jbas
  real(8) :: mass_factor 
  integer :: q,i
  
  mass_factor = 1.d0-1.d0/dfloat(rr%Aprot + rr%Aneut) 
  
  call calculate_h0_harm_osc(1.d0,jbas,rms,5)
  
  ! multiply by scale factors to make it into r^2 instead of u_ho 
  rms%fhh = rms%fhh * hbarc2_over_mc2 * 2.d0 * mass_factor / rms%hospace
  rms%fpp = rms%fpp * hbarc2_over_mc2 * 2.d0 * mass_factor / rms%hospace
  rms%fph = rms%fph * hbarc2_over_mc2 * 2.d0 * mass_factor / rms%hospace
  
  
  do q = 1, rms%nblocks
     do i = 1,6
        
        rms%mat(q)%gam(i)%X = -2.d0 * hbarc2_over_mc2 &
             / rr%hospace**2 / dfloat(rr%Aneut + rr%Aprot) * &
             rr%mat(q)%gam(i)%X 
        
     end do 
  end do 


end subroutine 
!=========================================================
subroutine initialize_CM_radius(rms,rr,jbas)
  implicit none 
  
  type(sq_op) :: rr, rms 
  type(spd) :: jbas
  real(8) :: mass_factor 
  integer :: q,i
  
  mass_factor = 1.d0/dfloat(rr%Aprot + rr%Aneut) 
  
  call calculate_h0_harm_osc(1.d0,jbas,rms,5)
  
  ! multiply by scale factors to make it into r^2 instead of u_ho 
  rms%fhh = rms%fhh * hbarc2_over_mc2 * 2.d0 * mass_factor / rms%hospace
  rms%fpp = rms%fpp * hbarc2_over_mc2 * 2.d0 * mass_factor / rms%hospace
  rms%fph = rms%fph * hbarc2_over_mc2 * 2.d0 * mass_factor / rms%hospace
  
  
  do q = 1, rms%nblocks
     do i = 1,6
        
        rms%mat(q)%gam(i)%X = 2.d0 * hbarc2_over_mc2 &
             / rr%hospace**2 * mass_factor * &
             rr%mat(q)%gam(i)%X 
        
     end do 
  end do 


end subroutine
!====================================================================
!=========================================================
subroutine initialize_CM_radius_onebody(rms,rr,jbas)
  implicit none 
  
  type(sq_op) :: rr, rms 
  type(spd) :: jbas
  real(8) :: mass_factor, elem
  integer :: q,a,b,ak,bk,na,nb,la,lb
  
  mass_factor = 1.d0/dfloat(rr%Aprot + rr%Aneut) 
  
  
  ! multiply by scale factors to make it into r^2 instead of u_ho 
  
  do a = 1, rr%belowEF
     
     ak = jbas%holes(a)
     na = jbas%nn(ak)
     la = jbas%ll(ak) 
     
     do b = a, rr%belowEF
     
        bk = jbas%holes(b)
        nb = jbas%nn(bk)
        lb = jbas%ll(bk)

  
        if (la == lb - 1) then 
           
           if (na == nb) then 
              
              elem = sqrt(nb + lb + 0.5d0)
           
           else if (na == nb + 1)  then 
              
              elem = -1.d0*sqrt(nb + 1.d0)

           else 
              elem = 0.d0 
              
           end if 
        
        else if (la == lb + 1) then 
           
           if (na == nb) then 
              
              elem = sqrt(nb + lb + 1.5d0)
           
           else if (na == nb - 1) then 
              
              elem = -1.d0*sqrt(float(nb)) 
           
           else 
              
              elem = 0.d0 
           end if 
        else 
           elem = 0.d0 
        end if 
  
        rms%fhh(a,b) = sqrt(hbarc2_over_mc2 / rms%hospace) * elem * mass_factor
        rms%fhh(b,a) = rms%fhh(a,b) 
     end do  
  end do 

!fph  
  do a = 1, rr%Nsp - rr%belowEF
     
     ak = jbas%parts(a)
     na = jbas%nn(ak)
     la = jbas%ll(ak) 
     
     do b = 1, rr%belowEF
     
        bk = jbas%holes(b)
        nb = jbas%nn(bk)
        lb = jbas%ll(bk)

  
        if (la == lb - 1) then 
           
           if (na == nb) then 
              
              elem = sqrt(nb + lb + 0.5d0)
           
           else if (na == nb + 1)  then 
              
              elem = -1.d0*sqrt(nb + 1.d0)

           else 
              elem = 0.d0 
              
           end if 
        
        else if (la == lb + 1) then 
           
           if (na == nb) then 
              
              elem = sqrt(nb + lb + 1.5d0)
           
           else if (na == nb - 1) then 
              
              elem = -1.d0*sqrt(float(nb)) 
           
           else 
              
              elem = 0.d0
           end if 
        else 
           elem = 0.d0 
        end if 
  
        rms%fph(a,b) = sqrt(hbarc2_over_mc2 / rms%hospace) * elem * mass_factor
       
     end do  
  end do 
 
!fpp
   do a = 1,rr%Nsp-rr%belowEF
     
     ak = jbas%parts(a)
     na = jbas%nn(ak)
     la = jbas%ll(ak) 
     
     do b = a, rr%Nsp - rr%belowEF
     
        bk = jbas%parts(b)
        nb = jbas%nn(bk)
        lb = jbas%ll(bk)

  
        if (la == lb - 1) then 
           
           if (na == nb) then 
              
              elem = sqrt(nb + lb + 0.5d0)
           
           else if (na == nb + 1)  then 
              
              elem = -1.d0*sqrt(nb + 1.d0)

           else 
              elem = 0.d0 
              
           end if 
        
        else if (la == lb + 1) then 
           
           if (na == nb) then 
              
              elem = sqrt(nb + lb + 1.5d0)
           
           else if (na == nb - 1) then 
              
              elem = -1.d0*sqrt(float(nb)) 
           
           else 
              
              elem = 0.d0 
           end if
        else 
           elem = 0.d0 
        end if 
  
        rms%fpp(a,b) = sqrt(hbarc2_over_mc2 / rms%hospace) * elem * mass_factor
        rms%fpp(b,a) = rms%fpp(a,b) 
     end do  
  end do 
 
  

end subroutine
!==================================================================== 
subroutine calculate_CM_energy(pp,rr,hw) 
  implicit none 
  
  type(sq_op) :: pp,rr,Hcm
  real(8) :: hw,wTs(2),Ecm(3) 
  integer :: i
  character(200) ::  spfile,intfile,prefix
  common /files/ spfile,intfile,prefix
 
  call duplicate_sq_op(pp,Hcm)
  call add_sq_op(pp,1.d0,rr,1.d0,Hcm)
  Ecm(1) = Hcm%E0 - 1.5d0*hw ! store Ecm for this Hcm frequency 
    
  ! new frequencies
  wTs = optimum_omega_for_CM_hamiltonian(hw,Ecm(1)) 
  
  do i = 1, 2
     call add_sq_op(pp,1.d0,rr,wTs(i)**2/hw**2,Hcm)
     Ecm(i+1) = Hcm%E0 - 1.5d0 * wTs(i) 
  end do   

  !write results to file
  open(unit=42,file='../../output/'//&
         trim(adjustl(prefix))//'_Ecm.dat')
  write(42,'(6(e14.6))') hw, wTs, Ecm 
  close(42)
  
end subroutine
!=====================================================
!==================================================================== 
subroutine calculate_CM_energy_TDA(TDA,rr,pp,ppTDA,rrTDA,hw) 
  implicit none 
  
  type(sq_op) :: rr,pp,Hcm 
  type(full_sp_block_mat) :: TDA,HcmTDA,ppTDA,rrTDA
  real(8) :: hw,wTs(2)
  real(8),allocatable,dimension(:) :: energies,omegas 
  integer :: i,q
  character(200) ::  spfile,intfile,prefix
  character(3) :: args
  common /files/ spfile,intfile,prefix
 
  
  call duplicate_sp_mat(TDA,HcmTDA) 
  call duplicate_sq_op(rr,Hcm)
  allocate(energies(3*ppTDA%map(1))) 
  allocate(omegas(2*TDA%map(1)))
  
  call add_sp_mat(ppTDA,1.d0,rrTDA,1.d0,HcmTDA)
  call add_sq_op(pp,1.d0,rr,1.d0,Hcm) 
  call TDA_expectation_value(TDA,HcmTDA) 

  energies(1:TDA%map(1)) = HcmTDA%blkM(1)%eigval + Hcm%E0 - 1.5d0*hw
  
  do q = 1, TDA%map(1) 
     ! new frequencies
     wTs = optimum_omega_for_CM_hamiltonian(hw,energies(q)) 
     omegas(q) = wTs(1) 
     omegas(q+TDA%map(1)) = wTs(2)
     
     do i = 1, 2
        call add_sq_op(pp,1.d0,rr,wTs(i)**2/hw**2,Hcm) 
        call add_sp_mat(ppTDA,1.d0,rrTDA,wTs(i)**2/hw**2,HcmTDA)
        call TDA_expectation_value(TDA,HcmTDA) 
        energies(i*TDA%map(1)+q) = HcmTDA%blkM(1)%eigval(q) + Hcm%E0 - 1.5d0*wTs(i) 
     end do
  end do 
  
  ! for formatting
   i = 1+5*TDA%map(1) 
   write(args,'(I3)') i 
   args = adjustl(args) 

  !write results to file
  open(unit=42,file='../../output/'//&
         trim(adjustl(prefix))//'_Ecm_excited.dat')
  write(42,'('//trim(args)//'(e14.6))') hw, omegas, energies 
  close(42)
  
end subroutine
end module


  
  
  
  
