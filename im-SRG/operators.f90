module operators
  use cross_coupled
  implicit none 
  

contains
subroutine initialize_transition_operator(trs_type,rank,Op,zr,jbas,tcalc)
  implicit none
  
  type(spd) :: jbas
  type(sq_op) :: Op,zr
  character(1) :: trs_type,ranklab
  integer :: rank 
  logical :: tcalc
  
  tcalc = .true. 
  Select case (trs_type) 
     case ('E') ! ELECTRIC TRANSITION
         Op%rank = 2*rank
         Op%herm = 1
         Op%dpar = abs((-1)**rank-1) 
         write(ranklab,'(I1)') rank 
         Op%trans_label = trs_type//ranklab         
         if (rank==0) then 
            print*, 'not implemented'
         else
            call allocate_tensor(jbas,Op,zr)
            Op%hospace=zr%hospace
            call calculate_EX(Op,jbas)
         end if 
     case ('M') ! MAGNETIC TRANSITION
        Op%rank = 2*rank
        Op%herm = 1
        Op%dpar = abs((-1)**(rank-1)-1) 
        write(ranklab,'(I1)') rank 
        Op%trans_label = trs_type//ranklab
        if (rank==0) then 
           print*, 'not implemented'
        else
           call allocate_tensor(jbas,Op,zr)
           Op%hospace=zr%hospace
           call calculate_MX(Op,jbas)
        end if
     case ('F') ! FERMI TRANSITION
        stop 'not implemented'
     case ('G') ! GAMOW-TELLER TRANSITION
        stop 'not implemented'
     case default
        tcalc = .false. 
  end select
end subroutine initialize_transition_operator
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
  TDA%blocks = 1   
  
  if (allocated(TDA%blkM)) then 
     deallocate(TDA%blkM,TDA%map) 
  end if 
  
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
              if (i < 10) cycle
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
              if (i < 10) cycle
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

       print*, TDA%blkM(1)%labels(:,1)
       print*, TDA%blkM(1)%labels(:,2)
end subroutine 
!==========================================
!==========================================
subroutine calc_TDA(TDA,HS,HSCC,jbas) 
  implicit none 
  
  type(spd) :: jbas
  type(full_sp_block_mat) :: TDA
  type(sq_op) :: HS 
  type(cc_mat) :: HSCC
  integer :: q,JT,r1,r2,a,b,i,j,x,g,NBindx,Rindx,q1,Tz,PAR,JTM
  integer :: ji,jj,ja,jb,phase
  
  JTM = jbas%Jtotal_max
  do q = 1, TDA%blocks
     JT = TDA%blkM(q)%lmda(1) 
     Tz = 0
     PAR = TDA%blkM(q)%lmda(2)
  
     q1 = JT/2+1 + Tz*(JTM+1) + 2*PAR*(JTM+1)
     
     do r1 = 1, TDA%map(q) 
        i = TDA%blkM(q)%labels(r1,1)
        a = TDA%blkM(q)%labels(r1,2)
        
        ji = jbas%jj(i)
        ja = jbas%jj(a) 
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
           
           jj = jbas%jj(j)
           jb = jbas%jj(b)
           phase = (-1)**((jj+jb)/2)
           ! get CCindex for this pair
           x = CCindex(b,j,HS%Nsp) 
           g = 1
           do while (HSCC%qmap(x)%Z(g) .ne. q1 )
              g = g + 1
           end do
              
           Rindx = HSCC%rmap(x)%Z(g)  
           
           ! one body piece
           TDA%blkM(q)%matrix(r1,r2) = f_elem(a,b,HS,jbas) * &
                kron_del(i,j) - f_elem(j,i,HS,jbas) * kron_del(a,b)
           
           ! two body piece
!           TDA%blkM(q)%matrix(r1,r2) = TDA%blkM(q)%matrix(r1,r2) - &
 !               Vcc(i,a,b,j,JT,HS,jbas)*(-1)**((ja+ji)/2) * &
  !              (-1) ** ((JT)/2) /sqrt(JT+1.d0) 
           
           TDA%blkM(q)%matrix(r1,r2) = TDA%blkM(q)%matrix(r1,r2) - &
                HSCC%CCX(q1)%X(Rindx,NBindx) * HSCC%herm * & !need scaling.
                (-1) ** (JT/2) /sqrt(JT+1.d0) * phase  
         
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
  open(unit=42,file=trim(OUTPUT_DIR)//&
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
  character(3) :: args
 
  
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
  open(unit=42,file=trim(OUTPUT_DIR)//&
         trim(adjustl(prefix))//'_Ecm_excited.dat')
  write(42,'('//trim(args)//'(e14.6))') hw, omegas, energies 
  close(42)
  
end subroutine calculate_CM_energy_TDA

subroutine calculate_EX(op,jbas) 
  ! calculates electromagnetic transition operator 
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: op   
  integer :: a,b,ax,bx,rank,na,nb
  integer :: ta,tb,tc,td,la,lb,lc,ld,ja,jb,jc,jd
  real(8) :: pol,charge(2),dcgi,dcgi00,x,dcg,hw,holength

  pol = 0.0d0   ! between -1 and 0 
  hw = op%hospace
  holength = sqrt(hbarc2_over_mc2/hw) 
  x = dcgi00()
  ! right now the units of this operator are: e fm^(X)
  ! so the charge is not multiplied by the fundamental charge. 
  charge(1) = (1+pol) 
  charge(2) = pol 
  
  ! access the correct charge with: 
  ! charge( (jbas%tz(i) + 1)/2 + 1 )  
  
  rank = op%rank
  
  do ax = 1, op%belowEF 
     do bx = 1, op%belowEF 
 
        a = jbas%holes(ax)
        b = jbas%holes(bx) 
        
        ta = jbas%itzp(a) 
        tb = jbas%itzp(b) 
       
        if (ta .ne. tb) cycle
        
        la = jbas%ll(a) 
        lb = jbas%ll(b) 
       
        if ( mod(la+lb+rank/2,2) == 1) cycle 
        
        ja = jbas%jj(a) 
        jb = jbas%jj(b) 
               
        if (.not.(triangle(ja,jb,rank))) cycle
        
        na = jbas%nn(a)
        nb = jbas%nn(b) 

        op%fhh(ax,bx) = charge((ta + 1)/2+1)/sqrt(4.d0*Pi_const)* &
             (-1) **((ja + rank - 1)/2) * &
             sqrt((ja +1.d0) * (jb+1.d0)) * dcgi(ja,1,jb,-1,rank,0) *&
             RabLAM(na,la,nb,lb,RANK/2)*holength**(rank/2)  
     end do 
  end do 
        
  do ax = 1, op%nsp-op%belowEF 
     do bx = 1, op%nsp-op%belowEF 
  
        a = jbas%parts(ax)
        b = jbas%parts(bx) 
        
        ta = jbas%itzp(a) 
        tb = jbas%itzp(b) 
        
        if (ta .ne. tb) cycle
        
        la = jbas%ll(a) 
        lb = jbas%ll(b) 
        
        if ( mod(la+lb+rank/2,2) == 1) cycle 
        
        ja = jbas%jj(a) 
        jb = jbas%jj(b) 
        
        if (.not.(triangle(ja,jb,rank))) cycle
        
        na = jbas%nn(a)
        nb = jbas%nn(b) 
        
        op%fpp(ax,bx) = charge((ta + 1)/2+1)/sqrt(4.d0*Pi_const)* &
             (-1) **((ja + rank - 1)/2) * &
             sqrt((ja +1.d0) * (jb+1.d0)) * dcgi(ja,1,jb,-1,rank,0) *&
             RabLAM(na,la,nb,lb,RANK/2)*holength**(rank/2)  
     end do 
  end do     

  do ax = 1, op%nsp-op%belowEF 
     do bx = 1, op%belowEF 
  
        a = jbas%parts(ax)
        b = jbas%holes(bx) 
        
        ta = jbas%itzp(a) 
        tb = jbas%itzp(b) 
        
        if (ta .ne. tb) cycle
        
        la = jbas%ll(a) 
        lb = jbas%ll(b) 
        
        if ( mod(la+lb+rank/2,2) == 1) cycle 
        
        ja = jbas%jj(a) 
        jb = jbas%jj(b) 
        
        if (.not.(triangle(ja,jb,rank))) cycle
        
        na = jbas%nn(a)
        nb = jbas%nn(b) 
        
        op%fph(ax,bx) = charge((ta + 1)/2+1)/sqrt(4.d0*Pi_const)* &
             (-1) **((ja + rank - 1)/2) * &
             sqrt((ja +1.d0) * (jb+1.d0)) * dcgi(ja,1,jb,-1,rank,0) *&
             RabLAM(na,la,nb,lb,RANK/2)*holength**(rank/2)  
     end do
  end do     
        
end subroutine calculate_EX
!==================================================================================
!==================================================================================
subroutine calculate_MX(op,jbas) 
  ! calculates electromagnetic transition operator 
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: op   
  integer :: a,b,ax,bx,rank,na,nb,kappa
  integer :: ta,tb,tc,td,la,lb,lc,ld,ja,jb,jc,jd
  real(8) :: pol,charge(2),dcgi,dcgi00,x,dcg,hw,holength
  real(8) :: mu_N_overC,gl(2),gs(2)
  

  pol = 0.0d0   ! between -1 and 0 
  hw = op%hospace
  holength = sqrt(hbarc2_over_mc2/hw) 
  mu_N_overc = 1.d0!sqrt(hbarc2_over_mc2/4.d0) 
  x = dcgi00()
  ! right now the units of this operator are: e fm^(X)
  ! so the charge is not multiplied by the fundamental charge. 
  charge(1) = (1+pol) 
  charge(2) = pol 
  gl(1) = 1
  gl(2) = 0
  gs(1) = 5.586
  gs(2) = -3.826
  ! access the correct charge with: 
  ! charge( (jbas%tz(i) + 1)/2 + 1 )  
  
  rank = op%rank
  
  do ax = 1, op%belowEF 
     do bx = 1, op%belowEF 
 
        a = jbas%holes(ax)
        b = jbas%holes(bx) 
        
        ta = jbas%itzp(a) 
        tb = jbas%itzp(b) 
       
        if (ta .ne. tb) cycle
        
        la = jbas%ll(a) 
        lb = jbas%ll(b) 

        if ( mod(la+lb+rank/2,2) == 0) cycle 
        
        ja = jbas%jj(a) 
        jb = jbas%jj(b) 
               
        if (.not.(triangle(ja,jb,rank))) cycle
        
        na = jbas%nn(a)
        nb = jbas%nn(b) 
        
        kappa = (-1) ** ((ja+1)/2+la)*(ja+1) &
             +(-1) ** ((jb+1)/2+lb)*(jb+1) 
        
        op%fhh(ax,bx) = mu_N_overC*charge(1) &
             /sqrt(4.d0*Pi_const)*(-1) **((ja + rank - 1)/2) * &
             sqrt((ja +1.d0) * (jb+1.d0)) * dcgi(ja,1,jb,-1,rank,0) *&
             RabLAM(na,la,nb,lb,RANK/2-1)*holength**(rank/2-1)*&
             (rank-kappa)/2*(gl((ta+1)/2+1)*(1+kappa/(rank+2.d0))-gs((ta+1)/2+1)/2.d0) 
     end do 
  end do 
        
  do ax = 1, op%nsp-op%belowEF 
     do bx = 1, op%nsp-op%belowEF 
  
        a = jbas%parts(ax)
        b = jbas%parts(bx) 
        
        ta = jbas%itzp(a) 
        tb = jbas%itzp(b) 
        
        if (ta .ne. tb) cycle
        
        la = jbas%ll(a) 
        lb = jbas%ll(b) 

        if ( mod(la+lb+rank/2,2) == 0) cycle 
        
        ja = jbas%jj(a) 
        jb = jbas%jj(b) 
        
        if (.not.(triangle(ja,jb,rank))) cycle
        
        na = jbas%nn(a)
        nb = jbas%nn(b) 

        kappa = (-1) ** ((ja+1)/2+la)*(ja+1) &
             +(-1) ** ((jb+1)/2+lb)*(jb+1) 
        
        op%fpp(ax,bx) = mu_N_overC*charge(1) &
             /sqrt(4.d0*Pi_const)*(-1) **((ja + rank - 1)/2) * &
             sqrt((ja +1.d0) * (jb+1.d0)) * dcgi(ja,1,jb,-1,rank,0) *&
             RabLAM(na,la,nb,lb,RANK/2-1)*holength**(rank/2-1)*&
             (rank-kappa)/2*(gl((ta+1)/2+1)*(1+kappa/(rank+2.d0))-gs((ta+1)/2+1)/2.d0) 
       
     end do 
  end do     

  do ax = 1, op%nsp-op%belowEF 
     do bx = 1, op%belowEF 
  
        a = jbas%parts(ax)
        b = jbas%holes(bx) 
        
        ta = jbas%itzp(a) 
        tb = jbas%itzp(b) 
        
        if (ta .ne. tb) cycle
        
        la = jbas%ll(a) 
        lb = jbas%ll(b) 
!        if (la .ne. lb) cycle
        if ( mod(la+lb+rank/2,2) == 0) cycle 
        
        ja = jbas%jj(a) 
        jb = jbas%jj(b) 
        
        if (.not.(triangle(ja,jb,rank))) cycle
        
        na = jbas%nn(a)
        nb = jbas%nn(b) 

        kappa = (-1) ** ((ja+1)/2+la)*(ja+1) &
             +(-1) ** ((jb+1)/2+lb)*(jb+1) 

        op%fph(ax,bx) = mu_N_overC*charge(1) &
             /sqrt(4.d0*Pi_const)*(-1) **((ja + rank - 1)/2) * &
             sqrt((ja +1.d0) * (jb+1.d0)) * dcgi(ja,1,jb,-1,rank,0) *&
             RabLAM(na,la,nb,lb,RANK/2-1)*holength**(rank/2-1)*&
             (rank-kappa)/2*(gl((ta+1)/2+1)*(1+kappa/(rank+2.d0))-gs((ta+1)/2+1)/2.d0) 
        
     end do
  end do     
        
end subroutine calculate_MX
  
real(8) function transition_to_ground_ME( Trans_op , Qdag,jbas ) 
  implicit none
  
  type(spd) :: jbas
  type(sq_op) :: Trans_op,Qdag 
  integer :: ja,jb,ji,jj,rank,Nsp,Abody
  integer :: a,b,i,j,ax,ix,jx,bx,J1,J2
  real(8) :: sm , phase
  
  sm = 0.d0 
  
  Nsp = jbas%total_orbits 
  Abody = sum(jbas%con) 
   
  rank = Trans_op%rank 
  
  if (rank .ne. Qdag%rank)  then 
     transition_to_ground_ME = 0.d0 
     return
  end if 
  
  do ax = 1,Nsp-Abody
     a = jbas%parts(ax)
     ja = jbas%jj(a) 

     do ix = 1,Abody 
        i = jbas%holes(ix)
        ji = jbas%jj(i) 
        
        phase = (-1) ** ((ja-ji)/2)
        sm = sm + f_tensor_elem(i,a,Trans_op,jbas)*&
             f_tensor_elem(a,i,Qdag,jbas)*phase
     end do 
  end do 

  do ax = 1,Nsp-Abody
     a = jbas%parts(ax)
     ja = jbas%jj(a) 

     do bx = 1,Nsp-Abody
        b = jbas%parts(bx)
        jb = jbas%jj(b) 
        
        do ix = 1,Abody 
           i = jbas%holes(ix)
           ji = jbas%jj(i) 
           
           do jx = 1,Abody 
              j = jbas%holes(jx)
              jj = jbas%jj(j) 
  
              do J1 = abs(ji-jj),ji+jj,2
                 do J2 = abs(ja-jb),ja+jb,2
                    
                    phase = (-1) **((J1+J2)/2) 
                    sm = sm + phase*0.25d0*&
                         tensor_elem(i,j,a,b,J1,J2,Trans_op,jbas)*&
                         tensor_elem(a,b,i,j,J2,J1,Qdag,jbas) 
                    
                 end do 
              end do 
           end do
        end do
     end do
  end do
        
  transition_to_ground_ME = sm * (-1.d0)**(rank/2)
  !  BY SUHONEN'S DEFINITION, I SHOULD BE DEVIDING BY Sqrt(2J+1) 
  !  BUT THE LANCZOS ALGORITHM DID THAT FOR US ALREADY. 
  
end function transition_to_ground_ME

subroutine construct_number_operator(op,H,jbas) 
  implicit none
  
  type(spd) :: jbas
  type(sq_op) :: op,H
  integer :: i,ix,ji
  
  call duplicate_sq_op(H,op) 
  
  do ix = 1, H%belowEF
     i= jbas%holes(ix) 
     ji = jbas%jj(i)
     op%fhh(ix,ix) = (ji+1.d0)
  end do 
  
  do ix = 1,jbas%total_orbits-H%belowEF
     i= jbas%parts(ix) 
     ji = jbas%jj(i)
     op%fpp(ix,ix) = (ji+1.d0)
  end do 

  op%E0 = sum(op%fhh)
end subroutine

real(8) function RabLAM(na,la,nb,lb,LAM) 
  ! radial integrals in suhonen chapter 6. 
  implicit none 
  
  integer,intent(in) :: na,la,nb,lb,LAM
  integer :: tau_A , tau_b,sig,sig_min,sig_max
  real(8) :: coef,sm,xxx
  
  if (mod(la+lb+LAM,2).ne.0) then 
     RabLAM = 0.d0 
     return 
  end if

  tau_a = max((lb - la + LAM)/2,0)
  tau_b = max((la - lb + LAM)/2,0)
  sig_min = max(0,na-tau_a,nb-tau_b)
  sig_max = min(na,nb) 

  coef = sqrt(factorial(na)/half_int_gamma(na+la+1))
  coef = coef*sqrt(factorial(nb)/half_int_gamma(nb+lb+1))
  coef = coef * ( -1) ** (na+nb)*factorial(tau_a)*factorial(tau_b) 

  sm = 0.d0
  do sig = sig_min,sig_max      
     xxx = (half_int_gamma((la+lb+lam)/2+sig+1)/factorial(sig))
     xxx = xxx / factorial(na-sig)
     xxx = xxx / factorial(nb-sig)
     xxx = xxx / factorial(sig+tau_a-na)
     xxx = xxx / factorial(sig+tau_b-nb)
     sm = sm + xxx
  end do

  RabLAM = coef*sm
end function

! These functions are often an overflow threat, but 
! here I don't think I need to worry, because emax= 2n+l doesn't get much
! bigger than 15, and these functions definitely work well into the 20s. 

real(8) function factorial(N) 
  ! N!
  implicit none
  integer,intent(in) :: N
  integer :: i
  real(8) :: sm 
  
  sm = 1
  
  do i = 1,N
     sm = sm * i 
  end do 
  
  factorial = sm
end function 

real(8) function half_int_gamma(N)
  !GAMMA(N+1/2) 
  implicit none
  
  integer,intent(in) :: N
  integer :: i
  real(8) :: sm

  sm = sqrt(acos(-1.d0))
  
  do i = 1,N 
     sm = (i-0.5)*sm
  end do 

  half_int_gamma = sm
end function   
  
  


end module
  
  
  
  
