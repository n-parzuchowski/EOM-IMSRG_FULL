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
              
        TDA%map(q) = r
        allocate(TDA%blkM(q)%matrix(r,r)) 
        allocate(TDA%blkM(q)%extra(10*r))
        allocate(TDA%blkM(q)%eigval(r)) 
        allocate(TDA%blkM(q)%labels(r,2)) 
        allocate(TDA%blkM(q)%states(r))
        TDA%blkM(q)%states = 0.d0
        
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
end module


  
  
  
  
