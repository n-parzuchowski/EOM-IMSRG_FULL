module HF_mod
  use basic_IMSRG 
  ! J-scheme Hartree-Fockery. 
  implicit none 
  
contains
!====================================================
subroutine calc_HF( H ,jbas )
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: H 
  type(full_sp_block_mat) :: T,F,Vgam,rho,D,Dx 
  integer :: A,q,r,i
  real(8) :: crit

  ! allocate the workspace  
  call allocate_sp_mat(jbas,T) 
  call duplicate_sp_mat(T,F) 
  call duplicate_sp_mat(T,Vgam) 
  call duplicate_sp_mat(T,rho)
  call duplicate_sp_mat(T,D)
  call duplicate_sp_mat(T,Dx)
  
  A = H%Abody

  call write_kin_matrix(T,H,jbas) 

  !initial eigenvectors
  do q = 1, D%blocks
     D%blkM(q)%matrix = 0.d0 
     do i =1,D%map(q)
        D%blkM(q)%matrix(i,i) = 1.d0
     end do 
  end do 
  
  crit = 10.d0 
  r = 0
  !!! HARTREE FOCK MAIN LOOP 
  do while (crit > 1e-6) 
     
     call density_matrix(rho,D,DX,jbas)      
     call gamma_matrix(Vgam,H,rho,jbas)

     ! fock matrix
     do q = 1,T%blocks
        F%blkM(q)%matrix = T%blkM(q)%matrix + Vgam%blkM(Q)%matrix
     end do
       
     call diagonalize_blocks(F) 
     
     ! new eigenvectors
     ! calculate conv. criteria
     ! store eigenvalues
 
     crit = 0.d0 
     do q = 1,T%blocks
        D%blkM(q)%matrix = F%blkM(q)%matrix
        crit = crit + sqrt(sum((D%blkM(q)%eigval-F%blkM(q)%eigval)**2))
        D%blkM(q)%eigval = F%blkM(q)%eigval       
     end do
      
 end do 
  
 print*, 'Hartree Fock Energy: '
 print*, e_HF(rho,Vgam,F,jbas)

end subroutine  
!====================================================
subroutine write_kin_matrix(T,H,jbas)
  ! the kinetic energy matrix has already been calculated
  ! in the setup of the problem, here we are just writing it
  ! to a more convenient block matrix for the HF calculation
  implicit none 
  
  type(full_sp_block_mat) :: T
  type(sq_op) :: H
  type(spd) :: jbas
  integer :: q,i,j,n1,n2,c1,c2,cx,AX
  
  AX = H%belowEF
  
  do q = 1, T%blocks
     
     do i = 1, T%map(q)
        do j = i,T%map(q) 
           
           n1 = T%blkM(q)%states(i) 
           n2 = T%blkM(q)%states(j) 
           
           c1 = jbas%con(n1) 
           c2 = jbas%con(n2) 
           
           ! ph nature
           cx = c1 + c2 
           
           select case (cx) 
           ! map from pp,ph,hh strategy to block strategy 
           case(0) 
              T%blkM(q)%matrix(i,j) = H%fpp(n1-AX,n2-AX)
              T%blkM(q)%matrix(j,i) = H%fpp(n1-AX,n2-AX)
           case(1) 
              if (c2 > c1) then 
                 T%blkM(q)%matrix(i,j) = H%fph(n1-AX,n2) 
                 T%blkM(q)%matrix(j,i) = H%fph(n1-AX,n2) 
              else
                 T%blkM(q)%matrix(i,j) = H%fph(n2-AX,n1)
                 T%blkM(q)%matrix(j,i) = H%fph(n2-AX,n1)
              end if 
           case(2) 
              T%blkM(q)%matrix(i,j) = H%fhh(n1,n2)
              T%blkM(q)%matrix(j,i) = H%fhh(n1,n2)
           end select
           
           end do 
        end do 
     end do 
     
end subroutine   
!==================================================== 
subroutine gamma_matrix(gam,int,rho,jbas) 
  ! hartree fock potential matrix
  implicit none
  
  type(full_sp_block_mat) :: gam,rho
  type(sq_op) :: int
  type(spd) :: jbas
  integer :: q,r,i,jmax,lmax,n1,n2,j,JJ,n3,n4,tzrho,tzfoc
  integer :: grho,hrho,qrho,jrho,lrho,PAR,jfoc,lfoc,TZ
  real(8) :: sm

  do q = 1, gam%blocks

     gam%blkM(q)%matrix = 0.d0 
     jfoc = rho%blkM(q)%lmda(2)
     lfoc = rho%blkM(Q)%lmda(1)
     tzfoc = rho%blkM(q)%lmda(3) 
     
     ! loop over states in this block
     do i = 1, gam%map(q)
        do j = i,gam%map(q) 
        
           n1 = rho%blkM(q)%states(i) 
           n3 = rho%blkM(q)%states(j) 
           
           ! sum over blocks of the density matrix
           sm = 0.d0 
           do qrho =  1, rho%blocks
           
              jrho = rho%blkM(qrho)%lmda(2) 
              lrho = rho%blkM(qrho)%lmda(1) 
              tzrho = rho%blkM(qrho)%lmda(3) 
      
              PAR =mod(lfoc + lrho,2) 
              TZ = (tzrho + tzfoc)/2          

              ! sum over elements of the block
              do grho = 1,rho%map(qrho) 
                 do hrho = 1,rho%map(qrho) 
                       
                    n2 = rho%blkM(qrho)%states(grho) 
                    n4 = rho%blkM(qrho)%states(hrho) 
                       
                    ! sum over allowed JJ values
                    
                    do JJ = abs(jrho - jfoc),jrho+jfoc,2
                       
                       
                       sm = sm + rho%blkM(qrho)%matrix(hrho,grho) &
                            * v_elem(n1,n2,n3,n4,JJ,TZ,PAR,int,jbas) &
                            * (JJ + 1.d0)/(jrho + 1.d0) * &
                            sqrt( 1.d0 + kron_del(n1,n2)*(-1)**(JJ/2)) * &
                            sqrt( 1.d0 + kron_del(n3,n4)*(-1)**(JJ/2))
 
                       ! in J-scheme we have these weird normalization factors
                       ! that assure we aren't violating the pauli principle. 
                       ! those two square roots are them. 
                       
                    end do 
                    
                 end do 
              end do 
           end do 
   
           ! the matrix elements are multiplied by (2J+1)/(2j+1)/(2j'+1) 
           ! which amounts to averaging over J 
           
           gam%blkM(q)%matrix(i,j) = sm/(jfoc + 1.d0)
           gam%blkM(q)%matrix(j,i) = sm/(jfoc + 1.d0)

        end do 
     end do 

  end do 
  
end subroutine
!===========================================================
!===========================================================
subroutine density_matrix(rho,D,DX,jbas)
  ! density matrix is scaled by (2j+1)
  implicit none 
  
  type(full_sp_block_mat) :: rho,D,DX 
  type(spd) :: jbas
  integer :: i,j,q,N
  
  
  do q = 1, rho%blocks
     N = rho%map(q)
     if (N == 0) cycle 
                       
    ! zero out eigenvectors which aren't holes.
    do i = 1, N
       DX%blkM(q)%matrix(:,i) = D%blkM(q)%matrix(:,i) *&
            sqrt(float(rho%blkM(q)%lmda(2))+1.d0) * &
            jbas%con(rho%blkM(q)%states(i)) 
    end do 
    
    ! rho = Dx*Dx^T 
    call dgemm('N','T',N,N,N,al,Dx%blkM(q)%matrix,N,&
         Dx%blkM(q)%matrix,N,bet,rho%blkM(q)%matrix,N)
  end do 

end subroutine
!===========================================================
!===========================================================
real(8) function e_HF(rho,Vgam,F,jbas)
  ! calculate the hartree fock energy
  implicit none 
  
  real(8) :: sm
  integer :: q,i,j
  type(spd) :: jbas
  type(full_sp_block_mat) :: rho,F,Vgam 
  
 
  sm = 0.d0

  do q = 1, F%blocks
     ! sum over eigenvalues scaled by degeneracy
     do i = 1, F%map(q) 
        sm = sm + F%blkM(q)%eigval(i) &
             * jbas%con(F%blkm(q)%states(i)) * &
             (F%blkM(q)%lmda(2) + 1.d0)
           
     end do 
 
     ! this is also scaled by the degeneracy, but it's 
     ! contained within rho 
     do i = 1, F%map(q)
        do j = 1,F%map(q) 
     
           sm = sm -0.5* rho%blkM(q)%matrix(i,j) * &
             Vgam%blkM(q)%matrix(j,i) 
  
        end do 
     end do

 end do 
 
 e_HF = sm
     
end function 
!===========================================================
!===========================================================              
end module         
                
              
                 
                 
           
            
