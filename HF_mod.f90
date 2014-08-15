module HF_mod
  use basic_IMSRG 
  ! J-scheme Hartree-Fockery. 
  implicit none 
  
contains
!====================================================
subroutine calc_HF( H ,jbas )
  ! returns H in the normal orderd Hartree Fock basis
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: H 
  type(full_sp_block_mat) :: T,F,Vgam,rho,D,Dx 
  integer :: q,r,i
  real(8) :: crit,sm

  ! allocate the workspace  

  call allocate_sp_mat(jbas,T) 
  call duplicate_sp_mat(T,F) 
  call duplicate_sp_mat(T,Vgam) 
  call duplicate_sp_mat(T,rho)
  call duplicate_sp_mat(T,D)
  call duplicate_sp_mat(T,Dx)
  
  call write_kin_matrix(T,H,jbas) 

  !initial eigenvectors
  do q = 1, D%blocks
     D%blkM(q)%eigval = 0.d0
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
 
 do q = 1,T%blocks
    F%blkM(q)%matrix = T%blkM(q)%matrix + Vgam%blkM(Q)%matrix
    T%blkM(q)%eigval = F%blkM(q)%eigval    
 end do

 call transform_1b_to_HF(D,Dx,F,T,H,jbas) 
  
 ! this needs to come after the transformation
 ! e_HF is calculated in the hartree fock basis
 H%E0 = e_HF(T,jbas)

 call transform_2b_to_HF(D,H,jbas) 

end subroutine calc_HF
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
           
           T%blkM(q)%matrix(i,j) = f_elem(n1,n2,H,jbas) 
           T%blkM(q)%matrix(j,i) = T%blkM(q)%matrix(i,j) 
           
           
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
     
     ! loop over states in this block
     do i = 1, gam%map(q)
        do j = i,gam%map(q) 
        
           n1 = rho%blkM(q)%states(i) 
           n3 = rho%blkM(q)%states(j) 
           
           ! sum over blocks of the density matrix
           sm = 0.d0 
           do qrho =  1, rho%blocks
           
              jrho = rho%blkM(qrho)%lmda(2) 
                  
              ! sum over elements of the block
              do grho = 1,rho%map(qrho) 
                 do hrho = 1,rho%map(qrho) 
                    
                    n2 = rho%blkM(qrho)%states(grho) 
                    n4 = rho%blkM(qrho)%states(hrho) 
                       
                    ! sum over allowed JJ values
                    
                    do JJ = abs(jrho - jfoc),jrho+jfoc,2
                       
                       sm = sm + rho%blkM(qrho)%matrix(hrho,grho) &
                            * v_elem(n1,n2,n3,n4,JJ,int,jbas) &
                            * (JJ + 1.d0)/(jrho + 1.d0)
                             
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
            jbas%con(rho%blkM(q)%states(i)) 
    end do 
    
    ! rho = Dx*Dx^T 
    call dgemm('N','T',N,N,N,al,Dx%blkM(q)%matrix,N,&
         Dx%blkM(q)%matrix,N,bet,rho%blkM(q)%matrix,N)
    rho%blkM(q)%matrix = rho%blkM(q)%matrix * (rho%blkM(q)%lmda(2)+1)
    
    
  end do 
 
end subroutine
!===========================================================
!===========================================================
real(8) function e_HF(T,jbas)
  ! calculate the hartree fock energy
  implicit none 
  
  real(8) :: sm
  integer :: q,i
  type(spd) :: jbas
  type(full_sp_block_mat) :: T
  
 
  sm = 0.d0
  
  do q = 1, T%blocks
     ! sum over eigenvalues scaled by degeneracy
     do i = 1, T%map(q) 
        sm = sm + (T%blkM(q)%eigval(i) + &
             T%blkM(q)%matrix(i,i))&
             * jbas%con(T%blkm(q)%states(i)) * &
             (T%blkM(q)%lmda(2) + 1.d0) 
                  
     end do   
 end do 
 
 e_HF = sm*0.5d0
 
 
     
end function 
!===========================================================
!===========================================================
subroutine transform_1b_to_HF(D,Dx,F,T,H,jbas) 
  ! typical transformation, remap to fancy array
  implicit none 
  
  type(sq_op) :: H
  type(spd) :: jbas
  type(full_sp_block_mat) :: D,F,Dx,T
  integer :: q,dm,i,j,a,b,c1,c2,cx
  
  do q = 1, F%blocks
     
   
     dm = F%map(q)
     if (dm == 0)  cycle
     
     ! transform the Fock matrix
     call dgemm('N','N',dm,dm,dm,al,F%blkM(q)%matrix&
          ,dm,D%blkM(q)%matrix,dm,bet,Dx%blkM(q)%matrix,dm) 
     call dgemm('T','N',dm,dm,dm,al,D%blkM(q)%matrix&
          ,dm,Dx%blkM(q)%matrix,dm,bet,F%blkM(q)%matrix,dm) 
     
     ! transform the KE matrix
     call dgemm('N','N',dm,dm,dm,al,T%blkM(q)%matrix&
          ,dm,D%blkM(q)%matrix,dm,bet,Dx%blkM(q)%matrix,dm) 
     call dgemm('T','N',dm,dm,dm,al,D%blkM(q)%matrix&
          ,dm,Dx%blkM(q)%matrix,dm,bet,T%blkM(q)%matrix,dm) 
     
     do a=1,dm
        do b=a,dm
           
           i = F%blkM(q)%states(a) 
           j = F%blkM(q)%states(b) 
           
           c1 = jbas%con(i) 
           c2 = jbas%con(j) 
           cx = c1 + c2 
           
           ! fancy array remap ( normal ordered now ) 
     select case (cx) 
        case(0) 
           H%fpp(i-jbas%holesb4(i),j-jbas%holesb4(j)) = &
                F%blkM(q)%matrix(a,b) 
           H%fpp(j-jbas%holesb4(j),i-jbas%holesb4(i)) = &
                F%blkM(q)%matrix(a,b) 
        case(1) 
           if (c2 > c1) then 
              H%fph(i-jbas%holesb4(i),j-jbas%partsb4(j)) = &
                F%blkM(q)%matrix(a,b)   
           else
              H%fph(j-jbas%holesb4(j),i-jbas%partsb4(i)) = &
                F%blkM(q)%matrix(b,a)
           end if 
        case(2) 
           H%fhh(i-jbas%partsb4(i),j-jbas%partsb4(j)) = &
                F%blkM(q)%matrix(a,b)
           H%fhh(j-jbas%partsb4(j),i-jbas%partsb4(i)) = &
                F%blkM(q)%matrix(a,b)
     end select
         
         end do 
     end do 
     
  end do 
    
end subroutine 
!===========================================================
!===========================================================  
subroutine transform_2b_to_HF(D,H,jbas)
  ! map TBME into HF basis.DFWT
  ! huge mess. It's way too much of a hastle if you don't
  ! put everything into large square matrices and work from there
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: H
  type(full_sp_block_mat) :: D 
  integer :: q,i,j,a,b,nh,np,nb,nt,II,JJ,int1,int2
  integer :: j_i,j_j,j_a,j_b,li,lj,la,lb,ti,tj,ta,tb
  integer :: q1,q2,q3,q4
  integer,allocatable,dimension(:,:) :: qnbig
  real(8),allocatable,dimension(:,:) :: Vfull,Cfull,temp,Dsmall
  real(8),allocatable,dimension(:,:) :: V1,V2,V3,V4,crevfull
  
  ! construct full D matrix (with zeros) 
  allocate(Dsmall(H%Nsp,H%Nsp)) 
  Dsmall = 0.d0 
  do q = 1, D%blocks
     if (D%map(q) == 0) cycle
     
     do i = 1, D%map(q) 
        do j = 1,D%map(q) 
           
           Dsmall(D%blkM(q)%states(i),D%blkM(q)%states(j)) = &
            D%blkM(q)%matrix(i,j) 
        end do 
     end do 
  end do 

  ! cycle over major blocks of the tp hamiltonian
  do q = 1, H%nblocks
     
     np = H%mat(q)%npp
     nh = H%mat(q)%nhh
     nb = H%mat(q)%nph
     nt = np+nh+nb
     if (nt == 0) cycle
     
     allocate(Vfull(nt,nt)) 
     allocate(V1(nt,nt),V2(nt,nt),V3(nt,nt),V4(nt,nt)) 
     allocate(Cfull(nt,nt)) 
     allocate(Crevfull(nt,nt))
     allocate(temp(nt,nt)) 
     allocate(qnbig(nt,2)) 
     
     ! mapping things out to a square matrix 
     
     Vfull(1:nh,1:nh) = H%mat(q)%gam(5)%X 
     
     Vfull(nh+1:nb+nh,1:nh) = H%mat(q)%gam(6)%X
     Vfull(1:nh,nh+1:nb+nh) = Transpose(H%mat(q)%gam(6)%X) 
     
     Vfull(nh+nb+1:nt,1:nh) = H%mat(q)%gam(3)%X 
     Vfull(1:nh,nh+nb+1:nt) = Transpose(H%mat(q)%gam(3)%X) 
     
     Vfull(nh+1:nb+nh,nh+1:nb+nh) = H%mat(q)%gam(4)%X
     
     Vfull(nh+nb+1:nt,nh+1:nh+nb) = H%mat(q)%gam(2)%X
     Vfull(nh+1:nh+nb,nh+nb+1:nt) = Transpose(H%mat(q)%gam(2)%X)
           
     Vfull(nh+nb+1:nt,nh+nb+1:nt) = H%mat(q)%gam(1)%X 
     
     qnbig(1:nh,:) = H%mat(q)%qn(3)%Y
     qnbig(nh+1:nh+nb,:) = H%mat(q)%qn(2)%Y
     qnbig(nh+nb+1:nt,:) = H%mat(q)%qn(1)%Y
     
     ! filling the C matrices
   
     do II = 1,nt
           
        i = qnbig(II,1)
        j = qnbig(II,2)
          
        do JJ = 1, nt
           
           a = qnbig(JJ,1)
           b = qnbig(JJ,2)
                              
           Cfull(JJ,II) = Dsmall(a,i)*Dsmall(b,j) 
           
           Crevfull(JJ,II) = Dsmall(b,i)*Dsmall(a,j) *&
                (1 - kron_del(a,b)) * &
           (-1)**( (jbas%jj(a) + jbas%jj(b) ) /2 ) 
         
           
        end do  
     end do 
        
     ! do the transformation
     call dgemm('N','N',nt,nt,nt,al,Vfull,nt,Cfull,nt,bet,temp,nt) 
     call dgemm('T','N',nt,nt,nt,al,Cfull,nt,temp,nt,bet,V1,nt) 
     
     call dgemm('N','N',nt,nt,nt,al,Vfull,nt,Crevfull,nt,bet,temp,nt) 
     call dgemm('T','N',nt,nt,nt,al,Cfull,nt,temp,nt,bet,V2,nt) 
     
     call dgemm('N','N',nt,nt,nt,al,Vfull,nt,Cfull,nt,bet,temp,nt) 
     call dgemm('T','N',nt,nt,nt,al,Crevfull,nt,temp,nt,bet,V3,nt) 
     
     call dgemm('N','N',nt,nt,nt,al,Vfull,nt,Crevfull,nt,bet,temp,nt) 
     call dgemm('T','N',nt,nt,nt,al,Crevfull,nt,temp,nt,bet,V4,nt)
     
     Vfull = V1 - (-1)**(H%mat(q)%lam(1)/2)*(V2+V3) + V4  
   
     ! put back into the fancy array
     H%mat(q)%gam(5)%X = Vfull(1:nh,1:nh) 
     H%mat(q)%gam(6)%X = Vfull(nh+1:nb+nh,1:nh) 
     H%mat(q)%gam(3)%X = Vfull(nh+nb+1:nt,1:nh) 
     H%mat(q)%gam(4)%X = Vfull(nh+1:nb+nh,nh+1:nb+nh)
     H%mat(q)%gam(2)%X = Vfull(nh+nb+1:nt,nh+1:nh+nb)      
     H%mat(q)%gam(1)%X = Vfull(nh+nb+1:nt,nh+nb+1:nt)
 
     ! deallocate everything
     deallocate(Vfull)
     deallocate(V1,V2,V3,V4)
     deallocate(Cfull)
     deallocate(Crevfull)
     deallocate(qnbig) 
     deallocate(temp)
     
  end do 
end subroutine     
!==============================================================
!==============================================================
end module         
                
              
                 
                 
           
            
