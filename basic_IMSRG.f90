module basic_IMSRG 
  ! contains basic type declarations, constants, and adminstrative routines
  real(8),public,parameter :: al = 1.d0 , bet = 0.d0  
  integer,public,dimension(9) :: adjust_index = (/0,0,0,0,0,0,0,0,-4/) 

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  TYPE :: spd    ! single particle discriptor
     INTEGER :: total_orbits,Jtotal_max,lmax
     INTEGER, ALLOCATABLE,DIMENSION(:) :: nn, ll, jj, itzp, nshell, mvalue
     INTEGER, ALLOCATABLE,DIMENSION(:) :: con,holesb4,partsb4 
     ! for clarity:  nn, ll, nshell are all the true value
     ! jj is j+1/2 (so it's an integer) 
     ! likewise itzp is 2*tz 
     CHARACTER (LEN=10), ALLOCATABLE,DIMENSION(:) :: orbit_status, model_space,basis
     REAL(8), ALLOCATABLE,DIMENSION(:) :: e, evalence, e_original
  END TYPE spd
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TYPE :: real_mat !element of an array of matrices
     real(8),allocatable,dimension(:,:) :: X 
  END TYPE real_mat
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TYPE :: int_mat 
     integer,allocatable,dimension(:,:) :: Y
  END TYPE int_mat
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  TYPE :: int_vec
     integer,allocatable,dimension(:) :: Z
  END TYPE int_vec
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  TYPE ::  six_index_store
     type(real_mat),allocatable,dimension(:,:) :: tp_mat
     integer :: nf,nb,nhalf
  END TYPE six_index_store
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TYPE :: sq_block
     integer :: lam(3) ! specifices J,Par,Tz of the block 
     type(real_mat),dimension(6) :: gam !Vpppp,Vhhhh,Vphph,Vpphh,Vphhh,Vppph
     integer :: npp, nph , nhh ! dimensions
     type(int_mat),dimension(3) :: qn !tp to sp map 
     type(int_vec),dimension(3) :: pnt ! for mapping from sp to tp (confusing) 
  END TYPE sq_block
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TYPE :: sq_op !second quantized operator 
     type(sq_block),allocatable,dimension(:) :: mat
     real(8),allocatable,dimension(:,:) :: fph,fpp,fhh
     integer :: nblock,Aprot,Aneut,Nsp,herm,belowEF
     real(8) :: E0 
  END TYPE sq_op
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
type sp_block_mat
   real(8),allocatable,dimension(:,:) :: matrix!,eigvec
   real(8),allocatable,dimension(:) :: Eigval,extra
   integer,allocatable,dimension(:,:) :: labels
   integer,allocatable,dimension(:) :: states
   integer,dimension(3) :: lmda
end type sp_block_mat
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
type full_sp_block_mat
   type(sp_block_mat),allocatable,dimension(:) :: blkM
   integer,allocatable,dimension(:) :: map
   integer :: blocks
end type full_sp_block_mat
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   type(six_index_store),public :: store6j

contains 
!====================================================
!====================================================
subroutine read_sp_basis(jbas,sp_input_file,hp,hn)
  ! fills jscheme_basis array with single particle data from file
  ! file format must be: 5 integers and one real 
  ! for each state we have:   | label |  n  | l | 2 * j | 2*tz | E_sp |  
  implicit none 
  
  type(spd) :: jbas
  character(50) :: sp_input_file
  integer :: ist,i,label,ni,li,ji,tzi,ix,hp,hn
  real(8) :: e
  
  open(unit=39,file='../sp_inputs/'//trim(adjustl(sp_input_file)))
  
  ix = 0 
  
  ! count the number of states in the file. 
  do 
  
     read(39,*,iostat=ist) label,ni,li,ji,tzi,e
     
     if (ist>0) stop 'input error in .sps file. See &
          subroutine read_sp_basis in basic_IMSRG.f90' 
     if (ist<0) exit
  
     ix = ix + 1 
  
  end do
  
  ! build the jscheme descriptor
  jbas%total_orbits=ix
  allocate(jbas%nn(ix)) ! n 
  allocate(jbas%ll(ix)) ! l
  allocate(jbas%jj(ix)) ! j*2
  allocate(jbas%itzp(ix)) !isospin*2
  allocate(jbas%nshell(ix)) !shell number
  allocate(jbas%e(ix)) ! sp energies
  allocate(jbas%con(ix)) ! hole or particle (1 or 0) 
  allocate(jbas%holesb4(ix)) !number of holes beneath this index
  allocate(jbas%partsb4(ix)) !number of particles beneath this index
  
  ! go back to the start and read them in. 
  rewind(39) 
  
  do i = 1,ix
     
     read(39,*) label,ni,li,ji,tzi,e
     
     jbas%nn(i) = ni 
     jbas%ll(i) = li
     jbas%jj(i) = ji
     jbas%itzp(i) = tzi 
     jbas%nshell(i) = 2*ni + li 
     jbas%e(i) =  e 

  end do 
  
  call find_holes(jbas,hp,hn) 
  
  jbas%Jtotal_max = maxval(jbas%jj) 
  jbas%lmax = maxval(jbas%ll) 
  ! this is the maximum value of J, this code always uses 2*J 
  call store_6j(jbas) 
  close(39) 
  
end subroutine 
!==============================================
!==============================================
subroutine find_holes(jbas,pholes,nholes) 
  implicit none 
  
  type(spd) :: jbas
  integer :: pholes,nholes,i,minpos(1),rn,rp
  real(8),dimension(jbas%total_orbits) :: temp
  
  temp = jbas%e
  jbas%con = 0 ! zero if particle, one if hole
  
  rn = 0
  rp = 0 
  do while ((rn < nholes) .or. (rp < pholes))
     minpos = minloc(temp) !find lowest energy
     jbas%con(minpos(1)) = 1 ! thats a hole
     temp(minpos(1)) = 9.e9  ! replace it with big number
     if ( jbas%itzp(minpos(1)) == -1 ) then
        if (rn < nholes ) then 
           rn = rn + (jbas%jj(minpos(1))+1)  
        else 
           jbas%con(minpos(1)) = 0 
        end if 
     else
        if (rp < pholes ) then 
           rp = rp + (jbas%jj(minpos(1))+1)  
        else 
           jbas%con(minpos(1)) = 0 
        end if 
     end if 
        
  end do 
  
  ! these arrays help us later in f_elem (in this module) 
  do i = 1, jbas%total_orbits
     jbas%holesb4(i) = sum(jbas%con(1:i-1)) 
     jbas%partsb4(i) = i - jbas%holesb4(i) - 1 
  end do 
end subroutine  
!==============================================
!==============================================
subroutine allocate_blocks(jbas,op) 
  ! allocate the sq_op arrays
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: op
  integer :: N,AX,q,i,j,j1,j2,tz1,tz2,l1,l2
  integer :: Jtot,Tz,Par,nph,npp,nhh,CX
  
  AX = sum(jbas%con) !number of shells under eF 
  op%belowEF = AX  
  op%Nsp = jbas%total_orbits
  N = op%Nsp  !number of sp shells
  
  op%nblock =  (jbas%Jtotal_max + 1) * 6    ! 6 possible values of par * Tz    
  
  allocate(op%mat(op%nblock)) 
  allocate(op%fpp(N-AX,N-AX))
  allocate(op%fph(N-AX,AX))
  allocate(op%fhh(AX,AX)) 
  
  q = 1  ! block index
  
  do Jtot = 0,2*jbas%Jtotal_max,2 !looping over blocks
     do Tz = -1,1
        do Par = 0,1
     
     op%mat(q)%lam(1)=Jtot
     op%mat(q)%lam(2)=Par
     op%mat(q)%lam(3)=Tz
           
     nph = 0 
     npp = 0 
     nhh = 0 
     
     do i = 1,N   !looping over sp states
        do j = i,N
 
           ! check if the two sp states can exist in this block 
           tz1 = jbas%itzp(i)
           tz2 = jbas%itzp(j) 
           if ( Tz .ne. (tz1 + tz2)/2 ) cycle
           
           l1 = jbas%ll(i) 
           l2 = jbas%ll(j)
           if ( Par .ne. (1 - (-1)**(l1 + l2))/2 ) cycle
           
           j1 = jbas%jj(i) 
           j2 = jbas%jj(j) 
           if (.not. triangle(j1,j2,Jtot) ) cycle
           
           cX = jbas%con(i) + jbas%con(j)
 
           select case (CX) ! pp, ph or hh ? 
              case (0) 
                 npp = npp + 1
              case (1) 
                 nph = nph + 1
              case (2) 
                 nhh = nhh + 1
           end select
                     
        end do 
     end do      

     op%mat(q)%npp = npp
     op%mat(q)%nph = nph
     op%mat(q)%nhh = nhh
     
     allocate(op%mat(q)%qn(1)%Y(npp,2)) !qnpp
     allocate(op%mat(q)%qn(2)%Y(nph,2)) !qnph
     allocate(op%mat(q)%qn(3)%Y(nhh,2)) !qnhh
     allocate(op%mat(q)%pnt(1)%Z(npp)) !qnpp
     allocate(op%mat(q)%pnt(2)%Z(nph)) !qnph
     allocate(op%mat(q)%pnt(3)%Z(nhh)) !qnhh
    
     allocate(op%mat(q)%gam(1)%X(npp,npp)) !Vpppp
     allocate(op%mat(q)%gam(5)%X(nhh,nhh)) !Vhhhh
     allocate(op%mat(q)%gam(3)%X(npp,nhh)) !Vpphh
     allocate(op%mat(q)%gam(4)%X(nph,nph)) !Vphph
     allocate(op%mat(q)%gam(2)%X(npp,nph)) !Vppph
     allocate(op%mat(q)%gam(6)%X(nph,nhh)) !Vphhh
     do i = 1,6
        op%mat(q)%gam(i)%X = 0.0
     end do 
     
     nph = 0 ; npp = 0 ; nhh = 0
     do i = 1,N   !looping over sp states
        do j = i,N
           
           ! check if the two sp states can exist in this block 
           tz1 = jbas%itzp(i)
           tz2 = jbas%itzp(j) 
           if ( Tz .ne. (tz1 + tz2)/2 ) cycle
           
           l1 = jbas%ll(i) 
           l2 = jbas%ll(j)
           if ( Par .ne. (1 - (-1)**(l1 + l2))/2 ) cycle
           
           j1 = jbas%jj(i) 
           j2 = jbas%jj(j) 
           if (.not. triangle(j1,j2,Jtot) ) cycle
           
           cX = jbas%con(i) + jbas%con(j)
           
           select case (CX)
              case (0) 
                 npp = npp + 1
                 op%mat(q)%qn(1)%Y(npp,1) = i
                 op%mat(q)%qn(1)%Y(npp,2) = j
                 op%mat(q)%pnt(1)%Z(npp) = N*(i-1) + j !search for an integer 
              case (1)                                 ! instead of a pair of them  
                 nph = nph + 1
                 op%mat(q)%qn(2)%Y(nph,1) = i
                 op%mat(q)%qn(2)%Y(nph,2) = j 
                 op%mat(q)%pnt(2)%Z(nph) = N*(i-1) + j 
              case (2) 
                 nhh = nhh + 1
                 op%mat(q)%qn(3)%Y(nhh,1) = i
                 op%mat(q)%qn(3)%Y(nhh,2) = j 
                 op%mat(q)%pnt(3)%Z(nhh) = N*(i-1) + j 
           end select
                     
        end do 
     end do    
     q = q + 1
        end do
     end do
  end do
end subroutine      
!==================================================================  
!==================================================================
subroutine read_interaction(H,intfile,jbas) 
  ! read interaction from ASCII file produced by Scott_to_Morten.f90 
  implicit none
  
  type(sq_op) :: H
  type(spd) :: jbas
  character(50) :: intfile
  integer :: ist,J,Tz,Par,a,b,c,d,q,qx,N 
  real(8) :: V,g1,g2,g3
  integer :: C1,C2,int1,int2,i1,i2,pre
  
  open(unit=39,file = '../TBME_input/'//trim(adjustl(intfile))) 
  
  read(39,*);read(39,*);read(39,*);read(39,*)
  read(39,*);read(39,*);read(39,*);read(39,*) !skip all of the garbage 
  
  N = jbas%total_orbits 
  do 
     read(39,*,iostat=ist) Tz,Par,J,a,b,c,d,V,g1,g2,g3
     !read(39,*) Tz,Par,J,a,b,c,d,V,g1,g2,g3
     if (ist > 0) STOP 'interaction file error' 
     if (ist < 0) exit
     
     q = block_index(J,Tz,Par) 
     
     C1 = jbas%con(a)+jbas%con(b) + 1 !ph nature
     C2 = jbas%con(c)+jbas%con(d) + 1
    
     qx = C1*C2
     qx = qx + adjust_index(qx)   !Vpppp nature  

     ! get the indeces in the correct order
     pre = 1
     if ( a > b )  then 
        int1 = N*(b-1) + a 
        pre = (-1)**( 1 + (jbas%jj(a) + jbas%jj(b) -J)/2 ) 
     else
        int1 = N*(a-1) + b
     end if
  
     if (c > d)  then     
        int2 = N*(d-1) + c
        pre = pre * (-1)**( 1 + (jbas%jj(c) + jbas%jj(d) -J)/2 ) 
     else 
        int2 = N*(c-1) + d
     end if
  
     i1 = 1
     i2 = 1
      
     ! find the indeces
    
     do 
        if (int1 == H%mat(q)%pnt(C1)%Z(i1) ) exit
        i1 = i1 + 1
     end do 
        
     do 
        if (int2 == H%mat(q)%pnt(C2)%Z(i2) ) exit
        i2 = i2 + 1
     end do 
 
     if ((qx == 1) .or. (qx == 5) .or. (qx == 4)) then 
        H%mat(q)%gam(qx)%X(i2,i1)  = V *pre
        H%mat(q)%gam(qx)%X(i1,i2)  = V *pre
     else if (C1>C2) then
        H%mat(q)%gam(qx)%X(i2,i1)  = V *pre
     else
        H%mat(q)%gam(qx)%X(i1,i2) = V * pre
     end if 
     ! I shouldn't have to worry about hermiticity here, it is assumed to be hermitian
     
  end do    
     
end subroutine
!==================================================================  
!==================================================================
logical function triangle(j1,j2,J) 
  ! all entries are multiplied by 2 (to eliminate half integers) 
  implicit none 
  
  integer :: j1,j2,J
  

  if ( ( J .le. j1+j2 ) .and. (J .ge. abs(j1-j2) ) )then
     triangle = .true. 
  else 
     triangle = .false.
  end if 
end function 
!=================================================================     
!=================================================================
integer function block_index(J,T,P) 
  ! input 2*J,Tz,and Parity to get block index
  integer :: J,T,P
  
  block_index = 3*J + 2*(T+1) + P + 1
  
end function 
!=================================================================     
!=================================================================
real(8) function f_elem(a,b,op,jbas) 
  implicit none 
  
  integer :: a,b,x1,x2,c1,c2
  type(spd) :: jbas
  type(sq_op) :: op 
  
  ! are they holes or particles
  c1 = jbas%con(a)
  c2 = jbas%con(b) 
  
  select case(c1+c2) 
     case(0) 
        ! pp 
        f_elem = op%fpp(a-jbas%holesb4(a),b-jbas%holesb4(b)) 
  
     case(1) 
        ! ph 
        if (c1 > c2) then 
           f_elem = op%fph(b-jbas%holesb4(b),a-jbas%partsb4(a)) 
        else 
           f_elem = op%fph(a-jbas%holesb4(a),b-jbas%partsb4(b)) 
        end if
     case(2) 
        ! hh 
        f_elem = op%fhh(a-jbas%partsb4(a),b-jbas%partsb4(b)) 
  end select

end function
!==============================================================
!==============================================================
real(8) function v_elem(a,b,c,d,J,T,P,op,jbas) 
  ! grabs the matrix element you are looking for
  implicit none
  
  integer :: a,b,c,d,J,T,P,q,qx,c1,c2
  integer :: int1,int2,pre,i1,i2
  type(sq_op) :: op 
  type(spd) :: jbas
 
  q = block_index(J,T,P) 
     
  !! this is a ridiculous but efficient? way of 
  !! figuring out which Vpppp,Vhhhh,Vphph.... 
  !! array we are referencing. QX is a number between 
  !! 1 and 6 which refers to a unique array for the pphh characteristic
  
  C1 = jbas%con(a)+jbas%con(b) + 1 !ph nature
  C2 = jbas%con(c)+jbas%con(d) + 1
    
  qx = C1*C2
  qx = qx + adjust_index(qx)   !Vpppp nature  
  
  ! see subroutine "allocate_blocks" for mapping from qx to each 
  ! of the 6 storage arrays
  
  i1 = 1
  i2 = 1
  
  pre = 1 
  
  ! get the indeces in the correct order
  if ( a > b )  then 
     int1 = op%Nsp*(b-1) + a 
     pre = (-1)**( 1 + (jbas%jj(a) + jbas%jj(b) -J)/2 ) 
  else
     int1 = op%Nsp*(a-1) + b
  end if 
  
  if (c > d)  then     
     int2 = op%Nsp*(d-1) + c
     pre = pre * (-1)**( 1 + (jbas%jj(c) + jbas%jj(d) -J)/2 ) 
  else 
     int2 = op%Nsp*(c-1) + d
  end if 
     
! int1 and int2 are unique integers which refer to the 
! two-p pair in the bra and ket respectively 

! find the indeces by searching for unique integers
  do 
     if (int1 == op%mat(q)%pnt(C1)%Z(i1) ) exit
     i1 = i1 + 1
  end do 
        
  do 
     if (int2 == op%mat(q)%pnt(C2)%Z(i2) ) exit
     i2 = i2 + 1
  end do 

  ! grab the matrix element
   If (C1>C2) then 
      v_elem = op%mat(q)%gam(qx)%X(i2,i1) * op%herm * pre 
   else
      v_elem = op%mat(q)%gam(qx)%X(i1,i2) * pre
   end if 
  
end function 
!=====================================================
!=====================================================
subroutine calculate_h0_harm_osc(hw,jbas,H,Htype) 
  ! fills out the one body piece of the hamiltonian
  implicit none 
  
  integer :: i,j,mass,Htype,c1,c2,cx
  integer :: ni,li,ji,nj,lj,jj,tzi,tzj,AX
  real(8) :: hw,kij,T
  type(sq_op) :: H 
  type(spd) :: jbas

  
  !Htype =  |[1]: T - Tcm + V |[2]: T + Uho + V |[3]: T+V | 

  AX = H%belowEF
  mass = H%Aprot + H%Aneut

  
  H%fpp = 0.d0
  H%fhh = 0.d0
  H%fph = 0.d0
  do i = 1, H%Nsp
     
     ! extract info
     ni = jbas%nn(i) 
     li = jbas%ll(i) 
     ji = jbas%jj(i)
     tzi = jbas%itzp(i) 
     c1 = jbas%con(i) 
     
     do j = i,H%Nsp
     
        ! extract info, cycle if irrelevant
        lj = jbas%ll(j) 
        if (li .ne. lj) cycle 
        jj = jbas%jj(j) 
        if (ji .ne. jj) cycle
        tzj = jbas%itzp(j) 
        if (tzi .ne. tzj) cycle
        c2 = jbas%con(j)
        nj = jbas%nn(j) 

        if ( abs(ni - nj) > 1) cycle
        ! standard HO relations
        if (ni == nj) then 
           kij = 0.5*hw*(2*nj + lj +1.5)         
        else if (ni == nj - 1) then 
           kij = 0.5*hw * sqrt(nj*(nj+lj+0.5))            
        else if (ni == nj + 1) then 
           kij = 0.5*hw * sqrt((nj+1)*(nj+lj+1.5))
        end if 
        
        select case (Htype) 
           case(1) 
              T =  kij*(1.d0-1.d0/mass) 
           case(2) 
              T = 2*kij*kron_del(ni,nj)
           case(3)
              T = kij 
        end select 
        
        cx = c1+c2 
        
        
        select case (cx) 
           case(0) 
              H%fpp(i-jbas%holesb4(i),j-jbas%holesb4(j)) = T
              H%fpp(j-jbas%holesb4(j),i-jbas%holesb4(i)) = T
           case(1) 
              if (c2 > c1) then 
                 H%fph(i-jbas%holesb4(i),j-jbas%partsb4(j)) = T
              else
                 H%fph(j-jbas%holesb4(j),i-jbas%partsb4(i)) = T
              end if 
           case(2) 
              H%fhh(i-jbas%partsb4(i),j-jbas%partsb4(j)) = T
              H%fhh(j-jbas%partsb4(j),i-jbas%partsb4(i)) = T
        end select
     
     end do
  end do 

end subroutine         
!====================================================
!====================================================  
integer function kron_del(i,j) 
  implicit none 
  
  integer :: i,j 
  
  if (i == j) then
     kron_del = 1
  else
     kron_del = 0
  end if 
end function
!====================================================
!==================================================== 
subroutine diagonalize_blocks(R)
  implicit none 
  
  type(full_sp_block_mat) :: R
  integer :: a,q,info
  
  do q=1,R%blocks
     
     a=R%map(q)
     if (a == 0)  cycle 
     
     call dsyev('V','U',a,R%blkM(q)%matrix,a, &
          R%blkM(q)%eigval,R%blkM(q)%extra,10*a,info)
   
  end do 
end subroutine
!======================================================
!======================================================
subroutine store_6j(jbas) 
  ! run this once. use sixj in code to call elements 
  ! kind of complicated, but if fills a public array (store6j) 
  ! which theoretically doesn't need to be touched when 
  ! making changes to the code ( this is all background stuff ) 
  ! so DFWT subroutine because it will be a huge mess
  implicit none 
  
  type(spd) :: jbas
  integer :: j1,j2,j3,j4,j5,j6,num_half,num_whole,r1,r2
  integer :: nbos,nferm,X12,X45,j3min,j3max,j6min,j6max
  real(8) :: d6ji
  
  call dfact0() ! prime the anglib 6j calculator 
  
  num_half = (jbas%Jtotal_max + 1)/2
  store6j%nhalf = num_half
  nbos = num_half*(num_half+1)/2
  nferm = num_half*(num_half-1)/2
  
  store6j%nb = nbos
  store6j%nf = nferm 
  
  allocate(store6j%tp_mat(nbos,nbos+nferm)) 
  ! The first index refers to j1,j2 
  ! which are forcibly ordered
  
  ! the second index refers to j4,j5 which cannot be ordered
  ! so we need an extra set of states for the reverse ordering (nferm) 
  
  do j1 = 1, num_half
     do j2 = j1, num_half 

        X12 = bosonic_tp_index(j1,j2,num_half)
        
        do j4 = 1, num_half
           do j5 = 1, j4-1
           
              ! if j5<j4 use fermionic index
              X45 = fermionic_tp_index(j5,j4,num_half) 
              
              ! find the allowed values of j3,j6
            j3min = max(  abs( 2*j1-1 -  (2*j2-1) ) , abs( 2*j4-1 -  (2*j5-1) )   )               
            j3max = min(  2*j1-1 +  2*j2-1 , 2*j4-1 +  2*j5-1 ) 

            j6min = max(  abs( 2*j1-1 -  (2*j5-1) ) , abs( 2*j4-1 -  (2*j2-1) )   )               
            j6max = min(  2*j1-1 +  2*j5-1 , 2*j4-1 +  2*j2-1 ) 
            
            ! allocate second part of storage array to hold them 
              allocate(store6j%tp_mat(X12,X45)%X((j3max-j3min)/2+1,(j6max-j6min)/2+1 ) ) 
              
              ! fill the array with 6j symbols
              r1 = 1 
             
              do j3 = j3min,j3max,2
                 r2 = 1
                 do j6 = j6min,j6max,2
                    store6j%tp_mat(X12,X45)%X(r1,r2) =  d6ji( 2*j1-1 , 2*j2 - 1 , j3, &
                         2*j4 - 1 , 2*j5-1 , j6 )  
                    ! d6ji is the anglib sixj calculator (takes 2*j as arguments) 
                    r2 = r2 + 1
                 end do 
                 r1 = r1 + 1
              end do 
                 
           end do 
           
           do j5 = j4,num_half
              ! same deal, but with bosonic index, offset by fermionic
              X45 = bosonic_tp_index(j4,j5,num_half) + nferm  
            !  print*, j4,j5,X45,nferm,bosonic_tp_index(j4,j5,nbos),nbos
            j3min = max(  abs( 2*j1-1 -  (2*j2-1) ) , abs( 2*j4-1 -  (2*j5-1) )   )               
            j3max = min(  2*j1-1 +  2*j2-1 , 2*j4-1 +  2*j5-1 ) 

            j6min = max(  abs( 2*j1-1 -  (2*j5-1) ) , abs( 2*j4-1 -  (2*j2-1) )   )               
            j6max = min(  2*j1-1 +  2*j5-1 , 2*j4-1 +  2*j2-1 ) 
            
              allocate(store6j%tp_mat(X12,X45)%X((j3max-j3min)/2+1,(j6max-j6min)/2+1 ) )   
 
              !print*, j1,j2,j4,j5,(j3max-j3min)/2+1,(j6max-j6min)/2+1 ,j3min,j3max,j6min,j6max
              r1 = 1
              do j3 = j3min,j3max,2
                 r2 = 1
                 do j6 = j6min,j6max,2
                    store6j%tp_mat(X12,X45)%X(r1,r2) =  d6ji( 2*j1-1 , 2*j2 - 1 , j3, &
                         2*j4 - 1 , 2*j5-1 , j6 )  
                    
                    r2 = r2 + 1
                 end do 
                 
                 r1 = r1 + 1
              end do 
                     
           end do 
        end do 
        
     end do 
  end do 
   
end subroutine    
!=========================================================
!=========================================================
real(8) function sixj(j1,j2,j3,j4,j5,j6)
  ! twice the angular momentum
  ! so you don't have to carry that obnoxious array around
  implicit none 
  
  integer :: j1,j2,j3,j4,j5,j6
  integer :: l1,l2,l3,l4,l5,l6
  integer :: x1,x2,j3min,j6min
  
  
  ! check triangle inequalities
  if ( (triangle(j1,j2,j3)) .and. (triangle(j4,j5,j3)) &
  .and. (triangle(j1,j5,j6)) .and. (triangle(j4,j2,j6)) ) then 
     
     l1 = j1; l2 = j2; l3 = j3; l4 = j4; l5 = j5; l6 = j6 
          
     ! find lowest possible momentums
     j3min = max(  abs( j1 -  j2)  , abs( j4 -  j5)    )/2               
     
     j6min = max(  abs( j1 -  j5)  , abs( j4 -  j2 )   )/2               
    

     ! make sure j1 < = j2
     if  ( j1 > j2 ) then 
        l1 = j2 ; l2 = j1 
        l4 = j5 ; l5 = j4
     end if 
     ! allowed by 6j symmetry
     
     ! find indeces
     x1 = bosonic_tp_index((l1+1)/2,(l2+1)/2,store6j%nhalf) 
     
     if (l4 > l5 )  then 
        x2 = fermionic_tp_index((l5+1)/2,(l4+1)/2,store6j%nhalf) 
     else
        x2 = bosonic_tp_index((l4+1)/2,(l5+1)/2,store6j%nhalf) + store6j%nf 
     end if 
     
     ! figure out indeces for j3,j6 based on lowest possible
     sixj = store6j%tp_mat(x1,x2)%X(l3/2 - j3min + 1 , l6/2 - j6min + 1 ) 
     
     
  else 
     sixj = 0.d0 
  end if 
end function
!=======================================================
!=======================================================
subroutine allocate_sp_mat(jbas,H) 
  ! allocate block sizes for sp matrix
  implicit none 
  
  type(spd) :: jbas
  type(full_sp_block_mat) :: H 
  integer :: jmax,lmax,m,blcks,d
  integer :: q,r,i,j ,l,tz
  
  ! how many blocks are there 
  jmax = jbas%Jtotal_max
  lmax = jbas%lmax
    
  blcks = (jmax+1)*(lmax+1) ! there is a factor of 2 for tz , which is canceled by 
                       ! the factor of 1/2 introduced by the fact that j is half-integer
  H%blocks = blcks
  
  allocate(H%blkM(blcks))
  allocate(H%map(blcks)) 
  
  ! find how many states in each block, allocate
  q = 1
  do tz = -1,1,2
     do l = 0,lmax
        do j = 1,jmax,2
           H%blkM(q)%lmda(1) = l
           H%blkM(q)%lmda(2) = j 
           H%blkM(q)%lmda(3) = tz 

           
           r = 0
           do i = 1, jbas%total_orbits
              if (l == jbas%ll(i)) then
                 if (j == jbas%jj(i)) then 
                    if (tz == jbas%itzp(i)) then 
                       r = r + 1
                    end if 
                 end if
              end if
           end do
        
           H%map(q) = r
           allocate(H%blkM(q)%matrix(r,r)) 
           allocate(H%blkM(q)%eigval(r))
           allocate(H%blkM(q)%extra(10*r)) 
           allocate(H%blkM(q)%states(r)) 
        
           r = 1
           do i = 1, jbas%total_orbits
              if (l == jbas%ll(i)) then
                 if (j == jbas%jj(i)) then 
                    if (tz == jbas%itzp(i)) then 
                       H%blkM(q)%states(r) = i 
                       r = r + 1
                    end if 
                 end if
              end if
           end do
        
           
           q = q + 1
        end do
     end do
  end do 
end subroutine  
!==========================================
!==========================================  
subroutine duplicate_sp_mat(H,T) 
  ! use all of the information about H to make 
  ! an empty block matrix of the same size
  implicit none 
  
  type(full_sp_block_mat) :: H,T 
  integer :: q,r
  
  T%blocks = H%blocks
  
  allocate(T%map(H%blocks)) 
  allocate(T%blkM(H%blocks)) 
  
  T%map = H%map 
  do q = 1, H%blocks
     T%blkM(q)%lmda = H%blkM(q)%lmda
     T%map(q) = H%map(q) 
     r = H%map(q) 
     allocate(T%blkM(q)%matrix(r,r)) 
     allocate(T%blkM(q)%eigval(r)) 
     allocate(T%blkM(q)%extra(10*r)) 
     allocate(T%blkM(q)%states(r))
     T%blkM(q)%states = H%blkM(q)%states
  end do 

end subroutine
!=====================================================
integer function bosonic_tp_index(i,j,n)
  ! n is total number of sp states
  ! assume i <= j 
  implicit none 
  
  integer :: i,j,n
  
  bosonic_tp_index = n*(i-1) + (3*i-i*i)/2 + j - i 
  
end function 
!=====================================================  
integer function fermionic_tp_index(i,j,n)
  ! n is total number of sp states
  ! assume i <= j 
  implicit none 
  
  integer :: i,j,n
  
  fermionic_tp_index = n*(i-1) + (i-i*i)/2 + j - i 
  
end function 
!=====================================================  
end module
       
  
