module basic_IMSRG 
  ! contains basic type declarations, constants, and adminstrative routines
  integer,public,dimension(9) :: adjust_index = (/0,0,0,0,0,0,0,0,-4/) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  TYPE :: spd    ! single particle discriptor
     INTEGER :: total_orbits,Jtotal_max
     INTEGER, ALLOCATABLE,DIMENSION(:) :: nn, ll, jj, itzp, nshell, mvalue,con
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
  TYPE :: sq_block
     integer :: lam(3) ! specifices J,Par,Tz of the block 
     type(real_mat),dimension(6) :: gam !Vpppp,Vhhhh,Vphph,Vpphh,Vphhh,Vpppph
     integer :: npp, nph , nhh ! dimensions
     type(int_mat),dimension(3) :: qn !tp to sp map 
     type(int_vec),dimension(3) :: pnt ! for mapping from sp to tp (confusing) 
  END TYPE sq_block
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TYPE :: sq_op !second quantized operator 
     type(sq_block),allocatable,dimension(:) :: mat
     real(8),allocatable,dimension(:,:) :: fph,fpp,fhh
     integer :: nblock,nbody,msp,herm
     real(8) :: E0 
  END TYPE sq_op
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

contains 
!====================================================
!====================================================
subroutine read_sp_basis(jbas,sp_input_file,holes)
  ! fills jscheme_basis array with single particle data from file
  ! file format must be: 5 integers and one real 
  ! for each state we have:   | label |  n  | l | 2 * j | 2*tz | E_sp |  
  implicit none 
  
  type(spd) :: jbas
  character(50) :: sp_input_file
  integer :: ist,i,label,ni,li,ji,tzi,ix,holes
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
  
  call find_holes(jbas,holes) 
  
  jbas%Jtotal_max = maxval(jbas%jj) 
  ! this is the maximum value of J, this code always uses 2*J 

  close(39) 
  
end subroutine 
!==============================================
!==============================================
subroutine find_holes(jbas,holes) 
  implicit none 
  
  type(spd) :: jbas
  integer :: holes,i,minpos(1)
  real(8),dimension(jbas%total_orbits) :: temp
  
  temp = jbas%e
  jbas%con = 0 ! zero if particle, one if hole
  
  do i = 1,holes
     minpos = minloc(temp) !find lowest energy 
     jbas%con(minpos(1)) = 1 ! thats a hole
     temp(minpos(1)) = 9.e9  ! replace it with big number
  end do 
     
end subroutine  
!==============================================
!==============================================
subroutine allocate_blocks(jbas,op) 
  ! allocate the sq_op arrays
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: op
  integer :: n,m,q,i,j,j1,j2,tz1,tz2,l1,l2
  integer :: Jtot,Tz,Par,nph,npp,nhh,CX
  
 
  n = op%nbody !number of holes
  op%msp = jbas%total_orbits
  m = op%msp  !number of sp states
  
  op%nblock =  (jbas%Jtotal_max + 1) * 6    ! 6 possible values of par * Tz    
  
  allocate(op%mat(op%nblock)) 
  allocate(op%fpp(m-n,m-n))
  allocate(op%fph(m-n,n))
  allocate(op%fhh(n,n)) 
  
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
     
     do i = 1,m   !looping over sp states
        do j = i,m
 
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
     
     nph = 0 ; npp = 0 ; nhh = 0
     do i = 1,m   !looping over sp states
        do j = i,m
           
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
                 op%mat(q)%pnt(1)%Z(npp) = m*(i-1) + j !search for an integer 
              case (1)                                 ! instead of a pair of them  
                 nph = nph + 1
                 op%mat(q)%qn(2)%Y(nph,1) = i
                 op%mat(q)%qn(2)%Y(nph,2) = j 
                 op%mat(q)%pnt(2)%Z(nph) = m*(i-1) + j 
              case (2) 
                 nhh = nhh + 1
                 op%mat(q)%qn(3)%Y(nhh,1) = i
                 op%mat(q)%qn(3)%Y(nhh,2) = j 
                 op%mat(q)%pnt(3)%Z(nhh) = m*(i-1) + j 
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
  integer :: ist,J,Tz,Par,a,b,c,d,q,qx
  real(8) :: V,g1,g2,g3
  integer :: C1,C2,int1,int2,i1,i2
  
  open(unit=39,file = '../TBME_input/'//trim(adjustl(intfile))) 
  
  read(39,*);read(39,*);read(39,*);read(39,*)
  read(39,*);read(39,*);read(39,*);read(39,*) !skip all of the garbage 
  
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
     
     i1 = 1
     i2 = 1
     int1 = H%msp*(a-1) + b
     int2 = H%msp*(c-1) + d
     
     ! find the indeces
     print*, c1,c2
     do 
        if (int1 == H%mat(q)%pnt(C1)%Z(i1) ) exit
        i1 = i1 + 1
     end do 
        
     do 
        if (int2 == H%mat(q)%pnt(C2)%Z(i2) ) exit
        i2 = i2 + 1
     end do 
     
     
     print*, q,qx,i1,i2,J,Tz,Par,a,b,c,d
     ! always make sure the input file has everything ordered so that a<=b, c<=d 
     If (C1>C2) then 
        H%mat(q)%gam(qx)%X(i2,i1)  = V 
     else
        H%mat(q)%gam(qx)%X(i1,i2) = V
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
end module
       
  
