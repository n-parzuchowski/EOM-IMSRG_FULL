module basic_IMSRG 
  ! contains basic type declarations, constants, and adminstrative routines  
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
  TYPE :: spd    ! single particle discriptor
     INTEGER :: total_orbits,Jtotal_max,lmax,spblocks
     INTEGER, ALLOCATABLE,DIMENSION(:) :: nn, ll, jj, itzp, nshell, mvalue
     INTEGER, ALLOCATABLE,DIMENSION(:) :: con,holesb4,partsb4,holes,parts 
     type(int_vec), allocatable,dimension(:) :: states
     ! for clarity:  nn, ll, nshell are all the true value
     ! jj is j+1/2 (so it's an integer) 
     ! likewise itzp is 2*tz  
     REAL(8), ALLOCATABLE,DIMENSION(:) :: e
  END TYPE spd
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TYPE :: TPD 
     integer,allocatable,dimension(:) :: direct_omp
     integer,dimension(3) :: chan
     integer,allocatable,dimension(:,:) :: ppp,hhh 
  END TYPE TPD
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TYPE ::  six_index_store
     type(real_mat),allocatable,dimension(:,:) :: tp_mat
     integer :: nf,nb,nhalf
  END TYPE six_index_store
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TYPE :: sq_block
     integer :: lam(3) ! specifices J,Par,Tz of the block
     type(real_mat),dimension(6) :: gam !Vpppp,Vhhhh,Vphph,Vpphh,Vphhh,Vppph
     integer :: npp, nph , nhh,ntot! dimensions
     type(int_mat),dimension(3) :: qn !tp to sp map 
  END TYPE sq_block
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TYPE :: tensor_block
     integer :: lam(3) ! specifices J,Par,Tz of the block
     integer :: Jpair(2) ! for tensor operators that have rank > 0  
     type(real_mat),dimension(9) :: tgam
     integer :: npp1, nph1 , nhh1 ,npp2,nph2,nhh2,ntot! dimensions
     type(int_mat),dimension(3,2) :: tensor_qn
  END TYPE tensor_block
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TYPE :: sq_op !second quantized operator 
     type(sq_block),allocatable,dimension(:) :: mat
     type(tensor_block),allocatable,dimension(:) :: tblck
     type(int_vec),allocatable,dimension(:) :: xmap 
     real(8),allocatable,dimension(:,:) :: fph,fpp,fhh
     integer,allocatable,dimension(:,:) :: exlabels
     integer,allocatable,dimension(:) :: direct_omp 
     integer :: nblocks,Aprot,Aneut,Nsp,herm,belowEF,neq
     integer :: Jtarg,Ptarg,valcut,Rank,dpar
     real(8) :: E0,hospace,lawson_beta,com_hw 
     logical :: pphh_ph
  END TYPE sq_op
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  TYPE :: mscheme_3body !second quantized operator 
     type(real_mat),allocatable,dimension(:) :: Wmat
     type(int_vec),allocatable,dimension(:) :: lam
     type(int_mat),allocatable,dimension(:) :: qn_h,qn_p
     integer,allocatable,dimension(:,:) :: qnm_holes,qnm_parts
     integer,allocatable,dimension(:) :: direct_omp 
     integer :: nblocks,Aprot,Aneut,Nsp,herm,belowEF,neq
     integer :: Jtarg,Ptarg,valcut,oc_m,un_m
     real(8) :: E0,hospace
  END TYPE mscheme_3body 
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
type cross_coupled_31_mat
   type(int_vec),allocatable,dimension(:) :: rmap,qmap,nbmap
   type(real_mat),allocatable,dimension(:) :: CCX,CCR
   integer,allocatable,dimension(:) :: Jval,Jval2,nph,rlen
   integer :: nblocks,Nsp,rank,herm,dpar
end type cross_coupled_31_mat
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

   ! THIS STUFF IS IMPORTANT, DON'T CHANGE IT.
  ! I try to keep public stuff to a minimum, and only use it where absolutely necessary 
   real(8),public,parameter :: al = 1.d0 , bet = 0.d0  ! parameters for dgemm that NEVER change
   real(8),allocatable,dimension(:,:) :: phase_pp,phase_hh 
   integer,public,dimension(9) :: adjust_index = (/0,0,0,0,0,0,0,0,-4/) ! for finding which of the 6 Vpppp arrays
   integer,public,dimension(6) :: tensor_adjust = (/0,6,4,0,0,3/)
   type(six_index_store),public :: store6j , half6j  ! holds 6j-symbols
   !The following public arrays give info about the 6 different categories
   ! of matrix elements: Vpppp, Vppph , Vpphh , Vphph , Vhhhh, Vphhh 
   ! holds the c values for qn and pn arrays
   integer,public,dimension(9) :: sea1 = (/1,1,1,2,3,2,3,2,3/), sea2 = (/1,2,3,2,3,3,1,1,2/) 
   ! true if square matrix
   logical,public,dimension(9) :: sqs = (/.true.,.false.,.false.,.true.,.true.,.false.,.false.,.false.,.false./)
   ! 100000 if square matrix, 1 if not. 
   integer,public,dimension(9) :: jst = (/10000000,1,1,10000000,10000000,1,1,1,1/)
   integer,public :: global_counter1=0,global_counter2=0,global_counter3=0
   character(500) :: TBME_DIR,SP_DIR,INI_DIR
   


   ! CHANGE THESE AS NEEDED. ==========================================================
   real(8),public,parameter :: hbarc = 197.326968d0, m_nuc = 938.918725 !2006 values 
   real(8),public,parameter :: hbarc2_over_mc2 = hbarc*hbarc/m_nuc
   real(8),public,parameter :: Pi_const = acos(-1.d0) 

   character(500) :: OUTPUT_DIR = '/home/nathan/nuclear_IMSRG/output/'   
   character(500),dimension(1) :: TBME_DIRECTORY_LIST=(/ '/home/nathan/nuclear_IMSRG/TBME_input/' /)  
   character(500),dimension(1) :: SP_DIRECTORY_LIST=(/ '/home/nathan/nuclear_IMSRG/sp_inputs/' /)  
   character(500),dimension(1) :: INI_DIRECTORY_LIST=(/ '/home/nathan/nuclear_IMSRG/inifiles/' /) 
   
contains
!====================================================
!====================================================
subroutine read_sp_basis(jbas,hp,hn,method)
  ! fills jscheme_basis array with single particle data from file
  ! file format must be: 5 integers and one real 
  ! for each state we have:   | label |  n  | l | 2 * j | 2*tz | E_sp |  
  implicit none 
  
  type(spd) :: jbas
  character(200) :: spfile,intfile,prefix
  character(2) :: hk
  character(200) :: interm
  integer :: ist,i,label,ni,li,ji,tzi,ix,hp,hn
  integer :: q,r,Jtarget,PARtarget,method
  real(8) :: e
  common /files/ spfile,intfile,prefix
  
  interm= adjustl(spfile)
  open(unit=39,file=trim(SP_DIR)//trim(interm))
  hk=interm(1:2)

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
 
  call find_holes(jbas,hp,hn,hk) 
  
  jbas%Jtotal_max = maxval(jbas%jj) 
  jbas%lmax = maxval(jbas%ll) 
  ! this is the maximum value of J, this code always uses 2*J 
  call store_6j(jbas,method) ! save six-j symbols to an array
  close(39) 
  
  jbas%spblocks =  (jbas%Jtotal_max + 1) * (jbas%lmax + 1)
  allocate(jbas%states(jbas%spblocks)) 

  ! this set of loops fills jbas%states
  ! which tells us which sp states 
  ! are contained in each block
  q = 1
  do tzi = -1,1,2
     do li = 0,jbas%lmax
        do ji = 1,jbas%Jtotal_max,2
           
           ! count states
           r = 0
           do i = 1, jbas%total_orbits
              
              if (jbas%jj(i) .ne. ji) cycle
              if (jbas%ll(i) .ne. li) cycle
              if (jbas%itzp(i) .ne. tzi) cycle 
              
              r = r + 1
           end do 
           
           allocate(jbas%states(q)%Z(r)) 
           
           ! fill array
           r = 0
           do i = 1, jbas%total_orbits
              
              if (jbas%jj(i) .ne. ji) cycle
              if (jbas%ll(i) .ne. li) cycle
              if (jbas%itzp(i) .ne. tzi) cycle 
              
              r = r + 1
              jbas%states(q)%Z(r) = i 
           end do 
           q = q + 1
        end do
     end do 
  end do 

end subroutine  
!==============================================
!==============================================
subroutine find_holes(jbas,pholes,nholes,hk) 
  implicit none 
  
  type(spd) :: jbas
  integer :: pholes,nholes,i,minpos(1),rn,rp,r1,r2
  real(8),dimension(jbas%total_orbits) :: temp
  character(2) :: hk
  

  temp = jbas%e
  jbas%con = 0 ! zero if particle, one if hole
  
  rn = 0
  rp = 0 
  ! protons have isospin -1/2 
 
  if (pholes==nholes) then 

     if (pholes == 2) then 
        jbas%con(1:2) = 1
     else if (pholes == 8) then 
        jbas%con(1:6) = 1
     else if (pholes == 20) then 
        jbas%con(1:12) = 1
     else if (pholes == 28) then 
        if (hk == 'hk') then 
           jbas%con(1:12) = 1
           jbas%con(19:20) = 1
        else
           jbas%con(1:14) = 1 
        end if    
     else 
        STOP 'this nucleus is not available' 
     end if 

   else
     if (pholes == 2) then 
        jbas%con(1:2) = 1
     else if (pholes == 8) then 
        jbas%con(1:6) = 1
        
        if (nholes == 14) then 
           if (hk == 'hk') then 
              jbas%con(12) = 1
           else 
              jbas%con(8) = 1 
           end if
        else if (nholes == 16) then
           jbas%con(12) = 1
           jbas%con(8) = 1
        else 
           STOP 'this nucleus is not available' 
        end if
     
     else if (pholes == 20) then 
        jbas%con(1:12) = 1
        if (nholes == 28) then 
           if (hk=='hk') then 
              jbas%con(20) = 1
           else 
              jbas%con(14) = 1
           end if 
         else 
            STOP 'this nucleus is not available' 
         end if 
     else if (pholes == 20) then 
        jbas%con(1:12) = 1
        if (nholes == 28) then 
           if (hk=='hk') then 
              jbas%con(20) = 1
           else 
              jbas%con(14) = 1
           end if 
         else 
            STOP 'this nucleus is not available' 
         end if 
      else if (pholes == 28) then 
        jbas%con(1:12) = 1
        if (nholes == 20) then 
           if (hk=='hk') then 
              jbas%con(19) = 1
           else 
              jbas%con(13) = 1
           end if  
        else 
            STOP 'this nucleus is not available' 
         end if 
         
      else if (pholes == 14) then 
         if (nholes == 8) then 
            jbas%con(1:6) = 1
            if (hk=='hk')then 
               jbas%con(11) = 1
            else
               jbas%con(7) = 1
            end if 
         end if
      else 
        STOP 'this nucleus is not available' 
     end if 
  end if 
    
  allocate(jbas%holes(sum(jbas%con)))
  allocate(jbas%parts(jbas%total_orbits - sum(jbas%con))) 
  ! these arrays help us later in f_elem (in this module) 
  r1 = 1
  r2 = 1
  do i = 1, jbas%total_orbits
     jbas%holesb4(i) = sum(jbas%con(1:i-1)) 
     jbas%partsb4(i) = i - jbas%holesb4(i) - 1 
     
     ! write down the position of the holes and parts
     if (jbas%con(i) == 1) then 
        jbas%holes(r1) = i 
        r1 = r1 + 1
     else
        jbas%parts(r2) = i 
        r2 = r2 + 1
     end if 
     
  end do 

end subroutine  
!==============================================
!==============================================
subroutine allocate_blocks(jbas,op) 
  ! allocate the sq_op arrays
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: op
  integer :: N,AX,q,i,j,j1,j2,tz1,tz2,l1,l2,x
  integer :: Jtot,Tz,Par,nph,npp,nhh,CX,j_min,j_max,numJ
  
  AX = sum(jbas%con) !number of shells under eF 
  op%belowEF = AX
  op%Nsp = jbas%total_orbits
  op%rank = 0 !scalar operator
  op%dPar = 0 
  op%pphh_ph = .false. 
  N = op%Nsp  !number of sp shells
  
  op%nblocks =  (jbas%Jtotal_max + 1) * 6 ! -3  ! 6 possible values of par * Tz  
                                               ! except for Jtot=Jmax has 3
  allocate(op%mat(op%nblocks)) 
  allocate(op%fpp(N-AX,N-AX))
  allocate(op%fph(N-AX,AX))
  allocate(op%fhh(AX,AX)) 
  
  op%fpp=0.d0
  op%fhh=0.d0
  op%fph=0.d0
  q = 1  ! block index

  allocate(op%xmap(N*(N+1)/2)) 
 
  ! allocate the map array, which is used by 
  ! v_elem to find matrix elements
  do i = 1,N
     do j = i,N
        
        j_min = abs(jbas%jj(i)-jbas%jj(j))
        j_max = jbas%jj(i) + jbas%jj(j) 
  
        numJ = (j_max - j_min)/2 + 2
  
        x = bosonic_tp_index(i,j,N) 
        allocate(op%xmap(x)%Z(numJ)) 
        op%xmap(x)%Z = 0
        op%xmap(x)%Z(1) = j_min
        
     end do 
  end do 
        
  do Jtot = 0,2*jbas%Jtotal_max,2 !looping over blocks
     do Tz = -1,1
        do Par = 0,1
           !if ((Jtot == 2*jbas%Jtotal_max) .and. (Par==1)) cycle
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
     
     allocate(op%mat(q)%gam(1)%X(npp,npp)) !Vpppp
     allocate(op%mat(q)%gam(5)%X(nhh,nhh)) !Vhhhh
     allocate(op%mat(q)%gam(3)%X(npp,nhh)) !Vpphh
     allocate(op%mat(q)%gam(4)%X(nph,nph)) !Vphph
     allocate(op%mat(q)%gam(2)%X(npp,nph)) !Vppph
     allocate(op%mat(q)%gam(6)%X(nph,nhh)) !Vphhh
     do i = 1,6
        op%mat(q)%gam(i)%X = 0.d0
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
           
           x = bosonic_tp_index(i,j,N) 
           j_min = op%xmap(x)%Z(1) 
           select case (CX)
              case (0) 
                 npp = npp + 1
                 op%mat(q)%qn(1)%Y(npp,1) = i
                 op%mat(q)%qn(1)%Y(npp,2) = j
                 op%xmap(x)%Z((Jtot-j_min)/2+2) = npp 
              case (1)                      
                 nph = nph + 1
                 op%mat(q)%qn(2)%Y(nph,1) = i
                 op%mat(q)%qn(2)%Y(nph,2) = j 
                 op%xmap(x)%Z((Jtot-j_min)/2+2) = nph
              case (2) 
                 nhh = nhh + 1
                 op%mat(q)%qn(3)%Y(nhh,1) = i
                 op%mat(q)%qn(3)%Y(nhh,2) = j 
                 op%xmap(x)%Z((Jtot-j_min)/2+2) = nhh
           end select
                     
        end do 
     end do    
     q = q + 1
        end do
     end do
  end do
  
  call divide_work(op) 
end subroutine
!==================================================================
!==================================================================
subroutine allocate_tensor(jbas,op,zerorank) 
  ! allocate the sq_op arrays
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: op,zerorank 
  integer :: N,AX,q,i,j,j1,j2,tz1,tz2,l1,l2,x,rank
  integer :: Jtot1,Jtot2,Jtot,Tz,Par1,Par2,nph1,npp1,nhh1,nph2,npp2,nhh2
  integer :: CX,j_min,j_max,numJ,q1,q2
  real(8) :: d6ji,lwake

  ! rank is multiplied by 2 as well
  rank = op%rank 
  op%hospace = zerorank%hospace
  if ( mod(rank,2) == 1 ) STOP 'RANK MUST BE MULTIPLIED BY 2'

  AX = sum(jbas%con) !number of shells under eF 
  op%belowEF = AX
  op%Nsp = jbas%total_orbits
  N = op%Nsp  !number of sp shells

  ! allocate the map array, which is used by 
  ! v_elem to find matrix elements
  allocate(op%xmap(N*(N+1)/2)) 
  do i = 1,N
     do j = i,N
        
        j_min = abs(jbas%jj(i)-jbas%jj(j))
        j_max = jbas%jj(i) + jbas%jj(j) 
  
        numJ = (j_max - j_min)/2 + 2
  
        x = bosonic_tp_index(i,j,N) 
        allocate(op%xmap(x)%Z(numJ)) 
        op%xmap(x)%Z = 0
        op%xmap(x)%Z(1) = j_min
        
     end do 
  end do 
  
  ! quantum numbers of the last block 
  Tz = 1 
  Par1 = 1
  Jtot1 = jbas%Jtotal_max*2 
  Jtot2 = Jtot1+ RANK 
  op%nblocks =  tensor_block_index(Jtot1,Jtot2,RANK,Tz,Par1)               
 
  allocate(op%tblck(op%nblocks)) 
  allocate(op%fpp(N-AX,N-AX))
  allocate(op%fph(N-AX,AX))
  allocate(op%fhh(AX,AX)) 
 
  op%fpp = 0.d0
  op%fph = 0.d0
  op%fhh = 0.d0 

  q = 1
  do Jtot1 = 0,2*jbas%Jtotal_max,2 
     do Jtot2 = max(abs(Jtot1 - rank),Jtot1), Jtot1+rank, 2   
       
        do Tz = -1,1
           do Par1 = 0,1
                    
              if (mod(op%dpar/2,2) == 1) then ! this only works for EX transitions
                 Par2 = abs(Par1-1)  
              else
                 Par2 = Par1 
              end if 
              
              op%tblck(q)%Jpair(1) = Jtot1
              op%tblck(q)%Jpair(2) = Jtot2
        
              op%tblck(q)%lam(1) = (-1)**((Jtot1+Jtot2)/2) ! phase instead of J 
              op%tblck(q)%lam(2) = Par1 !just remember that they change if the operator has odd parity.
              op%tblck(q)%lam(3) = Tz
              
       !       q = tensor_block_index(Jtot1,Jtot2,RANK,Tz,Par1)              
              q1 = block_index(Jtot1,Tz,Par1) 
              q2 = block_index(Jtot2,Tz,Par2) 
              
              ! we already figured this stuff out for 
              ! the rank 0 case, so lets re-use it
              
              if ( Jtot2 > Jbas%Jtotal_max*2 ) then 
                 npp1=0;nph1=0;nhh1=0
                 npp2=0;nph2=0;nhh2=0
              else
                 npp1 = zerorank%mat(q1)%npp
                 nph1 = zerorank%mat(q1)%nph
                 nhh1 = zerorank%mat(q1)%nhh
                 npp2 = zerorank%mat(q2)%npp
                 nph2 = zerorank%mat(q2)%nph
                 nhh2 = zerorank%mat(q2)%nhh
              end if 
              
              op%tblck(q)%npp1 = npp1
              op%tblck(q)%nph1 = nph1
              op%tblck(q)%nhh1 = nhh1              
              op%tblck(q)%npp2 = npp2
              op%tblck(q)%nph2 = nph2
              op%tblck(q)%nhh2 = nhh2              
                            
              allocate(op%tblck(q)%tensor_qn(1,1)%Y(npp1,2)) !qnpp1
              allocate(op%tblck(q)%tensor_qn(2,1)%Y(nph1,2)) !qnph1
              allocate(op%tblck(q)%tensor_qn(3,1)%Y(nhh1,2)) !qnhh1
              
              allocate(op%tblck(q)%tensor_qn(1,2)%Y(npp2,2)) !qnpp2
              allocate(op%tblck(q)%tensor_qn(2,2)%Y(nph2,2)) !qnph2
              allocate(op%tblck(q)%tensor_qn(3,2)%Y(nhh2,2)) !qnhh2

              if ( Jtot2 .le. Jbas%Jtotal_max*2 ) then 

              ! yeah this is a mess but it really facilitates the 
              ! commutators

              ! blocks such as pphh ALWAYS have pp first then hh
              ! looks scary but i'm just copying everything from the
              ! rank zero operators we already have
              op%tblck(q)%tensor_qn(1,1)%Y = zerorank%mat(q1)%qn(1)%Y
              op%tblck(q)%tensor_qn(2,1)%Y = zerorank%mat(q1)%qn(2)%Y
              op%tblck(q)%tensor_qn(3,1)%Y = zerorank%mat(q1)%qn(3)%Y
         
              op%tblck(q)%tensor_qn(1,2)%Y = zerorank%mat(q2)%qn(1)%Y
              op%tblck(q)%tensor_qn(2,2)%Y = zerorank%mat(q2)%qn(2)%Y
              op%tblck(q)%tensor_qn(3,2)%Y = zerorank%mat(q2)%qn(3)%Y
              
              end if 
              
              allocate(op%tblck(q)%tgam(1)%X(npp1,npp2)) !Vpppp
              allocate(op%tblck(q)%tgam(5)%X(nhh1,nhh2)) !Vhhhh
              allocate(op%tblck(q)%tgam(3)%X(npp1,nhh2)) !Vpphh
              allocate(op%tblck(q)%tgam(4)%X(nph1,nph2)) !Vphph
              allocate(op%tblck(q)%tgam(2)%X(npp1,nph2)) !Vppph
              allocate(op%tblck(q)%tgam(6)%X(nph1,nhh2)) !Vphhh
             
              allocate(op%tblck(q)%tgam(7)%X(nhh1,npp2)) !Vhhpp
              allocate(op%tblck(q)%tgam(8)%X(nph1,npp2)) !Vphpp
              allocate(op%tblck(q)%tgam(9)%X(nhh1,nph2)) !Vhhph
        
               do i = 1,9
                  op%tblck(q)%tgam(i)%X = 0.d0
               end do
  
               q = q + 1
            
            end do
         end do
         
      end do
   end do
 
  
   ! i need to fill this last array for tensor_elem
   do q = 1, zerorank%nblocks
   nph1 = 0 ; npp1 = 0 ; nhh1 = 0
   do i = 1,N   !looping over sp states
      do j = i,N

         Jtot = zerorank%mat(q)%lam(1)
         Tz = zerorank%mat(q)%lam(3)
         Par1 = zerorank%mat(q)%lam(2)
         ! check if the two sp states can exist in this block 
         tz1 = jbas%itzp(i)
         tz2 = jbas%itzp(j) 
         if ( Tz .ne. (tz1 + tz2)/2 ) cycle
         
         l1 = jbas%ll(i) 
         l2 = jbas%ll(j)
         if ( Par1 .ne. (1 - (-1)**(l1 + l2))/2 ) cycle
         
         j1 = jbas%jj(i) 
         j2 = jbas%jj(j) 
         if (.not. triangle(j1,j2,Jtot) ) cycle
                                                      
         cX = jbas%con(i) + jbas%con(j)
         
         x = bosonic_tp_index(i,j,N) 
         j_min = op%xmap(x)%Z(1) 
       
      
         select case (CX)
         case (0) 
            npp1 = npp1 + 1
            op%xmap(x)%Z((Jtot-j_min)/2+2) = npp1 
         case (1)                      
            nph1 = nph1 + 1
            op%xmap(x)%Z((Jtot-j_min)/2+2) = nph1
         case (2) 
            nhh1 = nhh1 + 1
            op%xmap(x)%Z((Jtot-j_min)/2+2) = nhh1
         end select
        
      end do
   end do
   end do 
              
   if (allocated( phase_hh) ) return 
   ! THESE ARE PUBLIC ARRAYS NEEDED FOR THE STORAGE OF 
   ! THE PHASE OF EACH TENSOR COMPONENT
   
   allocate(phase_hh(op%belowEF,op%belowEF))
   allocate(phase_pp(op%nsp-op%belowEF,op%nsp-op%belowEF))

   do i = 1, op%belowEF
      do j = i , op%belowEF 
         
         j1 = jbas%jj(jbas%holes(i))
         j2 = jbas%jj(jbas%holes(j))
         
         phase_hh(i,j) = (-1)**((j1-j2)/2) 
         phase_hh(j,i) = phase_hh(i,j) 
      end do 
   end do 
   
   do i = 1,op%nsp- op%belowEF
      do j = i ,op%nsp-op%belowEF                   
         j1 = jbas%jj(jbas%parts(i))
         j2 = jbas%jj(jbas%parts(j))
         
         phase_pp(i,j) = (-1)**((j1-j2)/2) 
         phase_pp(j,i) = phase_pp(i,j) 
      end do 
   end do 
   
   ! these are six j symbols that I don't already have,
   ! which the commutators need for this tensor.
   ! access with XXXsixj
   call store_6j_3halfint(jbas,rank)
  
   
   call divide_work_tensor(op) 

 end subroutine allocate_tensor
!==================================================================  
!==================================================================
subroutine divide_work(r1) 
  implicit none 
  
  type(sq_op) :: r1
  integer :: A,N,threads,omp_get_num_threads
  integer :: i ,g,q,k,b,j
  
!$omp parallel
  threads=omp_get_num_threads() 
!$omp end parallel
!threads = 1
  b = 0.d0
  do q = 1, r1%nblocks
     b = b + r1%mat(q)%nhh +r1%mat(q)%npp + r1%mat(q)%nph 
  end do
  
  allocate( r1%direct_omp(threads+1) )
  
  g = 0
  k = 0
  q = 1
  do i = 1,threads
     
     g = g+k 
     j = (b-g)/(threads-i+1) 
     
     k = 0
     
     do while ( q .le. r1%nblocks) 
        k = k + r1%mat(q)%nhh +r1%mat(q)%npp + r1%mat(q)%nph 
        q = q + 1
        
        if (k .ge. j) then 
           r1%direct_omp(i+1) = q-1
           exit
        end if
     end do 
  end do 
  r1%direct_omp(threads+1) = r1%nblocks
  r1%direct_omp(1) = 0 

end subroutine     
!==================================================================  
!==================================================================
subroutine divide_work_tensor(r1) 
  implicit none 
  
  type(sq_op) :: r1
  integer :: A,N,threads,omp_get_num_threads
  integer :: i ,g,q,k,b,j,spot
  
!$omp parallel
  threads=omp_get_num_threads() 
!$omp end parallel
!threads = 1
  b = 0.d0
  do q = 1, r1%nblocks
     do spot = 1, 9
        b = b + size(r1%tblck(q)%tgam(spot)%X) 
     end do 
  end do
  
  allocate( r1%direct_omp(threads+1) )
  
  g = 0
  k = 0
  q = 1
  do i = 1,threads
     
     g = g+k 
     j = (b-g)/(threads-i+1) 
     
     k = 0
     
     do while ( q .le. r1%nblocks) 
        do spot = 1, 9
           k = k + size(r1%tblck(q)%tgam(spot)%X) 
        end do
        q = q + 1
        
        if (k .ge. j) then 
           r1%direct_omp(i+1) = q-1
           exit
        end if
     end do 
  end do 
  r1%direct_omp(threads+1) = r1%nblocks
  r1%direct_omp(1) = 0 

end subroutine     
!==================================================================  
!==================================================================
subroutine divide_work_tpd(threebas) 
  implicit none 
  
  type(tpd),dimension(:) :: threebas
  integer :: A,N,threads,omp_get_num_threads
  integer :: i ,g,q,k,b,j,blocks
  
!$omp parallel
  threads=omp_get_num_threads() 
!$omp end parallel
!threads = 1

  blocks =size(threebas) 
  b = 0.d0
  do q = 1,blocks 
     b = b + size(threebas(q)%hhh(:,1))*size(threebas(q)%ppp(:,1)) 
  end do
  
  allocate( threebas(1)%direct_omp(threads+1) )
  
  g = 0
  k = 0
  q = 1
  do i = 1,threads
     
     g = g+k 
     j = (b-g)/(threads-i+1) 
     
     k = 0
     
     do while ( q .le. blocks) 
        k = k + size(threebas(q)%hhh(:,1))*size(threebas(q)%ppp(:,1)) 
        q = q + 1
        
        if (k .ge. j) then 
           threebas(1)%direct_omp(i+1) = q-1
           exit
        end if
     end do 
  end do 
  threebas(1)%direct_omp(threads+1) = blocks
  threebas(1)%direct_omp(1) = 0 

end subroutine     
!==================================================================  
!==================================================================
subroutine read_interaction(H,jbas,htype,hw,rr,pp) 
  ! read interaction from ASCII file produced by Scott_to_Morten.f90 
  implicit none
  
  type(sq_op) :: H
  type(sq_op),optional :: rr,pp
  type(spd) :: jbas
  character(200) :: spfile,intfile,prefix
  integer :: ist,J,Tz,Par,a,b,c,d,q,qx,N,j_min,x
  real(8) :: V,Vcm,g1,g2,g3,pre,hw,V1
  integer :: C1,C2,int1,int2,i1,i2,htype,COM
  logical :: rr_calc,pp_calc
  common /files/ spfile,intfile,prefix
  
  open(unit=39,file = trim(TBME_DIR)//trim(adjustl(intfile))) 
  
  read(39,*);read(39,*);read(39,*);read(39,*)
  read(39,*);read(39,*);read(39,*);read(39,*) !skip all of the garbage 
  
  COM = 0
  if (htype == 1) COM = 1 ! remove center of mass hamiltonian? 
  
  N = jbas%total_orbits
  
  ! check if we are concerned with other operators
  rr_calc = .false.
  pp_calc = .false. 
  if (present(rr)) rr_calc= .true. 
  if (present(pp)) pp_calc= .true. 

  do 
     read(39,*,iostat=ist) Tz,Par,J,a,b,c,d,V,g1,g2,g3
     !read(39,*) Tz,Par,J,a,b,c,d,V,g1,g2,g3
     
     ! g1 is COM expectation value, NOT CALCULATED WITH SCOTT'S CODE
     ! g2 is r1*r2 ME 
     ! g3 is p1*p2 ME 
     
     if (ist > 0) STOP 'interaction file error' 
     if (ist < 0) exit

     if (rr_calc) then 
        g2 = r1_r2( a, b, c, d, J ,jbas ) ! morten and Koshiroh are inconsistent with their definitions
     end if    
     g3 = p1_p2( a, b, c, d, J ,jbas ) ! morten and Koshiroh are inconsistent with their definitions
     ! so I am doing it explicitly. 
     
     g1 = (g3 + H%com_hw**2 /hw**2 *g2)*H%lawson_beta ! lawson term    

     V = V + (g1 - g3*COM) *hw/(H%Aneut + H%Aprot) ! center of mass correction
 
     
     q = block_index(J,Tz,Par) 
     
     C1 = jbas%con(a)+jbas%con(b) + 1 !ph nature
     C2 = jbas%con(c)+jbas%con(d) + 1
    
     qx = C1*C2
     qx = qx + adjust_index(qx)   !Vpppp nature  

     ! get the indeces in the correct order
     pre = 1

     if ( a > b )  then 
        
        x = bosonic_tp_index(b,a,N) 
        j_min = H%xmap(x)%Z(1)  
        i1 = H%xmap(x)%Z( (J-j_min)/2 + 2) 
        pre = (-1)**( 1 + (jbas%jj(a) + jbas%jj(b) -J)/2 ) 
     else
       ! if (a == b) pre = pre * sqrt( 2.d0 )
       
        x = bosonic_tp_index(a,b,N) 
        j_min = H%xmap(x)%Z(1)  
        i1 = H%xmap(x)%Z( (J-j_min)/2 + 2) 
     end if
  
     if (c > d)  then     
        
        x = bosonic_tp_index(d,c,N) 
        j_min = H%xmap(x)%Z(1)  
        i2 = H%xmap(x)%Z( (J-j_min)/2 + 2) 
        
        pre = pre * (-1)**( 1 + (jbas%jj(c) + jbas%jj(d) -J)/2 ) 
     else 
       ! if (c == d) pre = pre * sqrt( 2.d0 )
      
        x = bosonic_tp_index(c,d,N) 
        j_min = H%xmap(x)%Z(1)  
        i2 = H%xmap(x)%Z( (J-j_min)/2 + 2) 
     end if
     ! kets/bras are pre-scaled by sqrt(2) if they 
     ! have two particles in the same sp-shell
        
     ! get the units right. I hope 
   
     
     if ((qx == 1) .or. (qx == 5) .or. (qx == 4)) then 
        H%mat(q)%gam(qx)%X(i2,i1)  = V *pre
        H%mat(q)%gam(qx)%X(i1,i2)  = V *pre
        
        if (rr_calc) then 
           rr%mat(q)%gam(qx)%X(i2,i1)  = hw*g2*pre/(H%Aneut + H%Aprot)
           rr%mat(q)%gam(qx)%X(i1,i2)  = hw*g2*pre/(H%Aneut + H%Aprot)
        end if 

        if (pp_calc) then 
           pp%mat(q)%gam(qx)%X(i2,i1)  = hw*g3*pre/(H%Aneut + H%Aprot)
           pp%mat(q)%gam(qx)%X(i1,i2)  = hw*g3*pre/(H%Aneut + H%Aprot)
        end if 

        
     else if (C1>C2) then
        H%mat(q)%gam(qx)%X(i2,i1)  = V *pre
        
        if (rr_calc) then 
           rr%mat(q)%gam(qx)%X(i2,i1)  = hw*g2*pre/(H%Aneut + H%Aprot) 
        end if
        
        if (pp_calc) then 
           pp%mat(q)%gam(qx)%X(i2,i1)  = hw*g3*pre/(H%Aneut + H%Aprot) 
        end if

     else
        H%mat(q)%gam(qx)%X(i1,i2) = V * pre
        
        if (rr_calc) then 
           rr%mat(q)%gam(qx)%X(i1,i2)  = hw*g2*pre/(H%Aneut + H%Aprot) 
        end if

        if (pp_calc) then 
           pp%mat(q)%gam(qx)%X(i1,i2)  = hw*g3*pre/(H%Aneut + H%Aprot) 
        end if

     end if 
     ! I shouldn't have to worry about hermiticity here, input is assumed to be hermitian
     
  end do   
  close(39)
      
end subroutine
!==================================================================  
!==================================================================
subroutine read_binary(H,jbas,htype,hw,rr,pp) 
  ! read interaction from binary file produced by VRenormalize 
  implicit none
  
  type(sq_op) :: H
  type(sq_op),optional :: rr,pp
  type(spd) :: jbas
  character(200) :: spfile,intfile,prefix
  integer :: ist,J,Tz,Par,a,b,c,d,q,qx,N,j_min,x
  real(8) :: V,Vcm,g1,g2,g3,pre,hw
  integer :: C1,C2,int1,int2,i1,i2,htype,COM
  integer :: ntot,npp,npn,nnn,count
  logical :: rr_calc,pp_calc
  common /files/ spfile,intfile,prefix
  
  open(unit=39,file = trim(TBME_DIR)//trim(adjustl(intfile)),form='unformatted') 
  
  COM = 0
  if (htype == 1) COM = 1 ! remove center of mass hamiltonian? 
  
  N = jbas%total_orbits
  
  ! check if we are concerned with other operators
  rr_calc = .false.
  pp_calc = .false. 
  if (present(rr)) rr_calc= .true. 
  if (present(pp)) pp_calc= .true. 

  read(39) ntot,npp,npn,nnn 

  do count = 1, ntot

     read(39) Tz,Par,J,a,b,c,d,V
     !read(39) Tz,Par,J,a,b,c,d,V,g1,g2,g3

     if (rr_calc) then 
        g3 = r1_r2( a, b, c, d, J ,jbas ) ! morten and Koshiroh are inconsistent with their definitions
     end if    
     g3 = p1_p2( a, b, c, d, J ,jbas ) ! morten and Koshiroh are inconsistent with their definitions
     
     ! g1 is COM expectation value, NOT CALCULATED WITH SCOTT'S CODE
     ! g2 is r1*r2 ME 
     ! g3 is p1*p2 ME 
     g1 = (g3 + H%com_hw**2 /hw**2 *g2)*H%lawson_beta ! lawson term    

     V = V + (g1 - g3*COM) *hw/(H%Aneut + H%Aprot) ! center of mass correction
          
     q = block_index(J,Tz,Par) 
     
     C1 = jbas%con(a)+jbas%con(b) + 1 !ph nature
     C2 = jbas%con(c)+jbas%con(d) + 1
    
     qx = C1*C2
     qx = qx + adjust_index(qx)   !Vpppp nature  

     ! get the indeces in the correct order
     pre = 1
     if ( a > b )  then 
        
        x = bosonic_tp_index(b,a,N) 
        j_min = H%xmap(x)%Z(1)  
        i1 = H%xmap(x)%Z( (J-j_min)/2 + 2) 
        pre = (-1)**( 1 + (jbas%jj(a) + jbas%jj(b) -J)/2 ) 
     else
       ! if (a == b) pre = pre * sqrt( 2.d0 )
       
        x = bosonic_tp_index(a,b,N) 
        j_min = H%xmap(x)%Z(1)  
        i1 = H%xmap(x)%Z( (J-j_min)/2 + 2) 
     end if
  
     if (c > d)  then     
        
        x = bosonic_tp_index(d,c,N) 
        j_min = H%xmap(x)%Z(1)  
        i2 = H%xmap(x)%Z( (J-j_min)/2 + 2) 
        
        pre = pre * (-1)**( 1 + (jbas%jj(c) + jbas%jj(d) -J)/2 ) 
     else 
       ! if (c == d) pre = pre * sqrt( 2.d0 )
      
        x = bosonic_tp_index(c,d,N) 
        j_min = H%xmap(x)%Z(1)  
        i2 = H%xmap(x)%Z( (J-j_min)/2 + 2) 
     end if
     ! kets/bras are pre-scaled by sqrt(2) if they 
     ! have two particles in the same sp-shell
        
     ! get the units right. I hope 
   
     
     if ((qx == 1) .or. (qx == 5) .or. (qx == 4)) then 
        H%mat(q)%gam(qx)%X(i2,i1)  = V *pre
        H%mat(q)%gam(qx)%X(i1,i2)  = V *pre
        
        if (rr_calc) then 
           rr%mat(q)%gam(qx)%X(i2,i1)  = hw*g2*pre/(H%Aneut + H%Aprot)
           rr%mat(q)%gam(qx)%X(i1,i2)  = hw*g2*pre/(H%Aneut + H%Aprot)
        end if 

        if (pp_calc) then 
           pp%mat(q)%gam(qx)%X(i2,i1)  = hw*g3*pre/(H%Aneut + H%Aprot)
           pp%mat(q)%gam(qx)%X(i1,i2)  = hw*g3*pre/(H%Aneut + H%Aprot)
        end if 

        
     else if (C1>C2) then
        H%mat(q)%gam(qx)%X(i2,i1)  = V *pre
        
        if (rr_calc) then 
           rr%mat(q)%gam(qx)%X(i2,i1)  = hw*g2*pre/(H%Aneut + H%Aprot) 
        end if
        
        if (pp_calc) then 
           pp%mat(q)%gam(qx)%X(i2,i1)  = hw*g3*pre/(H%Aneut + H%Aprot) 
        end if

     else
        H%mat(q)%gam(qx)%X(i1,i2) = V * pre
        
        if (rr_calc) then 
           rr%mat(q)%gam(qx)%X(i1,i2)  = hw*g2*pre/(H%Aneut + H%Aprot) 
        end if

        if (pp_calc) then 
           pp%mat(q)%gam(qx)%X(i1,i2)  = hw*g3*pre/(H%Aneut + H%Aprot) 
        end if

     end if 
     ! I shouldn't have to worry about hermiticity here, input is assumed to be hermitian
     
  end do   
  close(39)
      
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
integer function halfint_index(j1,j2,RANK) 
  ! input 2*J1 2*J2,rank,Tz,and Parity to get block index
  integer :: J1,J2,T,P,RANK,MORE,n
  
  n = (j1+1)/2
  
  IF ( j1 < RANK ) THEN
     n = (j1-1)/2
     halfint_index = n*(n+1) + (j2 - abs(j1-rank))/2 + 1           
  ELSE
     n = rank/2
     halfint_index = n*(n+1) + (j2 - abs(j1-rank))/2 + 1  + &
          (j1-rank-1)/2*(rank+1.d0) 
  END IF

end function 
!=================================================================     
!=================================================================
integer function tensor_block_index(J1,J2,RANK,T,P) 
  ! input 2*J1 2*J2,rank,Tz,and Parity to get block index
  integer :: J1,J2,T,P,RANK,MORE,JX
  
  IF ( J1 < RANK/2 ) THEN
     tensor_block_index = 6*(J1*J1/4+(J2-rank+J1)/2) &
          + 2*(T+1) + P + 1      
  ELSE
     tensor_block_index = 6*(((rank-2)/4+1)**2 +&
          (J1/2-((rank-2)/2)/2-1)*(rank+2)/2+(J2-J1)/2 ) &
          + 2*(T+1) + P + 1 
  END IF
     
end function 
!=================================================================     
!=================================================================
integer function CCtensor_block_index(J1,J2,RANK,T,P) 
  ! input 2*J1 2*J2,rank,Tz,and Parity to get block index
  integer :: J1,J2,T,P,RANK
    
  
  IF ( J1 < RANK/2 ) THEN
     CCtensor_block_index = 4*(J1*J1/4+(J2-rank+J1)/2) &
          + 2*T + P + 1      
  ELSE
     CCtensor_block_index = 4*(((rank-2)/4+1)**2 + &
          (J1/2-((rank-2)/2)/2-1)*(rank+2)/2+(J2-J1)/2 ) &
          + 2*T + P + 1 
  END IF
       
  !CCtensor_block_index = 4*( (J1/2 - 1) * ( rank/2+1) +(J2 - J1)/2 + 1) &
   !                + 2*T + P + 1 

end function 
!=================================================================     
!=================================================================     
integer function sp_block_index(j,l,t,jbas) 
  implicit none 
  
  type(spd) :: jbas
  integer :: j,l,t
  
  sp_block_index = ((t+1) *  (jbas%Jtotal_max + 1) * (jbas%lmax + 1)) / 4 + &
       l * (jbas%Jtotal_max + 1) / 2 + (j+1) / 2
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
           f_elem = op%fph(b-jbas%holesb4(b),a-jbas%partsb4(a)) * &
                op%herm
        else 
           f_elem = op%fph(a-jbas%holesb4(a),b-jbas%partsb4(b)) 
        end if
     case(2) 
        ! hh 
        f_elem = op%fhh(a-jbas%partsb4(a),b-jbas%partsb4(b)) 
  end select

end function
!=================================================================     
!=================================================================
real(8) function ph_elem(a,b,op,jbas) 
  implicit none 
  
  integer :: a,b,x1,x2,c1,c2
  type(spd) :: jbas
  type(sq_op) :: op 
  
  ! are they holes or particles
  c1 = jbas%con(a)
  c2 = jbas%con(b) 
  
  if ( (c1 == 1) .or. (c2 == 0) ) then 
     ph_elem = 0.d0 
     return
  end if 
  
  select case(c1+c2) 
     case(0) 
        ! pp 
        ph_elem = op%fpp(a-jbas%holesb4(a),b-jbas%holesb4(b)) 
  
     case(1) 
        ! ph 
        if (c1 > c2) then 
           ph_elem = op%fph(b-jbas%holesb4(b),a-jbas%partsb4(a)) * &
                op%herm
        else 
           ph_elem = op%fph(a-jbas%holesb4(a),b-jbas%partsb4(b)) 
        end if
     case(2) 
        ! hh 
        ph_elem = op%fhh(a-jbas%partsb4(a),b-jbas%partsb4(b)) 
  end select

end function
!=================================================================     
!=================================================================
real(8) function f_tensor_elem(a,b,op,jbas) 
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
        f_tensor_elem = op%fpp(a-jbas%holesb4(a),b-jbas%holesb4(b)) 
  
     case(1) 
        ! ph 
        if (c1 > c2) then 
           f_tensor_elem = op%fph(b-jbas%holesb4(b),a-jbas%partsb4(a)) * &
              op%herm * (-1)**( (jbas%jj(a) - jbas%jj(b))/2 ) 
        else 
           f_tensor_elem = op%fph(a-jbas%holesb4(a),b-jbas%partsb4(b)) 
        end if
     case(2) 
        ! hh 
        f_tensor_elem = op%fhh(a-jbas%partsb4(a),b-jbas%partsb4(b)) 
  end select

end function
!=================================================================     
!=================================================================
real(8) function ph_tensor_elem(a,b,op,jbas) 
  implicit none 
  
  integer :: a,b,x1,x2,c1,c2
  type(spd) :: jbas
  type(sq_op) :: op 
  
  ! are they holes or particles
  c1 = jbas%con(a)
  c2 = jbas%con(b) 
  
  if ( (c1 == 1) .or. (c2 == 0) ) then 
     ph_tensor_elem = 0.d0 
     return
  end if 
  
  ph_tensor_elem = op%fph(a-jbas%holesb4(a),b-jbas%partsb4(b)) 
           
end function
!==============================================================
!==============================================================
real(8) function v_elem(a,b,c,d,J,op,jbas) 
  ! grabs the matrix element you are looking for
  implicit none
  
  integer :: a,b,c,d,J,T,P,q,qx,c1,c2,N,pphh_int
  integer :: int1,int2,i1,i2,j_min,x
  integer :: ja,jb,jc,jd,la,lb,lc,ld,ta,tb,tc,td
  integer :: c1_c,c2_c,q_c,qx_c,i1_c,i2_c  
  logical :: fail_c
  type(sq_op) :: op 
  type(spd) :: jbas
  real(8) :: pre,pre_c
  common /TBMEinfo/ c1_c,c2_c,q_c,qx_c,i1_c,i2_c,pre_c,fail_c  
  
 ! make sure the matrix element exists first
  if (op%pphh_ph) then 
     ! is this a generator with a bunch of zeros? 
     pphh_int = jbas%con(a) + jbas%con(b)
     pphh_int = pphh_int - jbas%con(c) - jbas%con(d) 
     if (abs(pphh_int) .ne. 2) then 
        v_elem = 0.d0 
        return
     end if
  end if
 
 fail_c = .true. 
 ja = jbas%jj(a)
 jb = jbas%jj(b)
 jc = jbas%jj(c)
 jd = jbas%jj(d)
  
 if ( .not. ((triangle(ja,jb,J)) .and. (triangle (jc,jd,J))) ) then 
    v_elem = 0.d0
    return
 end if 
     
  la = jbas%ll(a)
  lb = jbas%ll(b)
  lc = jbas%ll(c)
  ld = jbas%ll(d)
     
  P = mod(la + lb,2) 
     
  if ( mod(lc + ld,2) .ne. P ) then
    v_elem = 0.d0 
    return
  end if 
        
  ta = jbas%itzp(a)
  tb = jbas%itzp(b)
  tc = jbas%itzp(c)
  td = jbas%itzp(d)
     
  T = (ta + tb)/2
     
  
  if ((tc+td) .ne. 2*T) then     
    v_elem = 0.d0
    return
  end if 

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
    
  pre = 1 
  
  N = op%Nsp
  ! get the indeces in the correct order
  if ( a > b )  then 
     x = bosonic_tp_index(b,a,N) 
     j_min = op%xmap(x)%Z(1)  
     i1 = op%xmap(x)%Z( (J-j_min)/2 + 2) 
     pre = (-1)**( 1 + (jbas%jj(a) + jbas%jj(b) -J)/2 ) 
  else
     if (a == b) pre = pre * sqrt( 2.d0 )
     x = bosonic_tp_index(a,b,N)
     j_min = op%xmap(x)%Z(1)  
     i1 = op%xmap(x)%Z( (J-j_min)/2 + 2) 
  end if 
  
  if (c > d)  then     
     x = bosonic_tp_index(d,c,N) 
     j_min = op%xmap(x)%Z(1)  
     i2 = op%xmap(x)%Z( (J-j_min)/2 + 2) 
     pre = pre * (-1)**( 1 + (jbas%jj(c) + jbas%jj(d) -J)/2 ) 
  else 
     if (c == d) pre = pre * sqrt( 2.d0 )
     x = bosonic_tp_index(c,d,N) 
     j_min = op%xmap(x)%Z(1)  
     i2 = op%xmap(x)%Z( (J-j_min)/2 + 2)  
  end if 
 
  ! grab the matrix element
   If (C1>C2) then 
      v_elem = op%mat(q)%gam(qx)%X(i2,i1) * op%herm * pre 
   else
      v_elem = op%mat(q)%gam(qx)%X(i1,i2) * pre
   end if 
  
   ! stored info if we are looking for this same ME next time. 
   c1_c=C1;c2_c=C2;q_c=q;qx_c=qx
   i1_c=i1;i2_c=i2;pre_c=pre;fail_c=.false.   
end function
!==============================================================
!==============================================================
real(8) function pphh_elem(a,b,c,d,J,op,jbas) 
  ! grabs the matrix element you are looking for
  implicit none
  
  integer :: a,b,c,d,J,T,P,q,qx,c1,c2,N
  integer :: int1,int2,i1,i2,j_min,x
  integer :: ja,jb,jc,jd,la,lb,lc,ld,ta,tb,tc,td
  integer :: c1_c,c2_c,q_c,qx_c,i1_c,i2_c  
  logical :: fail_c,pphh
  type(sq_op) :: op 
  type(spd) :: jbas
  real(8) :: pre,pre_c
  
  !make sure the matrix element exists first
  
  pphh = .false. 
  if ( jbas%con(a) == 1) pphh = .true.  
  if ( jbas%con(b) == 1) pphh = .true. 
  if ( jbas%con(c) == 0) pphh = .true. 
  if ( jbas%con(d) == 0) pphh = .true. 
     
  if (pphh) then 
     pphh_elem = 0.d0 
     return
  end if 
 
  fail_c = .true. 
  ja = jbas%jj(a)
  jb = jbas%jj(b)
  jc = jbas%jj(c)
  jd = jbas%jj(d)
  
 if ( .not. ((triangle(ja,jb,J)) .and. (triangle (jc,jd,J))) ) then 
    pphh_elem = 0.d0
    return
 end if 
     
  la = jbas%ll(a)
  lb = jbas%ll(b)
  lc = jbas%ll(c)
  ld = jbas%ll(d)
     
  P = mod(la + lb,2) 
     
  if ( mod(lc + ld,2) .ne. P ) then
    pphh_elem = 0.d0 
    return
  end if 
        
  ta = jbas%itzp(a)
  tb = jbas%itzp(b)
  tc = jbas%itzp(c)
  td = jbas%itzp(d)
     
  T = (ta + tb)/2
     
  
  if ((tc+td) .ne. 2*T) then     
    pphh_elem = 0.d0
    return
  end if 

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
    
  pre = 1 
  
  N = op%Nsp
  ! get the indeces in the correct order
  if ( a > b )  then 
     x = bosonic_tp_index(b,a,N) 
     j_min = op%xmap(x)%Z(1)  
     i1 = op%xmap(x)%Z( (J-j_min)/2 + 2) 
     pre = (-1)**( 1 + (jbas%jj(a) + jbas%jj(b) -J)/2 ) 
  else
     if (a == b) pre = pre * sqrt( 2.d0 )
     x = bosonic_tp_index(a,b,N)
     j_min = op%xmap(x)%Z(1)  
     i1 = op%xmap(x)%Z( (J-j_min)/2 + 2) 
  end if 
  
  if (c > d)  then     
     x = bosonic_tp_index(d,c,N) 
     j_min = op%xmap(x)%Z(1)  
     i2 = op%xmap(x)%Z( (J-j_min)/2 + 2) 
     pre = pre * (-1)**( 1 + (jbas%jj(c) + jbas%jj(d) -J)/2 ) 
  else 
     if (c == d) pre = pre * sqrt( 2.d0 )
     x = bosonic_tp_index(c,d,N) 
     j_min = op%xmap(x)%Z(1)  
     i2 = op%xmap(x)%Z( (J-j_min)/2 + 2)  
  end if 
 
  ! grab the matrix element
   If (C1>C2) then 
      pphh_elem = op%mat(q)%gam(qx)%X(i2,i1) * op%herm * pre 
   else
      pphh_elem = op%mat(q)%gam(qx)%X(i1,i2) * pre
   end if 
   
end function
!==============================================================
!==============================================================
real(8) function tensor_elem(ax,bx,cx,dx,J1x,J2x,op,jbas) 
  ! grabs the matrix element you are looking for
  ! right now it's not as versatile as v_elem because 
  ! it doesn't know ahead of time what the symmetries are
  implicit none
  
  integer :: a,b,c,d,J1,J2,rank,T,P,q,qx,c1,c2,N,J1x,J2x
  integer :: int1,int2,i1,i2,j_min,x,k1,k2,ax,bx,cx,dx
  integer :: ja,jb,jc,jd,la,lb,lc,ld,ta,tb,tc,td
  integer :: c1_c,c2_c,q_c,qx_c,i1_c,i2_c  ,phase
  logical :: fail_c,switch
  type(sq_op) :: op 
  type(spd) :: jbas
  real(8) :: pre,pre_c
  
  J1 = min(J1x,J2x) 
  
  if (J1 == J1x) then 
     switch = .false. 
     J2 = J2x 
     phase = 1
     a = ax
     b = bx
     c = cx
     d = dx
  else
     J2 = J1x 
     switch = .true. 
     a = cx
     b = dx
     c = ax
     d = bx 
  end if
  
  !make sure the matrix element exists first
 
  rank = op%rank

  if ( .not. (triangle ( J1,J2,rank ))) then 
     tensor_elem = 0.d0 
     return
  end if 
 
  
  fail_c = .true. 
  ja = jbas%jj(a)
  jb = jbas%jj(b)
  jc = jbas%jj(c)
  jd = jbas%jj(d)
  
  if ( .not. ((triangle(ja,jb,J1)) .and. (triangle (jc,jd,J2))) ) then 
     tensor_elem = 0.d0
     return
  end if
     
  la = jbas%ll(a)
  lb = jbas%ll(b)
  lc = jbas%ll(c)
  ld = jbas%ll(d)

  P = mod(la + lb,2) 
     
  if ( mod(lc + ld,2) .ne. abs(P - ((-1)**(op%dpar/2+1)+1)/2) ) then
    tensor_elem = 0.d0 
    return
  end if 
        
  ta = jbas%itzp(a)
  tb = jbas%itzp(b)
  tc = jbas%itzp(c)
  td = jbas%itzp(d)
     
  T = (ta + tb)/2
     
  if ((tc+td) .ne. 2*T) then     
    tensor_elem = 0.d0
    return
  end if 

  q = tensor_block_index(J1,J2,rank,T,P) 

  if (switch) then 
     phase =  OP%tblck(q)%lam(1)*op%herm 
  end if 
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
    
  pre = 1 
  
  N = op%Nsp
  ! get the indeces in the correct order
  if ( a > b )  then 
     x = bosonic_tp_index(b,a,N) 
     j_min = op%xmap(x)%Z(1)  
     i1 = op%xmap(x)%Z( (J1-j_min)/2 + 2) 
     pre = (-1)**( 1 + (jbas%jj(a) + jbas%jj(b) -J1)/2 ) 
  else
     if (a == b) pre = pre * sqrt( 2.d0 )
     x = bosonic_tp_index(a,b,N)
     j_min = op%xmap(x)%Z(1)  
     i1 = op%xmap(x)%Z( (J1-j_min)/2 + 2) 
  end if 
  
  if (c > d)  then     
     x = bosonic_tp_index(d,c,N) 
     j_min = op%xmap(x)%Z(1)  
     i2 = op%xmap(x)%Z( (J2-j_min)/2 + 2) 
     pre = pre * (-1)**( 1 + (jbas%jj(c) + jbas%jj(d) -J2)/2 ) 
  else 
     if (c == d) pre = pre * sqrt( 2.d0 )
     x = bosonic_tp_index(c,d,N) 
     j_min = op%xmap(x)%Z(1)  
     i2 = op%xmap(x)%Z( (J2-j_min)/2 + 2)  
  end if 
 
  ! grab the matrix element

  ! right now i1 and i2 still refer to where the pair is located
  ! in the rank zero qn storage

  
  
   If (C1>C2) qx = qx + tensor_adjust(qx)       
   
   tensor_elem = op%tblck(q)%tgam(qx)%X(i1,i2) * pre *phase
    
end function
!==============================================================
!==============================================================
real(8) function pphh_tensor_elem(ax,bx,cx,dx,J1x,J2x,op,jbas) 
  ! grabs the matrix element you are looking for
  ! right now it's not as versatile as v_elem because 
  ! it doesn't know ahead of time what the symmetries are
  implicit none
  
  integer :: a,b,c,d,J1,J2,rank,T,P,q,qx,c1,c2,N,J1x,J2x
  integer :: int1,int2,i1,i2,j_min,x,k1,k2,ax,bx,cx,dx
  integer :: ja,jb,jc,jd,la,lb,lc,ld,ta,tb,tc,td
  integer :: c1_c,c2_c,q_c,qx_c,i1_c,i2_c  ,phase
  logical :: fail_c,switch,pphh
  type(sq_op) :: op 
  type(spd) :: jbas
  real(8) :: pre,pre_c
  
  pphh = .false. 
  if ( jbas%con(ax) == 1) pphh = .true.  
  if ( jbas%con(bx) == 1) pphh = .true. 
  if ( jbas%con(cx) == 0) pphh = .true. 
  if ( jbas%con(dx) == 0) pphh = .true. 
     
  if (pphh) then 
     pphh_tensor_elem = 0.d0 
     return
  end if 
  
  J1 = min(J1x,J2x) 
  
  if (J1 == J1x) then 
     switch = .false. 
     J2 = J2x 
     phase = 1
     a = ax
     b = bx
     c = cx
     d = dx
  else
     J2 = J1x 
     switch = .true. 
     a = cx
     b = dx
     c = ax
     d = bx 
  end if
  
  !make sure the matrix element exists first
 
  rank = op%rank

  if ( .not. (triangle ( J1,J2,rank ))) then 
     pphh_tensor_elem = 0.d0 
     return
  end if 
 
  
  fail_c = .true. 
  ja = jbas%jj(a)
  jb = jbas%jj(b)
  jc = jbas%jj(c)
  jd = jbas%jj(d)
  
  if ( .not. ((triangle(ja,jb,J1)) .and. (triangle (jc,jd,J2))) ) then 
     pphh_tensor_elem = 0.d0
     return
  end if
     
  la = jbas%ll(a)
  lb = jbas%ll(b)
  lc = jbas%ll(c)
  ld = jbas%ll(d)

  P = mod(la + lb,2) 
     
  if ( mod(lc + ld,2) .ne. abs(P - ((-1)**(op%dpar/2+1)+1)/2) ) then
    pphh_tensor_elem = 0.d0 
    return
  end if 
        
  ta = jbas%itzp(a)
  tb = jbas%itzp(b)
  tc = jbas%itzp(c)
  td = jbas%itzp(d)
     
  T = (ta + tb)/2
     
  if ((tc+td) .ne. 2*T) then     
    pphh_tensor_elem = 0.d0
    return
  end if 

  q = tensor_block_index(J1,J2,rank,T,P) 

  if (switch) then 
     phase =  OP%tblck(q)%lam(1)*op%herm 
  end if 
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
    
  pre = 1 
  
  N = op%Nsp
  ! get the indeces in the correct order
  if ( a > b )  then 
     x = bosonic_tp_index(b,a,N) 
     j_min = op%xmap(x)%Z(1)  
     i1 = op%xmap(x)%Z( (J1-j_min)/2 + 2) 
     pre = (-1)**( 1 + (jbas%jj(a) + jbas%jj(b) -J1)/2 ) 
  else
     if (a == b) pre = pre * sqrt( 2.d0 )
     x = bosonic_tp_index(a,b,N)
     j_min = op%xmap(x)%Z(1)  
     i1 = op%xmap(x)%Z( (J1-j_min)/2 + 2) 
  end if 
  
  if (c > d)  then     
     x = bosonic_tp_index(d,c,N) 
     j_min = op%xmap(x)%Z(1)  
     i2 = op%xmap(x)%Z( (J2-j_min)/2 + 2) 
     pre = pre * (-1)**( 1 + (jbas%jj(c) + jbas%jj(d) -J2)/2 ) 
  else 
     if (c == d) pre = pre * sqrt( 2.d0 )
     x = bosonic_tp_index(c,d,N) 
     j_min = op%xmap(x)%Z(1)  
     i2 = op%xmap(x)%Z( (J2-j_min)/2 + 2)  
  end if 
 
  ! grab the matrix element

  ! right now i1 and i2 still refer to where the pair is located
  ! in the rank zero qn storage

  
  
   If (C1>C2) qx = qx + tensor_adjust(qx)       
   
   pphh_tensor_elem = op%tblck(q)%tgam(qx)%X(i1,i2) * pre *phase
    
end function
!==============================================================
!==============================================================
real(8) function T_twobody(a,b,c,d,J,T,op,jbas) 
  ! grabs the matrix element you are looking for
  implicit none
  
  integer :: a,b,c,d,J,T,P,q,qx,c1,c2,N,i,JT,ax,bx,cx,dx
  integer :: int1,int2,i1,i2,j_min,x,ji,x1,x2,x3,x4
  integer :: ja,jb,jc,jd,la,lb,lc,ld,na,nb,nc,nd
  integer :: c1_c,c2_c,q_c,qx_c,i1_c,i2_c  
  logical :: fail_c
  type(sq_op) :: op 
  type(spd) :: jbas
  real(8) :: pre,pre_c,sm1,sm2,sm3,sm4
 
  !make sure the matrix element exists first
  
 ja = jbas%jj(a)
 jb = jbas%jj(b)
 jc = jbas%jj(c)
 jd = jbas%jj(d)
  
 la = jbas%ll(a)
 lb = jbas%ll(b)
 lc = jbas%ll(c)
 ld = jbas%ll(d)
 
 na = jbas%nn(a)
 nb = jbas%nn(b)
 nc = jbas%nn(c)
 nd = jbas%nn(d) 
 
 if ( .not. ((triangle(ja,jb,J)) .and. (triangle (jc,jd,J))) ) then 
    T_twobody = 0.d0
    return
 end if 
 
 pre = 1.d0 
 if (a == b) pre = pre * sqrt( 2.d0 )
 if (c == d) pre = pre * sqrt( 2.d0 )
 
 sm1=0.d0
 sm2=0.d0
 sm3=0.d0
 sm4=0.d0
 
 do i1 = 1,op%belowEF
    i = jbas%holes(i1) 
    ji = jbas%jj(i) 
    
    do JT = abs(ja-ji), ja +ji,2 
       sm1 = sm1 + v_elem(a,i,c,i,JT,op,jbas) * (JT+1.d0)  
    end do 
   
    do JT = abs(jb-ji), jb +ji,2 
       sm2 = sm2 + v_elem(b,i,d,i,JT,op,jbas) * (JT+1.d0)  
    end do 

    do JT = abs(jb-ji), jb +ji,2 
       sm3 = sm3 + v_elem(b,i,c,i,JT,op,jbas) * (JT+1.d0)  
    end do 

    do JT = abs(ja-ji), ja +ji,2 
       sm4 = sm4 + v_elem(a,i,d,i,JT,op,jbas) * (JT+1.d0)  
    end do
   
    sm1 = sm1/(ja +1.d0)
    sm4 = sm4/(ja +1.d0)
    sm2 = sm2/(jb +1.d0)
    sm3 = sm3/(jb +1.d0)
 end do 

 x1 = 1.d0
 x2 = 1.d0 
 x3 = 1.d0 
 x4 = 1.d0 
 if ((ja == jc).and.(la==lc).and.(na==nc)) x1 = 0.d0
 if ((ja == jd).and.(la==ld).and.(na==nd)) x3 = 0.d0
 if ((jb == jd).and.(lb==ld).and.(nb==nd)) x2 = 0.d0
 if ((jb == jc).and.(lb==lc).and.(nb==nc)) x4 = 0.d0
 
 T_twobody = pre*(kron_del(b,d)*x1* (f_elem(a,c,op,jbas)-sm1) + &
      kron_del(a,c)*x2 * (f_elem(b,d,op,jbas)-sm2) - &
      (-1)** ( (ja + jb - J)/2 ) * x3* kron_del(a,d) * (f_elem(b,c,op,jbas)-sm3) - &
      (-1)** ( (jc + jd - J)/2 ) * x4* kron_del(b,c) * (f_elem(a,d,op,jbas)-sm4))
 
end function 
!=====================================================
!=====================================================
real(8) function T_elem(a,b,op,jbas) 
  ! grabs the matrix element you are looking for
  implicit none
  
  integer :: a,b,c,d,J,T,P,q,qx,c1,c2,N,i,JT
  integer :: int1,int2,i1,i2,j_min,x,ji
  integer :: ja,jb,jc,jd,la,lb,lc,ld,ta,tb,tc,td
  integer :: c1_c,c2_c,q_c,qx_c,i1_c,i2_c  
  logical :: fail_c
  type(sq_op) :: op 
  type(spd) :: jbas
  real(8) :: pre,pre_c,sm1,sm2,sm3,sm4
  
 ja = jbas%jj(a)
 jb = jbas%jj(b)
 
 ta = jbas%itzp(a)
 tb = jbas%itzp(b)
 
 la = jbas%ll(a)
 lb = jbas%ll(b)
 
 if (( ja .ne. jb ) .or. (la .ne. lb) .or. (ta .ne. tb)) then 
    T_elem = 0.d0
    return
 end if 
 
 sm1 = 0.d0  
  do i1 = 1,op%belowEF
    i = jbas%holes(i1) 
    ji = jbas%jj(i) 
    
    do JT = abs(ja-ji), ja +ji,2 
       sm1 = sm1 + v_elem(a,i,b,i,JT,op,jbas) * (JT+1.d0)  
    end do 
  end do 
  
  T_elem = f_elem(a,b,op,jbas) - sm1/(ja+1.d0) 

end function
!=====================================================
!=====================================================
real(8) function v_same(op) 
  !produces the matrix element of op which is in the exact position 
  !as the one obtained from the previous call of v_elem 
  
  integer :: c1_c,c2_c,q_c,qx_c,i1_c,i2_c
  logical :: fail_c
  type(sq_op) :: op 
  real(8) :: pre_c
  common /TBMEinfo/ c1_c,c2_c,q_c,qx_c,i1_c,i2_c,pre_c,fail_c
  
  if (fail_c) then 
     v_same = 0.d0 
     return 
  end if 
  
   ! grab the matrix element
  If (C1_c>C2_c) then 
      v_same = op%mat(q_c)%gam(qx_c)%X(i2_c,i1_c) * op%herm * pre_c 
  else
      v_same = op%mat(q_c)%gam(qx_c)%X(i1_c,i2_c) * pre_c
  end if 
end function
!=====================================================
!=====================================================
subroutine calculate_h0_harm_osc(hw,jbas,H,Htype) 
  ! fills out the one body piece of the hamiltonian
  implicit none 
  
  integer,intent(in) :: Htype
  real(8),intent(in) :: hw
  integer :: i,j,mass,c1,c2,cx
  integer :: ni,li,ji,nj,lj,jj,tzi,tzj,AX
  real(8) :: kij,T,beta,cmhw
  type(sq_op) :: H 
  type(spd) :: jbas
 
  !Htype =  |[1]: T - Tcm + V |[2]: T + Uho + V |[3]: T+V | 

  AX = H%belowEF
  mass = H%Aprot + H%Aneut
  
  H%fpp = 0.d0
  H%fhh = 0.d0
  H%fph = 0.d0
  beta = H%lawson_beta
  cmhw = H%com_hw

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
              if (abs(beta) < 1e-5) then
                 T =  kij*(1.d0-1.d0/mass) 
              else
                 T =  kij*(1.d0+(beta-1.d0)/mass) & !lawson modified
                      + beta*cmhw*cmhw/(hw*hw) * &
                      kij*(kron_del(ni,nj) - kron_del(ni,nj+1) - &
                   kron_del(ni,nj-1) )/mass
              end if 
           case(2) 
              T = 2*kij*kron_del(ni,nj)
           case(3)
              T = kij
           case(4)
              T = kij/mass
              ! dividing by mass for Hcm
           case(5)
              ! Hcm, lacking (omegaT/hw)**2 factor
              T = kij*(kron_del(ni,nj) - kron_del(ni,nj+1) - &
                   kron_del(ni,nj-1) )/mass
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
subroutine store_6j(jbas,method) 
  ! run this once. use sixj in code to call elements 
  ! kind of complicated, but if fills a public array (store6j) 
  ! which theoretically doesn't need to be touched when 
  ! making changes to the code ( this is all background stuff ) 
  ! so DFWT subroutine because it will be a huge mess
  implicit none 
  
  type(spd) :: jbas
  integer :: j1,j2,j3,j4,j5,j6,num_half,num_whole,r1,r2
  integer :: nbos,nferm,X12,X45,j3min,j3max,j6min,j6max,method
  real(8) :: d6ji
  
  call dfact0() ! prime the anglib 6j calculator 
  
  if (method == 5) then  
     ! this is necessary for three-body calculations
     num_half = (3*jbas%Jtotal_max + 1)/2
     ! however, it clearly increases the memory requirements
  else
     num_half = (jbas%Jtotal_max+6 + 1)/2
  end if
  
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
          
            j3min = max(  abs( 2*j1-1 -  (2*j2-1) ) , abs( 2*j4-1 -  (2*j5-1) )   )               
            j3max = min(  2*j1-1 +  2*j2-1 , 2*j4-1 +  2*j5-1 ) 

            j6min = max(  abs( 2*j1-1 -  (2*j5-1) ) , abs( 2*j4-1 -  (2*j2-1) )   )               
            j6max = min(  2*j1-1 +  2*j5-1 , 2*j4-1 +  2*j2-1 ) 
            
              allocate(store6j%tp_mat(X12,X45)%X((j3max-j3min)/2+1,(j6max-j6min)/2+1 ) )   
 
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
!======================================================
!======================================================
subroutine store_6j_3halfint(jbas,rank) 
  ! THIS CODE ASSUMES THAT SIX-J SYMBOLS WITH 3 half-INTEGER j take the form:
  ! { J1 , J2  , X }
  ! { a  ,  b  , c }    WHERE X IS THE RANK OF SOME TENSOR    
  !
  ! run this once. use sixj in code to call elements 
  ! kind of complicated, but if fills a public array (half6j) 
  ! which theoretically doesn't need to be touched when 
  ! making changes to the code ( this is all background stuff ) 
  ! so DFWT subroutine because it will be a huge mess
  implicit none 
  
  type(spd) :: jbas
  integer :: j1,j2,j3,j4,j5,j6,halfmax,num_whole,r1,r2,rank,JTM
  integer :: nbos,nferm,X12,X45,j3min,j3max,j6min,j6max,method,TZ,PAR
  real(8) :: d6ji
  
  call dfact0() ! prime the anglib 6j calculator 
  
  JTM = jbas%jtotal_max*2
 
  halfmax = jbas%Jtotal_max + JTM
  !end if                        
  
  !these mean nothing physical
  ! they are just set to one so I can use tensor_block_index
  TZ = 1
  PAR = 1
  
  !half6j%nhalf = num_half
  nbos = tensor_block_index(JTM,JTM+RANK,RANK,TZ,PAR)/6 
  ! divide by six because TZ and PAR are irrelevant
  nferm = halfint_index(halfmax,halfmax+rank,rank) 
  
  half6j%nb = nbos
  half6j%nf = nferm 
  
  
  allocate(half6j%tp_mat(nbos,nferm)) 
  ! The first index refers to J1,J2 
  ! which are forcibly ordered
  
  ! the second index refers to j4,j5 which cannot be ordered
  ! so we need an extra set of states for the reverse ordering (nferm) 
  
  do J1 = 0,JTM,2 
     do J2 = max(abs(J1 - rank),J1), J1+rank , 2   

        X12 = tensor_block_index(J1,J2,RANK,TZ,PAR)/6 
        
        do j4 = 1, halfmax,2
           do j5 = abs(j4-rank),j4+rank,2
           
              ! if j5<j4 use fermionic index
              X45 = halfint_index(j4,j5,rank) 
              
              ! find the allowed values of j3
              j3min = max(  abs( J1 - j5 ) , abs( J2-j4 )   )               
              j3max = min(  J1+j5 , J2+j4 ) 

            
            ! allocate second part of storage array to hold them 
              allocate(half6j%tp_mat(X12,X45)%X((j3max-j3min)/2+1,1 ) ) 
              
              ! fill the array with 6j symbols
              r1 = 1 
             
              do j3 = j3min,j3max,2
                 r2 = 1
                    half6j%tp_mat(X12,X45)%X(r1,1) =  d6ji(J1,J2,rank,j4,j5,j3)   
                    ! d6ji is the anglib sixj calculator (takes 2*j as arguments) 
                 r1 = r1 + 1
              end do 
                 
           end do 
        end do
     end do 
  end do 
   
end subroutine    
!=========================================================
!=========================================================
subroutine normal_order(H,jbas) 
  ! puts H in normal ordering
  implicit none 
  
  type(sq_op) :: H,TEMP
  type(spd) :: jbas
  integer :: i,j,q,a,b,JT,T,P
  integer :: ax,bx,ix
  real(8) :: sm
  
  
  call duplicate_sq_op(H,TEMP) 
  call copy_sq_op(H,TEMP) 
  
  ! calculate vacuum expectation value
  
  sm = 0.d0 
  do i = 1,TEMP%Nsp
     sm = sm + f_elem(i,i,TEMP,jbas)*jbas%con(i)* &
          (jbas%jj(i) + 1)
  end do 
  
  sm = sm * 2.d0 
  
  do i = 1,TEMP%Nsp
     do j = 1,TEMP%Nsp
        
        if (jbas%con(i)*jbas%con(j) == 0) cycle
        
        do JT = abs(jbas%jj(i) - jbas%jj(j)) ,jbas%jj(i)+jbas%jj(j),2 
           sm = sm + v_elem(i,j,i,j,JT,TEMP,jbas)*(JT + 1) 
        end do
     end do 
  end do 
  
  
  H%E0 = sm *0.5d0 
  
  ! calculate one body part
  !fhh
  do a = 1,TEMP%belowEF
     do b = a,TEMP%belowEF 
        
        ax = jbas%holes(a) 
        bx = jbas%holes(b) 
        
        
        if (jbas%jj(ax) .ne. jbas%jj(bx) ) cycle
        if (jbas%ll(ax) .ne. jbas%ll(bx) ) cycle
        if (jbas%itzp(ax) .ne. jbas%itzp(bx) ) cycle
        
        sm = 0.d0 
        
        do i = 1, TEMP%belowEF
           ix = jbas%holes(i) 
           do JT = abs(jbas%jj(ix) - jbas%jj(ax)), &
                jbas%jj(ix) + jbas%jj(ax), 2
                  
              sm = sm + (JT+1.d0)/(jbas%jj(ax) + 1) *&
                   v_elem(ax,ix,bx,ix,JT,TEMP,jbas) 
              
           end do 
        end do 
           
        H%fhh(a,b) = f_elem(ax,bx,TEMP,jbas) + sm 
        H%fhh(b,a) = H%fhh(a,b) 
     end do 
  end do
 
  !fph          
  do a = 1,TEMP%Nsp-TEMP%belowEF
     do b = 1,TEMP%belowEF 
        
        ax = jbas%parts(a) 
        bx = jbas%holes(b) 
        
        
        if (jbas%jj(ax) .ne. jbas%jj(bx) ) cycle
        if (jbas%ll(ax) .ne. jbas%ll(bx) ) cycle
        if (jbas%itzp(ax) .ne. jbas%itzp(bx) ) cycle
        
        sm = 0.d0 
        
        do i = 1, TEMP%belowEF
           ix = jbas%holes(i) 
           do JT = abs(jbas%jj(ix) - jbas%jj(ax)), &
                jbas%jj(ix) + jbas%jj(ax), 2
              
              sm = sm + (JT+1.d0)/(jbas%jj(ax) + 1) *&
                   v_elem(ax,ix,bx,ix,JT,TEMP,jbas) 
              
           end do 
        end do 
           
        H%fph(a,b) = f_elem(ax,bx,TEMP,jbas) + sm 
        
     end do 
  end do 

  !fpp
  do a = 1,TEMP%Nsp-TEMP%belowEF
     do b = a,TEMP%Nsp - TEMP%belowEF 
        
        ax = jbas%parts(a) 
        bx = jbas%parts(b) 
       
        if (jbas%jj(ax) .ne. jbas%jj(bx) ) cycle
        if (jbas%ll(ax) .ne. jbas%ll(bx) ) cycle
        if (jbas%itzp(ax) .ne. jbas%itzp(bx) ) cycle
        
        sm = 0.d0 
        
        do i = 1, TEMP%belowEF
           ix = jbas%holes(i) 
           do JT = abs(jbas%jj(ix) - jbas%jj(ax)), &
                jbas%jj(ix) + jbas%jj(ax), 2
              
              sm = sm + (JT+1.d0)/(jbas%jj(ax) + 1) *&
                   v_elem(ax,ix,bx,ix,JT,TEMP,jbas) 
              
           end do 
        end do 
           
        H%fpp(a,b) = f_elem(ax,bx,TEMP,jbas) + sm 
        H%fpp(b,a) = H%fpp(a,b) 
     end do 
  end do 
  
  ! two body part is already normal ordered 
  ! UNLESS WE HAVE THREE BODY FORCES
  
end subroutine   
!=========================================================
!=========================================================
real(8) function sixj(j1,j2,j3,j4,j5,j6)
  ! twice the angular momentum
  ! j3 and j6 are INTEGERS, the rest are HALF-INTEGERS
  ! you should be able to re-write the six-j symbol to look like this.
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
!=========================================================
!=========================================================
real(8) function xxxsixj(j1,j2,RANK,j4,j5,j6)
  ! twice the angular momentum
  ! J1, J2, RANK  are INTEGERS, the rest are HALF-INTEGERS
  ! you should be able to re-write the six-j symbol to look like this.
  ! so you don't have to carry that obnoxious array around
  implicit none 
 
  integer :: j1,j2,j3,j4,j5,j6,RANK
  integer :: l1,l2,l3,l4,l5,l6
  integer :: x1,x2,j3min,j6min
    
  ! check triangle inequalities
  if ( (triangle(j1,j2,RANK)) .and. (triangle(j4,j5,RANK)) &
  .and. (triangle(j1,j5,j6)) .and. (triangle(j4,j2,j6)) ) then 
     
     l1 = j1; l2 = j2; l3 = RANK; l4 = j4; l5 = j5; l6 = j6 
          
     ! find lowest possible momentum
     j6min = max(  abs( j1 -  j5)  , abs( j4 -  j2 )   )/2               
    

     ! make sure j1 < = j2
     if  ( j1 > j2 ) then 
        l1 = j2 ; l2 = j1 
        l4 = j5 ; l5 = j4
     end if 
     ! allowed by 6j symmetry
    
     ! find indeces
     x1 = tensor_block_index(l1,l2,RANK,1,1)/6
     x2 = halfint_index(l4,l5,rank) 

     ! figure out indeces for j3,j6 based on lowest possible
     xxxsixj = half6j%tp_mat(x1,x2)%X(l6/2 - j6min + 1 ,1) 
     
     
  else 
     xxxsixj = 0.d0 
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
!==========================================
!==========================================  
subroutine add_sp_mat(A,ax,B,bx,C) 
  ! use all of the information about H to make 
  ! an empty block matrix of the same size
  implicit none 
  
  type(full_sp_block_mat) :: A,B,C 
  real(8),intent(in) :: ax,bx
  integer :: q
  
  do q = 1, A%blocks
     if (A%map(q) == 0) cycle
     C%blkM(q)%matrix = ax*A%blkM(q)%matrix + bx*B%blkM(q)%matrix
  end do 

end subroutine
!=====================================================
subroutine duplicate_sq_op(H,op) 
  ! make a copy of the shape of H onto op
  implicit none 
  
  type(sq_op) :: H,op
  integer :: q,i,j,holes,parts,nh,np,nb,N,nh2,np2,nb2
  
  op%herm = 1 ! default, change it in the calling program
             ! if you want anti-herm (-1) 
  op%Aprot = H%Aprot
  op%Aneut = H%Aneut
  op%Nsp = H%Nsp
  op%belowEF = H%belowEF
  op%nblocks = H%nblocks
  op%neq = H%neq
  op%hospace = H%hospace
  op%Jtarg = H%Jtarg
  op%Ptarg = H%Ptarg
  op%valcut = H%valcut 
  op%rank = H%rank 
  op%dpar = H%dpar
  op%pphh_ph = .false. 
  
  holes = op%belowEF ! number of shells below eF 
  parts = op%Nsp - holes 
  
  allocate(op%fhh(holes,holes)) 
  allocate(op%fpp(parts,parts)) 
  allocate(op%fph(parts,holes)) 
  
  ! everything is initialized to zero
  op%E0 = 0.d0 
  op%fhh = 0.d0
  op%fpp = 0.d0
  op%fph = 0.d0 
  
  if (allocated(H%mat)) then ! scalar operator. 
     
     allocate(op%mat(op%nblocks)) 
     allocate(op%direct_omp(size(H%direct_omp)))
     op%direct_omp = H%direct_omp

     N = op%nsp
     allocate(op%xmap(N*(N+1)/2))
     do q = 1, N*(N+1)/2 
        allocate(op%xmap(q)%Z(size(H%xmap(q)%Z)))
        op%xmap(q)%Z = H%xmap(q)%Z
     end do
  
     if ( allocated( H%exlabels ) ) then 
        allocate(op%exlabels(size(H%exlabels(:,1)),2)) 
     else 
        allocate(op%exlabels(1,2),H%exlabels(1,2)) 
        H%exlabels = 0
     end if
     op%exlabels = H%exlabels
  
     do q = 1, op%nblocks
     
        op%mat(q)%lam = H%mat(q)%lam

        nh = H%mat(q)%nhh
        np = H%mat(q)%npp
        nb = H%mat(q)%nph
         
        op%mat(q)%npp = np
        op%mat(q)%nph = nb
        op%mat(q)%nhh = nh
        op%mat(q)%ntot = H%mat(q)%ntot

        allocate(op%mat(q)%qn(1)%Y(np,2)) !qnpp
        allocate(op%mat(q)%qn(2)%Y(nb,2)) !qnph
        allocate(op%mat(q)%qn(3)%Y(nh,2)) !qnhh
     
        do i = 1,3
           op%mat(q)%qn(i)%Y = H%mat(q)%qn(i)%Y
        end do
    
        allocate(op%mat(q)%gam(1)%X(np,np)) !Vpppp
        allocate(op%mat(q)%gam(5)%X(nh,nh)) !Vhhhh
        allocate(op%mat(q)%gam(3)%X(np,nh)) !Vpphh
        allocate(op%mat(q)%gam(4)%X(nb,nb)) !Vphph
        allocate(op%mat(q)%gam(2)%X(np,nb)) !Vppph
        allocate(op%mat(q)%gam(6)%X(nb,nh)) !Vphhh
     
        do i = 1,6
           op%mat(q)%gam(i)%X = 0.0
        end do
    
     end do
  
  else  ! tensor operator


     allocate(op%tblck(op%nblocks)) 
    
     allocate(op%direct_omp(size(H%direct_omp)))
     op%direct_omp = H%direct_omp

     N = op%nsp
     allocate(op%xmap(N*(N+1)/2))
     do q = 1, N*(N+1)/2 
        allocate(op%xmap(q)%Z(size(H%xmap(q)%Z)))
        op%xmap(q)%Z = H%xmap(q)%Z
     end do
  
     if ( allocated( H%exlabels ) ) then 
        allocate(op%exlabels(size(H%exlabels(:,1)),2)) 
     else 
        allocate(op%exlabels(1,2),H%exlabels(1,2)) 
        H%exlabels = 0
     end if
     op%exlabels = H%exlabels
  
     do q = 1, op%nblocks
     
        op%tblck(q)%lam = H%tblck(q)%lam
        op%tblck(q)%Jpair = H%tblck(q)%Jpair
        nh = H%tblck(q)%nhh1
        np = H%tblck(q)%npp1
        nb = H%tblck(q)%nph1
         
        op%tblck(q)%npp1 = np
        op%tblck(q)%nph1 = nb
        op%tblck(q)%nhh1 = nh
        op%tblck(q)%ntot = H%tblck(q)%ntot
        op%tblck(q)%npp2 = H%tblck(q)%npp2
        op%tblck(q)%nph2 = H%tblck(q)%nph2
        op%tblck(q)%nhh2 = H%tblck(q)%nhh2

 
        allocate(op%tblck(q)%tensor_qn(1,1)%Y(np,2)) !qnpp
        allocate(op%tblck(q)%tensor_qn(2,1)%Y(nb,2)) !qnph
        allocate(op%tblck(q)%tensor_qn(3,1)%Y(nh,2)) !qnhh

        nh2 = H%tblck(q)%nhh2
        np2 = H%tblck(q)%npp2
        nb2 = H%tblck(q)%nph2
        
        allocate(op%tblck(q)%tensor_qn(1,2)%Y(np2,2)) !qnpp
        allocate(op%tblck(q)%tensor_qn(2,2)%Y(nb2,2)) !qnph
        allocate(op%tblck(q)%tensor_qn(3,2)%Y(nh2,2)) !qnhh
     
        do i = 1,2
           do j = 1, 3
              op%tblck(q)%tensor_qn(j,i)%Y = H%tblck(q)%tensor_qn(j,i)%Y
           end do 
        end do
    
        allocate(op%tblck(q)%tgam(1)%X(np,np2)) !Vpppp
        allocate(op%tblck(q)%tgam(5)%X(nh,nh2)) !Vhhhh
        allocate(op%tblck(q)%tgam(3)%X(np,nh2)) !Vpphh
        allocate(op%tblck(q)%tgam(4)%X(nb,nb2)) !Vphph
        allocate(op%tblck(q)%tgam(2)%X(np,nb2)) !Vppph
        allocate(op%tblck(q)%tgam(6)%X(nb,nh2)) !Vphhh
        allocate(op%tblck(q)%tgam(7)%X(nh,np2)) !Vhhpp
        allocate(op%tblck(q)%tgam(8)%X(nb,np2)) !Vphpp
        allocate(op%tblck(q)%tgam(9)%X(nh,nb2)) !Vhhph
        
        do i = 1,9
           op%tblck(q)%tgam(i)%X = 0.0
        end do
   
     end do

  end if 
end subroutine 
!=====================================================
!=====================================================
subroutine copy_sq_op(H,op) 
  ! make a copy of H onto op
  implicit none 
  
  type(sq_op) :: H,op
  integer :: q,i,j,holes,parts,nh,np,nb
     
  op%herm = H%herm
  op%E0 = H%E0
  op%fhh = H%fhh
  op%fpp = H%fpp
  op%fph = H%fph
  
  do q = 1, op%nblocks
              
     do i = 1,6
        op%mat(q)%gam(i)%X = H%mat(q)%gam(i)%X
     end do 
     
  end do 
  
end subroutine 
!======================================================
!======================================================
subroutine write_binary_operator(H,stage) 
  ! first figure out how many equations there are:
 
  type(sq_op) :: H 
  integer Atot,Ntot,q
  real(8),allocatable,dimension(:):: outvec 
  character(*),intent(in) :: stage 
  character(200) :: spfile,intfile,input,prefix
  common /files/ spfile,intfile,prefix 
  
  if (prefix(1:8) == 'testcase') return  
  Atot = H%belowEF
  Ntot = H%Nsp
  
  neq = 1 + Atot*Atot + Atot*(Ntot-Atot) + (Ntot - Atot)**2 
  
  do q = 1, H%nblocks
     
     nh = H%mat(q)%nhh
     np = H%mat(q)%npp
     nb = H%mat(q)%nph 
     
     neq = neq + (nh*nh+nh +  nb*nb+nb + np*np+np)/2 + nb*np + nh*np + nh*nb 
  end do 
  
  H%neq = neq
  allocate(outvec(neq)) 
  
  call vectorize(H,outvec) 
  
  open(unit=55,file='../../hamiltonians/'//&
       trim(adjustl(prefix))//'_'//trim(adjustl(stage)),&
       form='unformatted')
  
  do q =1,neq
     write(55) outvec(q)
  end do
  close(55) 
  
end subroutine
!======================================================
!======================================================
subroutine read_binary_operator(H,stage) 
  ! first figure out how many equations there are:
 
  type(sq_op) :: H 
  integer Atot,Ntot,q
  real(8),allocatable,dimension(:):: outvec 
  character(*),intent(in) :: stage 
  character(200) :: spfile,intfile,input,prefix
  common /files/ spfile,intfile,prefix 
  
  Atot = H%belowEF
  Ntot = H%Nsp
  
  neq = 1 + Atot*Atot + Atot*(Ntot-Atot) + (Ntot - Atot)**2 
  
  do q = 1, H%nblocks
     
     nh = H%mat(q)%nhh
     np = H%mat(q)%npp
     nb = H%mat(q)%nph 
     
     neq = neq + (nh*nh+nh +  nb*nb+nb + np*np+np)/2 + nb*np + nh*np + nh*nb 
  end do 
  H%neq = neq
  allocate(outvec(neq)) 
   
  open(unit=55,file='../../hamiltonians/'//&
       trim(adjustl(prefix))//'_'//trim(adjustl(stage)),&
       form='unformatted')
  
  do q =1,neq
     read(55) outvec(q) 
  end do
  close(55) 
  call repackage(H,outvec) 
end subroutine
!======================================================
!======================================================
subroutine split_1b_2b(Op,onebd,twobd) 
  ! make a copy of H onto op
  implicit none 
  
  type(sq_op) :: Op,onebd,twobd
  integer :: q,i,j,holes,parts,nh,np,nb

  onebd%herm = Op%herm
  twobd%herm = Op%herm 
  
  onebd%E0 = Op%E0
  onebd%fhh = Op%fhh
  onebd%fpp = Op%fpp
  onebd%fph = Op%fph
  
  do q = 1, op%nblocks
     i=3         
     do i = 1,6
        twobd%mat(q)%gam(i)%X = Op%mat(q)%gam(i)%X
     end do 
     
  end do 
  
end subroutine 

!=====================================================
!=====================================================
subroutine add_sq_op(A,ax,B,bx,C) 
  ! A*ax + B*bx = C 
  implicit none 
  
  type(sq_op) :: A,B,C
  real(8) :: ax,bx
  integer :: q,i,j,holes,parts,nh,np,nb
     
  C%E0 = A%E0*ax + B%E0*bx
  C%fhh = A%fhh*ax + B%fhh*bx
  C%fpp = A%fpp*ax + B%fpp*bx
  C%fph = A%fph*ax + B%fph*bx
  
  do q = 1, A%nblocks
              
     do i = 1,6
        C%mat(q)%gam(i)%X = A%mat(q)%gam(i)%X*ax + B%mat(q)%gam(i)%X*bx
     end do 
     
  end do 
  
end subroutine 
!=====================================================
!=====================================================
!=====================================================
!=====================================================
subroutine scale_sq_op(A,ax) 
  ! A = A*ax 
  implicit none 
  
  type(sq_op) :: A
  real(8) :: ax
  integer :: q,i,j,holes,parts,nh,np,nb
     
  A%E0 = A%E0*ax 
  A%fhh = A%fhh*ax
  A%fpp = A%fpp*ax
  A%fph = A%fph*ax 
  
  do q = 1, A%nblocks
              
     do i = 1,6
        A%mat(q)%gam(i)%X = A%mat(q)%gam(i)%X*ax 
     end do 
     
  end do 
  
end subroutine 
!=====================================================
!=====================================================
subroutine clear_sq_op(C) 
  ! A*ax + B*bx = C 
  implicit none 
  
  type(sq_op) :: C
  integer :: q,i,j,holes,parts,nh,np,nb
     
  C%E0 = 0.d0
  C%fhh = 0.d0
  C%fpp = 0.d0
  C%fph = 0.d0
  
  do q = 1, C%nblocks
              
     do i = 1,6
        C%mat(q)%gam(i)%X = 0.d0
     end do 
     
  end do 
  
end subroutine 
!=====================================================
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
  ! assume i < j 
  implicit none 
  
  integer :: i,j,n
  
  fermionic_tp_index = n*(i-1) + (i-i*i)/2 + j - i 
  
end function 
!===================================================== 
!=====================================================  
subroutine print_matrix(matrix)
	implicit none 
	
	integer :: i,m
	real(8),dimension(:,:) :: matrix
	character(1) :: y
	character(10) :: fmt2

    m=size(matrix(1,:))

    write(y,'(i1)') m
  
    fmt2= '('//y//'(f14.8))'	
	
	print*
	do i=1,m
	   write(*,fmt2) matrix(i,:)
	end do
	print* 
	
end subroutine 
!===============================================  
subroutine read_main_input_file(input,H,htype,HF,method,EXcalc,COM,R2RMS,&
     ME2J,ME2b,MORTBIN,hw,skip_setup,skip_gs)
  !read inputs from file
  implicit none 
  
  character(200) :: spfile,intfile,input,prefix
  character(50) :: valence
  type(sq_op) :: H 
  integer :: htype,jx,jy,Jtarg,Ptarg,excalc,com_int,rrms_int
  integer :: method,Exint,ISTAT ,i
  logical :: HF,COM,R2RMS,ME2J,ME2B,skip_setup,skip_gs,MORTBIN, found
  real(8) :: hw
  common /files/ spfile,intfile,prefix 
    
  input = adjustl(input) 
  if (trim(input) == '') then 
     print*, 'RUNNING TEST CASE: testcase.ini'
     i = 1 
     found = .false. 
     do while (.not. (found))   
        INI_DIR = INI_DIRECTORY_LIST(i) 
        inquire(file=trim(INI_DIR)//'testcase.ini',exist=found) 
        i = i + 1
     end do
     open(unit=22,file=trim(INI_DIR)//'testcase.ini')
  else       
     i = 1 
     found = .false. 
     do while (.not. (found))   
        INI_DIR = INI_DIRECTORY_LIST(i) 
        inquire(file=trim(INI_DIR)//trim(input),exist=found)      
        i = i + 1
     end do
     open(unit=22,file=trim(INI_DIR)//trim(input)) 
  end if 
  
  read(22,*);read(22,*);read(22,*)
  read(22,*);read(22,*);read(22,*)

  read(22,*) prefix
  read(22,*) 
  read(22,*) intfile
  read(22,*)
  read(22,*) spfile
  read(22,*);read(22,*)
  read(22,*) htype
  read(22,*)
  read(22,*) hw
  read(22,*)
  read(22,*) H%Aprot
  read(22,*)
  read(22,*) H%Aneut
  read(22,*);read(22,*)
  read(22,*) jx
  read(22,*);read(22,*)
  read(22,*) method
  read(22,*);read(22,*)
  read(22,*) excalc
  read(22,*)
  read(22,*) Jtarg
  read(22,*)
  read(22,*) Ptarg
  read(22,*)
  read(22,*) valence
  read(22,*) 
  read(22,*) com_int
  read(22,*)
  read(22,*) rrms_int 
  read(22,*)
  read(22,*) H%lawson_beta,H%com_hw
 
  valence = adjustl(valence)
  H%Jtarg = Jtarg
  H%Ptarg = Ptarg
 
  HF = .false. 
  if (jx == 1) HF = .true. 
  COM = .false.
  if (com_int == 1) COM = .true.
  R2RMS = .false.
  if (rrms_int == 1) R2RMS = .true.
  
  select case (trim(valence)) 
     case ('0p') 
        H%valcut = 6 
     case ('1s0d')
        H%valcut = 12
     case ('1p0f')
        H%valcut = 20 
     case ( '2s1d0g') ! probably not needed 
        H%valcut = 30
     case default
        stop 'valence space not listed in database, SEE basic_IMSRG.f90' 
  end select 
        
  me2j = .true.
  me2b = .false. 
  MORTBIN=.false.
  if( intfile(len(trim(intfile))-3:len(trim(intfile))) == '.int') then 
     me2j =.false.
     ! Morten format
  else if( intfile(len(trim(intfile))-3:len(trim(intfile))) == '.bin') then 
     me2j =.false.
     mortbin = .true. 
     ! binary version of Morten format
  else if( intfile(len(trim(intfile))-6:len(trim(intfile))) == 'me2b.gz') then 
     me2j =.false.
     me2b = .true.
     ! normal ordered Heiko format
  else 
     ! Heiko format
     if (spfile(1:2) .ne. 'hk') then 
        STOP 'inconsistent interaction and sps files' 
     end if 
  end if 
  
  inquire( file='../../hamiltonians/'//&
       trim(adjustl(prefix))//'_bare' , exist = skip_setup )

  if (EXcalc .ne. 0) then 
  inquire( file='../../hamiltonians/'//&
       trim(adjustl(prefix))//'_gs_decoup' , exist = skip_gs )
  else 
     skip_gs = .false.
  end if

  ! figure out where the TBME and SP files are....
  i = 1 
  found = .false. 
  do while (.not. (found))   
     TBME_DIR = TBME_DIRECTORY_LIST(i) 
     inquire(file=trim(TBME_DIR)//trim(intfile),exist=found)      
     i = i + 1
  end do
     
  i = 1 
  found = .false. 
  do while (.not. (found))   
     SP_DIR = SP_DIRECTORY_LIST(i) 
     inquire(file=trim(SP_DIR)//trim(spfile),exist=found)      
     i = i + 1
  end do

end subroutine
!=======================================================  
!=======================================================
subroutine allocate_CCMAT(HS,CCME,jbas) 
  ! allocates a cross-coupled ME storage structure
  ! currently the only CCME of interest are phab terms    |---<---| 
  ! coupling in the 3-1 channel                        <(pa)J|V|(hb)J>
  !                                                      |---<---|
  implicit none 
  
  
  type(spd) :: jbas
  type(sq_op) :: HS
  type(cross_coupled_31_mat) :: CCME
  integer :: JT,ji,jp,jj,jh,JC,q1,q2,g,li,lj,ti,tj
  integer :: a,b,p,h,i,j,r,Jmin,Jmax,NX,TZ,PAR,x,JTM
  integer :: int1,int2,IX,JX,i1,i2,nb,nh,np,numJ
  real(8) :: sm,sm2
  
  NX = HS%Nsp
  CCME%Nsp = NX
  CCME%rank = 0
  CCME%dpar = 0
  CCME%herm = HS%herm
  JTM = jbas%Jtotal_max
  CCME%nblocks = (JTM + 1) * 2 * 2
  ! 2 dof for parity, 2 for Tz (i'm only worried about abs(Tz) ) 
  allocate(CCME%CCX(CCME%nblocks)) ! < h a |v | p b> 
  allocate(CCME%CCR(CCME%nblocks)) ! < p a |v | h b> 
  allocate(CCME%nph(CCME%nblocks)) ! number of ph pairs in block 
  allocate(CCME%rlen(CCME%nblocks)) ! number of total pairs in block
  allocate(CCME%Jval(CCME%nblocks)) ! J value for block
  allocate(CCME%rmap(NX*NX))  ! map of pair index to r indeces
  allocate(CCME%qmap(NX*NX))  ! map of pair index to q indeces
  allocate(CCME%nbmap(NX*NX)) 
  
  do i = 1, NX
     do j = 1,NX
        
        ji = jbas%jj(i) 
        jj = jbas%jj(j) 
        
        numJ = (ji + jj - abs(ji-jj))/2 + 1
        
        x = CCindex(i,j,NX) 
        allocate(CCME%rmap(x)%Z(numJ)) 
        allocate(CCME%qmap(x)%Z(numJ))
        allocate(CCME%nbmap(x)%Z(numJ))
        CCME%qmap(x)%Z = 0
        CCME%nbmap(x)%Z = 0
     end do 
  end do 
  
  
  do q1 = 1, CCME%nblocks
     
     JC = mod(q1-1,JTM+1) * 2 
     PAR = (q1 - 1) / (2*JTM + 2) 
     TZ = mod(q1-1,(2*JTM+2))/(JTM+1)  
     ! fastest changing quantity : JC
     ! slowest: PAR 
     
     nb = 0 
     r = 0 
     do i = 1, NX
        do j = 1,NX 
           
           ji = jbas%jj(i) 
           jj = jbas%jj(j) 
           
           if (.not. (triangle(ji,jj,JC))) cycle 
           if ( mod(jbas%ll(i) + jbas%ll(j),2) .ne. PAR ) cycle
           if (abs(jbas%itzp(i) - jbas%itzp(j))/2 .ne. Tz ) cycle 
           
           if ( (jbas%con(i) == 0 ).and. (jbas%con(j) == 1)) then 
              nb = nb + 1 
           end if 
           
           r = r+1
        end do 
     end do 
     
     allocate( CCME%CCX(q1)%X(r,nb) ) 
     allocate( CCME%CCR(q1)%X(nb,r) ) 
     
     nb = 0 
     r = 0 
     do i = 1, NX
        do j = 1,NX 
           
           ji = jbas%jj(i) 
           jj = jbas%jj(j) 
           if (.not. (triangle(ji,jj,JC))) cycle 
           if ( mod(jbas%ll(i) + jbas%ll(j),2) .ne. PAR ) cycle
           if (abs(jbas%itzp(i) - jbas%itzp(j))/2 .ne. Tz ) cycle 
          
           x = CCindex(i,j,NX) 
           
           g = 1
           do while (CCME%qmap(x)%Z(g) .ne. 0) 
              g = g + 1
           end do 
         
           if ( (jbas%con(i) == 0 ).and. (jbas%con(j) == 1)) then
              nb = nb + 1
              CCME%nbmap(x)%Z(g) = nb 
           end if
           
           r = r+1
           
           CCME%qmap(x)%Z(g) = q1
           CCME%rmap(x)%Z(g) = r

        end do 
     end do 
               
     CCME%nph(q1) = nb
     CCME%rlen(q1) = r
     CCME%Jval(q1) = JC 
  end do 
end subroutine  
!=======================================================  
!=======================================================
subroutine allocate_tensor_CCMAT(OP,CCME,jbas) 
  ! allocates a cross-coupled ME storage structure
  ! currently the only CCME of interest are phab terms    |---<---| 
  ! coupling in the 3-1 channel                        <(pa)J|V|(hb)J>
  !                                                      |---<---|
  implicit none 
  
  
  type(spd) :: jbas
  type(sq_op) :: OP
  type(cross_coupled_31_mat) :: CCME
  integer :: Jtot1,Jtot2,ji,jp,jj,jh,JC,q1,q2,g,li,lj,ti,tj,q
  integer :: a,b,p,h,i,j,Jmin,Jmax,NX,TZ,PAR,x,JTM,RANK,Jold
  integer :: int1,int2,IX,JX,i1,i2,nb1,nb2,r1,r2,nh,np,numJ
  real(8) :: sm,sm2
  
  NX = OP%Nsp
  RANK = OP%rank
  CCME%rank = OP%rank
  CCME%dpar = OP%dpar
  CCME%Nsp = NX
  CCME%herm = OP%herm  
  JTM = jbas%Jtotal_max
  Jold = 1
! quantum numbers of the last block 

  Tz = 1 
  Par = 1
  Jtot1 = jbas%Jtotal_max*2 
  Jtot2 = Jtot1+ RANK 
  CCME%nblocks  =  tensor_block_index(Jtot1,Jtot2,RANK,Tz,Par)/3*2 
  ! 2 dof for parity, 2 for Tz (i'm only worried about abs(Tz) ) 

  allocate(CCME%CCX(CCME%nblocks)) ! < h a |v | p b>    h(p)b(a)
  allocate(CCME%CCR(CCME%nblocks)) ! < p a |v | h b>    p(h)a(b)
  allocate(CCME%nph(CCME%nblocks)) ! number of ph pairs in block 
  allocate(CCME%rlen(CCME%nblocks)) ! number of total pairs in block
  allocate(CCME%Jval(CCME%nblocks)) ! J value for block
  allocate(CCME%Jval2(CCME%nblocks)) ! J2 value for the block
  allocate(CCME%rmap(NX*NX))  ! map of pair index to r indeces
  allocate(CCME%qmap(NX*NX))  ! map of pair index to q indeces
  allocate(CCME%nbmap(NX*NX)) 
  
  do i = 1, NX
     do j = 1,NX
        
        ji = jbas%jj(i) 
        jj = jbas%jj(j) 
        
        numJ = (ji + jj - abs(ji-jj))/2 + 1
        
        x = CCindex(i,j,NX) 
        allocate(CCME%rmap(x)%Z(numJ)) 
        allocate(CCME%qmap(x)%Z(numJ))
        allocate(CCME%nbmap(x)%Z(numJ))
        CCME%qmap(x)%Z = 0
        CCME%nbmap(x)%Z = 0
     end do 
  end do 
  
  q1 = 0
  do Jtot1 = 0,2*jbas%jtotal_max,2 
     do Jtot2 = max(abs(Jtot1 - rank),Jtot1),Jtot1+rank,2
        do Tz = 0, 1    
           do PAR = 0,1 
              
              q1 = q1 + 1
              CCME%Jval(q1) = Jtot1
              CCME%Jval2(q1) = Jtot2
              
              if (jtot2 > 2*jbas%jtotal_max) cycle
              
              !PAR = mod(q1-1,2)
              !Tz = mod((q1-1)/2,2) 

              ! fastest changing quantity : JC
              ! slowest: PAR 

              nb1 = 0 
              r1 = 0 
              nb2 = 0 
              r2 = 0

              do i = 1, NX
                 do j = 1,NX 

                    ji = jbas%jj(i) 
                    jj = jbas%jj(j) 


                    if (abs(jbas%itzp(i) - jbas%itzp(j))/2 .ne. Tz ) cycle 

                    if ( mod(jbas%ll(i) + jbas%ll(j),2) == PAR ) then
                       if (triangle(ji,jj,Jtot1)) then 
                          if ( (jbas%con(i) == 0 ).and. (jbas%con(j) == 1)) then 
                             nb1 = nb1 + 1 
                          end if
                       
                          r1 = r1+1
                       end if
                    end if
                    
                    if ( mod(jbas%ll(i) + jbas%ll(j),2) == mod(PAR+op%dpar/2,2) ) then
                       if (triangle(ji,jj,Jtot2)) then 
                          if ( (jbas%con(i) == 0 ).and. (jbas%con(j) == 1)) then 
                             nb2 = nb2 + 1 
                          end if
                          
                          r2 = r2+1
                       end if
                    end if
                 end do
              end do

              allocate( CCME%CCX(q1)%X(r1,nb2) ) 
              allocate( CCME%CCR(q1)%X(nb1,r2) ) 


              CCME%nph(q1) = nb1
              CCME%rlen(q1) = r1



              nb1 = 0 
              r1 = 0 
              nb2 = 0 
              r2 = 0

              if (max(abs(Jtot1 - rank),Jtot1) .ne. Jtot2) cycle
              ! I only need one of these arrays per J, so I use the same shape as before. 
              do i = 1, NX
                 do j = 1,NX 
                   

                    ji = jbas%jj(i) 
                    jj = jbas%jj(j) 

                    if (.not. (triangle(ji,jj,Jtot1))) cycle 
                    if ( mod(jbas%ll(i) + jbas%ll(j),2) .ne. PAR ) cycle
                    if (abs(jbas%itzp(i) - jbas%itzp(j))/2 .ne. Tz ) cycle 

                    x = CCindex(i,j,NX) 

                    g = 1
                    do while (CCME%qmap(x)%Z(g) .ne. 0) 
                       g = g + 1
                    end do

                    if ( (jbas%con(i) == 0 ).and. (jbas%con(j) == 1)) then
                       nb1 = nb1 + 1
                       CCME%nbmap(x)%Z(g) = nb1 
                    end if

                    r1 = r1+1
                    q = block_index(Jtot1,Tz,Par)          
                    CCME%qmap(x)%Z(g) = q
                    CCME%rmap(x)%Z(g) = r1
                    
                    
                 end do
              end do
              
             ! Jold = Jtot1
              

           end do
        end do
     end do
  end do

end subroutine  
!=======================================================  
!=======================================================
subroutine duplicate_CCMAT(C1,CCME) 
  ! makes a copy of C1 onto CCME
  implicit none 
  
  type(cross_coupled_31_mat) :: CCME,C1
  integer :: JT,ji,jp,jj,jh,JC,q1,q2
  integer :: a,b,p,h,i,j,r,Jmin,Jmax,NX,TZ,PAR
  integer :: int1,int2,IX,JX,i1,i2,nb,nh,np,numJ
  real(8) :: sm,sm2
  
  NX = C1%Nsp
  CCME%Nsp = NX
  CCME%rank = C1%rank
  CCME%dpar = C1%dpar
  CCME%nblocks = C1%nblocks
  CCME%herm = C1%herm
  allocate(CCME%CCX(C1%nblocks))
  allocate(CCME%CCR(C1%nblocks))
  allocate(CCME%nph(C1%nblocks))
  allocate(CCME%rlen(C1%nblocks)) 
  allocate(CCME%Jval(C1%nblocks)) 
  allocate(CCME%rmap(NX*NX))  ! map of pair index to r indeces
  allocate(CCME%qmap(NX*NX))  ! map of pair index to q indeces
  allocate(CCME%nbmap(NX*NX)) ! map of ph pair index to nb indeces
  
  do i = 1, NX
     do j = 1,NX
        
        r = CCindex(i,j,NX)
        
        numJ = size(C1%rmap(r)%Z) 
        
        allocate(CCME%rmap(r)%Z(numJ)) 
        allocate(CCME%qmap(r)%Z(numJ))
        allocate(CCME%nbmap(r)%Z(numJ)) 
        CCME%rmap(r)%Z = c1%rmap(r)%Z
        CCME%qmap(r)%Z = c1%qmap(r)%Z
        CCME%nbmap(r)%Z = c1%nbmap(r)%Z
       
     end do 
  end do 
   
  CCME%nph = C1%nph
  CCME%rlen = C1%rlen
  CCME%Jval = C1%Jval
  if (allocated(C1%Jval2)) then
     allocate(CCME%Jval2(C1%nblocks)) 
     CCME%Jval2 = C1%Jval2
  end if 
  
  do q1 = 1, CCME%nblocks
     
     r = CCME%rlen(q1)
     nb = CCME%nph(q1) 
     allocate( CCME%CCX(q1)%X(r,nb) ) 
     allocate( CCME%CCR(q1)%X(nb,r) ) 
     
  end do 
end subroutine              
!=======================================================  
!=======================================================          
subroutine calculate_cross_coupled(HS,CCME,jbas,phase) 
  ! currently the only CCME of interest are phab terms    |---<---| 
  ! coupling in the 3-1 channel                        <(pa)J|V|(hb)J>
  !                                                      |---<---|
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: HS
  type(cross_coupled_31_mat) :: CCME
  integer :: JT,ja,jp,jb,jh,JC,q1,q2,TZ,PAR,la,lb,Ntot,th,tp,lh,lp
  integer :: a,b,p,h,i,j,Jmin,Jmax,Rindx,Gindx,g,ta,tb,Atot,hg,pg
  integer :: int1,int2,IX,JX,i1,i2,nb,nh,np,gnb,NBindx,x,JTM
  real(8) :: sm,sm2,pre,horse
  logical :: phase

  Atot = HS%belowEF
  Ntot = HS%Nsp
  JTM = jbas%Jtotal_max 
  pre = 1.d0 
  CCME%herm = HS%Herm
  
!$omp parallel do default(firstprivate),shared(CCME,HS,jbas) 
  do q1 = 1, CCME%nblocks
      
     JC = mod(q1-1,JTM+1) * 2 
     PAR = (q1 - 1) / (2*JTM + 2) 
     TZ = mod(q1-1,(2*JTM+2))/(JTM+1)  
     
     CCME%CCR(q1)%X = 0.d0
     CCME%CCX(q1)%X = 0.d0
         
     ! ab = ph 
     do hg = 1, Atot
        do pg = 1, Ntot - Atot 
           
           h = jbas%holes(hg) 
           p = jbas%parts(pg) 
           
           jp = jbas%jj(p) 
           jh = jbas%jj(h)
           lp = jbas%ll(p) 
           lh = jbas%ll(h)
           tp = jbas%itzp(p) 
           th = jbas%itzp(h)
        
           if (.not. triangle(jp,jh,JC) )  cycle
           if ( mod(lp + lh,2) .ne. PAR ) cycle
           if (abs(tp - th)/2 .ne. Tz ) cycle 
           
           x = CCindex(p,h,HS%Nsp)
           gnb = 1
           do while (CCME%qmap(x)%Z(gnb) .ne. q1 )
              gnb = gnb + 1
           end do
              
           NBindx = CCME%nbmap(x)%Z(gnb) 

           if (phase) pre = (-1)**((jp +jh)/2) !convenient to have this 
           ! for the ph  channel 2body derivative 
        
           do a = 1, HS%nsp
              do b = 1, HS%nsp
  
                 ja = jbas%jj(a) 
                 jb = jbas%jj(b)
                 la = jbas%ll(a) 
                 lb = jbas%ll(b)
                 ta = jbas%itzp(a) 
                 tb = jbas%itzp(b)
                 
                 if (.not. triangle(ja,jb,JC) )  cycle
                 if ( mod(la + lb,2) .ne. PAR ) cycle
                 if (abs(ta - tb)/2 .ne. Tz ) cycle 
                 
                 x = CCindex(a,b,HS%Nsp) 
                 g = 1
                 do while (CCME%qmap(x)%Z(g) .ne. q1 )
                    g = g + 1
                 end do
              
                 Rindx = CCME%rmap(x)%Z(g)
                 
                 x = CCindex(b,a,HS%Nsp) 
                 g = 1
                 do while (CCME%qmap(x)%Z(g) .ne. q1 )
                    g = g + 1
                 end do
              
                 Gindx = CCME%rmap(x)%Z(g)
                 
                 sm = 0.d0 
               
!                 horse = 0.d0 
                 if ( (mod(la + lh,2) == mod(lb + lp,2)) .and. &
                      ( (ta + th) == (tb + tp) ) ) then  
               
                    ! hapb 
                    Jmin = max(abs(jp - jb),abs(ja - jh)) 
                    Jmax = min(jp+jb,ja+jh) 
                    
                    sm = 0.d0 
                    do JT = Jmin,Jmax,2
                       sm = sm + (-1)**(JT/2) * (JT + 1) * &
                            sixj(jp,jh,JC,ja,jb,JT)  * &
                            v_elem(h,a,p,b,JT,HS,jbas) 
                    end do
                 
                    ! store  < h a | v | p b>    Pandya ( V )_h(p)b(a)
                    CCME%CCX(q1)%X(Rindx,NBindx) = sm * &
                         (-1) **( (jh + jb + JC) / 2) * pre * sqrt(JC + 1.d0)
                    ! scaled by sqrt(JC + 1) for convience in ph derivative

                    ! store < p b | v | h a>     Pandya ( V )_p(h)a(b)
                    CCME%CCR(q1)%X(NBindx,Gindx) = sm * HS%herm * &
                         (-1) **( (jp + ja + JC) / 2) * pre * sqrt(JC + 1.d0)
                 
                 end if
              end do
           end do
        end do
     end do

  end do 
!$omp end parallel do

end subroutine 
!=======================================================  
!=======================================================          
subroutine calculate_cross_coupled_pphh(HS,CCME,jbas,phase) 
  ! currently the only CCME of interest are phab terms    |---<---| 
  ! coupling in the 3-1 channel                        <(pa)J|V|(hb)J>
  !                                                      |---<---|
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: HS
  type(cross_coupled_31_mat) :: CCME
  integer :: JT,ja,jp,jb,jh,JC,q1,q2,TZ,PAR,la,lb,Ntot,th,tp,lh,lp
  integer :: a,b,p,h,i,j,Jmin,Jmax,Rindx,Gindx,g,ta,tb,Atot,hg,pg
  integer :: int1,int2,IX,JX,i1,i2,nb,nh,np,gnb,NBindx,x,JTM,ax,bx
  real(8) :: sm,sm2,pre,horse
  logical :: phase

  Atot = HS%belowEF
  Ntot = HS%Nsp
  JTM = jbas%Jtotal_max 
  pre = 1.d0 
  CCME%herm = HS%Herm
  
!$omp parallel do default(firstprivate),shared(CCME,HS,jbas) 
  do q1 = 1, CCME%nblocks
      
     JC = mod(q1-1,JTM+1) * 2 
     PAR = (q1 - 1) / (2*JTM + 2) 
     TZ = mod(q1-1,(2*JTM+2))/(JTM+1)  
     
     CCME%CCR(q1)%X = 0.d0
     CCME%CCX(q1)%X = 0.d0
         
     ! ab = ph 
     do hg = 1, Atot
        do pg = 1, Ntot - Atot 
           
           h = jbas%holes(hg) 
           p = jbas%parts(pg) 
           
           jp = jbas%jj(p) 
           jh = jbas%jj(h)
           lp = jbas%ll(p) 
           lh = jbas%ll(h)
           tp = jbas%itzp(p) 
           th = jbas%itzp(h)
        
           if (.not. triangle(jp,jh,JC) )  cycle
           if ( mod(lp + lh,2) .ne. PAR ) cycle
           if (abs(tp - th)/2 .ne. Tz ) cycle 
           
           x = CCindex(p,h,HS%Nsp)
           gnb = 1
           do while (CCME%qmap(x)%Z(gnb) .ne. q1 )
              gnb = gnb + 1
           end do
              
           NBindx = CCME%nbmap(x)%Z(gnb) 

           if (phase) pre = (-1)**((jp +jh)/2) !convenient to have this 
           ! for the ph  channel 2body derivative 
        
           do ax = 1, HS%belowEF
              a = jbas%holes(ax) 
              do bx = 1, HS%nsp-HS%belowEF
                 b = jbas%parts(bx) 
                 
                 ja = jbas%jj(a) 
                 jb = jbas%jj(b)
                 la = jbas%ll(a) 
                 lb = jbas%ll(b)
                 ta = jbas%itzp(a) 
                 tb = jbas%itzp(b)
                 
                 if (.not. triangle(ja,jb,JC) )  cycle
                 if ( mod(la + lb,2) .ne. PAR ) cycle
                 if (abs(ta - tb)/2 .ne. Tz ) cycle 
                 
                 x = CCindex(a,b,HS%Nsp) 
                 g = 1
                 do while (CCME%qmap(x)%Z(g) .ne. q1 )
                    g = g + 1
                 end do
              
                 Rindx = CCME%rmap(x)%Z(g)
                 
                 x = CCindex(b,a,HS%Nsp) 
                 g = 1
                 do while (CCME%qmap(x)%Z(g) .ne. q1 )
                    g = g + 1
                 end do
              
                 Gindx = CCME%rmap(x)%Z(g)
                 
                 sm = 0.d0 
               
!                 horse = 0.d0 
                 if ( (mod(la + lh,2) == mod(lb + lp,2)) .and. &
                      ( (ta + th) == (tb + tp) ) ) then  
               
                    ! hapb 
                    Jmin = max(abs(jp - jb),abs(ja - jh)) 
                    Jmax = min(jp+jb,ja+jh) 
                    
                    sm = 0.d0 
                    do JT = Jmin,Jmax,2
                       sm = sm + (-1)**(JT/2) * (JT + 1) * &
                            sixj(jp,jh,JC,ja,jb,JT)  * &
                            v_elem(h,a,p,b,JT,HS,jbas) 
                    end do
                 
                    ! store  < h a | v | p b>    Pandya ( V )_h(p)b(a)
                    CCME%CCX(q1)%X(Rindx,NBindx) = sm * &
                         (-1) **( (jh + jb + JC) / 2) * pre * sqrt(JC + 1.d0)
                    ! scaled by sqrt(JC + 1) for convience in ph derivative

                    ! store < p b | v | h a>     Pandya ( V )_p(h)a(b)
                    CCME%CCR(q1)%X(NBindx,Gindx) = sm * HS%herm * &
                         (-1) **( (jp + ja + JC) / 2) * pre * sqrt(JC + 1.d0)
                 
                 end if
              end do
           end do
        end do
     end do

  end do 
!$omp end parallel do

end subroutine 
!=======================================================  
!=======================================================          
subroutine calculate_generalized_pandya(OP,CCME,jbas,phase) 
  ! currently the only CCME of interest are phab terms    |---<--|  J1 
  ! coupling in the 3-1 channel                        <(pa)|V|(hb)> rank
  !                                                      |---<--| J2 
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: OP 
  type(cross_coupled_31_mat) :: CCME
  integer :: Jtot1,Jtot2,ja,jp,jb,jh,JC,q1,q2,q,TZ,PAR,la,lb,Ntot,th,tp,lh,lp
  integer :: a,b,p,h,i,j,Jmin1,Jmax1,Rindx,Gindx,g,ta,tb,Atot,hg,pg,J3,J4,NBindx2,qONE,qTWO
  integer :: int1,int2,IX,JX,i1,i2,nb,nh,np,gnb,NBindx1,x,JTM,rank,Jmin2,Jmax2
  real(8) :: sm,sm2,pre,horse
  logical :: phase,parflip

  Atot = OP%belowEF
  Ntot = OP%Nsp
  JTM = jbas%Jtotal_max 
  pre = 1.d0 
  rank = OP%rank

  parflip = .false. 
  if ( mod(op%dpar/2,2) == 1) parflip = .true. 

  CCME%herm = OP%Herm
!$omp parallel do default(firstprivate),shared(CCME,OP,jbas) 
  do q1 = 1, CCME%nblocks
      
     Jtot1 = CCME%Jval(q1)
     Jtot2 = CCME%Jval2(q1)
     if (Jtot2 > 2*JTM) cycle
     
     PAR = mod(q1-1,2)
     Tz = mod((q1-1)/2,2) 
     
     CCME%CCR(q1)%X = 0.d0
     CCME%CCX(q1)%X = 0.d0
         
     qONE = block_index(Jtot1,Tz,Par)
     qTWO = block_index(Jtot2,Tz,mod(Par+op%dpar/2,2))
     
     ! ab = ph 
     do hg = 1, Atot
        do pg = 1, Ntot - Atot 
           
           h = jbas%holes(hg) 
           p = jbas%parts(pg) 
           
           jp = jbas%jj(p) 
           jh = jbas%jj(h)
           lp = jbas%ll(p) 
           lh = jbas%ll(h)
           tp = jbas%itzp(p) 
           th = jbas%itzp(h)
        
           if (.not. (parflip)) then 
              if ( mod(lp + lh,2) .ne. PAR ) cycle
           end if 
           
           if (abs(tp - th)/2 .ne. Tz ) cycle 
        
           NBindx1 = 0
           NBindx2 = 0 
           if ( triangle(jp,jh,Jtot1) )  then 
              
              
              if ( mod(lp + lh,2) ==  PAR ) then   
                 x = CCindex(p,h,OP%Nsp)
                 gnb = 1
              
                 do while (CCME%qmap(x)%Z(gnb) .ne. qONE )
                    gnb = gnb + 1 
                 end do
              
                 NBindx1 = CCME%nbmap(x)%Z(gnb) 
             
              end if 
              
              if  ( triangle(jp,jh,Jtot2) )  then 
              
                 if ( mod(lp+lh+op%dpar/2,2) == PAR) then 
        
                    x = CCindex(p,h,OP%Nsp)
                    gnb = 1
                    
                    do while (CCME%qmap(x)%Z(gnb) .ne. qTWO )
                       gnb = gnb + 1 
                    end do
                 
                    NBindx2 = CCME%nbmap(x)%Z(gnb)
                 end if 
              end if
 
           else if  ( triangle(jp,jh,Jtot2) )  then 
              
              if ( mod(lp+lh+op%dpar/2,2) == PAR) then 
        
                 x = CCindex(p,h,OP%Nsp)
                 gnb = 1
              
                 do while (CCME%qmap(x)%Z(gnb) .ne. qTWO )
                    gnb = gnb + 1 
                 end do
              
                 NBindx2 = CCME%nbmap(x)%Z(gnb)
              end if
           else 
              cycle
           end if 
          
           
           if (phase) pre = (-1)**((jp +jh)/2) !convenient to have this 
           ! for the ph  channel 2body derivative 
        
           do a = 1, OP%nsp
              do b = 1, OP%nsp
           
                 ja = jbas%jj(a) 
                 jb = jbas%jj(b)
                 la = jbas%ll(a) 
                 lb = jbas%ll(b)
                 ta = jbas%itzp(a) 
                 tb = jbas%itzp(b)
                 

                 if (.not. (parflip)) then 
                    if ( mod(la + lb,2) .ne. PAR ) cycle
                 end if 

                 if (abs(ta - tb)/2 .ne. Tz ) cycle 

       
                 if ( (triangle(ja,jb,Jtot1)) .and. (NBindx2 .ne. 0) ) then 
       
                    if ( mod(la+lb,2) == PAR ) then 
                       x = CCindex(a,b,OP%Nsp) 
                       
                       g = 1
                       do while (CCME%qmap(x)%Z(g) .ne. qONE )
                          g = g + 1
                       end do
              
                       Rindx = CCME%rmap(x)%Z(g)
                 
                       if ( (mod(la + lh,2) == mod(lb + lp + op%dpar/2,2)) .and. &
                            ( (ta + th) == (tb + tp) ) ) then  
                         
                          ! hapb 
                          Jmin1 = abs(ja - jh) 
                          Jmax1 = ja+jh 
                          Jmin2 = abs(jp - jb)
                          Jmax2 = jp+jb 
                          
                          sm = 0.d0 
                          do J3 = Jmin1,Jmax1,2
                             do J4 = Jmin2,Jmax2,2
                                sm = sm - (-1)**(J4/2) * sqrt((J3 + 1.d0) * (J4+1.d0)) * &
                                     ninej(ja,jb,Jtot1,jh,jp,Jtot2,J3,J4,rank) * &
                                     tensor_elem(a,h,p,b,J3,J4,OP,jbas) 
                                
                             end do
                          end do
                          
                     
                          CCME%CCX(q1)%X(Rindx,NBindx2) = sm * & 
                               (-1) **( (jh+jb+Jtot2) / 2) * pre * sqrt((Jtot1 + 1.d0)*(Jtot2 + 1.d0))

                       end if
                    end if
                 end if 
                 
                 if ((triangle(ja,jb,Jtot2)) .and. (NBindx1 .ne. 0) ) then 
                    
                    if ( mod(la+lb+op%dpar/2,2) == PAR ) then 
                       x = CCindex(b,a,OP%Nsp) 
                       g = 1
                       do while (CCME%qmap(x)%Z(g) .ne. qTWO )
                          g = g + 1
                       end do
                    
                       Gindx = CCME%rmap(x)%Z(g)
                 

                       if ( (mod(la + lh,2) == mod(lb + lp + op%dpar/2,2)) .and. &
                            ( (ta + th) == (tb + tp) ) ) then  
               
                          ! hapb 
                          Jmin1 = abs(ja - jh) 
                          Jmax1 = ja+jh 
                          Jmin2 = abs(jp - jb)
                          Jmax2 = jp+jb 
                    
                          sm = 0.d0 
                          do J3 = Jmin1,Jmax1,2
                             do J4 = Jmin2,Jmax2,2
                                sm = sm - (-1)**(J4/2) * sqrt((J3 + 1.d0) * (J4+1.d0)) * &
                                     ninej(jh,jp,Jtot1,ja,jb,Jtot2,J3,J4,rank) * &
                                     tensor_elem(h,a,b,p,J3,J4,OP,jbas) 
                       
                             end do
                          end do

                          ! store  ( V )_h(p)b(a)
                          
                          CCME%CCR(q1)%X(NBindx1,Gindx) = sm * &
                               (-1) **( (jp+ja+Jtot2)/2) * pre * sqrt((Jtot1 + 1.d0)*(Jtot2 + 1.d0))
                       end if
                    end if
                 end if 
                 
              end do
           end do
        end do
     end do

  end do 
!$omp end parallel do

end subroutine 
!=======================================================  
!=======================================================          
subroutine EOM_generalized_pandya(OP,CCME,jbas,phase) 
  ! currently the only CCME of interest are phab terms    |---<--|  J1 
  ! coupling in the 3-1 channel                        <(pa)|V|(hb)> rank
  !                                                      |---<--| J2 
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: OP 
  type(cross_coupled_31_mat) :: CCME
  integer :: Jtot1,Jtot2,ja,jp,jb,jh,JC,q1,q2,q,TZ,PAR,la,lb,Ntot,th,tp,lh,lp
  integer :: a,b,p,h,i,j,Jmin1,Jmax1,Rindx,Gindx,g,ta,tb,Atot,hg,pg,J3,J4,NBindx2,qONE,qTWO
  integer :: int1,int2,IX,JX,i1,i2,nb,nh,np,gnb,NBindx1,x,JTM,rank,Jmin2,Jmax2,bx,ax
  real(8) :: sm,sm2,pre,horse
  logical :: phase,parflip

  Atot = OP%belowEF
  Ntot = OP%Nsp
  JTM = jbas%Jtotal_max 
  pre = 1.d0 
  rank = OP%rank

  parflip = .false. 
  if ( mod(op%dpar/2,2) == 1) parflip = .true. 

  CCME%herm = OP%Herm
!$omp parallel do default(firstprivate),shared(CCME,OP,jbas) 
  do q1 = 1, CCME%nblocks
      
     Jtot1 = CCME%Jval(q1)
     Jtot2 = CCME%Jval2(q1)
     if (Jtot2 > 2*JTM) cycle
     
     PAR = mod(q1-1,2)
     Tz = mod((q1-1)/2,2) 
     
     CCME%CCR(q1)%X = 0.d0
     CCME%CCX(q1)%X = 0.d0
         
     qONE = block_index(Jtot1,Tz,Par)
     qTWO = block_index(Jtot2,Tz,mod(Par+op%dpar/2,2))
     
     ! ab = ph 
     do hg = 1, Atot
        do pg = 1, Ntot - Atot 
           
           h = jbas%holes(hg) 
           p = jbas%parts(pg) 
           
           jp = jbas%jj(p) 
           jh = jbas%jj(h)
           lp = jbas%ll(p) 
           lh = jbas%ll(h)
           tp = jbas%itzp(p) 
           th = jbas%itzp(h)
        
           if (.not. (parflip)) then 
              if ( mod(lp + lh,2) .ne. PAR ) cycle
           end if 
           
           if (abs(tp - th)/2 .ne. Tz ) cycle 
        
           NBindx1 = 0
           NBindx2 = 0 
           if ( triangle(jp,jh,Jtot1) )  then 
              
              
              if ( mod(lp + lh,2) ==  PAR ) then   
                 x = CCindex(p,h,OP%Nsp)
                 gnb = 1
              
                 do while (CCME%qmap(x)%Z(gnb) .ne. qONE )
                    gnb = gnb + 1 
                 end do
              
                 NBindx1 = CCME%nbmap(x)%Z(gnb) 
             
              end if 
              
              if  ( triangle(jp,jh,Jtot2) )  then 
              
                 if ( mod(lp+lh+op%dpar/2,2) == PAR) then 
        
                    x = CCindex(p,h,OP%Nsp)
                    gnb = 1
                    
                    do while (CCME%qmap(x)%Z(gnb) .ne. qTWO )
                       gnb = gnb + 1 
                    end do
                 
                    NBindx2 = CCME%nbmap(x)%Z(gnb)
                 end if 
              end if
 
           else if  ( triangle(jp,jh,Jtot2) )  then 
              
              if ( mod(lp+lh+op%dpar/2,2) == PAR) then 
        
                 x = CCindex(p,h,OP%Nsp)
                 gnb = 1
              
                 do while (CCME%qmap(x)%Z(gnb) .ne. qTWO )
                    gnb = gnb + 1 
                 end do
              
                 NBindx2 = CCME%nbmap(x)%Z(gnb)
              end if
           else 
              cycle
           end if 
          
           
           if (phase) pre = (-1)**((jp +jh)/2) !convenient to have this 
           ! for the ph  channel 2body derivative 
        
           do ax = 1, OP%belowEF
              a = jbas%holes(ax)
              do bx = 1, OP%nsp - OP%belowEF
                 b = jbas%parts(bx) 
                 
                 ja = jbas%jj(a) 
                 jb = jbas%jj(b)
                 la = jbas%ll(a) 
                 lb = jbas%ll(b)
                 ta = jbas%itzp(a) 
                 tb = jbas%itzp(b)
                 

                 if (.not. (parflip)) then 
                    if ( mod(la + lb,2) .ne. PAR ) cycle
                 end if 

                 if (abs(ta - tb)/2 .ne. Tz ) cycle 

       
                 if ( (triangle(ja,jb,Jtot1)) .and. (NBindx2 .ne. 0) ) then 
       
                    if ( mod(la+lb,2) == PAR ) then 
                       x = CCindex(a,b,OP%Nsp) 
                       
                       g = 1
                       do while (CCME%qmap(x)%Z(g) .ne. qONE )
                          g = g + 1
                       end do
              
                       Rindx = CCME%rmap(x)%Z(g)
                 
                       if ( (mod(la + lh,2) == mod(lb + lp + op%dpar/2,2)) .and. &
                            ( (ta + th) == (tb + tp) ) ) then  
                         
                          ! hapb 
                          Jmin2 = abs(ja - jh) 
                          Jmax2 = ja+jh 
                          Jmin1 = abs(jp - jb)
                          Jmax1 = jp+jb 
                          
                          sm = 0.d0 
                          do J3 = Jmin1,Jmax1,2
                             do J4 = Jmin2,Jmax2,2
                                sm = sm - (-1)**(J4/2) * sqrt((J3 + 1.d0) * (J4+1.d0)) * &
                                     ninej(jp,jh,Jtot2,jb,ja,Jtot1,J3,J4,rank) * &
                                     tensor_elem(p,b,a,h,J3,J4,OP,jbas) 
                                
                             end do
                          end do
                          
                          ! STORED SUCH THAT THE PH JTOT is GREATER. 
                          CCME%CCX(q1)%X(Rindx,NBindx2) = sm * & 
                               (-1) **( (jh+jb+Jtot1) / 2) * pre * sqrt((Jtot1 + 1.d0)*(Jtot2 + 1.d0))
                          ! stored backwards because. 

                       end if
                    end if
                 end if 
                 
                 if ((triangle(ja,jb,Jtot2)) .and. (NBindx1 .ne. 0) ) then 
                    
                    if ( mod(la+lb+op%dpar/2,2) == PAR ) then 
                       x = CCindex(a,b,OP%Nsp) 
                       g = 1
                       do while (CCME%qmap(x)%Z(g) .ne. qTWO )
                          g = g + 1
                       end do
                    
                       Gindx = CCME%rmap(x)%Z(g)
                 

                       if ( (mod(la + lh,2) == mod(lb + lp + op%dpar/2,2)) .and. &
                            ( (ta + th) == (tb + tp) ) ) then  
               
                          ! p(h)a(b)  
                          Jmin2 = abs(ja - jh) 
                          Jmax2 = ja+jh 
                          Jmin1 = abs(jp - jb)
                          Jmax1 = jp+jb 
                    
                          sm = 0.d0 
                          do J3 = Jmin1,Jmax1,2
                             do J4 = Jmin2,Jmax2,2
                                sm = sm - (-1)**(J4/2) * sqrt((J3 + 1.d0) * (J4+1.d0)) * &
                                     ninej(jp,jh,Jtot1,jb,ja,Jtot2,J3,J4,rank) * &
                                     tensor_elem(p,b,a,h,J3,J4,OP,jbas) 
                       
                             end do
                          end do

                          ! store  ( V )_p(h)a(b)
                          
                          ! PH JTOT IS LESS THAN OR EQUAL TO.
                          CCME%CCR(q1)%X(NBindx1,Gindx) = sm * &
                               (-1) **( (jb+jh+Jtot2)/2) * pre * sqrt((Jtot1 + 1.d0)*(Jtot2 + 1.d0))
                       end if
                    end if
                 end if 
                 
              end do
           end do
        end do
     end do

  end do 
!$omp end parallel do

end subroutine 
!=======================================================  
!=======================================================          
subroutine EOM_scalar_cross_coupled(HS,CCME,jbas,phase) 
  ! currently the only CCME of interest are phab terms    |---<---| 
  ! coupling in the 3-1 channel                        <(pa)J|V|(hb)J>
  !                                                      |---<---|
  implicit none 
  
  type(spd) :: jbas
  type(sq_op) :: HS
  type(cross_coupled_31_mat) :: CCME
  integer :: JT,ja,jp,jb,jh,JC,q1,q2,TZ,PAR,la,lb,Ntot,th,tp,lh,lp
  integer :: a,b,p,h,i,j,Jmin,Jmax,Rindx,Gindx,g,ta,tb,Atot,hg,pg
  integer :: int1,int2,IX,JX,i1,i2,nb,nh,np,gnb,NBindx,x,JTM
  real(8) :: sm,sm2,pre,horse
  logical :: phase

  Atot = HS%belowEF
  Ntot = HS%Nsp
  JTM = jbas%Jtotal_max 
  pre = 1.d0 

!$omp parallel do default(firstprivate),shared(CCME,HS,jbas) 
  do q1 = 1, CCME%nblocks
      
     JC = mod(q1-1,JTM+1) * 2 
     PAR = (q1 - 1) / (2*JTM + 2) 
     TZ = mod(q1-1,(2*JTM+2))/(JTM+1)  
     
     CCME%CCR(q1)%X = 0.d0
     CCME%CCX(q1)%X = 0.d0
         
     ! ab = ph 
     do hg = 1, Atot
        do pg = 1, Ntot - Atot 
           
           h = jbas%holes(hg) 
           p = jbas%parts(pg) 
           
           jp = jbas%jj(p) 
           jh = jbas%jj(h)
           lp = jbas%ll(p) 
           lh = jbas%ll(h)
           tp = jbas%itzp(p) 
           th = jbas%itzp(h)
        
           if (.not. triangle(jp,jh,JC) )  cycle
           if ( mod(lp + lh,2) .ne. PAR ) cycle
           if (abs(tp - th)/2 .ne. Tz ) cycle 
           
           x = CCindex(p,h,HS%Nsp)
           gnb = 1
           do while (CCME%qmap(x)%Z(gnb) .ne. q1 )
              gnb = gnb + 1
           end do
              
           NBindx = CCME%nbmap(x)%Z(gnb) 

           if (phase) pre = (-1)**((jp +jh)/2) !convenient to have this 
           ! for the ph  channel 2body derivative 
        
           do a = 1, HS%nsp
              if ( jbas%con(a) .ne. 1 )  cycle

              do b = 1, HS%nsp
                 if (jbas%con(b) .ne. 0) cycle

                 ja = jbas%jj(a) 
                 jb = jbas%jj(b)
                 la = jbas%ll(a) 
                 lb = jbas%ll(b)
                 ta = jbas%itzp(a) 
                 tb = jbas%itzp(b)
                 
                 if (.not. triangle(ja,jb,JC) )  cycle
                 if ( mod(la + lb,2) .ne. PAR ) cycle
                 if (abs(ta - tb)/2 .ne. Tz ) cycle 
                 
                 x = CCindex(a,b,HS%Nsp) 
                 g = 1
                 do while (CCME%qmap(x)%Z(g) .ne. q1 )
                    g = g + 1
                 end do
              
                 Rindx = CCME%rmap(x)%Z(g)
                 
                 x = CCindex(b,a,HS%Nsp) 
                 g = 1
                 do while (CCME%qmap(x)%Z(g) .ne. q1 )
                    g = g + 1
                 end do
              
                 Gindx = CCME%rmap(x)%Z(g)
                 
                 sm = 0.d0 
               
!                 horse = 0.d0 
                 if ( (mod(la + lh,2) == mod(lb + lp,2)) .and. &
                      ( (ta + th) == (tb + tp) ) ) then  
               
                    ! hapb 
                    Jmin = max(abs(jp - jb),abs(ja - jh)) 
                    Jmax = min(jp+jb,ja+jh) 
                    
                    sm = 0.d0 
                    do JT = Jmin,Jmax,2
                       sm = sm + (-1)**(JT/2) * (JT + 1) * &
                            sixj(jp,jh,JC,ja,jb,JT)  * &
                            v_elem(h,a,p,b,JT,HS,jbas) 
                    end do
                 
                    ! store  < p b | v | h a> 
                    ! scaled by sqrt(JC + 1) for convience in ph derivative
                             
                    CCME%CCR(q1)%X(NBindx,Gindx) = sm * HS%herm * &
                         (-1) **( (jp + ja + JC) / 2) * pre * sqrt(JC + 1.d0)
                 
                 end if
              end do
           end do
        end do
     end do

  end do 
!$omp end parallel do

end subroutine 
!===========================================================
!===========================================================     
integer function CCindex(a,b,N)
  implicit none 
  
  integer :: a,b,N 
  
  CCindex = N*(a-1) + b
end function 
!===========================================================
!===========================================================  
subroutine allocate_CC_wkspc(CCHS,WCC)
  implicit none 
  
  type(cross_coupled_31_mat) :: CCHS,WCC 
  integer :: q,r
  
  allocate(WCC%CCX(CCHS%nblocks))
  allocate(WCC%CCR(CCHS%nblocks))
  
  do q = 1,CCHS%nblocks
     
     r = CCHS%rlen(q) 
     
     allocate(WCC%CCX(q)%X(r,r)) 
     allocate(WCC%CCR(q)%X(r,r)) 
     WCC%CCX(q)%X = 0.d0
     WCC%CCR(q)%X = 0.d0
     
  end do
end subroutine 
!===========================================================
!===========================================================  
subroutine allocate_CCtensor_wkspc(CCOP,WCC)
  implicit none 
  
  type(cross_coupled_31_mat) :: CCOP,WCC 
  integer :: q,r1,r2
  
  allocate(WCC%CCX(CCOP%nblocks))
  allocate(WCC%CCR(CCOP%nblocks))
  
  do q = 1,CCOP%nblocks
     
     r1 = size(CCOP%CCX(q)%X(:,1)) 
     r2 = size(CCOP%CCR(q)%X(1,:))
     
     allocate(WCC%CCX(q)%X(r1,r2)) 
     allocate(WCC%CCR(q)%X(r1,r2)) 
     WCC%CCX(q)%X = 0.d0
     WCC%CCR(q)%X = 0.d0
     
  end do
end subroutine 
!===========================================================
!===========================================================
function optimum_omega_for_CM_hamiltonian(hw,Ew) result(epm) 
! G. Hagen, T. Papenbrock, D.J. Dean (2009)  COM diagonistic 
  real(8),intent(in) :: hw,Ew 
  real(8),dimension(2) :: epm 
  
  if (Ew < 0.d0) then 
     epm = hw 
     return
  end if 
  
  epm(1) = hw + 2.d0/3.d0*Ew + sqrt( 4.d0/9.d0*Ew*Ew + 4.d0/3.d0*hw*Ew )   
  epm(2) = hw + 2.d0/3.d0*Ew - sqrt( 4.d0/9.d0*Ew*Ew + 4.d0/3.d0*hw*Ew )
end function  
!=================================================
!=================================================
subroutine vectorize(rec,vout)
  !!! maps full_ham to vector
  implicit none 

  integer ::  i,j,k,l,gx,Atot,Ntot
  type(sq_op) :: rec
  real(8),dimension(rec%neq) :: vout
  
 
  Atot = rec%belowEF
  Ntot= rec%Nsp
  
  Ntot = Ntot - Atot 
  
  vout(1) = rec%E0
  
  k=2
  do i=1,Atot
     do j=1,Atot
        
        vout(k) =  rec%fhh(j,i) 
        k=k+1
       
     end do 
  end do 

  do i=1,Ntot
     do j=1,Ntot
        
        vout(k) = rec%fpp(j,i) 
        k=k+1
        
     end do 
  end do 

  do i=1,Atot
     do j=1,Ntot
        
        vout(k) = rec%fph(j,i) 
        k=k+1

     end do 
  end do 
  
  do l=1,rec%nblocks
    
     do gx = 1, 6
        
        do i=1,size(rec%mat(l)%gam(gx)%X(1,:) )
           do j=min(jst(gx),i),size(rec%mat(l)%gam(gx)%X(:,1))
              
              vout(k) = rec%mat(l)%gam(gx)%X(j,i) 
            
              k=k+1
              
           end do
        end do
     end do 
  end do    

end subroutine 
!=================================================
!=================================================
subroutine repackage(rec,vout)
  !!! maps full_ham to vector
  implicit none 

  integer :: i,j,k,l,gx,Atot,Ntot
  type(sq_op) :: rec
  real(8),dimension(rec%neq) :: vout
  

  Atot = rec%belowEF
  Ntot= rec%Nsp
  
  Ntot = Ntot - Atot 
  
  rec%E0 = vout(1)
  
  k=2
  do i=1,Atot
     do j=1,Atot
        
        rec%fhh(j,i) = vout(k) 
        k=k+1
       
     end do 
  end do 

  do i=1,Ntot
     do j=1,Ntot
        
        rec%fpp(j,i) = vout(k)
        k=k+1
        

     end do 
  end do 

  do i=1,Atot
     do j=1,Ntot
        
        rec%fph(j,i) = vout(k) 
        k=k+1

     end do 
  end do 
    
  do l=1,rec%nblocks
     
     do gx = 1, 6
       
        do i=1,size(rec%mat(l)%gam(gx)%X(1,:) )
           do j=min(jst(gx),i),size(rec%mat(l)%gam(gx)%X(:,1))
           
              rec%mat(l)%gam(gx)%X(j,i) = vout(k) 
        if (sqs(gx)) rec%mat(l)%gam(gx)%X(i,j) = vout(k) * rec%herm 
              k=k+1
        
           end do
        end do
     end do 
  end do    
   
end subroutine 
!===============================================================
!===============================================================  
real(8) function mat_frob_norm(op) 
  implicit none 
  
  type(sq_op) :: op
  integer :: q,g
  real(8) :: sm

  sm = sum(op%fhh**2) +  sum(op%fph**2) + sum(op%fpp**2)
  do q = 1, op%nblocks
     do g = 1,6
        sm = sm + sum(op%mat(q)%gam(g)%X**2)
     end do 
  end do 
  
  mat_frob_norm = sqrt(sm)
end function 
!===============================================================
!===============================================================
subroutine write_excited_states(steps,s,TDA,e0,un) 
  implicit none
  
  integer :: steps,sm,q,r,un
  real(8) :: s,e0
  type(full_sp_block_mat) :: TDA
  real(8),allocatable,dimension(:) :: vec
  character(5) :: num 
  
  
  sm = 0
  do q = 1,TDA%blocks
     sm = sm + TDA%map(q)
  end do 
  
  allocate(vec(sm)) 
  
  r = 1
  do q = 1,TDA%blocks 
     if ( TDA%map(q) > 0 ) then 
        vec(r:r+TDA%map(q)-1) = TDA%blkM(q)%eigval !+e0
        r = r + TDA%map(q)
     end if
  end do 
  
  sm = sm + 1
  write(num,'(I5)') sm 
  num = adjustl(num) 
  
  write(un,'(I6,'//trim(num)//'(e14.6))') steps,s,vec

end subroutine 
!=====================================
real(8) function mbpt2(H,jbas) 
  implicit none 
  
  integer :: i,j,k,l,II,JJ,q
  type(sq_op) :: H 
  type(spd) :: jbas
  real(8) :: sm , sm_singlej, eden,fi,fj,fk,fl
  
  sm = 0.d0 
  do q = 1, H%nblocks
     
     sm_singlej = 0.d0
     do II = 1,H%mat(q)%npp
        i = H%mat(q)%qn(1)%Y(II,1)
        j = H%mat(q)%qn(1)%Y(II,2)
        fi = f_elem(i,i,H,jbas)
        fj = f_elem(j,j,H,jbas)
        
        do JJ = 1,H%mat(q)%nhh            
           k = H%mat(q)%qn(3)%Y(JJ,1)
           l = H%mat(q)%qn(3)%Y(JJ,2)
           fk = f_elem(k,k,H,jbas)
           fl = f_elem(l,l,H,jbas)
           
           eden = fk + fl - fi - fj 
          
           sm_singlej = sm_singlej + H%mat(q)%gam(3)%X(II,JJ)**2/eden
        end do 
     end do 
     
     sm = sm + (H%mat(q)%lam(1)+1.d0)*sm_singlej
  end do 

  mbpt2 = sm
end function           
!===================================================================
!===================================================================  
subroutine write_tilde_from_Rcm(rr) 
  implicit none 
  
  type(sq_op) :: rr 
  character(200) :: spfile,intfile,input,prefix
  common /files/ spfile,intfile,prefix 

  open(unit=53, file=OUTPUT_DIR//&
       trim(adjustl(prefix))//'_omegatilde.dat') 
  write(53,'(2(e14.6))') rr%hospace,&
       hbarc2_over_mc2*1.5d0/float(rr%Aprot+rr%Aneut)/rr%E0 
  close(53) 
end subroutine 
!===================================================================
!===================================================================  
subroutine enumerate_three_body(threebas,jbas) 
  implicit none 
  
  type(spd) :: jbas
  type(tpd),allocatable,dimension(:) :: threebas
  integer :: a,b,c,i,j,k,JTot, PAR, Tz,q
  integer :: Jmaxx,blocks,holes,parts
  integer :: ix,jx,kx,ax,bx,cx,Jab_max,Jab_min
  integer :: ja,jb,jc,ji,jj,jk,Jij_max,Jij_min 
  integer :: ta,tb,tc,ti,tj,tk,num_hhh,Jab
  integer :: la,lb,lc,li,lj,lk,num_ppp,Jij
  
  holes = sum(jbas%con) 
  parts = size(jbas%con) - holes
  
  Jmaxx = maxval(jbas%jj)*3
  
  blocks= (Jmaxx + 1)/2 * 8 
  allocate(threebas(blocks))
  q = 1
  do Jtot = 1,Jmaxx,2 ! Jtot is odd 
     
     do Tz = -3,3,2 
        
        do PAR = 0, 1
       
        threebas(q)%chan(1) = Jtot
        threebas(q)%chan(2) = Tz
        threebas(q)%chan(3) = PAR 
        
   num_hhh = 0
! sum over all possible un-symmetrized hhh combinations 
   do ix = 1,holes
      i =jbas%holes(ix)
      ji = jbas%jj(i)
      li = jbas%ll(i)      
      ti = jbas%itzp(i)
      
      do jx = 1,holes
         j =jbas%holes(jx)
         jj = jbas%jj(j)
         lj = jbas%ll(j)      
         tj = jbas%itzp(j)
                  
         Jij_min = abs(ji - jj) 
         Jij_max = ji + jj 
         
         do kx= 1,holes
            
            k =jbas%holes(kx)
            jk = jbas%jj(k)
            lk = jbas%ll(k)      
            tk = jbas%itzp(k)
            
            if (ti + tj + tk .ne. Tz) cycle
            if (mod(li + lj + lk,2) .ne. PAR)  cycle
            do Jij = Jij_min, Jij_max 
               if (triangle(Jk,Jtot,Jij)) then 
                  exit
               end if
            end do 
            
            if (Jij > Jij_max) cycle 
            
            num_hhh = num_hhh + 1 
        end do 
     end do 
   end do 
           
   allocate( threebas(q)%hhh(num_hhh,3)) 
   
   num_hhh = 0
! sum over all possible un-symmetrized hhh combinations 
   do ix = 1,holes
      i =jbas%holes(ix)
      ji = jbas%jj(i)
      li = jbas%ll(i)      
      ti = jbas%itzp(i)
      
      do jx = 1,holes
         j =jbas%holes(jx)
         jj = jbas%jj(j)
         lj = jbas%ll(j)      
         tj = jbas%itzp(j)
                  
         Jij_min = abs(ji - jj) 
         Jij_max = ji + jj 
         
         do kx= 1,holes
            
            k =jbas%holes(kx)
            jk = jbas%jj(k)
            lk = jbas%ll(k)      
            tk = jbas%itzp(k)
           
            if (ti + tj + tk .ne. Tz) cycle
            if (mod(li + lj + lk,2) .ne. PAR)  cycle
           
            do Jij = Jij_min, Jij_max 
               if (triangle(Jk,Jtot,Jij)) then 
                  exit
               end if
            end do 
            
            if (Jij > Jij_max) cycle 
           
            num_hhh = num_hhh + 1
            threebas(q)%hhh(num_hhh,1) = i 
            threebas(q)%hhh(num_hhh,2) = j 
            threebas(q)%hhh(num_hhh,3) = k 

        end do 
     end do 
   end do      

  
   num_ppp = 0
! run over all possible un-symmetrized ppp combinations 
   do ax = 1,parts
      a =jbas%parts(ax)
      ja = jbas%jj(a)
      la = jbas%ll(a)      
      ta = jbas%itzp(a)
      
      do bx = 1,parts
         b =jbas%parts(bx)
         jb = jbas%jj(b)
         lb = jbas%ll(b)      
         tb = jbas%itzp(b)
                  
         Jab_min = abs(ja - jb) 
         Jab_max = ja + jb 
         
         do cx= 1,parts
            
            c =jbas%parts(cx)
            jc = jbas%jj(c)
            lc = jbas%ll(c)      
            tc = jbas%itzp(c)
            
            
            if (ta + tb + tc .ne. Tz) cycle
            if (mod(la + lb + lc,2) .ne. PAR)  cycle
            do Jab = Jab_min, Jab_max 
               if (triangle(Jc,Jtot,Jab)) then 
                  exit
               end if
            end do 
            
            if (Jab > Jab_max) cycle 
            
            num_ppp = num_ppp + 1 
        end do 
     end do 
   end do 

   allocate( threebas(q)%ppp(num_ppp,3)) 
  
   num_ppp = 0
! run over all possible un-symmetrized ppp combinations 
   do ax = 1,parts
      a =jbas%parts(ax)
      ja = jbas%jj(a)
      la = jbas%ll(a)      
      ta = jbas%itzp(a)
      
      do bx = 1,parts
         b =jbas%parts(bx)
         jb = jbas%jj(b)
         lb = jbas%ll(b)      
         tb = jbas%itzp(b)
                  
         Jab_min = abs(ja - jb) 
         Jab_max = ja + jb 
         
         do cx= 1,parts
            
            c =jbas%parts(cx)
            jc = jbas%jj(c)
            lc = jbas%ll(c)      
            tc = jbas%itzp(c)
            
            if (ta + tb + tc .ne. Tz) cycle
            if (mod(la + lb + lc,2) .ne. PAR)  cycle
              
            do Jab = Jab_min, Jab_max 
               if (triangle(Jc,Jtot,Jab)) then 
                  exit
               end if
            end do 
            
            if (Jab > Jab_max) cycle 
            
            num_ppp = num_ppp + 1 
            threebas(q)%ppp(num_ppp,1) = a 
            threebas(q)%ppp(num_ppp,2) = b 
            threebas(q)%ppp(num_ppp,3) = c
 
        end do 
     end do 
   end do 
   
            q = q + 1 
        end do
     end do
   end do 

   call divide_work_tpd(threebas) 
end subroutine 

real(8) function ninej(a,b,J1,c,d,J2,J3,J4,RANK) 
  ! fast ninej using stored 6j
  implicit none 

  integer ::a ,b, c,d,J1,J2,J3,J4,RANK,x,xmin,xmax
  real(8) :: sm 

  xmin = max(abs(J2-b),abs(a-rank),abs(J4-c))  
  xmax = min( J2+b , a+rank ,J4+c )
  
  sm = 0.d0 
  do x = xmin,xmax,2
     sm = sm - (x+1.d0) * xxxsixj(J1,J2,rank,x,a,b) * &
          sixj(c,d,J2,b,x,J4)* xxxsixj(J3,J4,rank,x,a,c) 
  end do 
  
  ninej = sm 
end function


real(8) function p1_p2( a, b, c, d, J ,jbas ) 
  !This is 1/b^2 times del1_del2
  implicit none 
  
  type(spd) :: jbas
  integer :: a,b,c,d,J,na,nb,nc,nd,ja,jb,jc,jd
  integer :: la,lb,lc,ld,ta,tb,tc,td
  real(8) :: d6ji, j_ac,j_bd,j_ad,j_bc
  
  ja = jbas%jj(a)
  jb = jbas%jj(b)
  jc = jbas%jj(c)
  jd = jbas%jj(d)  
  
  la = jbas%ll(a)
  lb = jbas%ll(b)
  lc = jbas%ll(c)
  ld = jbas%ll(d)  

  na = jbas%nn(a)
  nb = jbas%nn(b)
  nc = jbas%nn(c)
  nd = jbas%nn(d)  

  ta = jbas%itzp(a)
  tb = jbas%itzp(b)
  tc = jbas%itzp(c)
  td = jbas%itzp(d)  

  j_ac = 0.d0
  j_bd = 0.d0 
  j_bc = 0.d0 
  j_ad = 0.d0 
  
  if( ta == tc) then 
     j_ac = sqrt( (ja + 1.d0)*(jc+1.d0) ) * (-1) ** ((3+2*la +jc)/2) * &
          d6ji( ja, 2 , jc, 2*lc , 1 , 2*la ) * reduced_p( a, c, jbas) 
  end if 
  
  if (ta == td) then 
     j_ad = sqrt( (ja + 1.d0)*(jd+1.d0) ) * (-1) ** ((3+2*la +jd)/2) * &
          d6ji( ja, 2 , jd, 2*ld , 1 , 2*la ) * reduced_p( a, d, jbas) 
  end if 
  
  if (tb == tc) then
     j_bc = sqrt( (jb + 1.d0)*(jc+1.d0) ) * (-1) ** ((3+2*lb +jc)/2) * &
          d6ji( jb, 2 , jc, 2*lc , 1 , 2*lb ) * reduced_p( b, c, jbas) 
  end if 
  
  if (tb == td) then 
     j_bd = sqrt( (jb + 1.d0)*(jd+1.d0) ) * (-1) ** ((3+2*lb +jd)/2) * &
          d6ji( jb, 2 , jd, 2*ld , 1 , 2*lb ) * reduced_p( b, d, jbas) 
  end if 
  
  p1_p2 = ((-1) ** ((jb - jc - J ) /2 ) * d6ji(ja,jb,J,jd,jc,2) * &
       j_ac * j_bd  &
       -(-1)**((jc +jd-J)/2)* (-1) ** ((jb - jd - J ) /2 ) * d6ji(ja,jb,J,jc,jd,2) * &
       j_ad * j_bc) / sqrt( (1.d0+kron_del(a,b))*(1.d0+kron_del(c,d))) 

end function
  
real(8) function r1_r2( a, b, c, d, J ,jbas ) 
  !This is b^2 times r1_r2
  implicit none 
  
  type(spd) :: jbas
  integer :: a,b,c,d,J,na,nb,nc,nd,ja,jb,jc,jd
  integer :: la,lb,lc,ld,ta,tb,tc,td
  real(8) :: d6ji, j_ac,j_bd,j_ad,j_bc
  
  ja = jbas%jj(a)
  jb = jbas%jj(b)
  jc = jbas%jj(c)
  jd = jbas%jj(d)  
  
  la = jbas%ll(a)
  lb = jbas%ll(b)
  lc = jbas%ll(c)
  ld = jbas%ll(d)  

  na = jbas%nn(a)
  nb = jbas%nn(b)
  nc = jbas%nn(c)
  nd = jbas%nn(d)  

  ta = jbas%itzp(a)
  tb = jbas%itzp(b)
  tc = jbas%itzp(c)
  td = jbas%itzp(d)  

  j_ac = 0.d0
  j_bd = 0.d0 
  j_bc = 0.d0 
  j_ad = 0.d0 
  
  if( ta == tc) then 
     j_ac = sqrt( (ja + 1.d0)*(jc+1.d0) ) * (-1) ** ((3+2*la +jc)/2) * &
          d6ji( ja, 2 , jc, 2*lc , 1 , 2*la ) * reduced_r( a, c, jbas) 
  end if 
  
  if (ta == td) then 
     j_ad = sqrt( (ja + 1.d0)*(jd+1.d0) ) * (-1) ** ((3+2*la +jd)/2) * &
          d6ji( ja, 2 , jd, 2*ld , 1 , 2*la ) * reduced_r( a, d, jbas) 
  end if 
  
  if (tb == tc) then
     j_bc = sqrt( (jb + 1.d0)*(jc+1.d0) ) * (-1) ** ((3+2*lb +jc)/2) * &
          d6ji( jb, 2 , jc, 2*lc , 1 , 2*lb ) * reduced_r( b, c, jbas) 
  end if 
  
  if (tb == td) then 
     j_bd = sqrt( (jb + 1.d0)*(jd+1.d0) ) * (-1) ** ((3+2*lb +jd)/2) * &
          d6ji( jb, 2 , jd, 2*ld , 1 , 2*lb ) * reduced_r( b, d, jbas) 
  end if 
  
  r1_r2 = ((-1) ** ((jb + jc - J ) /2 ) * d6ji(ja,jb,J,jd,jc,2) * &
       j_ac * j_bd  &
       -(-1)**((jc +jd-J)/2)* (-1) ** ((jb + jd - J ) /2 ) * d6ji(ja,jb,J,jc,jd,2) * &
       j_ad * j_bc) / sqrt( (1.d0+kron_del(a,b))*(1.d0+kron_del(c,d))) 

end function

real(8) function reduced_p( a, b , jbas ) 
  ! < n l || vec{r} || n' l' > 
  implicit none 
  
  type(spd) :: jbas
  integer :: a,b,la,lb,na,nb
  real(8) :: c1, p1 
  
  la = jbas%ll(a) 
  lb = jbas%ll(b) 
  
  if ( la == lb+1) then 
     c1 = sqrt(lb +1.d0)
     
     na = jbas%nn(a) 
     nb = jbas%nn(b)
     if (na == nb) then         
        p1 = -1 * sqrt( nb + lb + 1.5d0) 
     else if (na == nb - 1) then 
        p1 = -1 * sqrt( nb*1.d0 )
     else 
        p1 = 0.d0 
     end if 
        
  else if ( la == lb-1) then 
     c1 = -1*sqrt(1.d0*lb)

     na = jbas%nn(a) 
     nb = jbas%nn(b)
     if (na == nb) then         
        p1 =  sqrt( nb + lb + 0.5d0) 
     else if (na == nb + 1) then 
        p1 =  sqrt( nb + 1.d0 )
     else 
        p1 = 0.d0 
     end if 

  else
     c1 = 0.d0 
     p1 = 0.d0
  end if 
  
  reduced_p = p1*c1
  
end function

real(8) function reduced_r( a, b , jbas ) 
  ! < n l || vec{r} || n' l' > 
  implicit none 
  
  type(spd) :: jbas
  integer :: a,b,la,lb,na,nb
  real(8) :: c1, p1 
  
  la = jbas%ll(a) 
  lb = jbas%ll(b) 
  
  if ( la == lb+1) then 
     c1 = sqrt(lb +1.d0)
     
     na = jbas%nn(a) 
     nb = jbas%nn(b)
     if (na == nb) then         
        p1 = sqrt( nb + lb + 1.5d0) 
     else if (na == nb - 1) then 
        p1 = -1 * sqrt( nb*1.d0 )
     else 
        p1 = 0.d0 
     end if 
        
  else if ( la == lb-1) then 
     c1 = -1*sqrt(1.d0*lb)

     na = jbas%nn(a) 
     nb = jbas%nn(b)
     if (na == nb) then         
        p1 =  sqrt( nb + lb + 0.5d0) 
     else if (na == nb + 1) then 
        p1 =  -1*sqrt( nb + 1.d0 )
     else 
        p1 = 0.d0 
     end if 

  else
     c1 = 0.d0 
     p1 = 0.d0
  end if 
  
  reduced_r = p1*c1
  
end function
  
end module       



