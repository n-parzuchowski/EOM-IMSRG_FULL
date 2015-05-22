
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       subroutine setupAB 
C======================================================================
C Calculates the RPA matrices and solves the RPA equations
C======================================================================
      use RPAmatrices
      use wfn
      use spspace
      use nuclide
      use HAMILTONIAN
      use phBASIS
      use flags
      implicit none
C======================================================================
C
C INPUT:
C
C OUTPUT:
C
C Subroutines called:
C         make_Gamma:         calculates Gamma matrix
C         make_hmatrix:       calculates hpp,hnn
C         eig/eigval:         calculates eigenvalues and/ eigenvectors of
C                             a symmetric matrix
C         eigsrt:             sorts the eigenvalues
C         map_ph:             maps the ph states in pairs and calculates the number
C                             of pairs
C         setAandB:           calculates diagonal and off-diagonal RPA matrices
C                             for protons and neutrons
C         get_AB:             adds to A and B the pn parts

C Called by: the main routine
C======================================================================
C      include 'gcb_dim.inc'

C-------INPUT-----------------
      real gammap(nsps(1),nsps(1)),gamman(nsps(2),nsps(2))
      real hpp(nsps(1),nsps(1)),hnn(nsps(2),nsps(2))

C===========INTERNAL VARIABLES======================================

C------ph Basis--------------------
      real,allocatable  :: work(:)         ! dummy working array for
                                           ! diagonalization
      real,dimension(nsps(1),nsps(1))    :: array1,array2
C      integer,allocatable :: occ(:),uocc(:)  ! occupied/unoccupied states


C-----Temporary variables
      integer i,j,k,m

C
!      call make_projectors(psd,rhop,qprot,np,nsps(1),overp)
!      call make_projectors(nsd,rhon,qneutr,nn,nsps(2),overn)

      call make_Gamma(gammap,gamman)
      call make_hmatrix(gammap,e_spepp,hpp,nsps(1))
      call make_hmatrix(gamman,e_spenn,hnn,nsps(2))

      nphp=np*(nsps(1)-np)
      nphn=nn*(nsps(2)-nn)

c        call hmult(nsps,nmatpp,nmatnn,nmatpn,
c     &  map_combpp,map_combnn,map_combpn,
c     &  hvecpp,hvecnn,hvecpn,e_spepp,e_spenn,
c     &  overp,overn,rhop,rhon,over,hme)
c      print*,Ehf,hme


      if(allocated(vecp))deallocate(vecp)
      if(allocated(vecn))deallocate(vecn)
      if(allocated(energp))deallocate(energp)
      if(allocated(energn))deallocate(energn)
      if(allocated(map_prot))deallocate(map_prot)
      if(allocated(map_neutr))deallocate(map_neutr)

      allocate(vecp(nsps(1),nsps(1)),vecn(nsps(2),nsps(2)))
      allocate(energp(nsps(1)),energn(nsps(2)))
      allocate(map_prot(nphp,2),map_neutr(nphn,2))


C      write(12,*)nsps
C.........Calculate ph basis.....
      allocate(work(nsps(1)))
      call eig(hpp,nsps(1),nsps(1),energp,vecp,work)     ! for protons
C      do k=1,nsps(1)
C         write(12,'(20(F10.5,1x))')(vecp(i,k),i=1,nsps(1))
C      enddo
      call eigsrt(energp,vecp,nsps(1),nsps(1))
      deallocate(work)
C      write(10,*)' '
C      write(10,'(12F12.6)')energp

C      do i=1,nsps(1)
C         do j=1,nsps(1)
C           array2(i,j)=rhop(i,j)
C         enddo
C      enddo
C      array1=matmul(transpose(vecp),array2)
C      array2=matmul(array1,vecp)

C      write(6,'(12F9.5)')array2
C      write(6,*)' '
C      write(6,'(12F9.5)')energp

C      stop




c      if(noRPA)return

      allocate(work(nsps(2)))
      call eig(hnn,nsps(2),nsps(2),energn,vecn,work)        ! for neutrons
      call eigsrt(energn,vecn,nsps(2),nsps(2))
      deallocate(work)

      open(unit=68,file='hf.spe',status='unknown')
      write(68,*)np,nn
      do i = 1,nsps(1)
        write(68,*)i,energp(i)
      enddo
      do i = 1,nsps(2)
        write(68,*)i,energn(i)
      enddo
      close(unit=68)

!.... ADDED IN 2014.... PRINTS OUT HF ENERGIES AND SINGLE-PARTICLE BASIS....

      open(unit=78,file='hf.basis',status='unknown')
!...... WRITE OUT SINGLE PARTICLE QUANTUM NUMBERS
      write(78,*)nsps(1)
      do i = 1,nsps(1)
            write(78,'(i3,2x,5i3)')i,(spsqn(1,i,j),j=1,5)
      end do

      write(78,*)np,nn
!....... WRITE OUT SINGLE PARTICLE ENERGIES
      do i = 1,nsps(1)
        write(78,*)i,energp(i)
      enddo
      do i = 1,nsps(2)
        write(78,*)i,energn(i)
      enddo
!........ WRITE OUT ASSOCIATED HF STATES.....
      do i = 1,nsps(1)
        write(78,*)i,energp(i)
        write(78,'(7f10.6)')(vecp(j,i),j=1,nsps(1))
      enddo      
      do i = 1,nsps(2)
        write(78,*)i,energn(i)
        write(78,'(7f10.6)')(vecn(j,i),j=1,nsps(2))
      enddo    
      close(unit=78)

      print*,' HF energies written to file hf.spe '
      print*,' full HF basis written to file hf.basis '
      print*,' '
!-----------------CORRECTION
!                 FIXES TR JUMPING BACK AND FORTH
!                CWJ June 2009 (@ UW INT )
!
      do i=1,np                                ! Fill SD ---- protons
         do m=1,nsps(1)
            psd(m,i)=vecp(m,i)
	 enddo
      enddo 

      do i=1,nn                                ! Fill SD ---- neutrons
         do m=1,nsps(2)
	    nsd(m,i)=vecn(m,i)
	 enddo
      enddo
      call make_projectors(psd,rhop,qprot,np,nsps(1),overp)
      call make_projectors(nsd,rhon,qneutr,nn,nsps(2),overn)
      call commuter(nsps(1),rhop,hpp)
      call commuter(nsps(2),rhon,hnn)
c      do k=1,nsps(1)
c         write(12,'(20(F10.5,1x))')(vecn(i,k),i=1,nsps(2))
c      enddo
      call map_ph(nsps(1),np,rhop,vecp,map_prot,nphp)

      call map_ph(nsps(2),nn,rhon,vecn,map_neutr,nphn)
      if(noRPA)return

C.... Start setting up the RPA matrices A and B
C     subroutine setAandB(nph,map,hvec,map_comb,nmat,energ,vec,ns,a,b)

      nph=nphp+nphn
      if(allocated(App))deallocate(App)
      if(allocated(Ann))deallocate(Ann)
      if(allocated(Apn))deallocate(Apn)
      if(allocated(Bpp))deallocate(Bpp)
      if(allocated(Bnn))deallocate(Bnn)
      if(allocated(Bpn))deallocate(Bpn)
      if(allocated(A))deallocate(A)
      if(allocated(B))deallocate(B)

      allocate(App(nphp,nphp),Bpp(nphp,nphp))
      allocate(Ann(nphn,nphn),Bnn(nphn,nphn))
      allocate(Apn(nphp,nphn),Bpn(nphp,nphn))
      allocate(A(nph,nph),B(nph,nph))
      call setAandB(nphp,map_prot,hvecpp,map_combpp,nmatpp,
     &                energp,vecp,nsps(1),app,bpp)                    ! PROTONS
      call setAandB(nphn,map_neutr,hvecnn,map_combnn,nmatnn,
     &                energn,vecn,nsps(2),ann,bnn)                    ! NEUTRONS

      call setAandB_pn(nphp,nphn,map_prot,map_neutr,hvecpn,map_combpn,
     &    nmatpn,vecp,vecn,nsps,apn,bpn)        ! PROTON-NEUTRON INTERACTION

C......SET UP THE full A & B MATRICES...................................
      call get_AorB(app,ann,apn,nphp,nphn,a)
      call get_AorB(bpp,bnn,bpn,nphp,nphn,b)


      nph = nphp + nphn

C      if(2*nph.gt.max_sps4)then
C         write(6,*)2*nph,' > ', max_sps4, '. Icrease max_sps4'
C         stop
C      endif
      print*,' A and B matrices formed '
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine solve_stability(nvalues,errRPA)
C
C  solves the stability matrix and checks to make sure okay
C  
C  THIS VERSION COMPARES DIRECTLY REAL, IMAGINARY PARTS
C  note: because this program uses real wfns, the stability 
C  matrix becomes A+B (real p-h amplitudes) and A-B (imaginary 
C  p-h amplitudes) 
C  
C INPUT:
C
C OUTPUT:
C    stabe:     eigenvalues of stability matrix A+B, A-B
C    nvalues:   # of eigenvalues = 2x # of p-h states
C    errRPA:    flag for err in stability matrix
C
      use RPAmatrices
      use RPAstab
      use flags
      implicit none

      integer,intent(OUT)      :: nvalues

C------Stability Matrix-----------------------------------------
      real stabr(nph),stabi(nph),stab(nph,nph)
      real stabv(nph,nph)               ! stability eigenvectors not necess.
      real stabw(nph)                   ! dummy working array for
                                        ! diagonalization
C..............MISC..........................
      integer i,j
      logical errRPA



      if(allocated(stabe))deallocate(stabe)
      allocate(stabe(2*nph))
C........Test if the stability matrix is positive......................
      errRPA = .false.

C............first REAL AMPLITUDES.............
C      allocate(stab(nph,nph))
      stab=a+b
      nvalues=2*nph

      call eigval(stab,nph,nph,stabr,stabv,stabw)
      call eigsrt(stabr,stabv,nph,nph)

      
      
      do i=1,nph
        stabe(i) = stabr(i)
        if(stabr(i).lt.-0.1)then                
           errRPA=.true.
        endif
      enddo

C..................IMAGINARY AMPLITUDES..................      

      stab=a-b

      call eigval(stab,nph,nph,stabi,stabv,stabw)
      call eigsrt(stabi,stabv,nph,nph)

      do i=1,nph
        stabe(i+nph) = stabi(i)
        if(stabi(i).lt.-0.1)then                
           errRPA=.true.
           print*,i,stabi(i)
        endif
      enddo

      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine solve_RPA
C======================================================================
C Calculates the RPA matrices and solves the RPA equations
C======================================================================
      use wfn
      use nuclide
      use spspace
      use RPAsolutions
      use RPAmatrices
      implicit none
C======================================================================
C
C INPUT:
C
C OUTPUT:
C
C Subroutines called:
C
C         get_rpa:            calculates the full RPA matrix, i.e.
C                              / A   B\
C                              \-B  -A/
C         RPA_solver:         solves the RPA equations by reducing the problem
C                             to a symmetric matrix (DOESN'T WORK YET)
C         rg:                 calculates the eigenvalues and eigenvectors
C                             of a general matrix
C         RPA_eigv:           extracts the positive RPA frequencies and
C                             corresponding eigenvectors if use rg routine
C
C Called by: the main routine
C======================================================================
C      include 'gcb_dim.inc'

C------OUTPUT--------
      real ecor       ! correlation energy
      integer errRPA  ! -1: stability matrix negatively defined
                      !  0: stability matrix positively defined
                      

C===========INTERNAL VARIABLES======================================


C-------RPA---------------------
      real*8,allocatable   :: rpa_matrix(:,:)  ! full RPA matrix
      real*8,allocatable   :: w(:),w_j(:)      ! real & imaginary part of
                                               ! RPA frequencies
                                               ! when using rg routine
      real*8,allocatable    :: z(:,:)

      integer rpadim
      real*8,allocatable   :: work3(:)  ! dummy working array for diagonalization
      integer,allocatable   :: work2(:) ! dummy working array for diagonalization
      integer ierr            ! flag for RPA routine not used for now
      real eigtol             ! tolerance for zero eigenvalues
      parameter(eigtol=5.0e-2)

      allocate( rpa_matrix(2*nph,2*nph) )
      allocate( w(2*nph),w_j(2*nph) )
      allocate( z(2*nph,2*nph) )
      allocate( work3(2*nph), work2(2*nph))
C.......SET THE FULL RPA MATRIX.................
      call get_rpa(a,b,nph,rpa_matrix)
C......Calculate the RPA frequencies.............................
      call rg(2*nph,2*nph,rpa_matrix,w,w_j,1,z,work2,work3,ierr)
C...and extract the positive ones and corresponding eigenvectors......
      call RPA_eigv(z,w,nph,2*nph,eigtol)
      w1=>w0
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine map_ph(n_states,n_part,rho,vec,map,n)
C===================================================================
C   Maps particle-hole state pairs
C===================================================================
      use flags
      implicit none
C===================================================================
C
C   INPUT:
C        n_states: # of states
C        n_part:   # of particles
C        rho:      s.p. density matrix in fundamental basis
C        vec:      the transformation matrix to ph basis
C
C   OUTPUT:
C        map:       map(k,1) the particle state corresponding to the k-th 
C                            ph pair
C                   map(k,2) the hole state corresponding to the k-th
C                            ph pair
C        n:         total # of ph pairs
C
C  Subroutines called: none
C  Called by:          Solve_RPA
C============================================================================
C----INPUT----------
      integer n_states,n_part
      real*8,intent(IN)   :: rho(n_states,n_states)
      real,intent(IN)     :: vec(n_states,n_states)
      integer n

C----OUTPUT---------
      integer,intent(OUT) :: map(n,2)

C----Dummies----
      real,dimension(n_states,n_states) :: array1,array2,array3
      integer occ(n_part),uocc(n_states-n_part)  ! occupied/unoccupied states

      real sum


      integer i,j,k,m

C      write(12,*)
C      do k=1,n_states
C         write(12,'(20(F5.3,1x))')(rho(i,k),i=1,n_states)
C      enddo

c      sum=0.0
      do i=1,n_states
         do j=1,n_states
         array2(i,j)=sngl(rho(i,j))
         enddo
      enddo
c      print*,'sum=',sum
C      print*,n_part,n_states

      array1=transpose(vec)
      array3=matmul(array1,array2)
      array1=0.0
      array1=matmul(array3,vec)

CCCCCCCCCCC      call matrix_mult(vec,vec,array3,n_states,n_states,n_states)
C      call matrix_mult(array1,array2,array3,n_states,n_states,n_states)
C      call matrix_mult(array3,vec,array1,n_states,n_states,n_states)

C      print*,' '
C      print20,array1
 20   format(12F10.6)

C... Array1 is now rho in particle-hole basis: ones and zeros on diagonal

      k=0
      m=0
      do i=1,n_states
        if(array1(i,i).gt.0.9)then
           k=k+1
           if(k>n_part)then
             print*,' Problem with this SD: it seems that [rho,h].neq.0'
             print*,k,n_part,' (1) '
C             stop
             noRPA=.true.
             return
           endif
           occ(k)=i
c           print*,'k=',k
        else
           m=m+1
           if(m>n_states-n_part)then
             print*,' Problem with this SD: it seems that [rho,h].neq.0'
             print*,m,n_states-n_part,' (2)'
C             stop
             noRPA=.true.
             return
           endif
C           print*,'m=',m
           uocc(m)=i
        endif
      enddo

      if(k/=n_part)then
         print*,'Error in mapping particle-hole states'
         stop
      endif

C      write(12,180)k,m
C      write(12,180)(occ(i),i=1,k)
C      write(12,180)(uocc(i),i=1,m)
  180 format(20I6)

      k=0
      do i=1,n_part                    ! Mapp particle-hole states
         do j=1,n_states-n_part
            k=k+1
            if(k.gt.n_part*(n_states-n_part))then
              write(6,*)'Problem with dimensions'
              stop
            endif
            map(k,2)=occ(i)            ! hole states
            map(k,1)=uocc(j)           ! particle states
	 enddo
      enddo
      if(n/=k)stop 'Something wrong 111'
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine setAandB(nph,map,hvec,map_comb,nmat,energ,vec,ns,a,b)
C=================================================================
C Calculate A and B for each species
C=================================================================
      implicit none
C=================================================================
C
C  INPUT:
C         nph:    # of ph paires
C         nmat:   # of nonzero matrix elements
C         map_comb(i): list of (packed) nonzero indices, i = 1 to nmat
C         energ:  s.p. energies
C         vec:    transformation matrix from fundamental to ph basis
C
C  OUTPUT:
C         A:      diagonal block of RPA matrix
C         B:      off-diagonal block of RPA matrix
C WARNING: in the case of two species of particles one has to add also
C          the pn part (this subroutine calculates only pp or nn part)
C          Still a problem the sing of B (3/28/01)
C
C Subroutines called: none
C
C Called by:          Solve_RPA
C=================================================================

C-------INPUT-----------------
      integer,intent(IN)    :: nmat,ns
      integer,intent(IN)    :: nph                  ! number of p-h paires
      integer,intent(IN)    :: map_comb(nmat)
      real,intent(IN)       :: hvec(nmat)
      real,intent(IN)       :: energ(ns),vec(ns,ns)
      integer,intent(IN)    :: map(nph,2)

C------OUTPUT--------
      real,intent(OUT)      :: a(nph,nph)
      real,intent(OUT)      :: b(nph,nph)

C-----Dummies--------
      integer i,j,k,l,m,i1,i2,i3,i4,ii,n

      real*8 tmp1,temp1
      real*8 tmp2,temp2
      integer deltakr


      do k=1,nph
         m=map(k,1)
	 i=map(k,2)
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP& PRIVATE(l,n,j,temp2,ii,tmp1,tmp2,i1,i2,i3,i4)
!$OMP DO
         do l=1,nph
	   n=map(L,1)
	   j=map(L,2)
           tmp2= 0.0
           temp2=0.0
!!!!$OMP PARALLEL DEFAULT(SHARED)
!!!!!$OMP& PRIVATE(ii,i1,i2,i3,i4) REDUCTION(+:temp2,tmp2)
!!!!!$OMP DO
           do ii=1,nmat
             call unpack0(map_comb(ii),i1,i2,i3,i4)
             tmp1=    -vec(i1,m)*vec(i2,j)*vec(i3,n)*vec(i4,i)
             tmp1=tmp1-vec(i2,m)*vec(i1,j)*vec(i4,n)*vec(i3,i)
             tmp1=tmp1+vec(i2,m)*vec(i1,j)*vec(i3,n)*vec(i4,i)
             tmp1=tmp1+vec(i1,m)*vec(i2,j)*vec(i4,n)*vec(i3,i)
             tmp2=tmp2+hvec(ii)*tmp1

             temp1=     -vec(i1,m)*vec(i2,n)*vec(i3,j)*vec(i4,i)
             temp1=temp1-vec(i2,m)*vec(i1,n)*vec(i4,j)*vec(i3,i)
             temp1=temp1+vec(i2,m)*vec(i1,n)*vec(i3,j)*vec(i4,i)
             temp1=temp1+vec(i1,m)*vec(i2,n)*vec(i4,j)*vec(i3,i)
             temp2=temp2+hvec(ii)*temp1
           enddo
!!!$OMP END DO
!!!$OMP END PARALLEL
           a(k,l)=deltakr(i,j)*deltakr(m,n)*(energ(m)-energ(i))
     &                  -sngl(tmp2) !/2.
           b(k,l)=sngl(-temp2)      !/2.
         enddo
!$OMP END DO
!$OMP END PARALLEL
      enddo


      return
      end

      subroutine get_AorB(app,ann,apn,nphp,nphn,a)
C====================================================================
C            Calculates the full A & B matrices
C====================================================================
      implicit none
C====================================================================
C
C INPUT:
C     app/ann/apn:      pp/nn/pn parts of the diagonal or off-diagonal
C                       RPA matrix
C     nphp/nphn:        proton/neutron ph pair number for p/n
C
C OUTPUT:
C     A                 diagonal or off-diagonal blocks of the RPA
C                       matrices
C
C Subroutines called:   none
C
C Called by:            Solve_RPA
C====================================================================

C-----INPUT---------
      integer nphp,nphn
      real app(nphp,nphp),apn(nphp,nphn),
     &     ann(nphn,nphn)

C-----OUTPUT----------
      real a(nphp+nphn,nphp+nphn)

C-----Dummies---------
      integer i,j,i1,j1

      do i=1,nphp
        do j=1,nphp
	  a(i,j)  =app(i,j)
        enddo
        do j=1,nphn
          j1=j+nphp
	  a(i,j1) =apn(i,j)
        enddo
      enddo
      do j=1,nphn
        j1=j+nphp
        do i=1,nphp
           a(j1,i)=apn(i,j)
        enddo
        do i=1,nphn
           i1=i+nphp
           a(i1,j1)=ann(i,j)
        enddo
      enddo

      return
      end


      subroutine get_stability(a,b,nph,stab)
C==================================================================
C  Calculates the stability matrix
C==================================================================
      implicit none
C==================================================================
C
C INPUT:
C        A:  diagonal block of the RPA matrix
C        B:  off-diagonal block of the RPA matrix
C
C OUTPUT:
C        STAB: stability matrix
C               / A   B \
C        STAB= |         |
C               \ B   A /
C
C  Subroutines called: none
C
C  CAlled by: SOLVE_RPA
C==================================================================
C      include 'gcb_dim.inc'

C-------INPUT-------------------
      integer,intent(IN)    :: nph
      real,intent(IN)       :: a(nph,nph),b(nph,nph)

C-----OUTPUT---------------------
      real,intent(OUT)      :: stab(2*nph,2*nph)         ! stability matrix

C-----Dummies
      integer i,j

      stab=0.0

      do i=1,nph
       do j=1,nph
         stab(i,j)=a(i,j)
         stab(i+nph,j+nph)=a(i,j)
       enddo
      enddo
      do i=nph+1,2*nph
         do j=1,nph
            stab(i,j)=b(i-nph,j)
         enddo
      enddo
      do i=1,nph
         do j=nph+1,2*nph
            stab(i,j)=b(i,j-nph)
         enddo
      enddo

      return
      end

      subroutine get_rpa(a,b,nph,rpa_matrix)
C===================================================================
C  Calculates the full RPA matrix
C===================================================================
      implicit none
C===================================================================
C
C INPUT:
C      A/B: diagonal/off-diagonal blocks of the RPA matrix
C      nph: the total number of ph pairs
C
C OUTPUT:
C      RPA_matrix: the full RPA matrix in double precision
C
C                    / A   B \
C      RPA_MATRIX = |         |
C                    \-B  -A /
C
C  Subroutines called: none
C
C  CAlled by: SOLVE_RPA
C===================================================================

C-------INPUT-------------------
      integer,intent(IN)  :: nph
      real,intent(IN)     :: a(nph,nph),b(nph,nph)

C-----OUTPUT---------------------
      real*8,intent(OUT)   :: rpa_matrix(2*nph,2*nph)         ! RPA matrix

C-----Dummies
      integer i,j
      do i=1,2*nph
         do j=1,2*nph
            rpa_matrix(i,j)=0.0
         enddo
      enddo
      do i=1,nph
       do j=1,nph
         rpa_matrix(i,j)=dble(a(i,j))
         rpa_matrix(i+nph,j+nph)=-dble(a(i,j))
       enddo
      enddo
      do i=nph+1,2*nph
         do j=1,nph
            rpa_matrix(i,j)=-dble(b(i-nph,j))
         enddo
      enddo
      do i=1,nph
         do j=nph+1,2*nph
            rpa_matrix(i,j)=dble(b(i,j-nph))
         enddo
      enddo

      return
      end

      subroutine RPA_eigv(z,w,n,np,eigtol)
C==============================================================
C Extracts the RPA eigenvectors when using rg routine
C==============================================================
      use RPAsolutions
      use RPAmatrices
      implicit none
C==============================================================

C-------INPUT-------------
      integer,intent(IN) :: np           ! physical dimension of the arrays
      integer,intent(IN) :: n            ! real dimension of the arrays
      real*8,intent(IN)  :: w(np)        ! real part of eigenfrequencies
      real,intent(IN)    :: eigtol

      real*8,intent(IN)  :: z(np,np)     ! eigenvectors returned by rg routine

C-------Miscelaneous------
      integer i,j,ii,k,kk
      real*8 norm
      real tol
      parameter(tol=1.0E-3)
      real*8 temp1,temp2

      real*8              :: x0(n),y0(n)
      real                :: ZTMP(np,np),WTMP(np)


      ii=0
      kk=0

      

      do i=1,2*n
         if(w(i).gt.eigtol)then
           ii=ii+1
         endif
      enddo
      n0=N-ii
      write(63,*)2*n,n0
      write(63,*)w
      allocate(x(N,ii),Y(N,ii),w0(ii))
      X=0.0
      Y=0.0
      W0=0.0
      W1=>W0

      ii=0
      do i=1,2*n
         if(w(i).gt.eigtol)then
           ii=ii+1
           w0(ii)=w(i)
           do j=1,n
             x(j,ii) =z(j,i)
             y(j,ii) =z(j+n,i)
           enddo
         endif
      enddo
      write(63,*)ii
      write(63,*)w0
c      write(6,*)'kk=',kk
C      write(6,*)'Number of nonzero frequencies=',ii
      n0=N-ii

C----Enforce the normalization-------
      do i=1,ii
        norm=0.0
        do j=1,n
          norm=norm+x(j,i)*x(j,i)-y(j,i)*y(j,i)
        enddo
        if(norm.lt.1.D-4)stop 'Zero norm'
        if(norm.lt.0.0)then
           write(6,*)'Error, Y > X'
           return
        endif
        norm=dsqrt(norm)

        do j=1,n
          x(j,i)=x(j,i)/norm
          y(j,i)=y(j,i)/norm
        enddo
      enddo

      do i =1,ii
        temp1= 0.0
        temp2= 0.0
        norm = 0.0
      	do j =1,N
      	    x0(j) = 0.0
      	    y0(j) = 0.0
      	    do k = 1,N
      	    	x0(j) = x0(j) + a(j,k)*x(k,i)+b(j,k)*y(k,i)
      	    	y0(j) = y0(j) - b(j,k)*x(k,i)-a(j,k)*y(k,i)
      	    enddo

      	    norm = norm+(x0(j)**2-y0(j)**2)/w0(i)**2
      	    temp1=temp1+(x0(j)/w0(i)-x(j,i))**2.0
      	    temp2=temp2+(y0(j)/w0(i)-y(j,i))**2.0
	enddo
	temp1=dsqrt(temp1)
	temp2=dsqrt(temp2)

c	write(10,'(5(F10.6,2x))')norm,temp1,temp2
	if(dabs(norm-1.0).gt.tol .or. temp1.gt.tol .or. temp2.gt.tol)then
	  write(6,100)i,w0(i),norm,temp1,temp2
c	  write(6,*)' vectors written to unit 66 '
c	  do j =1,N
C	      write(6,*)x0(j)/w(i),x(j,i),y0(j)/w(i),y(j,i)
c	      write(66,*)x0(j)/w(i),x(j,i),y0(j)/w(i),y(j,i)
c	  enddo
	
	endif
      enddo
100     format(1X,'error in vectors ',I3,2X,5(F10.5,2X))

      return
      end

C....


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCALVINCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine project1ph(op1,nsps,nph,vec,map,xph,yhp)
C-------------------------------------------------------------------------
C
C   this subroutine projects a one-body operator onto vectors 
C   in the p-h space
C
C
C    INPUT
C	op1(alpha,beta): one-body operator in fundamental basis
C
C       nsps = # of single-particle states in fundamental basis
C       nph = # of p-h states
C	vec = p-h vectors
C       map maps ph states to fundamental basis
C
C    OUTPUT
C	xph(i)	: ph (mi) elements
C       yhp(i)  : hp (im) elements
C		NOTE: if op1 is hermitian, real, xph = yhp
C			if op1 is anti-hermitian, xph= -yhp
C
      implicit none
C=================================================================

C-------INPUT-----------------
      integer,intent(IN)      :: nsps
      real,intent(IN)         :: op1(nsps,nsps)		! one-body operator
      integer,intent(IN)      :: nph			! # of particle-hole states
      real,intent(IN)         :: vec(nsps,nsps)
      integer,intent(IN)      :: map(nph,2)

C------OUTPUT--------
      real,intent(OUT)        :: xph(nph)		! particle-hole amplitudes
      real,intent(OUT)        :: yhp(nph)		! hole-particle amplitudes

C-----Dummies--------
      integer i,k,m,i1,i2
      real*8 size


      size=0.0
      do k=1,nph
         m=map(k,1)
	 i=map(k,2)
	 xph(k) = 0.0
	 yhp(k) = 0.0
	 do i1=1,nsps
	   do i2 = 1,nsps
c                xph(k) = xph(k)+vec(m,i1)*op1(i1,i2)*vec(i,i2)
c                yhp(k) = yhp(k)+vec(i,i1)*op1(i1,i2)*vec(m,i2)
                xph(k) = xph(k)+vec(i1,m)*op1(i1,i2)*vec(i2,i)
                yhp(k) = yhp(k)+vec(i1,i)*op1(i1,i2)*vec(i2,m)
                	   
	   enddo
	 enddo
	 size = size + xph(k)**2.0
c	 write(32,'(2(F10.6,1x))')xph(k),yhp(k)
      enddo
c      write(32,*)

c      do k=1,nph
c         write(64,'(2(F10.5))')xph(k),yhp(k)
c      enddo
C       write(64,*)'Size=',size
c      write(64,*)
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine ZeroCorr(A,B,nph,eigtol,Ecorr0,n0)
C========================================================================
C This subroutine calculates the corrections due to zero
C modes in the RPA frequencies. See Ring & Shuck pp. 311-313.
C========================================================================
      implicit none
C========================================================================
C     include 'gcb_dim.inc'

C-----INPUT--------
      integer,intent(IN)      :: nph ! number of particle-hole states
      integer,intent(IN)      :: n0
      real,intent(IN)         :: a(nph,nph),b(nph,nph)
      real,intent(IN)         :: eigtol ! smalness parameter

C----OUTPUT--------
      real ecorr0

C---INTERNAL VAriables---
      real s(nph,nph)   ! A+B
      real d(nph,nph)   ! A-B

      real e(nph)           ! eigenvalues for A+/-B
      real vec(nph,nph)     ! eigenvectors for A+/-B
      real w(nph)            ! dummy working array
c      real prod(max_sps4)         ! vec*((A+/-B)^-1)*vec

      real u(nph,nph)   ! SVD decompozition
      real v(nph,nph)
      real*8 M0(n0)  ! M0 parameter
      real*8 sum1,sum2,sum3
      real p0(nph,n0),q0(nph,n0)
      integer    :: n01

      integer sqp0(n0)  ! =-1 if P0= complex, =1 if P0=real
C     common/zeromodes/P0,Q0,M0,sqp0,n0

C Rotation test
c      real xjp(max_sps2),xjn(max_sps2)
c      real yjp(max_sps2),yjn(max_sps2)
c      integer nphp,nphn

      integer i,j,k
      real small

      small=eigtol/10.


      m0=0.0

      s=a+b
      d=a-b

      call eig(s,nph,nph,e,vec,w)

      sum2=0.0
      n01 = 0
      k=0
      do i=1,nph
         if(e(i).lt.small)then
            k  = k+1
            n01 = n01+1
c            write(64,*)'Imaginary P'
            if(k.eq.1) call prepareInv(d,nph,nph,u,v,w)
            call getProd(u,v,w,vec(1,i),nph,nph,q0(1,n0),M0(n0))

c            sum1=0.0
c            sum3=0.0
            sqp0(n01)=-1
            p0(:,n01)=vec(:,i)
C            do j=1,nphp
c               sum1=sum1+xjp(j)*vec(j,i)
c               sum3=sum3+vec(j,i)*vec(j,i)
C	       p0(j,n01)=vec(j,i)
C            enddo
C            do j=1,nphn
c               sum1=sum1+xjn(j)*vec(j+nphp,i)
C	       p0(j+nphp,n01)=vec(j+nphp,i)
c               sum3=sum3+vec(j,i)*vec(j,i)
C            enddo
c            write(64,*)'s: ',sum1,sum3
C            sum2=sum2+sum1**2
         endif
      enddo

      call eig(d,nph,nph,e,vec,w)
      k = 0
      do i=1,nph
         if(e(i).lt.small)then
            k = k+1
            n01=n01+1
	    sqp0(n0)=1
c            write(64,*)'Real P'
            if(k.eq.1) call prepareInv(s,nph,nph,u,v,w)
            call getProd(u,v,w,vec(1,i),nph,nph,q0(1,n0),M0(n0))
            p0(:,n01)=vec(:,i)

C            sum1=0.0
C            sum3=0.0
C            do j=1,nphp
C               sum1=sum1+xjp(j)*vec(j,i)
C               sum3=sum3+vec(j,i)*vec(j,i)
C	       p0(j,n0)=vec(j,i)
C            enddo
C            do j=1,nphn
C               sum1=sum1-xjn(j)*vec(j+nphp,i)
C	       p0(j+nphp,n0)=vec(j+nphp,i)
C               sum3=sum3+vec(j,i)*vec(j,i)
C            enddo
C            write(64,*)'d: ',sum1,sum3
C            sum2=sum2+sum1**2

         endif
      enddo
c      write(64,*)'***** ',sum2
c      write(64,*)
c      call eigsrt(e,vec,nph,max_sps4)
c      do i=1,nph
c         write(78,'(F10.5)')e(i)
c      enddo
c      write(78,*)

      ecorr0=0.0
      do i=1,n0
          ecorr0 = ecorr0 - 0.25/M0(i)
      enddo

c      write(60,100)n0,(m0(i),i=1,n0)

100   format(1x,I2,' M0= ',5(F10.5,2X))
      return
      end

      subroutine prepareInv(A,n,np,u,v,w)
C==========================================================================
C  Performs the SVD decomposition of the matrix A, and zeroes the small
C  eigenvalues
C==========================================================================
      implicit none

C-----INPUT---------
      integer np,n
      real a(np,np)

C-----OUTPUT------
      real u(np,np),v(np,np),w(np)

C----Dummies-------
      integer i,j
      real wmax,wmin
      real small
      parameter(small=1.0E-3)

      do i=1,np
        do j=1,n
           u(i,j)=0.0
        enddo
      enddo
      
      do i=1,n
          do j=i,n
             u(i,j)=a(i,j)
             if(j.ne.i)u(j,i)=u(i,j)
          enddo
      enddo

      call svdcmp(u,n,n,np,np,w,v)
      
      wmax=0.0
      do i=1,n
         if(w(i).gt.wmax)wmax=w(i)
      enddo
      
      wmin=small*wmax

      do i=1,n
         if(w(i).lt.wmin)w(i)=0.0
      enddo

      return
      end


      subroutine getProd(u,v,w,v1,n,np,v2,prod)
C==========================================================================
      implicit none

C------Input-------
      integer n,np
      real u(np,np),v(np,np),w(np)
      real v1(np)
      real v2(np) ! working vector

C------OUTPUT------
      real*8 prod

      integer i

      do i=1,np
         v2(i)=0.0
      enddo

      call svbksb(u,w,v,n,n,np,np,v1,v2)

      prod=0.0
      do i=1,n
          prod=prod+v1(i)*v2(i)
      enddo

      return
      end

      
      subroutine test_eig(xph,yhp,a,b,nph)
C============================================================================
C Verifies that xph and yph are eigenvectors corresponding to null eigenvalue
C of the RPA matrix
C============================================================================
      implicit none
C============================================================================
C    INPUT:
C	xph(i)	: ph (mi) elements
C       yhp(i)  : hp (im) elements
C		NOTE: if op1 is hermitian, real, xph = yhp
C			  if op1 is anti-hermitian, xph= -yhp
C       nsps = # of single-particle states in fundamental basis
C       nph = # of p-h states
C       A,B = RPA particle-hole matrices
C
C     OUTPUT:
C       none
C
C=============================================================================
C      include 'gcb_dim.inc'
      
C-------INPUT-----------------
      integer nph			    ! # of particle-hole states
      real xph(nph)		! particle-hole amplitudes
      real yhp(nph)		! hole-particle amplitudes
      real A(nph,nph),B(nph,nph)
      
C-----Internal variables-----------
      real v(nph)
      real*8 norm

      integer i,j
      
      norm=0.0
      do i=1,nph
         v(i)=0.0
         
         do j=1,nph
            v(i)=v(i)+a(i,j)*xph(j)+B(i,j)*yhp(j)
         enddo
         norm=norm+v(i)**2
      enddo

C      write(6,*)'Norm=',norm
C      write(65,'(100(D12.5,1X))')(v(i),i=1,nph)
C      write(65,*)

      return
      end

      subroutine RotGen(rot,np,n,rot_gen)
C===================================================================
C Calculates the generator for rotations around Oy axix
C Takes the antisymmetrical part of the rotation matrix
C===================================================================
      implicit none
C===================================================================

      integer np
      real rot(np,np)
      real rot_gen(np,np)
      integer n

      integer i,j

      do i=1,n
         do j=i+1,n
            rot_gen(i,j)=rot(i,j)-rot(j,i)
            rot_gen(j,i)=rot_gen(i,j)
         enddo
         rot_gen(i,i)=0.0
      enddo

      return
      end

      subroutine setAandB_pn(nphp,nphn,mapp,mapn,hvec,map_comb,
     &    nmat,vecp,vecn,nsps,a,b)

C=================================================================
C Calculate A and B for each the pn interaction
C=================================================================
      use nuclide
      implicit none
C=================================================================
C
C  INPUT:
C         nphp/nphn:    # of ph paires for protons/neutrons
C         nmat:   # of nonzero matrix elements
C         map_comb(i): list of (packed) nonzero indices, i = 1 to nmat
C         vecp/vecn:    transformation matrix from fundamental to ph
C                       basis, for protons/neutrons
C
C  OUTPUT:
C         A:      diagonal block of RPA matrix
C         B:      off-diagonal block of RPA matrix
C WARNING: in the case of two species of particles one has to add also
C          the pn part (this subroutine calculates only pp or nn part)
C          Still a problem the sing of B (3/28/01)
C
C Subroutines called: none
C
C Called by:          Solve_RPA
C=================================================================
C      include 'gcb_dim.inc'

C-----INPUT-----
      integer,intent(IN)    :: nsps(2)

C------ph Basis--------------------
      real vecp(nsps(1),nsps(1)),vecn(nsps(2),nsps(2))
      integer nphp,nphn     ! number of proton/neutron matrix elements
                            ! in ph basis
      integer mapp(np*(nsps(1)-np),2),mapn(nn*(nsps(2)-nn),2)
                            ! map of ph states for protons/neutrons

C------interaction----
      integer nmat
      integer map_comb(nmat)
      real hvec(nmat)

C-----OUTPUT----
      real a(nphp,nphn)
      real b(nphp,nphn)


C-----Dummies------
      real*8 tmp1,tmp2
      integer i,j,m,n
      integer k,l
      integer ii,i1,i2,i3,i4


      A=0.0
      B=0.0

      do k=1,nphp          ! PROTON-NEUTRON INTERACTION
         n=mapp(k,1)
	 j=mapp(k,2)
         do l=1,nphn
	   m=mapn(l,1)
	   i=mapn(l,2)
           tmp1=0.0
           tmp2=0.0
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP& PRIVATE(ii,i1,i2,i3,i4) REDUCTION(+:tmp1,tmp2)
!$OMP DO
           do ii=1,nmat
             call unpack0(map_comb(ii),i1,i2,i3,i4)
                tmp1=tmp1+hvec(ii)*vecn(i1,i)*vecn(i2,m)
     &                           *vecp(i4,j)*vecp(i3,n)
                tmp2=tmp2+hvec(ii)*vecn(i1,i)*vecn(i2,m)
     &                           *vecp(i4,n)*vecp(i3,j)
c                tmp2=tmp2+hvec(ii)*vecn(i1,m)*vecn(i2,i)  ! 2/5/2002
c     &                           *vecp(i3,j)*vecp(i4,n)   ! 2/5/2002
           enddo
!$OMP END DO NOWAIT
!$OMP END PARALLEL
           a(k,l)=sngl(tmp1)
           b(k,l)=sngl(tmp2)
         enddo
      enddo
      return
      end


      subroutine traceM2(a,b,n,np,trace)
C===================================================================
C Calculates the trace of a product of 2 matrices
C===================================================================
      implicit none

      integer n,np
      real a(np,np),b(np,np)
      real*8 trace

      integer i,j

      trace=0.0
      do i=1,n
         do j=1,n
           trace=trace+a(i,j)*b(j,i)
         enddo
      enddo


      return
      end


c      subroutine newCorr(ecor,X,Y,W,energp,energn,nphp,nphn,
c     1            mapp,mapn,app,apn,bpp,bpn)
C========================================================================
C Calculates higher order corrections as prescribed by NPA535(1991)1-22
C========================================================================
c      implicit none
C========================================================================
c      include 'gcb_dim.inc'
c
C-----INPUT--------
c      integer nphp,nphn ! number of particle-hole states
c      integer mapp(max_sps2,2),mapn(max_sps2,2)
c                            ! map of ph states for protons/neutrons
c
c      real x(max_sps4,max_sps4),y(max_sps4,max_sps4)
c                                            ! RPA eigenvectors
c      real app(max_sps2,max_sps2),apn(max_sps2,max_sps2),
c     &     ann(max_sps2,max_sps2)
c      real bpp(max_sps2,max_sps2),bpn(max_sps2,max_sps2),
c     &     bnn(max_sps2,max_sps2)
c

