
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       subroutine setupABpn
C======================================================================
C Calculates the pnRPA matrices and solves the RPA equations
C======================================================================
      use spspace
      use wfn
      use Hamiltonian
      use phBASIS
      use flags
      use nuclide
      use RPAmatrices
      use RPAsolutions
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

      real,dimension(nsps(1),nsps(1))  :: gammap,hpp
      real,dimension(nsps(2),nsps(2))  :: gamman,hnn

      real,allocatable     :: work(:)       ! dummy working array for
                                            ! diagonalization

                            ! map of ph states for protons/neutrons

      integer n1,n2
      integer n11,n22

      real efp,efn

      integer occp(np),uoccp(nsps(1)-np)
      integer occn(nn),uoccn(nsps(2)-nn)

C-----Temporary variables
      integer i,j,k,m
      character*1 RefSPE

      if(RPAflag)then ! eliminate the like-RPA
          RPAflag=.false.
          deallocate(A,B,X,Y,W0,App,Ann,Apn,Bpp,Bnn,Bpn)
      else
         call make_Gamma(gammap,gamman)

         call make_hmatrix(gammap,e_spepp,hpp,nsps(1))
         call make_hmatrix(gamman,e_spenn,hnn,nsps(2))


         allocate(vecp(nsps(1),nsps(1)),energp(nsps(1)))
         allocate(work(nsps(1)))
         call eig(hpp,nsps(1),nsps(1),energp,vecp,work)     ! for protons
         deallocate(work)


C         call map_ph(nsps(1),np,rhop,vecp,map_prot,nphp)

         allocate(vecn(nsps(2),nsps(2)),energn(nsps(2)))
         allocate(work(nsps(2)))
         call eig(hnn,nsps(2),nsps(2),energn,vecn,work)        ! for neutrons
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
!----------------- END OF CORRECTION

C         call map_ph(nsps(2),nn,rhon,vecn,map_neutr,nphn)
      endif

      nphpn=(nsps(1)-np)*nn
      nphnp=(nsps(2)-nn)*np
      nph_pn=nphpn+nphnp

      if(allocated(map_pn))deallocate(map_pn)
      if(allocated(map_np))deallocate(map_np)
      allocate(map_pn(nphpn,2),map_np(nphnp,2))

C      subroutine occupation(vec,uocc,occ,rho,nsps,np)
      call occupation(vecp,uoccp,occp,rhop,nsps(1),np)
      call occupation(vecn,uoccn,occn,rhon,nsps(2),nn)

      call map_phpn(uoccp,occn,map_pn,nsps(1),np,nn,nphpn)
      call Fen(occp,energp,nsps(1),np,efp)
      call map_phpn(uoccn,occp,map_np,nsps(2),nn,np,nphnp)
      call Fen(occn,energn,nsps(2),nn,efn)

      if(allocated(Apnnp))deallocate(Apnnp)
      if(allocated(Anppn))deallocate(Anppn)
      if(allocated(Bpnnp))deallocate(Bpnnp)
      if(allocated(Bnppn))deallocate(Bnppn)
      if(allocated(A))deallocate(A)
      if(allocated(B))deallocate(B)

      allocate(Apnnp(nphpn,nphpn),Anppn(nphnp,nphnp))
      allocate(Bpnnp(nphpn,nphnp),Bnppn(nphnp,nphpn))
      allocate(A(nph_pn,nph_pn),B(nph_pn,nph_pn))

 111  continue
      write(6,*)
C.... Start setting up the pnRPA matrices A and B
      write(6,*)'Introduce option for computing A matrix:'
      write(6,*)'   (F)S.P. energies relative to the Fermi levels'
      write(6,*)'   (B)S.P. energies relative to the bottom of each',
     1' well'
      read(5,*)RefSPE
      if(RefSPE.eq.'F' .or. RefSPE.eq.'f')then
        call setA_pn(nphpn,map_pn,hvecpn,map_combpn,
     &   nmatpn,energp,energn,vecp,vecn,apnnp,efp,efn,1,nsps(1),nsps(2))
        call setA_pn(nphnp,map_np,hvecpn,map_combpn,
     &   nmatpn,energn,energp,vecp,vecn,anppn,efn,efp,2,nsps(2),nsps(1))
       else if(RefSPE.eq.'B' .or. RefSPE.eq.'b')then
        call setA_pn(nphpn,map_pn,hvecpn,map_combpn,
     &   nmatpn,energp,energn,vecp,vecn,apnnp,0.0,0.0,1,nsps(1),nsps(2))
        call setA_pn(nphnp,map_np,hvecpn,map_combpn,
     &   nmatpn,energn,energp,vecp,vecn,anppn,0.0,0.0,2,nsps(2),nsps(1))
       else
         write(6,*)'This option is not available. Please try again.'
	 goto 111
c       endif
      endif
      call setB_pn(nphpn,nphnp,map_pn,map_np,hvecpn,map_combpn,nmatpn,
     &             vecp,vecn,Bpnnp,Bnppn,nsps(1),nsps(2))


C......SET UP THE full A & B MATRICES...................................


      call get_fullA(Apnnp,Anppn,nphpn,nphnp,A)
      call get_fullB(Bpnnp,Bnppn,nphpn,nphnp,B)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine solve_pnstability(errRPA)
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

      use RPAstab
      use RPAmatrices
      implicit none

C------Stability Matrix-----------------------------------------
      real stab(2*nph_pn,2*nph_pn)         ! stability matrix
      real stabv(2*nph_pn,2*nph_pn)        ! stability eigenvectors not necess.
      real stabw(2*nph_pn)                 ! dummy working array for
                                           ! diagonalization
C..............MISC..........................
      integer i
      logical,intent(OUT)        :: errRPA

C........Test if the stability matrix is positive......................
      errRPA = .false.
C      return

      if(allocated(stabe))deallocate(stabe)
      allocate(stabe(2*nph_pn))


      call get_stability(a,b,nph_pn,stab)
      call eigval(stab,2*nph_pn,2*nph_pn,stabe,stabv,stabw)
      call eigsrt(stabe,stabv,2*nph_pn,2*nph_pn)

      do i=1,2*nph_pn
        if(stabe(i).lt.-0.1)then
           errRPA=.true.
        endif
      enddo
      
C      print1000,stabe
 1000 format(10F12.6)

      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine solve_pnRPA
C======================================================================
C Calculates the RPA matrices and solves the RPA equations
C======================================================================
      use RPAmatrices
      use RPAsolutions
      use nuclide
C      use spspace
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


C------OUTPUT--------
      real ecor       ! correlation energy
      integer errRPA  ! -1: stability matrix negatively defined
                      !  0: stability matrix positively defined


C===========INTERNAL VARIABLES======================================


C-------RPA---------------------
c      real a(max_sps4,max_sps4),b(max_sps4,max_sps4)

      real*8,dimension(nph_pn,nph_pn)  :: rpa_matrix  ! full RPA matrix

      real*8,dimension(nph_pn)           :: w,w_j       ! real & imaginary part of
                                                          ! RPA frequencies
                                                          ! when using rg routine
                                                          ! RPA eigenvectors

      real*8 z(nph_pn,nph_pn)
      real*8 work3(nph_pn)  ! dummy working array for diagonalization
      integer work2(nph_pn) ! dummy working array for diagonalization
      integer ierr            ! flag for RPA routine not used for now
      real,parameter :: eigtol=0.05             ! tolerance for zero eigenvalues

      integer i

C      print*,'Nph_pn=',nph_pn

      if(allocated(Xpn))deallocate(Xpn)
      if(allocated(Xnp))deallocate(Xnp)
      if(allocated(Ypn))deallocate(Ypn)
      if(allocated(Ynp))deallocate(Ynp)
      if(allocated(wpn))deallocate(wpn)
      if(allocated(wnp))deallocate(wnp)

C      allocate(wpn(nphpn),wnp(nphnp))

      call get_RPApn(anppn,bnppn,apnnp,nphnp,nphpn,rpa_matrix)

C......Calculate the RPA frequencies.............................
      call rg(nph_pn,nph_pn,rpa_matrix,w,w_j,1,z,work2,work3,ierr)
      do i=1,nph_pn
       if(dabs(w_j(i)).gt.1D-3)write(6,*)'Warning! Imaginary W!',w_j(i)
      enddo
C      print1000,w
 1000 format(10F12.6)
C...and extract the positive ones and corresponding eigenvectors......
      call RPA_eigvpn(z,w,nphnp,nphpn,nph_pn)
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine map_phpn(uocc,occ,map12,nsps,n1,n2,nph12)
C======================================================================
C Maps paricle- one kind with hole - the other kind
C=======================================================================
      implicit none
C      include 'gcb_dim.inc'
      integer nsps,n1,n2
      integer nph12         ! the number of such elements
      integer uocc(nsps-n1),occ(n2)

C---OUTPUT
      integer map12(nph12,2)
                            ! map of ph states for protons and neutrons

      integer i,j,k

      k=0

      do i=1,nsps-n1
         do j=1,n2
	    k=k+1
	    map12(k,1)=uocc(i)
	    map12(k,2)=occ(j)
	 enddo
      enddo
      if(k/=nph12)stop 'Wrong ph dimension in pn mapping.'

      return
      end



      subroutine Fen(occ,e,np,n,ef)
C======================================================================
C Calculates the Fermi energy
C======================================================================
      implicit none

      integer np,n
      integer occ(np)
      real e(np),ef

      integer i

      ef=e(occ(1))
      do i=2,n
         if(e(occ(i)).gt.ef)ef=e(occ(i))
      enddo

      return
      end


      subroutine RPA_eigvpn(z,w,n1,n2,np)
C==============================================================
C Extracts the pnRPA eigenvectors when using rg routine
C==============================================================
      use RPAmatrices
      use RPAsolutions
      implicit none
C==============================================================

C-------INPUT-------------
      integer np          ! physical dimension of the arrays
      integer n1,n2       ! real dimensions of arrays
      real*8 w(np)        ! real part of eigenfrequencies

      real*8 z(np,np)     ! eigenvectors returned by rg routine

C-------OUTPUT------------
      real,pointer,dimension(:,:)  :: x1,y1               ! X & Y RPA vectors
      real,pointer,dimension(:,:)  :: x2,y2               ! X & Y RPA vectors
C      real,pointer,dimension(:)  :: w1,w2               ! positive RPA eigenfrequencies


C-------Miscelaneous------
      integer i,j,ii,k,kk
      real*8 norm
      integer size
      real*8 small

      small=1.D-4


C      if(n1+n2/=np)stop
      allocate(Xnp(nphnp,nphnp),Ynp(nphpn,nphnp))
      allocate(Xpn(nphpn,nphpn),Ypn(nphnp,nphpn))
      allocate(wnp(nphnp),wpn(nphpn))
      X2=>Xpn
      X1=>Xnp
      Y2=>Ypn
      Y1=>Ynp


 101  continue
      ii=0
      kk=0
      do i=1,n1+n2
         norm=0.0
         do j=1,n1
          norm=norm+z(j,i)**2
         enddo
	 do j=1,n2
	  norm=norm-Z(j+n1,i)**2
	 enddo
	 if(dabs(norm).lt.small)goto 100
         if(norm .gt. 0.0)then
	   norm=dsqrt(norm)
           ii=ii+1
C           if(.not.countonly)then
              wnp(ii)=w(i)
              do j=1,n1
                x1(j,ii)=z(j,i)/norm
              enddo
              do j=1,n2
               y1(j,ii)=z(j+n1,i)/norm
              enddo
C           endif
	 else
	   norm=dsqrt(-norm)
	   kk=kk+1
C           print*,kk,n2,nphpn
C           if(.not.countonly)then
              wpn(KK)=-w(i)
              do j=1,n1
                 y2(j,kk)=z(j,i)/norm
              enddo
	      do j=1,n2
                 x2(j,kk)=z(j+n1,i)/norm
	      enddo
           endif
C         endif
  100    continue
      enddo
C      if(countonly)then
C         w1=>wpn
C         w2=>wnp
C         countonly=.false.
         print*,'Renormalizable pnRPA vectors:',nphpn,nphnp
C         goto 101
C      endif

      return
      end

      
      
      subroutine get_RPApn(a1,b,a2,n1,n2,rpa)
C==============================================================================
      implicit none
C      include 'gcb_dim.inc'

      integer,intent(IN)      :: n1,n2
      real,intent(IN)         :: a1(n1,n1),a2(n2,n2)
      real,intent(IN)         :: b(n1,n2)

      real*8,intent(OUT)      :: rpa(n1+n2,n1+n2)

      integer i,j

      rpa=0.0

      do i=1,n1
         do j=1,n1
	    rpa(i,j)=a1(i,j)
	 enddo
      enddo

      do i=1,n2
         do j=1,n2
	    rpa(i+n1,j+n1)=-a2(i,j)
	 enddo
      enddo

      do j=1,n2
         do i=1,n1
	    rpa(i,j+n1)=b(i,j)
	    rpa(j+n1,i)=-b(i,j)
	 enddo
      enddo

      return
      end

      

      subroutine setA_pn(nph,mapph,hvec,map_comb,nmat,
     &             e1,e2,vecp,vecn,A,ef1,ef2,iwhich,ns1,ns2)
C=================================================================
C Calculate Apnpn and Bpnpn for each the pn interaction
C=================================================================
      implicit none
C=================================================================
C
C  INPUT:
C         nmat:   # of nonzero matrix elements
C         map_comb(i): list of (packed) nonzero indices, i = 1 to nmat
C         vecp/vecn:    transformation matrix from fundamental to ph
C                       basis, for protons/neutrons
C
C  OUTPUT:
C         A:      diagonal block of RPA matrix
C         B:      off-diagonal block of RPA matrix
C
C Subroutines called: none
C
C Called by:          Solve_RPA
C=================================================================
C      include 'gcb_dim.inc'

C-----INPUT-----

      integer,intent(IN)    :: ns1,ns2,nmat,nph
C------ph Basis--------------------
      real,intent(IN)       :: vecp(ns1,ns1),vecn(ns2,ns1)
      integer,intent(IN)    :: mapph(nph,2)
      real,intent(IN)       :: e1(ns1),e2(ns2)
      real,intent(IN)       :: ef1,ef2
      integer,intent(IN)    :: iwhich

C------interaction----
      integer,intent(IN)    :: map_comb(nmat)
      real,intent(IN)       :: hvec(nmat)

C-----OUTPUT----
      real,intent(OUT)      :: A(nph,nph)

C-----Dummies------
      real*8 tmp
      integer i,j,m,n
      integer k,l,k1,l1
      integer ii,i1,i2,i3,i4
      integer in,jp,mn,np
      integer alpha
      REAL SUM
C      integer occ(max_sps)


      A=0.0

      do k=1,nph
         n=mapph(k,1)
	 j=mapph(k,2)
	 a(k,k)=(e1(n)-ef1)-(e2(j)-ef2)
      enddo

C      return

      do k=1,nph
         if(iwhich.eq.1)then
           n=mapph(k,1)
  	   j=mapph(k,2)
         else
           j=mapph(k,1)
           n=mapph(k,2)
         endif
         do l=k,nph
           if(iwhich.eq.1)then
             m=mapph(l,1)
             i=mapph(l,2)
           else
             i=mapph(l,1)
             m=mapph(l,2)
           endif
	   sum=0.0
           do ii=1,nmat
            call unpack0(map_comb(ii),i1,i2,i3,i4)
            sum=sum+hvec(ii)*vecn(i1,j)*vecn(i2,i)
     &                        *vecp(i3,m)*vecp(i4,n)
	   enddo
	   a(k,l)=a(k,l)-sum
	   a(l,k)=a(k,l)
	enddo
      enddo

      return
      end


      subroutine setB_pn(nphpn,nphnp,mappn,mapnp,hvec,map_comb,nmat,
     &             vecp,vecn,Bpnnp,Bnppn,ns1,ns2)
C=================================================================
C Calculate Apnpn and Bpnpn for each the pn interaction
C=================================================================
      implicit none
C=================================================================
C
C  INPUT:
C         nmat:   # of nonzero matrix elements
C         map_comb(i): list of (packed) nonzero indices, i = 1 to nmat
C         vecp/vecn:    transformation matrix from fundamental to ph
C                       basis, for protons/neutrons
C
C  OUTPUT:
C         A:      diagonal block of RPA matrix
C         B:      off-diagonal block of RPA matrix
C
C Subroutines called: none
C
C Called by:          Solve_RPA
C=================================================================
C      include 'gcb_dim.inc'

C-----INPUT-----

      integer,intent(IN)    :: ns1,ns2,nmat,nphpn,nphnp
C------ph Basis--------------------
      real,intent(IN)       :: vecp(ns1,ns1),vecn(ns2,ns2)
      integer,intent(IN)    :: mappn(nphpn,2),mapnp(nphnp,2)

C------interaction----
      integer,intent(IN)    :: map_comb(nmat)
      real,intent(IN)       :: hvec(nmat)



C-----OUTPUT----
      real Bpnnp(nphpn,nphnp),Bnppn(nphnp,nphpn)

C-----Dummies------
      integer i,j,m,n
      integer k,l,k1,l1
      integer ii,i1,i2,i3,i4
      integer in,jp,mn,np
      REAL SUM1,SUM2


      do k=1,nphpn
         n=mappn(k,1)
         j=mappn(k,2)
         do l=1,nphnp
           m=mapnp(l,1)
           i=mapnp(l,2)
	   sum1=0.0
           sum2=0.0
           do ii=1,nmat
              call unpack0(map_comb(ii),i1,i2,i3,i4)
              sum1=sum1+hvec(ii)*vecn(i1,j)*vecn(i2,m)
     &                        *vecp(i3,i)*vecp(i4,n)
              sum2=sum2+hvec(ii)*vecn(i1,m)*vecn(i2,j)
     &                        *vecp(i3,n)*vecp(i4,i)
           enddo
           bnppn(l,k)=-sum1
           bpnnp(k,l)=-sum1
C           write(63,'(F8.3)')abs(sum1-sum2)
         enddo
c         stop
      enddo
c      write(62,*)
c      do k=1,nphpn
c        write(62,'(80F12.4)')(bpnnp(k,l),l=1,nphnp)
c      enddo
c      write(62,*)
C      stop


      return
      end


      subroutine get_fullA(apnnp,anppn,nph1,nph2,a)
C=========================================================================
C Calculates the full A matrix for pnRPA out of Apnnp and Anppn
C=========================================================================
      implicit none
C=========================================================================
C      include 'gcb_dim.inc'

      integer nph1,nph2
      real Apnnp(nph1,nph1),Anppn(nph2,nph2)
      real A(nph1+nph2,nph1+nph2)


      integer i,j

      A=0.0

      do i=1,nph2
        do j=1,nph2
           A(i,j)=Anppn(i,j)
        enddo
      enddo

      do i=1,nph1
        do j=1,nph1
           A(nph2+i,nph2+j)=Apnnp(i,j)
        enddo
      enddo

      return
      end


      subroutine get_fullB(bpnnp,bnppn,nph1,nph2,B)
C=========================================================================
C Calculates the full A matrix for pnRPA out of Apnnp and Anppn
C=========================================================================
      implicit none
C=========================================================================
C      include 'gcb_dim.inc'

      integer,intent(IN)        :: nph1,nph2
      real,intent(IN)           :: Bpnnp(nph1,nph2),Bnppn(nph2,nph1)
      real,intent(OUT)          :: B(nph1+nph2,nph1+nph2)


      integer i,j

      B=0.0
      do i=1,nph2
        do j=1,nph1
           B(i,j+nph2)=Bnppn(i,j)
        enddo
      enddo

      do i=1,nph1
        do j=1,nph2
           B(nph2+i,j)=Bpnnp(i,j)
        enddo
      enddo

c
c      do i=1,nph1+nph2
c         write(62,'(20(F10.3))')(B(i,j),j=1,nph1+nph2)
c      enddo
c
      return
      end


      subroutine occupation(vec,uocc,occ,rho,nsps,np)
C==============================================================================
C
C Lists the occupied and unoccupied states by transforming the density matrix
C to the particle-hole basis
C
C===============================================================================
      use flags
      implicit none

      integer,intent(IN)  :: nsps,np
      real,intent(IN)     :: vec(nsps,nsps)
      real*8,intent(IN)   :: rho(nsps,nsps)
      real*8              :: part

      integer,intent(OUT) :: uocc(nsps-np),occ(np)

      real,dimension(nsps,nsps) :: rhoph,tmp1,tmp2,tmp0

      integer i,j,m

      tmp2=0.0
      do i=1,nsps
         do j=i,nsps
            tmp1(i,j)=sngl(rho(i,j))
            tmp1(j,i)=sngl(rho(j,i))
            tmp0(i,j)=vec(j,i)
            tmp0(j,i)=vec(i,j)
         enddo
      enddo

      tmp2=matmul(tmp0,tmp1)
      rhoph=matmul(tmp2,vec)

      m=0
      j=0

C      print1,rhoph !(i,i),i=1,nsps)
      part=0.0
      do i=1,nsps
        part=part+rhoph(i,i)
      enddo
C      print*,part,' particles.'
      do i=1,nsps
         if(rhoph(i,i)>0.9)then
            j=j+1
            if(j>np)then
             print*,' Problem with this SD: it seems that [rho,h].neq.0'
             noRPA=.true.
             return
C             stop
            endif
            occ(j)=i
         else
            m=m+1
            if(m>nsps-np)then
             print*,' Problem with this SD: it seems that [rho,h].neq.0'
             noRPA=.true.
C             stop
             return
            endif
            uocc(m)=i
         endif
      enddo
      if(m/=nsps-np .or. j/=np)stop ' Problems with occupation'

 1    format(11F9.4)
      return
      end subroutine occupation
