C
C  Package OpCorr.f
C  Calculates the corrections to different operators read from
C  interaction files
C
C  modified April 2003 CWJ to conform to SHERPA
C
      subroutine OpCorr(X,Y,nphp,nphn,mapp,mapn,
     &      vecp,vecn,nsps,n0,np,nn,O_HF,O_corr)
C================================================================
C Calculates the corrections to the Operator O introduced
C as an interaction file.
C================================================================
C
C  INPUT
C     X,Y
C     nphp,nphn
C     mapp,mapn
C     vecp,vecn
C     nsps
C     n0 	= # of zeroes
C     np,nn
C
C  OUTPUT
C     O_HF	: HF expectation value
C     O_corr    : RPA corrections
C
C  SUBROUTINES CALLED
C
C
      use wfn
      use observables
      implicit none
      
C      include 'gcb_dim.inc'
      
      integer                :: nsps(2),np,nn
      integer,intent(IN)     :: n0   ! the mumber of zero frequencies
      integer,intent(IN)     :: nphp,nphn
      real,intent(IN),dimension(nphp+nphn,nphp+nphn-n0) :: x,y
                                            ! RPA eigenvectors

C---------------PARTICLE-HOLE BASIS----------------------------

      real vecp(nsps(1),nsps(1)),vecn(nsps(2),nsps(2))
      integer mapp(nphp,2),mapn(nphn,2)
      integer nph
     
      real over
      real O_HF   ! HF average
      real*8 O_corr   ! RPA correction
      real*8 O_X(nphp+nphn-n0)    ! mean values for the Operator in
                                  ! the excited states
      real O_corr0  ! correction due to the zero eigenmodes

C-----One-body Part-----
      real hpp(nsps(1),nsps(1)),hnn(nsps(1),nsps(1))
      real gammap(nsps(1),nsps(1)),gamman(nsps(2),nsps(2))

C-----INTERNAL VARIABLES------
      real zero(nsps(1))
C      real zminj(max_sps2,max_sps2)  !ph correlations in the ground state
      real*8 oc,nrm

      real app(nphp,nphp),apn(nphp,nphn),
     &     ann(nphn,nphn)
      real bpp(nphp,nphp),bpn(nphp,nphn),
     &     bnn(nphn,nphn)

      real Opp_ph(nsps(1),nsps(1)),Onn_ph(nsps(1),nsps(1))

      real,dimension(nphp+nphn,nphp+nphn) :: a,b
      integer ierr

C      real xjp(max_sps2),xjn(max_sps2)
C      real yjp(max_sps2),yjn(max_sps2)
C      integer n1,n2
C      common/rotTest/xjp,xjn,yjp,yjn,n1,n2


C------Dummies-----
      integer i,j,k

      nph=nphp+nphn

      zero=0.0
      Opp_ph=0.0
      Onn_ph=0.0

C------Hartree-Fock value---------
      call hmult(nsps,nmatpp_O,nmatnn_O,nmatpn_O,
     &  map_combpp_O,map_combnn_O,map_combpn_O,
     &  hvecpp_O,hvecnn_O,hvecpn_O,e_spepp_O,e_spenn_O,
     &  overp,overn,rhop,rhon,over,O_HF)

      call make_GammaO(gammap,gamman)
      call make_hmatrix(gammap,e_spepp_O,hpp,nsps(1))
      call make_hmatrix(gamman,e_spenn_O,hnn,nsps(2))

C------Calculate the 1-body part in ph basis-----
      call OB2ph(vecp,hpp,nsps(1),nsps(1),Opp_ph)
      call OB2ph(vecn,hnn,nsps(1),nsps(2),Onn_ph)

C------Calculate App, Ann, Apn, Bpp, Bnn, Bpn----------
      call setAandB(nphp,mapp,hvecpp_O,map_combpp_O,
     &    nmatpp_O,zero,vecp,nsps(1),app,bpp)                         ! PROTONS
      call add_OneB(app,nphp,Opp_ph,mapp,nsps(1))

      call setAandB(nphn,mapn,hvecnn_O,map_combnn_O,
     &    nmatnn_O,zero,vecn,nsps(1),ann,bnn)                         ! NEUTRONS
      call add_OneB(ann,nphn,Onn_ph,mapn,nsps(2))

      call setAandB_pn(nphp,nphn,mapp,mapn,hvecpn_O,map_combpn_O,
     &    nmatpn_O,vecp,vecn,nsps,apn,bpn)        ! PROTON-NEUTRON INTERACTION

C------Calculate full A & B matrices----------
      call get_AorB(app,ann,apn,nphp,nphn,a)
      call get_AorB(bpp,bnn,bpn,nphp,nphn,b)

      O_corr0=0.0
      call Ocorrections(X,Y,nph-n0,nph,A,B,O_corr,O_X)

C      call nonZeroCorr(X,Y,A,B,O_corr,nph-n0,nph)
c      call ZeroCorr(P,Q,isign,A,O_corr0,n0)
c      O_corr=O_corr-O_corr0
C      write(15,*)'Op. Corr:',O_HF,O_corr,O_HF+O_corr
C      do i=1,nph-n0
C         write(42,100)W(i),O_HF+O_corr+O_X(i),O_X(i)
C      enddo
C      write(42,*)
c      if(n0.eq.0)then
c         call sphZminj(nphp,nphn,X,Y,Zminj,nrm)
c         call getZminj(nph-n0,nph,X,Y,Zminj,nrm)
c         oc=0.0
c         do i=1,nph
c          do j=1,nph
c            oc=oc+Zminj(i,j)*b(j,i)
c          enddo
c         enddo
c         oc=O_HF+oc/2.0
C         oc=oc/nrm
c         write(15,*)'Altern. corrections: ',oc,' Norm=',nrm
c      endif


  100 format(100(F10.5,2x))
      return
      end


      subroutine add_OneB(a,nph,O_ph,map,nsps)
C============================================================
C Adds the one-body contribution to the A matrix
C=============================================================
      implicit none
C=============================================================
C  INPUT:
C     A:     the matrix to be completed
C     O_ph:  the one body part in ph basis
C     nph:   # of ph states
C
C  OUTPUT:
C     A:     completed matrix A with the one-body part
C
C=============================================================

C------INPUT-----
      integer nph,nsps
      real a(nph,nph)
      integer map(nph,2)
      real O_ph(nsps,nsps)

C------OUTPUT----
C        MATRIX A
C

C-----INTERNAL-----

      integer i,j,m,n,k,l
      integer deltaKR  ! Kroneker delta

      do k=1,nph
         m=map(k,1)
	 i=map(k,2)
         do l=1,nph
	    n=map(l,1)
	    j=map(l,2)
 	     a(k,l)=a(k,l)+deltakr(i,j)*O_ph(m,n)
 	     a(k,l)=a(k,l)-deltaKr(m,n)*O_ph(i,j)
	 enddo
      enddo

      return
      end

      subroutine nonZeroCorr(X,Y,A,B,O_corr,n,nph)
C===============================================================
C
C  Calculates the corrections due to nonzero modes
C
C===============================================================
      implicit none
C===============================================================
C      include 'gcb_dim.inc'

C-------INPUT------
      integer n  ! # of nonzero frequencies
      integer nph   ! # of ph states
      real x(nph,n),y(nph,n)
                                            ! RPA eigenvectors
      real a(nph,nph),b(nph,nph)

C---OUTPUT-----                           
      real*8 O_corr

C----INTERNAL---
      integer i,j,k,l

      O_corr = 0.0
      do i=1,n
         do j=1,nph
           do k=1,nph
             O_corr=O_corr+(A(j,k)*Y(k,i)+B(j,k)*X(k,i))*Y(j,i)
           enddo
         enddo
      enddo

      return
      end


      subroutine Ocorrections(X,Y,n,nph,A,B,Ocorr,O_X)
C====================================================================
C  Calculates the corrections to operator averages in RPA
C  Calculates Also the values for the operator in the excited
C  states
C====================================================================
      implicit none


C      include 'gcb_dim.inc'


C-----INPUT--------
      integer n      ! # of non-zero frequencies
      integer nph    ! # of ph states
      real,intent(IN),dimension(nph,n)   :: x,y
                                            ! RPA eigenvectors
      real,intent(IN),dimension(nph,nph) :: a,b

C----OUTPUT--------
      real*8 Ocorr            ! correction to the HF value for GS
      real*8 O_X(n)    ! mean values for the Operator in
                              ! the excited states
C---INTERNAL-------
      integer i,j,lambda,k,l

c      write(23,*)((b(k,l),k=1,nph),l=1,nph)
      do lambda=1,n
         O_X(lambda)=0.0
         do i=1,nph
           do j=1,nph
              O_X(lambda)=O_X(lambda)+A(i,j)*(X(i,lambda)*X(j,lambda)+
     &                                        Y(i,lambda)*Y(j,lambda))
              O_X(lambda)=O_X(lambda)+B(i,j)*(X(i,lambda)*Y(j,lambda)+
     &                                       Y(i,lambda)*X(j,lambda))
           enddo
         enddo
      enddo

      Ocorr=0.0
      do i=1,n
        Ocorr=Ocorr+O_X(i)-A(i,i)
c	write(23,*)O_X(i)
      enddo
      do i=n+1,nph
        Ocorr=Ocorr-A(i,i)
      enddo
      Ocorr=Ocorr/2.0

      return
      end


      subroutine getInv(A,n,np,Ainv,indx)
C=====================================================
C Calculates the inverse of a matrix
C=====================================================
      implicit none
      integer n,np
      real A(np,np),Ainv(np,np)
      integer indx(np)
      integer parity

      integer i,j

      do i=1,n
        do j=i+1,n
           Ainv(i,j)=0.0
           Ainv(j,i)=0.0
        enddo
        Ainv(i,i)=1.0
      enddo

      call ludcmp(A,n,np,indx,parity)
      do j=1,n
         call lubksb(A,n,np,indx,Ainv(1,j))
      enddo

      return
      end

 

      subroutine make_GammaO(gammap,gamman)
C==========================================================================
C  Calculates Gamma matrix (eq. (28) CWJ notes) for a observable, other than
C  the Hamiltonian
C==========================================================================
      use spspace
      use nuclide
      use observables
      use wfn
      implicit none
C==========================================================================
C
C INPUT:
C	nmatpp,nmatnn,nmatpn:
C   -- the routine decodes the hamiltonian into all m-scheme TBME's
C   nmatpp = # of nonzero pp matrix elements; similar for nn, pn
C      
C	map_combpp(i): list of (packed) nonzero indices, i = 1 to nmatpp
C   to get indices IJKL of i'th matrix element do UNPACK(map_combpp)
C	hvecpp(i) = ith TBME associated with indices unpacked from map_combpp(i)
C	similar for map_combnn, map_combpn, hvecnn, hvecpn
C	e_spepp(i,j) = one-body matrix element between i,j m-states
C	e_spenn(i,j)
C	nsps = # of m-states
C       rhopij/rhonij: density matices for protons/neutrons
C
C OUTPUT:
C       gammap/gamman: Gamma matrix for protons/neutrons
C
C CALLED BY: main program
C SUBROUTINES CALLED: unpack0
C==========================================================================
C      include 'gcb_dim.inc'
      real,intent(OUT)      :: gammap(nsps(1),nsps(1)),
     &                         gamman(nsps(2),nsps(2))

      integer i,j,k,l,ii
      real*8,allocatable    :: doublg(:,:)

      allocate(doublg(nsps(1),nsps(1)))
      doublg=0.0


      do ii=1,nmatpp_O
         call unpack0(map_combpp_O(ii),i,j,k,l)
	 doublg(l,i)=doublg(l,i)+hvecpp_O(ii)*rhop(j,k)
	 doublg(k,j)=doublg(k,j)+hvecpp_O(ii)*rhop(i,l)
	 doublg(k,i)=doublg(k,i)-hvecpp_O(ii)*rhop(j,l)
	 doublg(l,j)=doublg(l,j)-hvecpp_O(ii)*rhop(i,k)
      enddo

      do ii=1,nmatpn_O
         call unpack0(map_combpn_O(ii),i,j,k,l)
	 doublg(l,k)=doublg(l,k)+hvecpn_O(ii)*rhon(i,j)
      enddo

      do i=1,nsps(1)
         do j=1,nsps(1)
	    gammap(i,j)=sngl(doublG(i,j))
	 enddo
      enddo

      deallocate(doublg)
      allocate(doublg(nsps(2),nsps(2)))

      doublG=0.0

      do ii=1,nmatnn_O
         call unpack0(map_combnn_O(ii),i,j,k,l)
	 doublg(l,i)=doublg(l,i)+hvecnn_O(ii)*rhon(j,k)
	 doublg(k,j)=doublg(k,j)+hvecnn_O(ii)*rhon(i,l)
	 doublg(k,i)=doublg(k,i)-hvecnn_O(ii)*rhon(j,l)
	 doublg(l,j)=doublg(l,j)-hvecnn_O(ii)*rhon(i,k)
      enddo

      do ii=1,nmatpn_O
         call unpack0(map_combpn_O(ii),i,j,k,l)
	 doublg(j,i)=doublg(j,i)+hvecpn_O(ii)*rhop(k,l)
      enddo

      do i=1,nsps(2)
         do j=1,nsps(2)
	    gamman(i,j)=sngl(doublG(i,j))
	 enddo
      enddo
      deallocate(doublg)
      
      return
      end subroutine make_GammaO
C..................End Subroutine make_Gamma...........................
