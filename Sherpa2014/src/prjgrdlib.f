CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C This library uses gradient descent to compute the minimum
C on the energy surface.
C Uses routines from libgrd.f, which have been initially written
C for this.
C IS, June 2004, UA
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine prjgrd(pade,emin,iter,ierr)
C=================================================================
C Performs the minimization using the gradient descent
C=================================================================
      use wfn
      use spspace
      use Hamiltonian
      use observables
      use nuclide
      implicit none
C=================================================================
C INPUT: pade=integer which controls computation of exp(Z)
C            = 0 exact
C            /=0 Pade Approximation
C        Hamiltonian in module Hamiltonian
C
C OUTPUT:
C      psd/nsd the proton/neutron SDs  in module wfn
C      emin:  HF energy
C      iter:  # of iterations
C      ierr:  error flag
C             >=0 no error
C             <0  minimum not found
C
C SUBROUTINES USED
C      make_Gamma
C      make_hmatrix
C      hmult
C      findZ
C      padeApprpn
C      exp_Z
C      smallZmi
C
C==================================================================
      integer,intent(IN)                :: pade
      real,intent(OUT)                  :: emin
      integer,intent(OUT)               :: iter,ierr

      real,dimension(nsps(1),nsps(1))   :: zpp,hpp,gammap,gammap_old
      real,dimension(nsps(1),nsps(1))   :: znn,hnn,gamman,gamman_old
      real                              :: psdmin(nsps(1),np)
      real                              :: nsdmin(nsps(2),nn)
      real                              :: dlambda
      real,parameter                    :: tol=1.E-6
      real,parameter                    :: tole=1.E-10
      real                              :: e,e1min
      logical                           :: exitflag
      integer                           :: iter1
      integer,parameter                 :: ITERMAX=100,ITERMAX1=100
      real                              :: size,over
      real                              :: ee1,ee2,ee3
 
      exitflag=.true.
      ierr=0
      iter=0


      call hmult(nsps,nmatpp,nmatnn,nmatpn,
     &  map_combpp,map_combnn,map_combpn,
     &  hvecpp,hvecnn,hvecpn,e_spepp,e_spenn,
     &  overp,overn,rhop,rhon,over,e1min)
      emin=e1min
C      write(22,*)emin

      iter=0
      iter1=0
 1    continue
C....Compute the one-body Hamiltonian
      call make_Gamma(gammap,gamman)
      call make_hmatrix(gammap,e_spepp,hpp,nsps(1))
      call make_hmatrix(gamman,e_spenn,hnn,nsps(2))
      call findZ(hpp,Zpp,rhop,qprot,nsps(1),.true.)
      call findZ(hnn,Znn,rhon,qneutr,nsps(2),.true.)
      size=0.0
      call smallZmi(Zpp,nsps(1),size)
      call smallZmi(Znn,nsps(2),size)
C      size=sqrt(size)
      if(size<tol)return


      if(iter/=0)then
          gammap=0.9*gammap+0.1*gammap_old
          gamman=0.9*gamman+0.1*gamman_old
      endif
      gammap_old=gammap
      gamman_old=gamman
      call make_hmatrix(gammap,e_spepp,hpp,nsps(1))
      call make_hmatrix(gamman,e_spenn,hnn,nsps(2))
      
C.....Compute the Thouless particle-hole matrix Z=Q*h*P
      call findZ(hpp,Zpp,rhop,qprot,nsps(1),.false.)
      call findZ(hnn,Znn,rhon,qneutr,nsps(2),.false.)
c      if(size<0.001)then
c         dlambda=1.
c      else
         dlambda=0.05
c      endif
C      write(21,*)iter+1,size,e
      call findZ(hpp,Zpp,rhop,qprot,nsps(1),.true.)
      call findZ(hnn,Znn,rhon,qneutr,nsps(2),.true.)

      if(pade==0)then
         call exp_Z(Zpp,nsps(1),dlambda)
         call exp_Z(Znn,nsps(2),dlambda)
      else
         call padeApprpn(Zpp,nsps(1),dlambda)
         call padeApprpn(Znn,nsps(2),dlambda)
      endif

      iter1=0
      psdmin=psd
      nsdmin=nsd
      do
         call getnewSD(Zpp,psd,np,nsps(1))
         call getnewSD(Znn,nsd,nn,nsps(2))
c         if(pade/=0)then
            call orth_SVD(psd,nsps(1),np)
            call orth_SVD(nsd,nsps(2),nn)
c         endif
         call make_projectors(psd,rhop,qprot,np,nsps(1),overp)
         call make_projectors(nsd,rhon,qneutr,nn,nsps(2),overn)
         call hmult(nsps,nmatpp,nmatnn,nmatpn,
     &                map_combpp,map_combnn,map_combpn,
     &                hvecpp,hvecnn,hvecpn,e_spepp,e_spenn,
     &                overp,overn,rhop,rhon,over,e)
         iter1=iter1+1
         write(20,*)iter1,e
         if(abs(e1min-e)<tole)exit
         if(e1min<e)exit
         e1min=e
C         print*,e
         if(e1min<emin)then
            psdmin=psd
            nsdmin=nsd
            e1min=e
         else
            exit
         endif
         if(iter1>ITERMAX1)exit
      enddo

      iter=iter+1
      psd=psdmin
      nsd=nsdmin
      call make_projectors(psd,rhop,qprot,np,nsps(1),overp)
      call make_projectors(nsd,rhon,qneutr,nn,nsps(2),overn)
      if(iter>ITERMAX)then
         ierr=-1
         return
C         exitflag=.false.
      endif

      goto 1

      return
      end subroutine prjgrd

      subroutine smallZmi(Z,ns,size)
C=========================================================================
C Evaluates the size of the QhP MEs, which have to be zero
C=========================================================================
      implicit none

      integer,intent(IN)                 :: ns
      real,intent(IN)                    :: Z(ns,ns)
      real,intent(INOUT)                 :: size

      integer                            :: i,j
      real*8                             :: sum

      sum=0.0d0
      do i=1,ns
       do j=i+1,ns
          sum=sum+Z(i,j)**2
       enddo
      enddo

      size=size+sum

      return
      end subroutine smallZmi

















