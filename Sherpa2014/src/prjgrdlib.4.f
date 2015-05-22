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

      real,dimension(nsps(1),nsps(1))   :: zpp,hpp,gammap
      real,dimension(nsps(1),nsps(1))   :: znn,hnn,gamman
      real                              :: psdmin(nsps(1),np)
      real                              :: nsdmin(nsps(2),nn)
      real,parameter                    :: dlambda=0.001
      real,parameter                    :: tol=1.E-4
      real,parameter                    :: tole=1.E-6
      real                              :: e,e1min
      logical                           :: exitflag
      integer                           :: iter1
      integer,parameter                 :: ITERMAX=100,ITERMAX1=100
      real                              :: size,over
C      real                              :: Zpp(nsps(1),nsps(1)),Znn0(nsps(2),nsps(2))
      real                              :: bpp(nsps(1),nsps(1)),
     &                                     bnn(nsps(2),nsps(2))
      real,parameter                    :: lstep=0.001
 

      integer                            :: i,j
 
      exitflag=.true.
      ierr=0
      iter=0


      call hmult(nsps,nmatpp,nmatnn,nmatpn,
     &  map_combpp,map_combnn,map_combpn,
     &  hvecpp,hvecnn,hvecpn,e_spepp,e_spenn,
     &  overp,overn,rhop,rhon,over,e1min)
      emin=e1min

      iter=0
 1    continue
C....Compute the one-body Hamiltonian
      call make_Gamma(gammap,gamman)
      call make_hmatrix(gammap,e_spepp,hpp,nsps(1))
      call make_hmatrix(gamman,e_spenn,hnn,nsps(2))
      
C.....Compute the Thouless particle-hole matrix Z=Q*h*P
      call findZ(hpp,Zpp,rhop,qprot,nsps(1),.false.)
      call findZ(hnn,Znn,rhon,qneutr,nsps(2),.false.)
      size=0.0
      call smallZmi(Zpp,nsps(1),size)
      call smallZmi(Znn,nsps(2),size)
      size=sqrt(size)
      if(size<tol)return
      write(21,*)iter+1,size
C      bpp=Zpp
C      bnn=znn

      call findZ(hpp,bpp,rhop,qprot,nsps(1),.true.)
      call findZ(hnn,bnn,rhon,qneutr,nsps(2),.true.)

      do i=1,nsps(1)
        do j=i,nsps(1)
           Zpp(i,j)=dlambda
           Zpp(j,i)=-dlambda
        enddo
        Zpp(i,i)=0.0
      enddo
      do i=1,nsps(2)
        do j=i,nsps(2)
           Znn(i,j)=dlambda
           Znn(j,i)=-dlambda
        enddo
        Znn(i,i)=0.0
      enddo

C      call findZ(hpp,Zpp,rhop,qprot,nsps(1),.true.)
C      call findZ(hnn,Znn,rhon,qneutr,nsps(2),.true.)

      if(pade==0)then
         call exp_Z(Zpp,nsps(1),1.0)
         call exp_Z(Znn,nsps(2),1.0)
      else
         call padeApprpn(Zpp,nsps(1),1.0)
         call padeApprpn(Znn,nsps(2),1.0)
      endif

      psdmin=psd
      nsdmin=nsd

      call getnewSD(Zpp,psdmin,np,nsps(1))
      call getnewSD(Znn,nsdmin,nn,nsps(2))

      call orth_SVD(psdmin,nsps(1),np)
      call orth_SVD(nsdmin,nsps(2),nn)

      call make_projectors(psdmin,rhop,qprot,np,nsps(1),overp)
      call make_projectors(nsdmin,rhon,qneutr,nn,nsps(2),overn)

      call make_Gamma(gammap,gamman)
      call make_hmatrix(gammap,e_spepp,hpp,nsps(1))
      call make_hmatrix(gamman,e_spenn,hnn,nsps(2))

      call findZ(hpp,Zpp,rhop,qprot,nsps(1),.true.)
      call findZ(hnn,Znn,rhon,qneutr,nsps(2),.true.)

      do i=1,nsps(1)
        do j=i+1,nsps(1)
           Zpp(i,j)=-bpp(i,j)*dlambda/(Zpp(i,j)-bpp(i,j))
           Zpp(j,i)=-Zpp(i,j)
        enddo
        Zpp(i,i)=0.
      enddo
c      Zpp=Zpp-transpose(Zpp)

      do i=1,nsps(2)
        do j=i+1,nsps(2)
           Znn(i,j)=-bnn(i,j)*dlambda/(Znn(i,j)-bnn(i,j))
           Znn(j,i)=-Znn(i,j)
        enddo
        Znn(i,i)=0.0
      enddo
c      Znn=Znn-transpose(Znn)

      if(pade==0)then
         call exp_Z(Zpp,nsps(1),1.0)
         call exp_Z(Znn,nsps(2),1.0)
      else
         call padeApprpn(Zpp,nsps(1),lstep)
         call padeApprpn(Znn,nsps(2),lstep)
      endif

      call getnewSD(Zpp,psd,np,nsps(1))
      call getnewSD(Znn,nsd,nn,nsps(2))

      call orth_SVD(psd,nsps(1),np)
      call orth_SVD(nsd,nsps(2),nn)

      call make_projectors(psd,rhop,qprot,np,nsps(1),overp)
      call make_projectors(nsd,rhon,qneutr,nn,nsps(2),overn)

      call hmult(nsps,nmatpp,nmatnn,nmatpn,
     &                map_combpp,map_combnn,map_combpn,
     &                hvecpp,hvecnn,hvecpn,e_spepp,e_spenn,
     &                overp,overn,rhop,rhon,over,e)

      iter=iter+1
      write(20,*)iter,e
C      psd=psdmin
C      nsd=nsdmin
C      call make_projectors(psd,rhop,qprot,np,nsps(1),overp)
C      call make_projectors(nsd,rhon,qneutr,nn,nsps(2),overn)
C      if(iter>ITERMAX)then
C         ierr=-1
C         return
C         exitflag=.false.
C      endif

      goto 1
c      deallocate(hpp,hnn,zpp,znn,psdmin,nsdmin,gammap,gamman)

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

















