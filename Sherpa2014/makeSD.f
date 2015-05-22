C
      subroutine random_sd(iseed,psd,nsd)
C
C   G.C. Basis version 0.1
C
C   creates a random slater determinant
C   WARNING: not terribly sophisicated
C
C   INPUT:
C	nsps(x)		= # s.p. m-states for x = p, or n
C	spsqn		= quantum numbers of those m-states
C       np,nn		= # of protons, neutrons
C	Jz		= total Jz
C
C   OUTPUT:
C	pSD(m,i) 	= proton Slater determinant as matrix
C			i labels which particle, m labels m-state
C	nSD(m,i)	= neutron Slater determinant
C
C   SUBROUTINES CALLED:
C	eig:		diagonalizes mock Hamiltonian
C	eigsrt:		sorts eigenstates into ascending order
C
C   FUNCTIONS CALLED:
C       sdjz:           computes 2 x jz of p,n slater determinants         
C	gauss:		generates gaussian random variables
C

      use spspace
      use nuclide
      implicit none 

C-------------SLATER DETERMINANTS------------------------------
   
C      integer np,nn		! # of protons, neutrons
      
      real pSD(nsps(1),np) ! proton slater determinant as matrix
      				! columns = particles, rows = m-states
      real nSD(nsps(2),nn) ! neutron slater determinant       
      
C----------RANDOM NUMBERS--------------------------------
      
      integer iseed
      real gauss		! gaussian random number generator
      
C---------- MOCK S.P. HAMILTONIAN---------------------------

      real,allocatable   :: h(:,:)    ! single-particle hamiltonian  
      real,allocatable   :: e(:)      ! eigenvalues
      real,allocatable   :: vec(:,:)  ! eigenvectors
      real,allocatable   :: work(:)	! dummy needed for subroutine eig

C----------FIXING Jz------------------------------------------

      integer Jz		!  2 x total Jz of SD's
C      integer mu(max_sps)	!  2 x Jz of s.p. m-states
      real sum			! dummy
      
C      integer map(max_sps)
      integer sdjz		! function to compute jz of a SD
      integer jztest		! temp value of Jz
      real lambda

C---------MISC------------------------------------------------

      integer i,j,m		! dummy counters
      integer ns
C==============================================================
      
      ns=nsps(1)
      allocate(h(ns,ns),e(ns),work(ns),vec(ns,ns))
      if(nsps(1).ne.nsps(2))then 
      	write(6,*)' whoa, i am way too stupid to handle this '
      	write(6,*)' proton, neutron spaces different! '
        write(6,*)nsps(1),nsps(2)
      	stop
      endif
      
C-------------------GENERATE MOCK S.P. HAMILTONIAN
      do i = 1,nsps(1)
        do j = i,nsps(1)
	   h(i,j)=gauss(iseed)
	   h(j,i)=h(i,j)
      	enddo
      enddo

C------------DIAGONALIZE MOCK HAMILTONIAN---------------
      call eig(h,ns,ns,e,vec,work)
C      call eigsrt(e,vec,nsps(1),max_sps)


C      print*,nsps

      do i = 1,np
      	do m = 1,nsps(1)
      	    pSD(m,i) = vec(m,i)
      	enddo
      enddo
      	
C       next, neutrons
      h=0.0
      do i = 1,nsps(1)
        do j = i,nsps(1)
	   h(i,j)=gauss(iseed)
	   h(j,i)=h(i,j)
      	enddo
      enddo

C------------DIAGONALIZE MOCK HAMILTONIAN---------------
      call eig(h,nsps(1),nsps(1),e,vec,work)

      do i = 1,nn
      	do m = 1,nsps(1)
      	    nSD(m,i) = vec(m,i)
      	enddo
      enddo
      deallocate(h,e,work)
c      iseed=107*gauss(iseed)	            
      return
      
      end subroutine random_sd
      
