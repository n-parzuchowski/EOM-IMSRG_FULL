CCC Date:09/19/00
C Added subroutines to calculate J^2 on 11/28/2000

      subroutine make_projectors(sd,rhoij,q,n_part,n_states,overlap)
C==========================================================================
C   Calculates the density matrices, as well as the Q projection operator
C==========================================================================
      implicit none
C==========================================================================
C Subroutines Called: none
C Called by: the main program
C
C INPUT:
C          sd:      Slater determinant
C          n_part:  particle number
C          n_states: number of states
C
C OUTPUT:
C          rhopij/rhonij: density matices for protons/neutrons
C	   qprot/qneutr:  projectors in unoccupied space for
C	                  protons/neutrons
C  
C==========================================================================
C      include 'gcb_dim.inc'
      integer,intent(IN)                :: n_part,n_states
      real,intent(IN)                   :: sd(n_states,N_part)
      real*8,intent(OUT)                :: rhoij(n_states,n_states)
      real*8,intent(OUT)                :: q(n_states,n_states)

      real,allocatable                  :: Unity(:,:)
      
      real*8,allocatable,dimension(:,:) :: psir,psil
      real*8 overlap
      
      integer i,j,m                     !! dummies
      real*8 delta

      if(n_part ==0)then
            overlap = 1.d0
            return

      endif
      allocate(psir(n_states,n_part),psil(n_states,n_part))

C Assume psd and nsd have real coefficients.
      do i=1,n_part
	 do m=1,n_states
	    psir(m,i)=dble(sd(m,i))
	    psil(m,i)=dble(sd(m,i))
         enddo
      enddo

C..............Calculate rho...............      
        
      call make_rhoij(n_states,n_part,psil,psir,overlap,rhoij)

c      do i = 1,n_states
c        write(27,101)(rhoij(i,j),j = 1,n_states)
c101     format(6f8.4)
c      enddo
C..............Calculate q...............      
      do i=1,n_states
         do j=1,n_states
	    delta=0.0
	    if(i.eq.j)delta=1.0
	    q(i,j)=delta-rhoij(i,j)
	 enddo
      enddo
      deallocate(psir,psil)

      return
      end
C............End Subroutine make_projectors.................................

      subroutine make_Gamma(gammap,gamman)
C==========================================================================
C  Calculates Gamma matrix (eq. (28) CWJ notes)
C==========================================================================
      use spspace
      use nuclide
      use HAMILTONIAN
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


      do ii=1,nmatpp
         call unpack0(map_combpp(ii),i,j,k,l)
	 doublg(l,i)=doublg(l,i)+hvecpp(ii)*rhop(j,k)
	 doublg(k,j)=doublg(k,j)+hvecpp(ii)*rhop(i,l)
	 doublg(k,i)=doublg(k,i)-hvecpp(ii)*rhop(j,l)
	 doublg(l,j)=doublg(l,j)-hvecpp(ii)*rhop(i,k)
      enddo

      do ii=1,nmatpn
         call unpack0(map_combpn(ii),i,j,k,l)
	 doublg(l,k)=doublg(l,k)+hvecpn(ii)*rhon(i,j)
      enddo

      do i=1,nsps(1)
         do j=1,nsps(1)
	    gammap(i,j)=sngl(doublG(i,j))
	 enddo
      enddo

      deallocate(doublg)
      allocate(doublg(nsps(2),nsps(2)))

      doublG=0.0

      do ii=1,nmatnn
         call unpack0(map_combnn(ii),i,j,k,l)
	 doublg(l,i)=doublg(l,i)+hvecnn(ii)*rhon(j,k)
	 doublg(k,j)=doublg(k,j)+hvecnn(ii)*rhon(i,l)
	 doublg(k,i)=doublg(k,i)-hvecnn(ii)*rhon(j,l)
	 doublg(l,j)=doublg(l,j)-hvecnn(ii)*rhon(i,k)
      enddo

      do ii=1,nmatpn
         call unpack0(map_combpn(ii),i,j,k,l)
	 doublg(j,i)=doublg(j,i)+hvecpn(ii)*rhop(k,l)
      enddo

      do i=1,nsps(2)
         do j=1,nsps(2)
	    gamman(i,j)=sngl(doublG(i,j))
	 enddo
      enddo
      deallocate(doublg)
      
      return
      end subroutine make_Gamma
C..................End Subroutine make_Gamma...........................


      subroutine make_hmatrix(gamma,e_spe,h,n_states)
C=======================================================================
C    Calculates matrix h=T+Gamma
C=======================================================================
      implicit none
C=======================================================================
C
C INPUT: 
C         gamma      : Gamma matrix
C	  e_spe(i,j) = one-body matrix element between i,j m-states
C         n_states
C         max_sps
C
C OUTPUT:
C         h: 
C
C=======================================================================
c      include 'gcb_dim.inc'
      integer,intent(IN)  :: n_states
      real,intent(IN)     :: e_spe(n_states,n_states)
      real,intent(IN)     :: gamma(n_states,n_states)
      real,intent(OUT)    :: h(n_states,n_states)

      integer i,j
      
      h=0.0

C....Start calculating.................................

      h=e_spe+gamma

      return
      end subroutine make_hmatrix

C.................End subroutine make_hmatrix.........................

      subroutine findZ(h,Z,rho,q,n_states,flag)
C======================================================================
C Calculates Thouless matrix Z
C======================================================================
      implicit none
C=======================================================================
C
C INPUT: h:
C        n_states
C        q:
C        rhoij:
C
C OUTPUT: zpp/znn= Thouless particle-hole matrix for protons/neutrons
C Called by: main program
C Subroutines called: matrix_mult
C=======================================================================
C      include 'gcb_dim.inc'
      integer,intent(IN)     :: n_states
      real,intent(OUT)       :: z(n_states,n_states)
      real*8,intent(IN)      :: rho(n_states,n_states)
      real*8,intent(IN)      :: q(n_states,n_states)
      real,intent(IN)        :: h(n_states,n_states)
      logical,intent(IN)     :: flag  ! = true compute the antisymmetric Z
                                      ! = false compute just the gradient

      integer i,j

      real*8,dimension(n_states,n_states)::temp,z_temp,doubleH
      

      do i=1,n_states
         do j=1,n_states
	    doubleH(i,j)=-dble(h(i,j))
	 enddo
      enddo

      temp=matmul(q,doubleH)
      z_temp=matmul(temp,rho)

      if(flag)then
      temp=0.0
      do i=1,n_states
         do j=1,n_states
	    temp(i,j)=z_temp(i,j)-z_temp(j,i)
	 enddo
      enddo      
      
      do i=1,n_states
         do j=1,n_states
	    z(i,j)=sngl(temp(i,j))
	 enddo
      enddo     
      else
        do i=1,n_states
          do j=1,n_states
             z(i,j)=-sngl(z_temp(i,j))
          enddo
C          write(99,'(20F10.4)')z(i,1:n_states)
        enddo
C        write(99,*)
      endif

      return
      end
C.............End subroutine findZ....................................


      subroutine getnewSD(Z,sd,n_part,n_states)
C=======================================================================
C  Calculates the Slater determinants for proton and neutrons after
C  applying Thouless transformation.
C=======================================================================
      implicit none
C=======================================================================
C INPUT:
C       z      : Thouless matrix
C       sd:      Slater determinant
C       n_states: # of states   
C       n_part:   number of particles
C
C OUTPUT:
C       psdnew/nsdnew: result of applying Thouless transformation
C
C Subroutines Called: matrix_mult
C Called by: main program
C=======================================================================
C      include 'gcb_dim.inc'
      integer,intent(IN)  :: n_states,n_part
      real,intent(INOUT)     :: sd(n_states,n_part)
      real                   :: sdnew(n_states,n_part)
      real,intent(IN)     :: z(n_states,n_states)
      
C      call matrix_mult(z,sd,sdnew,n_states,n_states,n_part)
      sdnew=matmul(z,sd)
      sd=sdnew

      return
      end

C.................END Subroutine getnewSD.............................

      
      subroutine padeApprpn(Z,n_states,lambda)
C==========================================================================
C Calculates exp(Ztilde) by Pade approximant
C Warning: destrois the Ztilde matrix, and replace it by the approximation
C          of exp(Ztilde)
C
C Pade Apprx: EXP(Z)=(1-Z/2)^{-1} * (1+Z/2)
C==========================================================================
      implicit none
C==========================================================================
C Subroutiones called:
C              ludcmp: LU decomposition of a matrix
C              lubksb: back substitution routine for a LU decomposed matrix
C Called by:
C              gradesc_sd
C
C INPUT: 
C              Zpp: Thouless particle-hole matrix 
C              max_sps: the maximum number of s.p. states
C              n_states:the number of states
C
C OUTPUT:
C==========================================================================
C
C This subroutine calculates separately for protons and neutrons
C
C      include 'gcb_dim.inc'
      integer,intent(IN)     :: n_states
      real,intent(INOUT)     :: z(n_states,n_states)

      real,allocatable       :: tmp(:,:)
      integer,allocatable    :: indx(:)
      
      integer i,j                ! dummies
      real lambda,parity
      real delta
      
      allocate(tmp(n_states,n_states),indx(n_states))
      do i=1,n_states
         do j=1,n_states
            delta = 0.0
            if(i.eq.j)delta=1.0
	    tmp(i,j)=delta - lambda*z(i,j)/2.0
	    z(i,j)=delta + lambda*z(i,j)/2.0
	 enddo
      enddo    
      
      call ludcmp(tmp,n_states,n_states,indx,parity)
      do j=1,n_states
         call lubksb(tmp,n_states,n_states,indx,z(1,j))
      enddo
      deallocate(tmp,indx)

      return
      end subroutine padeApprpn
C...........End subroutine padeApprpn.......................................


      subroutine make_Jz(jz_matrix,spsqn,n_states)
C==========================================================================
C   Calculates the s.p. Jz matrix.
C==========================================================================
      implicit none
C==========================================================================
C Subroutines Called: none
C Called by: the main program
C
C INPUT:
C    spsqn: the s.p. quantum numbers
C    n_states:  the # of s.p. states
C
C OUTPUT:
C    jz_matrix: the matrix corresponding to Jz
C  
C==========================================================================

      integer n_states
      real,intent(OUT)   :: Jz_matrix(n_states,n_states)
      integer,intent(IN) :: spsqn(2,n_states,6)
      				! quantum #'s of s.p. m-states - from WEO
      				! spsqn(t,i,j)  
      				! t = 1 = p, t = 2 = n
      				! i = label of s.p. m-state
      				! j = 1 -> label of j-orbit (note p,n different)
      				! j = 2 -> radial quantum number N
      				! j = 3 -> L 
      				! j = 4 -> 2 x J of state
      				! j = 5 -> 2 x M of state
      				! j = 6 -> 2 x Tz of state ( p/n=+/-1 )
				
      integer m

      Jz_matrix=0.0

      do m=1,n_states
	Jz_matrix(m,m)=spsqn(1,m,5)/2.
      enddo

      return
      end
C.....End subroutine make_Jz

      subroutine make_Jxy(jx_matrix,jy_matrix,spsqn,n_states)
C==========================================================================
C   Calculates Jx and Jy matrices in the state space.
C==========================================================================
      implicit none
C==========================================================================
C Called by: the main program
C
C INPUT:
C    spsqn: the s.p. quantum numbers
C    n_states:  the # of s.p. states
C
C OUTPUT:
C    jx_matrix: the matrix corresponding to Jx
C    jy_matrix: the matrix corresponding to Jy
C
C USES:
C    function deltaKR:  Kroneker delta
C  
C==========================================================================
C      include 'gcb_dim.inc'

      integer,intent(IN)    :: n_states      
      real,intent(OUT)      :: Jx_matrix(n_states,n_states)
      real,intent(OUT)      :: Jy_matrix(n_states,n_states)
      integer,intent(IN)    :: spsqn(2,n_states,6)
      				! quantum #'s of s.p. m-states - from WEO
      				! spsqn(t,i,j)  
      				! t = 1 = p, t = 2 = n
      				! i = label of s.p. m-state
      				! j = 1 -> label of j-orbit (note p,n different)
                                ! j = 2 -> radial quantum number N
      				! j = 3 -> L 
      				! j = 4 -> 2 x J of state
      				! j = 5 -> 2 x M of state
      				! j = 6 -> 2 x Tz of state ( p/n=+/-1 )
				
      integer i,j,n1,n2,j1,j2,mi,mf,l1,l2
      real m1,m2
      integer deltaKr
      real t1,t2

      do i=1,n_states
         n1=spsqn(1,i,2)
         j1=spsqn(1,i,4)
         m1=spsqn(1,i,5)/2.0
         mi=spsqn(1,i,5)
         l1=spsqn(1,i,3)
         do j=1,n_states
            n2=spsqn(1,j,2)
            j2=spsqn(1,j,4)
            m2=spsqn(1,j,5)/2.0
            mf=spsqn(1,j,5)
            l2=spsqn(1,j,3)
            t1=0.25*j2*(j2+2)-m2*(m2+1)
            if(t1.gt.0.0)then
               t1=sqrt(t1)
            else
               t1=0.0
            endif
            t2=0.25*j2*(j2+2)-m2*(m2-1)
            if(t2.gt.0.0)then
               t2=sqrt(t2)
            else
               t2=0.0
            endif
	    Jx_matrix(i,j)=deltakr(n1,n2)*deltakr(j1,j2)*
     &        deltakr(l1,l2)*
     &        (t1*deltakr(mf,mi-2)+t2*deltakr(mf,mi+2))/2.
	    Jy_matrix(i,j)=deltakr(n1,n2)*deltakr(j1,j2)*     !*(-I)
     &         deltakr(l1,l2)*
     &        (t1*deltakr(mf,mi-2)-t2*deltakr(mf,mi+2))/2.
         enddo
      enddo
      return
      end subroutine make_Jxy
C end subroutine make_Jxy..............................................

      subroutine get_OB2(OB,rho,O2_av,n_states)
C===========================================================
C Calculate <psi|O^2|psi>, where O is an one-body operator
C===========================================================
      implicit none
C===========================================================
C INPUT:
C      OB: one-body operator matrix
C      rho: density matrix
C      n_states: the number of states
C OUTPUT:
C      O2_av: the average <psi|O^2|psi>
C===========================================================
C----INPUT-------
      integer,intent(IN) :: n_states
      real*8,intent(IN)  :: rho(n_states,n_states)
      real,intent(IN)    :: OB(n_states,n_states)

C----OUTPUT------
      real,intent(OUT)   :: O2_av

C----Temporary variables and dummies
      real,allocatable   :: OB2(:,:)
      real*8,allocatable :: tmp(:,:)
      real t1,t2
      real*8 t3
      integer i,j,k

      allocate(OB2(n_states,n_states),tmp(n_states,n_states))

C------Calculate <PSI| OB^2 | PSI>
      ob2=matmul(ob,ob)

C      call matrix_mult(OB,OB,OB2,n_states,n_states,n_states)
      call compute_Jz(OB2,rho,n_states,t1)

C------Calculate [Tr(O*rho)]^2-------------------------------
      call  compute_Jz(OB,rho,n_states,t2)
      t2=t2*t2

C------Calculate Tr(O*rho*O*rho)-------------------------------
      do i=1,n_states
         do j=1,n_states
            tmp(i,j)=0.0
            do k=1,n_states
               tmp(i,j)=tmp(i,j)+OB(i,k)*rho(k,j)
            enddo
         enddo
      enddo

      t3=0.0
      do i=1,n_states
         do j=1,n_states
            t3=t3+tmp(i,j)*tmp(j,i)
         enddo
      enddo

      o2_av=t1+t2-sngl(t3)

      deallocate(tmp,OB2)

      return
      end
C....End subroutine


      subroutine calculateJ2(Jx_matrix,Jy_matrix,Jz_matrix,nsps,
     &                        rhopij,rhonij,Jx2,jy2,jz2)
C======================================================================
C Calculate J^2 and also Jx^2,Jy^2,Jz^2
C======================================================================
      implicit none
C======================================================================
C INPUT: Jx_matrix,Jy_matrix,Jz_matrix... : the matrices corresponding to
C                                           Jx,Jy,Jz
C        nsps(2): the number of states
C        rhopij/rhonij: density matrix corresponding to protons/neutrons
C
C OUTPUT: J2: the value for J^2
C========================================================================
C      include 'gcb_dim.inc'
C----INPUT-----
      integer,intent(IN)  :: nsps(2)
      real Jx_matrix(nsps(1),nsps(1)),Jy_matrix(nsps(1),nsps(1)),
     &     Jz_matrix(nsps(1),nsps(1))
      real*8 rhopij(nsps(1),nsps(1)),rhonij(nsps(1),nsps(1))
C----OUTPUT----
      real Jx2,Jy2,Jz2
C----Temporary Variables 
      real jx2p_av,jx2n_av
      real jxp_av,jxn_av

      real jy2p_av,jy2n_av
      real jyp_av,jyn_av

      real jz2p_av,jz2n_av
      real jzp_av,jzn_av

C......Jx^2
      call get_OB2(Jx_matrix,rhopij,Jx2p_av,nsps(1))
      call get_OB2(Jx_matrix,rhonij,Jx2n_av,nsps(1))
      call compute_Jz(Jx_matrix,rhopij,nsps(1),jxp_av)
      call compute_Jz(Jx_matrix,rhonij,nsps(1),jxn_av)
      jx2=(Jx2p_av+Jx2n_av+2*jxp_av*jxn_av)

C......Jy^2
      call get_OB2(Jy_matrix,rhopij,Jy2p_av,nsps(1))
      call get_OB2(Jy_matrix,rhonij,Jy2n_av,nsps(1))
      call compute_Jz(jy_matrix,rhopij,nsps(1),jyp_av)
      call compute_Jz(jy_matrix,rhonij,nsps(1),jyn_av)
      jy2=-(Jy2p_av+Jy2n_av+2.*jyp_av*jyn_av)           ! -1=I*I

C......Jz^2
      call get_OB2(Jz_matrix,rhopij,Jz2p_av,nsps(1))
      call get_OB2(Jz_matrix,rhonij,Jz2n_av,nsps(1))
      call compute_Jz(jz_matrix,rhopij,nsps(1),jzp_av)
      call compute_Jz(jz_matrix,rhonij,nsps(1),jzn_av)
      jz2=Jz2p_av+Jz2n_av+2.*jzp_av*jzn_av
C      write(24,*)Jxp_av+Jyp_av,Jzp_av+Jzp_av

      return
      end
C....End subroutine calculateJ2

      subroutine compute_Jz(jz_matrix,rho,n_states,jz_av)
C====================================================================
C Calculates <Jz>=Tr(rho*Jz)
C===================================================================
      implicit none
C==========================================================================
C Subroutines Called: none
C Called by: the main program
C
C INPUT:
C
C OUTPUT:
C  
C==========================================================================
      integer,intent(IN)    :: n_states
      real,intent(IN)       :: Jz_matrix(n_states,n_states)
      real*8,intent(IN)     :: rho(n_states,n_states)
      real jz_av
      
      integer i,j
      real*8 sum
      
      sum=0.0
      do i=1,n_states
         do j=1,n_states
	    sum=sum+rho(i,j)*Jz_matrix(j,i)
	 enddo
      enddo
      
      jz_av=sngl(sum)
      
      return
      end

      subroutine compute_trace(a,rho,n_states,trace)
C====================================================================
C Calculates trace=Tr(a*rho)
C===================================================================
      implicit none
C==========================================================================
C Subroutines Called: none
C Called by: the main program
C
C INPUT:
C
C OUTPUT:
C  
C==========================================================================
      integer,intent(IN)     :: n_states
      real,intent(IN)        :: a(n_states,n_states)
      real*8,intent(IN)      :: rho(n_states,n_states)

      real,intent(OUT)       :: trace
      
      integer i,j
      real*8 sum
      
      sum=0.0d0
      do i=1,n_states
         do j=1,n_states
	    sum=sum+rho(i,j)*a(j,i)
	 enddo
      enddo
      
      trace=sngl(sum)
      
      return
      end
      
      subroutine exp_Z(Z,n_states,lambda)
C==========================================================================
C Calculates exp(Ztilde) 
C Warning: destrois the Z matrix, and replace it by exp(Ztilde)
C
C==========================================================================
      implicit none
C==========================================================================
C Subroutiones called:
C Called by:
C
C INPUT: 
C              Z: Thouless particle-hole matrix 
C              lambda: real parameter
C              n_states:the number of states
C
C OUTPUT:
C==========================================================================
C
C
C      include 'gcb_dim.inc'
      integer,intent(IN)      :: n_states
      real,intent(INOUT)      :: z(n_states,n_states)
      real,allocatable        :: tmp(:,:)
      real,allocatable        :: x(:,:)
                               !eigenvectors corresponding to Z^2
                               !x_i and x_(i+1) correspond to the same
                               !eigenvalue
      real,allocatable        :: work(:)
      real,allocatable        :: lambda_Z2(:)  ! eigenvalues corresponding to Z^2 
      real,intent(IN)         :: lambda
      
      integer i,j,m    ! dummies
      real temp
      real small
      real,allocatable        :: vec(:)
      integer,allocatable     :: map0(:)   ! maps the zero eigenvalues
      integer count0          ! # of zero eigenvalues
      integer,allocatable     :: map(:)    ! maps the non-zero eigenvalues
      integer count           ! # of nonzero eigenvalues
      integer,allocatable     :: map2(:,:) ! maps the degenerated eigenvectors
                              ! for the nonzero eigenvalues
      real,allocatable        :: w(:)         ! non-zero eigenvalues
      integer m1,m2
      real maxim
      integer maxim_i
      integer one
      
      small=1.0e-4
      one=1
      maxim=0.0
      

      allocate(map0(n_states),map(n_states),map2(n_states,2))
      allocate(tmp(n_states,n_states),lambda_Z2(n_states))
      allocate(x(n_states,n_states),work(n_states),w(n_states))      
      allocate(vec(n_states))
C Calculate first Z*Z
C      call matrix_mult(Z,Z,tmp,n_states,n_states,n_states)
      tmp=matmul(Z,Z)

C and find the eigenvectors
      call eig(tmp,n_states,n_states,lambda_Z2,x,work)

C Calculate sqrt of each eigenvalue
      count0=0
      count =0
      do i=1,n_states
         if(abs(lambda_Z2(i)).lt.small)then
            lambda_Z2(i)=0.0
            count0=count0+1
            map0(count0)=i
         else
            count=count+1
            map(count)=i
         endif
         if(lambda_Z2(i).gt.0.)then
            write(6,*)'eigenvalue=',lambda_Z2(i)
            stop 'Error in exp_Z, eigenvalue>0'
         endif
         lambda_Z2(i)=sqrt(-lambda_z2(i))
      enddo
      
      do i=1,count
         w(i)=lambda_z2(map(i))
      enddo
       
      call check_pairs(w,map2,count,n_states)
       
      do m=1,count/2
         do i=1,n_states
          vec(i)=0.0
          do j=1,n_states
            vec(i)=vec(i)+Z(i,j)*x(j,map((map2(m,1))))
          enddo
          if(abs(maxim).lt.abs(vec(i)))then
             maxim=vec(i)
             maxim_i=i
          endif
         enddo
         
         if(maxim*x(maxim_i,map(map2(m,2))).gt.0.0)then
            do i=1,n_states
               x(i,map(map2(m,2)))=-x(i,map(map2(m,2)))
            enddo
         endif 
2        continue
      enddo    
      
C      write(*,*)count,count0
C      write(16,*)(lambda_z2(i),i=1,n_states)
c      do i=1,count/2
c        write(16,*)map(map2(i,1)),map(map2(i,2))
c      enddo
c      write(16,*)
c      stop  
      
      do i=1,n_states
         do j=1,n_states
            z(i,j)=0.0
            do m=1,count/2
               m1=map(map2(m,1))
               m2=map(map2(m,2))
               temp=x(i,m1)*x(j,m1)+x(i,m2)*x(j,m2)
               z(i,j)=z(i,j)+temp*cos(lambda*lambda_Z2(m1))
               temp=x(i,m1)*x(j,m2)-x(i,m2)*x(j,m1)
               z(i,j)=z(i,j)+temp*sin(lambda*lambda_z2(m1))
            enddo
            do m=1,count0
               z(i,j)=z(i,j)+x(i,map0(m))*x(j,map0(m))
            enddo
         enddo
      enddo
      deallocate(tmp,w,lambda_Z2,x,map,map0,map2,work,vec)

      return
      end subroutine exp_Z
C...........End subroutine exp_Z.......................................

      subroutine check_pairs(w,map2,n,np)
      implicit none
      
C----INPUT---
      integer np,n
      real w(np)
      
C----OUTPUT----
      integer map2(np,2)
      
      integer i,j,k
      integer kk,flag
      integer count
      real small
      parameter(small=1.E-4)
      
      do i=1,n
        do j=1,2
          map2(i,j)=0
        enddo
      enddo
      
      count=0
      do k=1,n-1
         i=k
         do kk=1,count
            if(i.eq.map2(kk,1) .or. i.eq.map2(kk,2))goto 1
         enddo
         j=i+1
2        continue
         flag=1
         do kk=1,count
            if(flag.eq.-1) goto 3
            if(j.eq.map2(kk,1) .or. j.eq.map2(kk,2))then
               j=j+1
               flag=-1
            endif 
3           continue
         enddo
         if(flag.eq.-1)goto 2
         
         do while(j .le. n)
           if(abs(w(i)-w(j)) .lt. small)then
              count=count+1
              map2(count,1)=i
              map2(count,2)=j
              j=n+1
           else
              j=j+1
           endif
         enddo
1        continue
      enddo
c      write(16,*)count
      
      return
      end
      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      

      subroutine writeoutsd(nsps,n,ifile,psi) 

      implicit none
C      include 'gcb_dim.inc'

      
      integer,intent(IN)   :: nsps,n
      real,intent(IN)      :: psi(nsps,n) ! slater determinant
      integer i,j
      integer ifile
      
      do i = 1,n
        write(ifile)(psi(j,i),j=1,nsps)
      enddo
100   format(20(F12.6,1x))
      return
      end
      

      subroutine readinsd(nsps,n,ifile,psi,errflag)

      implicit none
c      include 'gcb_dim.inc'


      integer,intent(IN)      :: nsps,n
      real,intent(OUT)        :: psi(nsps,n) ! slater determinant
      integer i,j
      integer ifile
      logical errflag

      errflag=.false.

      do i = 1,n
        read(ifile,err=103,end=103)(psi(j,i),j=1,nsps)
C        write(23,*)(psi(j,i),j=1,nsps)
      enddo
      return
  103 continue
      errflag=.true.
      return

      end   
      	
       
     
      subroutine matrix_mult(a,b,ab,n1,n2,n3)
C======================================================================
C  Multiplies 2 matrices a(n1,n2)*b(n2,n3) with maximum dimesions max
C======================================================================
      implicit none
C======================================================================
C  Subroutine called: none
C  Called by: make_hmatrix
C             findZ
C  
C  INPUT: a,b matrices to calculate the product
C  OUTUP: ab  the product matrix
C======================================================================
C      include 'gcb_dim.inc'
      integer n1,n2,n3
      real a(n1,n2),b(n2,n3)
      real ab(n1,n3)
      
      integer i,j,k

C      if(n1.gt.max_sps)then
C         print *,'Dimension exceeded in matrix_mult (n1)'
C	 print *,'Program terminated'
C	 stop
C      endif
      
C      if(n2.gt.max_sps)then
C         print *,'Dimension exceeded in matrix_mult (n2)'
C	 print *,'Program terminated'
C	 stop
C      endif
      
C      if(n3.gt.max_sps)then
C         print *,'Dimension exceeded in matrix_mult (n3)'
C	 print *,'Program terminated'
C	 stop
C      endif
      
      do i=1,n1
         do j=1,n3
	   ab(i,j)=0.0
	   do k=1,n2
	      ab(i,j)=ab(i,j)+a(i,k)*b(k,j)
	   enddo
	 enddo
      enddo
      return
      end


















