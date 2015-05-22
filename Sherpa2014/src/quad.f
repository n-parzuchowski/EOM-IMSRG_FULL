C Modified 5/8/2002, INT--Seattle
C Change some numerical factors so that it takes into account
C the core nucleons.
C 5/9/2002 Changed again to the correct factors
C
C May 2004 Fortran 90 version
C=======================================================================

       subroutine deform(beta,gamma) !,qp,qn)
C==========================================================================
C   Calculates the angular part of the matrix element for the
C   quadrupole operator
C   Computes the deformation parameters beta and gamma
C==========================================================================
      use SPspace
      use wfn
      implicit none
C==========================================================================
C INPUT: 
C
C OUTPUT:
C       beta
C       gamma   deformation parameters
C
C==========================================================================
c      include 'gcb_dim.inc'



      real,intent(OUT)       :: beta,gamma

      real,allocatable       ::  q0(:,:,:)   ! there are 6 independent matrix elements
                                             ! since the matrix is symmetric
      integer n_states

      real qp(3,3),qn(3,3),q(3,3)
      real q_diag(3)
      real q1,q2
      real vecq(3,3),work(3)

      real trace,r2trace

      integer i,j,m,l1,l2,n1,n2,ij1,ij2,im1
      integer map2(6,2)
      real m1,m2,f1,f2,f3,f4,f5,f6,fact,sign,delta
      real j1,j2
      real,allocatable,dimension(:,:)  :: tmp,r2
      real two
      integer itwo
      real pi,eps

      real yredME,yredME1,y2

      real tj, HORadInt

      two = 2.0
      itwo=2
      eps=1.E-4
      pi=acos(-1.)
      n_states=nsps(1)

      allocate(r2(n_states,n_states),q0(n_states,n_states,6))

      do i=1,n_states
        l1=spsqn(1,i,3)
        ij1=spsqn(1,i,4)
        j1=ij1/2.0
        n1=spsqn(1,i,2)
        m1=spsqn(1,i,5)/2.0
        im1=spsqn(1,i,5)
        sign=(-1)**((ij1-im1)/2)
        do j=1,n_states
          l2=spsqn(1,j,3)
          n2=spsqn(1,j,2)
          ij2=spsqn(1,j,4)
          j2=ij2/2.0
          m2=spsqn(1,j,5)/2.0

          delta=0.0
          if((m1.eq.m2) .and. (ij1.eq.ij2) 
     &                 .and. (l1.eq.l2))delta=1.0
     
          f1=sign*tj(j1,2.0,j2,-m1,-2.,m2)
          f1=f1*yredme(2,l1,ij1,l2,ij2)/2.

          f2=sign*tj(j1,2.0,j2,-m1,-1.,m2)
          f2=f2*yredme(2,l1,ij1,l2,ij2)/2.

          f3=sign*tj(j1,2.0,j2,-m1, 0.,m2)
          f3=f3*yredme(2,l1,ij1,l2,ij2)/2.

          f4=sign*tj(j1,2.0,j2,-m1,+1.,m2)
          f4=f4*yredme(2,l1,ij1,l2,ij2)/2.

          f5=sign*tj(j1,2.0,j2,-m1,+2.,m2)
          f5=f5*yredme(2,l1,ij1,l2,ij2)/2.

          f6=HORadInt(2,n1,l1,n2,l2)
          r2(i,j)=f6*delta

          q0(i,j,1)= f6*((f1+f5)*sqrt(6./5.)-2.*f3/sqrt(5.))      ! xx
          q0(i,j,2)=-f6*((f1+f5)*sqrt(6./5.)+2.*f3/sqrt(5.))      ! yy
          q0(i,j,3)= 4.*f6*f3/sqrt(5.)                            ! zz
          q0(i,j,4)= f6*sqrt(6./5.)*(f1-f5) ! *I (purely imaginary) xy
          q0(i,j,5)=-f6*sqrt(6./5.)*(f4-f2)                       ! xz
          q0(i,j,6)= f6*sqrt(6./5.)*(f4+f2) ! *I (purely imaginary) yz
        enddo
      enddo
C      stop

      do m=1,3
        map2(m,1)=m
        map2(m,2)=m
      enddo
      map2(4,1)=1
      map2(4,2)=2
      map2(5,1)=1
      map2(5,2)=3
      map2(6,1)=2
      map2(6,2)=3

      allocate(tmp(n_states,n_states))
      do m=1,6
         do i=1,n_states
           do j=1,n_states
             tmp(i,j)=q0(i,j,m)
           enddo
         enddo
         call compute_trace(tmp,rhop,n_states,trace)
         qp(map2(m,1),map2(m,2))=trace
         qp(map2(m,2),map2(m,1))=trace
         call compute_trace(tmp,rhon,n_states,trace)
         qn(map2(m,1),map2(m,2))=trace
         qn(map2(m,2),map2(m,1))=trace
      enddo

      call compute_trace(r2,rhop,n_states,trace)
      call compute_trace(r2,rhon,n_states,r2trace)
      r2trace=r2trace+trace  !+r20

C Set the imaginary elements to zero

      if(abs(qp(1,2)).gt.eps .or. abs(qp(2,3)).gt.eps .or.
     &      abs(qn(1,2)).gt.eps .or. abs(qn(2,3)).gt.eps)then
          write(*,*)qp(1,2),qp(2,3),qn(1,2),qn(2,3)
          pause 'Imaginary ME for Q'
      endif

      qp(1,2)=0.0
      qp(2,1)=0.0
      qp(2,3)=0.0
      qp(3,2)=0.0

      qn(1,2)=0.0
      qn(2,1)=0.0
      qn(2,3)=0.0
      qn(3,2)=0.0

      do i=1,3
         do j=1,3
            q(i,j)=qp(i,j)+qn(i,j)
         enddo
      enddo
C      do i=1,3
C          write(82,192)(q(i,j),j=1,3)
C      enddo
C      write(82,*)
192   format(3(F12.6,1X))
      call eigval(q,3,3,q_diag,vecq,work)
      call eigsrt(q_diag,vecq,3,3)

      call intrinsic(q_diag,q1,q2)

      gamma=atan(q2*sqrt(2.)/q1)*180./pi
      if(gamma.lt.0.0)gamma=-gamma
c      beta=q1*cos(gamma)+q2*sqrt(2.)*sin(gamma)
c      beta=sqrt(5.*pi)*beta/3.

      beta=0.0
      do i=1,3
         beta=beta+q_diag(i)**2
      enddo
      beta=5.*sqrt(beta/6.)/6.  ! 5/8/2002 Why?
C      beta=5.*sqrt(beta)/4./pi   ! The original factor in Eq. 27
C                                 ! (PRC 49, 1442) has been changed to
C                                 ! 5/(4*Pi) to account for core nucleons

C      beta=beta*2*pi*sqrt(5./4./Pi)/3.
C      beta=4*Pi*beta/5.
      beta=beta/r2trace

      return
      end

      subroutine intrinsic(q_diag,q1,q2)
C===================================================================
      implicit none
C===================================================================
C Calculates the intrinsic values for the quadrupole operator
C as defined in PRC 61, 034303 (2000)
C===================================================================

      real q_diag(3),q1,q2

      real pi

      pi=acos(-1.)

      q1=q_diag(3)*sqrt(5./pi)/4.

      q2=sqrt(5./(2.*pi))*q_diag(1)+q1*sqrt(2.)
      q2=q2/(2.*sqrt(3.))

      return
      end

      real function YredME(L,l1,j1,l2,j2)
C==================================================================
C Calculates the reduced matrix element 
C              <(l1 1/2) j1 || Y_L || (l2 1/2) j2> * sqrt(4*pi)
C==================================================================

      integer j1,j2
      integer L,l1,l2

      real tj,sj,xj1,xj2

      real fact,ll1,ll2,LL

      LL=float(L)
      ll1=float(l1)
      ll2=float(l2)
      xj1=j1/2.
      xj2=j2/2.

      fact=(-1)**(L+(j2+1)/2)
      fact=fact*sqrt((2.*ll1+1)*(2.*ll2+1.)*(2.*LL+1.)
     &              *(j1+1)*(j2+1))

      yredme=sj(ll1,xj1,0.5,xj2,ll2,LL)*tj(ll1,LL,ll2,0.0,0.0,0.0)
      yredME=fact*yredME

C      write(*,*)sixj(ll1,j1,0.5,xj2,ll2,LL)
      return
      end

      real function YredME1(L,l1,j1,l2,j2)
C==================================================================
C Calculates the reduced matrix element
C              <(l1 1/2) j1 || Y_L || (l2 1/2) j2> * sqrt(4*pi)
C yet another method
C==================================================================

      real j1,j2
      integer L,l1,l2

      real tj  ! three J

      real fact,ll1,ll2,LL,parity

      LL=float(L)
      ll1=float(l1)
      ll2=float(l2)

      parity=(1.+(-1)**(LL+ll1+ll2))/2.0
      fact=(-1)**(l1+l2+0.5-j2)*parity
      fact=fact*sqrt((2.*j1+1)*(2.*j2+1.)*(2.*LL+1.))

      yredme1=tj(j1,LL,j2,-0.5,0.0,0.5)
      yredME1=fact*yredME1
C      write(*,*)tj(j1,LL,j2,-0.5,0.0,0.5)

      return
      end

CCCCCCCCCCCCCCCCALVINCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function HORadInt(L,n1,l1,n2,l2)
C
C   computes radial integral R(n1,l1, alpha*r) R(n2,l2, alpha*r) r^(L+2) 
C   where R's are h.o. radial functions 
C
C   uses Lawson 1.11a; alpha = 1  must scale by 1/alpha^2 to get final 
C   answer
C
      implicit none 
    
      integer L			! weighting of r ( = L in Ylm)
      integer n1,n2,l1,l2       ! q#'s of states 

      real sum			! result 
      integer q			! dummy for summation
      integer qmin,qmax		! limits of summation 

      real lnprefact,lnsum,lnfact0,ln2fact		

C............let's get going..................
      if( mod(l1+l2+L,2) .ne. 0)then 	  ! must be even overall
	HORadInt = 0.0    
	return
      endif

C............. l1,l2,L must satisfy triangle relation....
C      if( ( l1+l2 .lt. L) .or. (abs(l1-l2) .gt. L) )then
C	HORadInt = 0.0
C	return
C      endif
      lnprefact = (lnfact0(n1)+lnfact0(n2) + log(2.)*(n1+n2-L)
     & - ln2fact(2*n1+2*l1+1) - ln2fact(2*n2+2*l2+1))/2. 
     & + lnfact0( (l2-l1+L)/2) + lnfact0( (l1-l2+L)/2)

      qmax = min(n1,n2)
      qmin = max(0,max( n1-(l2-l1+L)/2,n2-(l1-l2+L)/2))
      sum = 0.0
      do q = qmin,qmax
	lnsum =  ln2fact(l1+l2+L+2*q+1) 
     & -q*log(2.) - lnfact0(q) - lnfact0(n1-q) -lnfact0(n2-q) 
     & - lnfact0(q+(l1-l2+L)/2-n2) - lnfact0(q+(l2-l1+L)/2-n1)

	sum = sum + exp( lnprefact+lnsum)
C       write(20,*)exp( lnprefact+lnsum) 
      enddo

      if( mod( abs(n1-n2),2) .ne.0) sum = -sum
      
C      write(*,*)L,n1,l1,n2,l2,sum
      HORadInt = sum
      return
      end


C==========================================================================
      real function lnfact0(n)

      implicit none

      integer n
      integer i

      lnfact0 = 0.0

      do i = 2,n
	lnfact0 = lnfact0 + log(float(i))
      enddo
      return
      end

C==========================================================================
      real function ln2fact(n)
C
C   log of double factorial ln n!!
C
      implicit none 

      integer n
      integer i

      ln2fact = 0.0
      if(n .eq. 0 .or. n .eq. 1) return 

      do i = n,1,-2
	ln2fact = ln2fact + log(float(i))
      enddo
      return
      end

C==========================================================================

