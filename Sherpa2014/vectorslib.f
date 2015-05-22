CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Package vectorslib
C Handles different vector operations:
C          orth_GH: orthogonalizes linearly independent vectors by
C                   Gram-Schmidt procedure
C          sc_prod: calculates the scalar product of 2 vectors
C          normOne: normalizes a vector
C
C IS, LSU 9/15/00
C
C IS, LSU 8/25/01
C * added subroutine orth_SVD which orthonormalizes through SDV decompozition
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine orth_GS(a,nvect,ncomp)
C============================================================================
C Orthogonalizes a set of nvect vectors, and normalizes them
C============================================================================
      implicit none
C============================================================================
C Subroutines called: sc_prod, normOne
C Called by:          none in this library
C
C INPUT:   a=matrix, the vectors to be orthogonalised are on the columns
C          max=maximum dimension of vectors
C          nvect= the number of vectors to be orthogonalized
C          ncomp=the number of components of each vector
C OUTPUT:  a=the matrix containing orthonormalized vectors
C============================================================================
      integer,intent(IN)    :: nvect,ncomp
C      PARAMETER(NCOMP1=20,nvect1=20)
      real a(ncomp,ncomp)

C----------Dummies and Temporar Variables-----------------------------
      integer i,j,n,change
      real lambda
      real,allocatable    :: solver(:,:),b(:,:),temp(:,:)

      allocate(solver(nvect,nvect),b(nvect,1),temp(ncomp,ncomp))

      temp(1:ncomp,1)=a(1:ncomp,1)

      call normOne(temp(1:ncomp,1),ncomp,ncomp)
      call normOne(a(1:ncomp,1),ncomp,ncomp)

      if(nvect.eq.1)return

      temp(1:ncomp,2)=a(1:ncomp,2)

c      do i=1,ncomp
c         temp(i,2)=a(i,2)
c      enddo

      call sc_prod(temp(1:ncomp,1),temp(1:ncomp,2),lambda,ncomp,ncomp)

      if(lambda.ne.0.0)then
         do i=1,ncomp
            temp(i,2)=temp(i,1)-temp(i,2)/lambda
         enddo
      endif

      call normOne(temp(1:ncomp,2),ncomp,ncomp)
      call normOne(a(1:ncomp,2),ncomp,ncomp)

      change=0
      do i=1,nvect-1
         do j=i+1,nvect
	    call sc_prod(a(1:ncomp,i),a(1:ncomp,j),lambda,ncomp,ncomp)
	    if(abs(lambda).gt.1.0E-5)change=change+1
c	    write(*,*)lambda,change
         enddo
      enddo
      if(change.eq.0)return


      call sc_prod(temp(1,1),temp(1,2),lambda,ncomp,ncomp)

      if(nvect.gt.2)then
      do n=3,nvect

	 do i=1,n-1
            do j=2,n
               call sc_prod(a(1:ncomp,j),temp(1:ncomp,i),lambda,ncomp,
     &                                                           ncomp)
               solver(i,j-1)=lambda
	    enddo
	 enddo

	 b(1,1)=-1.0
         do i=2,n-1
            b(i,1)=0.0
         enddo

c         call gaussp(solver,n-1,tmp,b)
 	 call gaussj(solver,n-1,nvect,b,1,1)

	 do i=1,ncomp
	    temp(i,n)=temp(i,1)
	    do j=2,n
	       temp(i,n)=temp(i,n)+b(j-1,1)*a(i,j)
	    enddo
	 enddo

         call normOne(temp(1:ncomp,n),ncomp,ncomp)
      enddo
      endif

      do i=1,ncomp
         do j=1,nvect
            a(i,j)=temp(i,j)
	 enddo
      enddo
      deallocate(temp,solver,b)
c      stop

1     format(3(1X,g12.5,1X))
      return
      end
C.....................

      subroutine sc_prod(vec1,vec2,lambda,ncomp,np)
C=====================================================================
C Calculates the scalar product of 2 vectors
C=====================================================================
      implicit none
C=====================================================================
C Subroutines Called: none
C Called by: normOne, orth_GS
C
C INPUT:    vect1, vect2 the to vectors
C           ncomp: the number of components of each vector
C OUTPUT:   lambda: the value of the scalar product
C
C=====================================================================
      integer ncomp
      integer np
      real lambda,vec1(np),vec2(np)

      integer i

      if(np.lt.ncomp)stop 'Problem in normOne'
      lambda=0.0
      do i=1,ncomp
         lambda=lambda+vec1(i)*vec2(i)
      enddo
      return
      end
C.......................................................................
      subroutine normOne(vec,ncomp,np)
C======================================================================
C Normalizes a vector vec
C======================================================================
      implicit none
C======================================================================
C  Subroutines Called: sc_prod
C  Called by:          orth_GS
C
C  INPUT:   vec=the vector to be normalized
C           ncomp=the numcer of components of the vector
C  OUTPUT:  vec=the orthonormalized vector
C
C======================================================================
      integer ncomp
      integer np
      real vec(np),norm

      integer i

      call sc_prod(vec,vec,norm,ncomp,np)
      norm=sqrt(norm)
      do i=1,ncomp
         vec(i)=vec(i)/norm
      enddo
      return
      end subroutine normOne
C...End subroutine normOne(vec,ncomp)


      subroutine orth_SVD(a,ncomp,nvect)
C============================================================================
C Orthogonalizes a set of nvect vectors, and normalizes them.
C One uses for this the SVD decomposition which is more stable numerically than
C the C Gram-Schmidt algorithm.
C============================================================================
      implicit none
C============================================================================
C Subroutines called: svdcmp, normOne
C Called by:          none in this library
C
C INPUT:   a=matrix, the vectors to be orthogonalised are on the columns
C          max=maximum dimension of vectors
C          nvect= the number of vectors to be orthogonalized
C          ncomp=the number of components of each vector
C OUTPUT:  a=the matrix containing orthonormalized vectors
C============================================================================
      integer max,nvect,ncomp
      real a(ncomp,nvect)

C----------Dummies and Temporar Variables-----------------------------
C      integer ncomp1,nvect1
      integer i
      real,allocatable  :: v(:,:),w(:)

      allocate(v(ncomp,ncomp),w(ncomp))

      call svdcmp(a,ncomp,nvect,ncomp,nvect,w,v)
      do i=1,nvect
         call normOne2(a(1,i),ncomp,ncomp)
      enddo

      deallocate(v,w)

      return
      end

      subroutine normOne2(v,n,np)

      implicit none
      integer n,np
      real v(np)

      real*8 norm
      integer i

      norm=0.0
      do i=1,n
         norm=norm+v(i)*v(i)
      enddo
c      write(65,*)norm

      norm=dsqrt(norm)

      do i=1,n
         v(i)=sngl(v(i)/norm)
      enddo

      return
      end





