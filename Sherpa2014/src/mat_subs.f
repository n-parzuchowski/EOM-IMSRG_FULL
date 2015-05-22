CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine make_rhoij(nsps,numpart,psil,psir,
     &                      overlap,rhoij)
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c             Subroutine to compute rho(i,j) for two SD's that are
c             not orthogonal using a matrix representation.
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C  INPUT:
C	nsps:	# of of m-states
C	numpart:# of particles (single species)
C  	psir:	right SD; psil:	left SD
C		NOTE: wavefunctions stored in REAL*8
C
C OUTPUT:
C	overlap:  < psil |  psir>
C	rhoij	:  density matrix < psil | a^\dagger_i a_j | psir >
C
C  SUBROUTINES CALLED:
C	dgetrf:  LAPACK double precision LU routine
C	dgetri:	 LAPACK double precision inversion routine using dgetrf
C
C  FUNCTIONS CALLED:
C
C---------------------------------------------------------------------------

      implicit none
      
C      include 'gcb_dim.inc'

c============ internal variables
      integer i,j,k,num1,num2,iswitch
      real*8 sum
c============ Incoming variables
      integer nsps,numpart
      real*8,intent(IN)    :: psir(nsps,numpart),psil(nsps,numpart)

c============ Resultant - density matrix
      real*8 overlap,phase
      real*8 rhoij(nsps,nsps)
c============ Variables passed onto LAPACK routines dgetrf and dgetri
      integer info
      real*8,allocatable       :: over_matrix(:,:),temp(:,:)
      integer,allocatable      :: ipiv(:)
      real*8,allocatable       :: work(:)
c============ Useful constants
      real*8 zero,one

C==========================================================================
c============ Start program - Initialize constants
      zero=0.0d0
      one=1.0d0

c==============   Add small random part

C      call transform_psi(max_sps,nsps,numpart,tran1,
C     &                         psi_r,psirtemp)
C      call ortho_norm(max_sps,nsps,numpart,psirtemp)
C      call transform_psi(max_sps,nsps,numpart,tran2,
C     &                         psi_l,psiltemp)
C      call ortho_norm(max_sps,nsps,numpart,psiltemp)

c============ Make over_matrix
      allocate(over_matrix(numpart,numpart),ipiv(numpart),work(numpart))
C      if(max_sps.gt.max_dim)stop 'dimensions in overlap exceeded'
      do num1=1,numpart
         do num2=1,numpart
            sum=zero
            do i=1,nsps
               sum=sum+psil(i,num1)*psir(i,num2)
            end do
            over_matrix(num1,num2)=sum
         end do
C         write(6,*)(over_matrix(num1,num2),num2=1,numpart)
      end do
c============ Use LAPACK routines to compute inverse of
c============ over_matrix
c============ First LU decomposition
      call dgetrf(numpart,numpart,over_matrix,numpart,ipiv,info)
c============ Compute overlap - Determinant of over_matrix
c============ Determinant is product of diagonals. Phase determined 
c============ interchanges of orbits in ipv(i).
      overlap=one
      iswitch=0
      do i=1,numpart
         overlap=overlap*over_matrix(i,i)
         if(ipiv(i).ne.i)iswitch=iswitch+1
      end do
      if(iand(iswitch,1).eq.1)overlap=-overlap    ! odd number of interchanges
c============ Get the inverse
      call dgetri(numpart,over_matrix,numpart,ipiv,work,numpart,info)
c============ Construct rhoij matrix
c============ zero rhoij
       do i=1,nsps
          do j=1,nsps
              rhoij(j,i)=zero
              sum=zero
              do num1=1,numpart
                 do num2=1,numpart
                    sum=sum+psir(i,num1)*over_matrix(num1,num2)*
     &                      psil(j,num2)
                 end do
              end do
              rhoij(j,i)=sum
           end do
       end do
      deallocate(over_matrix,ipiv,work)
      return
      end



