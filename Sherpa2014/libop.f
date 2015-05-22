C===================================================================
C Package libop.f
C    used to calculate the RPA corrections to the one- and two- body
C    operators
C===================================================================
      subroutine TwoBodyA11(O1,O2,map,isign,nph,nsps,vec,A,B)
C===================================================================
C
C Calculates the matrices A and B associated with the operator
C O = Sum_{ab}(O1_{ab}c^+_a c_b) + Sum_{abcd}(O2_{abcd}c^+_a c^+_b c_d c_a)
C
C====================================================================
      implicit none
C====================================================================
C
C  INPUT:
C     O1: the 'one'-body matrix
C     O2: the two-body matrix
C     map:       map(k,1) the particle state corresponding to the k-th 
C                         ph pair
C                map(k,2) the hole state corresponding to the k-th 
C                         ph pair
C     nph:  # of particle-hole states
C     nsps: # of s.p. states
C     vec:    transformation matrix from fundamental to ph basis
C
C  OUTPUT:
C     A,B: particle-hole matrices associated with O
C
C====================================================================
C      include 'gcb_dim.inc'

C---INPUT---
      integer nsps,nph
      real o1(nsps,nsps),o2(nsps,nsps)
C      real*8 rho(max_sps,max_sps),q(max_sps,max_sps)
      integer map(nph,2)
      integer isign    ! +1 if Jx or Jz
                       ! -1 if Jy
      real vec(nsps,nsps)

C---OUTPUT--
      real a(nph,nph),b(nph,nph)

C---INTERNAL VARIABLES---
      integer deltakr   ! Delta Kroneker

      integer i,j,k,m,l,n,ii,jj,kk,ll
      real*8 tmp1,tmp2,tmp3,tmp4


      do k=1,nph
        m=map(k,1)
        i=map(k,2)
        do l=1,nph
           n=map(l,1)
           j=map(l,2)
           tmp1=0.0
           tmp2=0.0
           tmp3=0.0
           tmp4=0.0
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP& PRIVATE(ii,jj,kk,ll) REDUCTION(+:tmp1,tmp2,tmp3,tmp4)
!$OMP DO
           do ii=1,nsps
             do jj=1,nsps
              tmp1=tmp1+vec(ii,n)*O1(ii,jj)*vec(jj,m)
              tmp2=tmp2+vec(ii,i)*O1(ii,jj)*vec(jj,j)
               do kk=1,nsps
                 do ll=1,nsps
                     tmp3=tmp3+vec(ii,m)*vec(jj,j)*vec(kk,i)*vec(ll,n)*
     &                         (O2(ii,ll)*O2(jj,kk)-O2(ii,kk)*O2(jj,ll))
                     tmp4=tmp4+vec(ii,m)*vec(jj,n)*vec(kk,i)*vec(ll,j)*
     &                         (O2(ii,ll)*O2(jj,kk)-O2(ii,kk)*O2(jj,ll))
                 enddo
               enddo
             enddo
           enddo
!$OMP END DO
!$OMP END PARALLEL
c           write(*,*)k,l,'test'
           a(k,l)=a(k,l)+deltakr(i,j)*tmp1*0.-deltakr(m,n)*tmp2*0.
     &                                        -2.*isign*tmp3
           b(k,l)=b(k,l)-2.*isign*tmp4
        enddo
      enddo
c
c      do k=1,nph
c         write(89,'(100(F10.6,1X))')(a(k,l),l=1,nph)
c      enddo
c      write(89,*)

      return
      end

      subroutine make_gammaO(gamma,Op,On,rhop,rhon,nsps,isign)
C===================================================================
C
C Calculates the matrix Gamma associated with the operator O
C It has a pp(or nn) part and pn part.
C
C===================================================================
      implicit none
C      include 'gcb_dim.inc'

C---INPUT---
      integer nsps(2)
      real op(nsps(1),nsps(1)),on(nsps(2),nsps(2))
      real*8 rhop(nsps(1),nsps(1)),rhon(nsps(2),nsps(2))
      integer isign

C---OUTPUT
      real gamma(nsps(1),nsps(1))

C---Interbal Variables---

      integer i,j,ii,jj,kk
C      integer O2(nsps(1),nsps(1))
      real*8 tmp1

C      O2=matmul(Op,On)

C      call matrix_mult(Op,On,O2,max_sps,nsps(1),nsps(1),nsps(1))

c      do i=1,nsps(1)
c         write(90,'(20(F10.5,2x))')(o2(i,ii),ii=1,nsps(1))
c         write(91,'(20(F10.5,2x))')(op(i,j),j=1,nsps(1))
c      enddo
c      write(90,*)
c      write(91,*)

      do i=1,nsps(1)
         do j=1,nsps(1)
            tmp1=0.0
            do ii=1,nsps(1)
              do jj=1,nsps(1)
                tmp1=tmp1+rhop(ii,jj)*
     &                (Op(i,ii)*Op(jj,j)-Op(i,j)*Op(jj,ii))
C     &                  Op(i,j)*Op(jj,ii)*2.
                enddo
            enddo
            gamma(i,j)=gamma(i,j)-sngl(tmp1)*isign
         enddo
      enddo

      do i=1,nsps(1)
         do j=1,nsps(1)
            do ii=1,nsps(2)
            do kk=1,nsps(2)
               gamma(i,j)=gamma(i,j)+isign*Op(i,j)*2
     &                        *On(ii,kk)*rhon(kk,ii)
            enddo
            enddo
         enddo
      enddo

      return
      end
