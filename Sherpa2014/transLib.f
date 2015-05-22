C
C transLib.f
C   This package contains routines for calculating Sum rules
C   for transition operators.
C IS 11/01/2001 @ LSU
C
      subroutine main_trans
C====================================================================
C This is the main subroutine for calculating sum rules for transition
C operators.
C The matrix elements in the s.p. basis are brought through common
C statement. The quantities calculated by RPA are arguments for
C main_trans
C Modified 1/19/2002: isovector/isoscalar components
C          8/01/2002: isospin included in the one body matrix elements
C====================================================================
      use spspace
      use RPAsolutions
      use RPAmatrices
      use obop
      use phBasis
      use nuclide
      implicit none

      real,pointer                     :: w(:)

C-----The matrix for the transition operator--------------
      real opp_ph(nsps(1),nsps(1)),onn_ph(nsps(2),nsps(2))
                                 ! in the original ph basis

C.....Sum rule................
      real*8 sum1,sum2
      real*8 COM
      real*8 sf(nph-n0)  ! strength function
      real small
      parameter(small=1.E-3)
      real Op(nsps(1),nsps(1)),On(nsps(2),nsps(2))
      real xm1,xm2,xj1,xj2
      real Jtr
      real vs
      real charge
      integer M
      real XM
      real clebr,tj
      real*8 BB,B2
      integer nph_p,nph_n

      real*8 sr1,sr2
      real*8 SHF,Sy,Sz


      integer i,j  ! Dummies
      real tx(nph),ty(nph)

      nph=nphp+nphn
      Jtr=float(JJ)

      write(ITROUT,20)EHF
  20  format(2H# ,'HF energy: ',F12.6)


      w=>w0

      sf=0.0

      SHF=0.0
      Sy=0.0
      Sz=0.0
      B2=0.0
      COM=0.0

      nph_p=np*(nsps(1)-np)
      nph_n=nn*(nsps(2)-nn)

      do M=-JJ,JJ
      XM=float(M)


      do i=1,nsps(1)
         xj1=float(spsqn(1,i,4))/2.
         xm1=float(spsqn(1,i,5))/2.
         do j=1,nsps(1)
            xj2=float(spsqn(1,j,4))/2.
            xm2=float(spsqn(1,j,5))/2.
            Op(i,j)=Opp(i,j)*tj(xj1,Jtr,xj2,-xm1,XM,xm2)
         enddo
      enddo

      do i=1,nsps(2)
         xj1=float(spsqn(2,i,4))/2.
         xm1=float(spsqn(2,i,5))/2.
         do j=1,nsps(2)
            xj2=float(spsqn(2,j,4))/2.
            xm2=float(spsqn(2,j,5))/2.
            On(i,j)=Onn(i,j)*tj(xj1,Jtr,xj2,-xm1,XM,xm2)
         enddo
      enddo

C------Calculate the 1-body part in ph basis-----
      call OB2ph(vecp,op,nsps(1),nsps(1),Opp_ph)
      call OB2ph(vecn,on,nsps(2),nsps(2),Onn_ph)
C... and that's it, we work with 1-body operators for transitions.


C-----------USEFUL ROUTINES FOR CHECKING
C      call commZeroCorr(mapp,mapn,nphp,nphn,opp_ph,onn_ph,com)
C      call totalStr(mapp,mapn,nphp,nphn,X,Y,opp_ph,onn_ph,SHF,Sy,Sz)

      do i=1,nph-n0
C         if(w(i).gt.small)then
         do j=1,nph
            tx(j)=X(j,i)
            ty(j)=Y(j,i)
         enddo
c      subroutine sumOandXY(map,map_size,ni,nf,O_ph,sizeO,X,Y,sum,nXY)
         call sumOandXY(map_prot,nph_p,1,nphp,Opp_ph,nsps(1),
     &                                       X(1,i),Y(1,i),sum1,nph)
         call sumOandXY(map_neutr,nph_n,nphp+1, nph,Onn_ph,nsps(2),
     &                                               tx,ty,sum2,nph)
         sf(i)=sf(i)+(sum1+sum2)**2
         B2=B2+(sum1+sum2)**2
C         endif
      enddo

      enddo

      BB=0.0
      sr1=0.0
      sr2=0.0
      do i=1,nph-n0
c        sf(i)=sf(i)**2
        BB=BB+sf(i)
	sr1=sr1+w(i)*sf(i)
	sr2=sr2+w(i)*w(i)*sf(i)
      enddo

      do i=1,nph-n0
         write(ITROUT,1000)w(i),sf(i)
      enddo
C      write(ITROUT,1002)B2
      write(ITROUT,1011)BB,SHF,Sy,Sz,SHF+Sy+Sz
      write(ITROUT,1001)SR1,SR1+COM  !/2.0
      write(ITROUT,1002)SR2
      SR1=SR1/BB
      write(ITROUT,1003)SR1
      write(ITROUT,1004)dsqrt(SR2/BB - SR1*SR1)

      write(6,*)' total strength = ',bb
      write(6,*)' centroid of strength = ',sr1,' MeV '
      write(6,*)' width of strength = ', dsqrt(SR2/BB - SR1*SR1),' MeV'

 1000 FORMAT('  ',F8.4,2X,2F12.6)
 1001 FORMAT(2H# ,'SR1=  ',F12.6,2H  ,'COM= ',F12.6)
 1011 FORMAT(2H# ,'B  =  ',F12.6,'  SHF=',F7.3,' Ycorr=',F7.3,
     1' Zcorr= ',F7.3,
     1' S0=B= ',F7.3)
 1002 FORMAT(2H# ,'SR2=  ',F12.6)
 1003 FORMAT(2H# ,'CENTR=',F12.6)
 1004 FORMAT(2H# ,'WIDTH=',F12.6)

      return
      end

      subroutine sumOandXY(map,map_size,ni,nf,O_ph,sizeO,X,Y,sum,nXY)
C==============================================================
C Calculates the sum Sum_{mi}(X_{mi}*O_{mi}+Y_{mi}*O_{im})
C==============================================================
      implicit none
C      include 'gcb_dim.inc'
C---Input----
      integer  :: map_size,sizeO,nXY
      integer map(map_size,2)
      integer ni,nf                 ! Initial and final indexes
                                    ! (p: 1,nphp; n: nphp+1,nph)
      real o_ph(sizeO,sizeO)    ! Transition matrix in ph basis
      real X(nXY),Y(nXY)  ! RPA amplitudes

C---OUTPUT---
      real*8 sum

      integer n,k
      integer i,m

      sum=0.0
      do k=ni,nf
         n=k-ni+1
         m=map(n,1)
         i=map(n,2)
         sum=sum+O_ph(m,i)*X(k)+O_ph(i,m)*Y(k) !CWJ notes
      enddo
      return
      end

      subroutine readTransOp
C==================================================================================
C This routine reads the reduced matrix elements written in the format
C used by make_phonon (.ph files), and returns matrices for proton and
C neutron transitions in the s.p. space, each elemnt containing just the
C reduced ME and the sign.
C Modified 1/19/2002
C    Changes in order to be able to calculate the ISOSPIN/ISOVECTOR components
C    for electromagnetic transitions.
C Modified 8/1/2002
C    Included the isospin part in the one-body terms (IS/SDSU)
C==================================================================================
      use spspace
      use obop
      implicit none

C      common/TRANSOPpn/opn

      integer ilast
      character*15 filename
      character*60 flabel
      integer ns
      real jtr,Ttr
      integer :: nrorb1(nsps(1)),jorb1(nsps(1))
      real redME(nsps(1)**2)
      integer st1(nsps(1)**2),st2(nsps(1)**2)
      integer flag(nsps(1)**2)
      integer nME
      integer choice

      integer i
      integer n1,n2,j1,j2
      integer nn1,nn2
      real sgn,sgn0
      real xm1,xm2,xj1,xj2
      integer k,l
      real clebr,tj

C WARNING: same spaces assumed for both protons and neutrons
C      write(6,*)'Choose the type of transition:'
C      write(6,*)'    1)Electromagnetic'
C      write(6,*)'    2)Beta decay'
C      write(6,*)'These options not available yet'
c      read(5,*)choice
c      if(choice.ne.1)

      if(allocated(opp))deallocate(opp)
      if(allocated(onn))deallocate(onn)
      if(allocated(opn))deallocate(opn)

 1000 continue
      rewind(ITROUT)
      rewind(ITROUT-1)
      rewind(itrout+1)
      write(6,*)' Enter Transition operator filename (.ph) '
      read(5,'(a)')filename
      ilast=index(filename,' ')
      if(ilast.ne.0)then
         ilast=ilast-1
      else
         ilast=15
      endif

      open(unit=9,file=filename(1:ilast)//'.ph',status='old',
     & err=1007)
      goto 1008
 1007 continue
      write(6,*)' that file does not exist.'
      goto 1000
 1008 continue

      read(9,'(60a)')flabel
      write(6,*)'  '
      write(6,120)flabel
      write(ITROUT,120)flabel
      write(ITROUT+1,120)flabel
      write(ITROUT-1,120)flabel
 120  format('# Transition Strength for ',60a)
      read(9,*)ns
c      if(ns/=nsps(1))then
c          print*,' The phonon in file',filename(1:ilast),
c     &           ' seems to describe a different s.p. space'
c          close(9)
c          goto 1000
c      endif
      write(6,*)' '

      do i=1,ns
        read(9,*)nrorb1(i),jorb1(i)
        n1=nrorb1(i)
        j1=jorb1(i)
        nn1=1
        do while(n1.ne.spsqn(1,nn1,2) .and. j1.ne.spsqn(1,nn1,4))
           nn1=nn1+1
           if(nn1.gt.nsps(1))then
             write(6,*)'The s.p. states do not match in',
     &        ' readTransOp(transLib).'
             stop
           endif
        enddo
        flag(i)=nn1
      enddo

      read(9,*)jj,tt
      Jtr=real(JJ)
      TTr=real(TT)
      write(ITROUT,100)JJ,TT
      write(ITROUT-1,100)JJ,TT
      write(ITROUT+1,100)JJ,TT
  100 format('# J=',I2,2X,'T=',I2)
      nME=ns*ns  ! The number of reduced matrix elements

      allocate(opp(nsps(1),nsps(1)))
      allocate(onn(nsps(1),nsps(1)))
      allocate(opn(nsps(1),nsps(1)))

C      if(nME.gt.size)stop 'Increase size in readTransOp (transLib)'

      do i=1,nME
         read(9,*)st1(i),st2(i),redME(i)
      enddo
      close(unit=9)

      do i=1,nME
         n1=nrorb1(st1(i))
         j1=jorb1(st1(i))
         xj1=float(j1)/2.0
         n2=nrorb1(st2(i))
         j2=jorb1(st2(i))
         xj2=float(j2)/2.0

         nn1=1
         do while(n1.ne.spsqn(1,nn1,2) .and. j1.ne.spsqn(1,nn1,4))
            nn1=nn1+1
            if(nn1.gt.nsps(1))then
              write(6,*)'The s.p. states do not match in',
     &         ' readTransOp(transLib).'
              stop
            endif
         enddo

         nn2=1
         do while(n1.ne.spsqn(1,nn1,2) .and. j2.ne.spsqn(1,nn2,4))
            nn2=nn2+1
            if(nn2.gt.nsps(1))then
              write(6,*)'The s.p. states do not match in',
     &         ' readTransOp(transLib).'
              stop
            endif
         enddo

         do k=nn1,nn1+j1
           xm1=spsqn(1,k,5)/2.0
           sgn=(-1)**(xj1-xm1)
           do L=nn2,nn2+j2
               opp(k,l)= sgn*redME(i)*tj(.5,ttr,.5,-.5,.0,.5)
               onn(k,l)=-sgn*redME(i)*tj(.5,ttr,.5,.5,.0,-.5)
               opn(k,l)= sgn*redME(i)
           enddo
         enddo
      enddo

      return
      end

      subroutine OB2ph(vec,O_spe,n,np,O_ph)
C=======================================================
C Calculates the 1-body matrix part in ph basis
C========================================================
      implicit none
C========================================================
C  INPUT:
C     vec:  transformation matrix from fundamental to
C           ph basis
C     n:    # of states
C     np:   physical dimension
C
C  OUTPUT:
C     O_ph: transformed matrix
C
C  Subroutines called: none
C  Called by subroutine OppCorr in this package
C=========================================================

C-------INPUT-------
      integer n,np
      real vec(np,np),O_spe(np,np)

C-------OUTPUT------
      real O_ph(np,np)

      integer i,j,ii,jj

      do i=1,n
        do j=1,n
           O_ph(i,j)=0.0
           do ii=1,n
             do jj=1,n
                 O_ph(i,j)=O_ph(i,j)+
     &               vec(ii,i)*O_spe(ii,jj)*vec(jj,j)
             enddo
           enddo
        enddo
      enddo

      return
      end

