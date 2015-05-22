C
C pnTrans.f
C   This package contains routines for calculating Sum rules
C   for transition operators in pn RPA
C IS 6/5/2003 @ LSU
C
      subroutine rpa_transitionpn
C=============================================================================
C Compute the beta decay strengths
C=============================================================================
      use obop
      implicit none
C..............FILE HANDLING..........................

      character*1 ychar
      character filename*15  		! 
      integer ilast

      
      ITROUT=22

      if(allocated(opp))deallocate(opp)
      if(allocated(onn))deallocate(onn)
      if(allocated(opn))deallocate(opn)

      call readTransOp


C............. set up files for output......................

311    continue
      write(6,*)' Enter GT output filename ', 
     & '(.GTnp,.GTpn)'
      read(5,'(a)')filename
      ilast=index(filename,' ')
      if(ilast.ne.0)then
         ilast=ilast-1
      else
         ilast=15
      endif

      open(unit=(ITROUT-1),file=filename(1:ilast)//'.GTnp',
     &status='new',err=33)
      open(unit=(ITROUT+1),file=filename(1:ilast)//'.GTpn',
     &status='new',err=33)
      goto 44
33    continue
      write(6,*)' That file already exists;',
     &' do you want to overwrite (y/n)? '
      read(5,'(a)')ychar
 
      if(ychar.eq.'y' .or.ychar.eq.'Y')then
      open(unit=ITROUT-1,file=filename(1:ilast)//'.GTnp',
     & status='old')
      open(unit=ITROUT+1,file=filename(1:ilast)//'.GTpn',
     & status='old') 
      else
 	goto 311
      endif
44    continue

      call main_transPN
      return
      end subroutine rpa_transitionpn

      subroutine main_transPN
C====================================================================
C This is the main subroutine for calculating sum rules for transition
C operators. This routine calculates beta decays (GT)!!
C The matrix elements in the s.p. basis are brought through common
C statement. The quantities calculated by RPA are arguments for
C main_transPN
C====================================================================
      use phBASIS
      use RPAsolutions
      use RPAmatrices
      use obop
      use spspace
      implicit none

C-----The matrix for the transition operator--------------
      real opn_ph(nsps(1),nsps(1)),onp_ph(nsps(1),nsps(1))
                                 ! in the original ph basis

C.....Sum rule................
      real*8 sum
      real*8 COM
      real*8 sfpn(nphpn)  ! strength function
      real*8 sfnp(nphnp)  ! strength function
      real small
      parameter(small=1.E-5)
      real O_pn(nsps(1),nsps(1)),O_np(nsps(1),nsps(1))
      real xm1,xm2,xj1,xj2
      real Jtr
      real vs
      real charge
      integer M
      real XM
      real clebr,tj
      real*8 B2,BGTpn,BGTnp

      real*8 sr1pn,sr1np
      real*8 SHF,Sy,Sz
      real*8 SHFpn
      real*8 SHFnp
      real*8 sr2pn,sr2np


      integer i,j  ! Dummies
      real tx(nph_pn),ty(nph_pn)


      IF(TT.ne.1)then
         write(ITROUT+1,*)'No beta-decay possible.'
         write(ITROUT-1,*)'No beta-decay possible.'
         return
      endif

C      nph=nphp+nphn
      Jtr=float(JJ)

      write(ITROUT+1,20)EHF
      write(ITROUT-1,20)EHF
  20  format(2H# ,'HF energy: ',F12.6)


      sfpn=0.0
      sfnp=0.0
      SHF=0.0
      Sy=0.0
      Sz=0.0
      B2=0.0
      COM=0.0


      do M=-JJ,JJ
      XM=float(M)


      do i=1,nsps(1)
         xj1=float(spsqn(1,i,4))/2.
         xm1=float(spsqn(1,i,5))/2.
         do j=1,nsps(2)
            xj2=float(spsqn(2,j,4))/2.
            xm2=float(spsqn(2,j,5))/2.
            O_pn(i,j)= Opn(i,j)*tj(xj1,Jtr,xj2,-xm1,XM,xm2)
     &                        *tj(0.5,1.0,0.5,-0.5,1.0,-0.5)
         enddo
      enddo

      do i=1,nsps(2)
         xj1=float(spsqn(2,i,4))/2.
         xm1=float(spsqn(2,i,5))/2.
         do j=1,nsps(1)
            xj2=float(spsqn(1,j,4))/2.
            xm2=float(spsqn(1,j,5))/2.
            O_np(i,j)=-Opn(i,j)*tj(xj1,Jtr,xj2,-xm1,XM,xm2)
     &                        *tj(0.5,1.0,0.5,0.5,-1.0,0.5)
         enddo
      enddo

C------Calculate the 1-body part in ph basis-----
      call OB2phPN(vecp,vecn,o_pn,nsps(1),nsps(2),nsps(1),Opn_ph)
      call OB2phPN(vecn,vecp,o_np,nsps(2),nsps(1),nsps(1),Onp_ph)
C... and that's it, we work with 1-body operators for transitions.

      do i=1,nphnp
C         if(abs(wnp(i)).gt.small)then
         call sumOandXYpn(map_np,map_pn,nphnp,nphpn,Onp_ph,
     1                                   Xnp(1,i),Ynp(1,i),sum,nsps(1))
	 sfnp(i)=sfnp(i)+sum**2
C         endif
      enddo
      do i=1,nphpn
C         if(abs(wpn(i)).gt.small)then
         call sumOandXYpn(map_pn,map_np,nphpn,nphnp,Opn_ph,
     1                                   Xpn(1,i),Ypn(1,i),sum,nsps(2))
	 sfpn(i)=sfpn(i)+sum**2
C         endif
      enddo

      enddo  ! cycle over M

      BGTpn=0.0
      BGTnp=0.0
      sr1pn=0.0
      sr1np=0.0
      sr2pn=0.0
      sr2np=0.0

C............. need to sort here
      call srtval(wpn,sfpn,nphpn,nphpn)
      call srtval(wnp,sfnp,nphnp,nphnp)
      do i=1,nphpn
        BGTpn=BGTpn+sfpn(i)
	sr1pn=sr1pn+wpn(i)*sfpn(i)
	sr2pn=sr2pn+wpn(i)*wpn(i)*sfpn(i)
        if(abs(wpn(i)).gt.small .and. sfpn(i).gt.small)
     &         write(ITROUT+1,1000)wpn(i),sfpn(i)
      enddo
      do i=1,nphnp
        BGTnp=BGTnp+sfnp(i)
	sr1np=sr1np+wnp(i)*sfnp(i)   !**2
	sr2np=sr2np+wnp(i)*wnp(i)*sfnp(i)
        if(abs(wnp(i)).gt.small .and. sfnp(i).gt.small)
     & write(ITROUT-1,1000)wnp(i),sfnp(i)
      enddo
      write(ITROUT+1,1002)BGTpn,BGTpn-BGTnp
      write(ITROUT-1,1002)BGTnp,BGTpn-BGTnp
c      write(6,1002)BGTpn,BGTpn-BGTnp
c      write(6,1002)BGTnp,BGTpn-BGTnp
      write(ITROUT+1,1001)sr1pn
      write(ITROUT-1,1001)sr1np
      write(ITROUT+1,1003)sr2pn,sr1pn/BGTpn,
     1                   dsqrt(sr2pn/BGTpn-sr1pn*sr1pn/BGTpn/BGTpn)
      write(ITROUT-1,1003)sr2np,sr1np/BGTnp,
     1                   dsqrt(sr2np/BGTnp-sr1np*sr1np/BGTnp/BGTnp)
 1000 FORMAT('  ',F8.4,2X,2F12.6)
 1001 FORMAT(2H# ,'SR1=  ',F12.6,2H  ,'COM= ',F12.6)
 1002 FORMAT(2H# ,'BGT=  ',F12.6,2H  ,'DIF(+/-)= ',F12.6)
 1003 FORMAT(2H# ,'SR2=  ',F12.6,2H  ,'CENTR=',F10.2,' WIDTH=',F10.2)

      return
      end
      
      
      subroutine OB2phPN(v1,v2,O,ns1,ns2,np,O12)
C-------------------------------------------------------------------------
C Transforms a pn/np fundamental basis matrix to a deformed basis
C-------------------------------------------------------------------------
      implicit none
      
      integer np
      real v1(np,np),v2(np,np)
      integer ns1,ns2
      real O(np,np)
      
      real O12(np,np)
      
      integer i,j,k,l
      real*8 sum
      
      do i=1,ns1
         do j=1,ns2
	    sum=0.0
	    do k=1,ns1
	       do l=1,ns2
	          sum=sum+v1(k,i)*O(k,l)*v2(l,j)
	       enddo
	    enddo
	    O12(i,j)=sngl(sum)
	 enddo
      enddo
      
      return
      end
      
      

      subroutine sumOandXYpn(maph1,maph2,n1,n2,O,X,Y,sum,ns)
C==================================================================================
C Calculates the RPA matrix element for the transition operator O
C==================================================================================
C      implicit none
      
C      include 'gcb_dim.inc'
      
      integer n1,n2,ns
      integer maph1(n1,2),maph2(n2,2)
      real O(ns,ns)
      real X(n1),Y(n2)
      
      real*8 sum
      
      integer m,i
      integer k,kk
      
      sum=0.0
      kk=0
      do k=1,n1
         m=maph1(k,1)
	 i=maph1(k,2)
	 sum=sum+O(m,i)*X(k)
      enddo
      
      do k=1,n2
         m=maph2(k,1)
	 i=maph2(k,2)
	 sum=sum+O(i,m)*Y(k)
      enddo
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE srtval(e,f,n,np)
      INTEGER n,np
      REAL e(np),p
      real*8 f(np),pp
      INTEGER i,j,k

      do i=1,n-1
        k=i
        p=e(i)
        do  j=i+1,n
          if(e(j).le.p)then
            k=j
            p=e(j)
          endif
        enddo

        if(k.ne.i)then
          e(k)=e(i)
          e(i)=p

            pp=f(i)
            f(i)=f(k)
            f(k)=pp
        endif
      enddo

      return
      END






























