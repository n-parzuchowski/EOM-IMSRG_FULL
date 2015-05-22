CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   package GCLIB
C
C   library of routines needed for generator-coordinate program GCBasis
C
C   written by CWJ, modifications of routines by WEO
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      subroutine GetSPS
C
C   reads in .sps file information on j-orbits 
C
C  REVISION 8 April 2004: reads in REDSTICK-compatible .sps files
C  Revision 13 May 2004: use allocatable arrays
C
C
C   OUTPUT:	norb = # of j-orbits
C		jorb(i) = 2*j of ith j-orbit 
C		nrorb(i)=  radial quantum number of ith j-orbit
C		torb(i) = 2*Tz of j-orbit
C		nsps = # of m-states
C		orb_qn = quantum numbers of j-orbits, stored differently
C		spsqn  = quantum numbers of m-states
C
C   CALLED BY: main
C
      use spspace
      use chf_dim
      implicit none

C..............INTERMEDIATES........................

      integer,allocatable   :: ind(:)     
      real,allocatable      ::  xn(:),xl(:),xj(:)
      integer,allocatable   ::  lorb(:),ilabel(:,:)
      integer i,j,m
      integer ns

C...............FINDING CLEBSCH GORDON ARRAY.....................

      real clebr
      integer ia,ja,k
      real xji,xjj,xmi,xmj,xk
      real xxn,xxl,xxj
      integer iphase
      
C..............FILE HANDLING..........................
 
      character spsfil*15  		! name of .sps file 
      integer ilast 
      character isoread*3
      logical isoflag
      real yy				! dummy

C================================================================      


      if(allocated(spsqn))deallocate(spsqn)
      if(allocated(orb_qn))deallocate(orb_qn)
      if(allocated(jorb))deallocate(jorb)
      if(allocated(torb))deallocate(torb)
      if(allocated(nrorb))deallocate(nrorb)
      if(allocated(spsqn))deallocate(spsqn)
      if(allocated(kmax))deallocate(kmax)
      if(allocated(kmin))deallocate(kmin)
      if(allocated(indx))deallocate(indx)
      if(allocated(clb))deallocate(clb)

1     continue
      write(6,*)' Enter shell-model space file name (.sps)'
      read(5,'(a)')spsfil
      ilast=index(spsfil,' ')
      if(ilast.ne.0)then
         ilast=ilast-1
      else
         ilast=15
      endif
C----------  Now get orbit information for the calculation
C----------  njl is the number of single-particle states

      open(unit=1,file=spsfil(1:ilast)//'.sps',status='old',err=3)
      goto 4
3     continue
      write(6,*)' That file does not exist '
      goto 1

4     continue
      write(6,'('' Shell-Model space file name'',1x,a15)')spsfil

C.............................READS IN REDSTICK-COMPATIBLE FORMAT........

      read(1,'(a3)')isoread
      if(isoread.eq.'iso' .or. isoread.eq.'ISO')then
	isoflag=.true.

      elseif(isoread.eq.'pn' .or. isoread.eq.'PN')then
	isoflag=.false.
      else
	write(6,*)' You do not have a redstick compatible .sps file '
        write(6,*)' Nonetheless, I will try to read with old format '
        rewind(1)
        goto 1011
      endif

 
C....................ISOSPIN FORMALISM..........................

      if(isoflag)then
	read(1,*)norb(1)
	allocate(xn(norb(1)),xl(norb(1)),xj(norb(1)))
        allocate(nrorb(2*norb(1)),lorb(2*norb(1)),jorb(2*norb(1)),
     &                                            torb(2*norb(1)))
        do i = 1,norb(1)
C        ind(i)=i
      	read(1,*,err=44)xn(i),xl(i),xj(i),yy
C        write(6,*)i,xn(i),xl(i),xj(i),yy
      	nrorb(i)=int(xn(i))
      	jorb(i)=int(2*xj(i))
c        if(jorb(ind(i)).gt.maxk)maxk=jorb(ind(i))
      	torb(i)=1
      	lorb(i) = int(xl(i))
C................ERROR TRAPS ...............

        spinless = .false.
        if(xl(i) == xj(i))spinless = .true.
!        if(abs(xl(i)-xj(i)).gt.0.6)then
!	   write(6,*)' l, j : ' , xl(i),xj(i), ' not compatible '
!	   goto 3033
!        endif
!        if(abs(xl(i)-xj(i)).lt.0.4)then
!	   write(6,*)' l, j : ' , xl(i),xj(i), ' not compatible '
!	   goto 3033
!        endif         
!        if(abs( mod(jorb(i),2) ) .ne.1)then
!	   write(6,*)' j :',xj(i),jorb(i) ,' not half-integer '
!	   goto 3033
!        endif

C..............END ERROR TRAPS.............
     	
         enddo
         goto 35
   44    continue
	 write(6,*)' some error in reading file ',i
	 stop

      endif


C.......................PN FORMALISM.............................


 	if(.not.isoflag)then
		write(6,*)' This formalism not yet implemented '
		goto 1
	endif

      return
C.............................THIS SECTION IS FOR OLD FORMAT.........

 1011 continue

      write(6,*)' This version assumes equal proton and neutron orbits '
      norb(1) = 0
C      maxk = 0
      i=0
      do
        read(1,*,err=3033,end=34)i,xxn,xxl,xxj
        norb(1)=norb(1)+1
C................ERROR TRAPS ...............
        spinless = .false.
        if(xl(i) == xl(j))spinless = .true.


C..............END ERROR TRAPS.............
     	
      enddo

 34   continue


      allocate(ind(norb(1)),xn(norb(1)),xl(norb(1)),xj(norb(1)))
      allocate(nrorb(2*norb(1)),lorb(2*norb(1)),jorb(2*norb(1)),
     &                                          torb(2*norb(1)))

      rewind(1)


      do i=1,norb(1)
      	read(1,*,end=35,err=3033)ind(i),xn(i),xl(i),xj(i)
        write(6,*)i,ind(i),xn(i),xl(i),xj(i)
      	nrorb(ind(i))=int(xn(i))
      	jorb(ind(i))=int(2*xj(i))
c        if(jorb(ind(i)).gt.maxk)maxk=jorb(ind(i))
      	torb(ind(i))=1
      	
      	lorb(ind(i)) = int(xl(i))


C..............END ERROR TRAPS.............
     	
      enddo      

      deallocate(xn,xj,xl,ind)

   35 continue
      close(unit=1) 
c      deallocate(xl,xj,xn)
      norb(2) = norb(1)
      
      do i =1,norb(1)
      	write(6,*)i,nrorb(i),jorb(i)

C-----------FILL NEUTRON SPACES --------------

      	nrorb(i+norb(1))=nrorb(i)
      	jorb(i+norb(1)) = jorb(i)
      	lorb(i+norb(1)) = lorb(i)
      	torb(i+norb(1)) = -1
      	
      enddo

C-------- FILL IN WEO ARRAYS

      allocate(orb_qn(norb(1)+norb(2),3))

      do i =1,norb(1)+norb(2)
      	orb_qn(i,1) = float(nrorb(i))
      	orb_qn(i,2) = float(lorb(i))
      	orb_qn(i,3) = float(jorb(i))/2.
      enddo

      ns = 0
      do i=1,norb(1)
         do m=-jorb(i),jorb(i),2
           ns=ns+1
         enddo
      enddo

      allocate(spsqn(2,ns,6))

      j2max=0
      ns=0
      do i =1,norb(1)
        j2max=max(j2max,jorb(i))
      	do m = -jorb(i),jorb(i),2
      	    ns=ns+1
      	    spsqn(1,ns,1) = i
      	    spsqn(1,ns,2) = nrorb(i)
      	    spsqn(1,ns,3) = lorb(i)
      	    spsqn(1,ns,4) = jorb(i)
      	    spsqn(1,ns,5) = m
            spsqn(1,ns,6) = 1
      	enddo
      enddo

      ns=0
      do i=norb(1)+1,norb(1)+norb(2)
        do m=-jorb(i),jorb(i),2
           ns=ns+1
      	    spsqn(2,ns,1) = i
      	    spsqn(2,ns,2) = nrorb(i)
      	    spsqn(2,ns,3) = lorb(i)
      	    spsqn(2,ns,4) = jorb(i)
      	    spsqn(2,ns,5) = m
            spsqn(2,ns,6) = -1
        enddo
      enddo

      nsps(1) = ns
      nsps(2) = ns


      nindx=0
Cweo---  Set up all possible two particle combinations along with spin
Cweo---  Also find min and max angular momentum that can be coupled
Cweo---  Do all like particles first

C CWJ -- NOTE NOTE I assume proton neutron have identical orbits --

      allocate(kmax(norb(1),norb(1)),kmin(norb(1),norb(1)))

      do i=1,norb(1)
        do j=1,norb(1)
            if(torb(i).eq.torb(j))then
                  nindx=nindx+1
            end if
      	end do
      end do

      allocate(indx(2,nindx))

      nindx=0
      allocate(ilabel(norb(1),norb(1)))
c      print*,'norb',norb(1)
      do i=1,norb(1)
        do j=1,norb(1)
            if(torb(i).eq.torb(j))then
                  nindx=nindx+1
                  indx(1,nindx)=i
                  indx(2,nindx)=j
		  ilabel(i,j) = nindx
c                  print*,i,j
                  kmax(i,j)=int(xj(i)+xj(j))
                  kmin(i,j)=int(abs(xj(i)-xj(j)))
            end if
      	end do
      end do
Cweo---- Now different types  ---- Note xtz(i)=0. in isospin formalism
C      do i=1,norb(1)+norb(2)
C        do j=1,norb(1)+norb(2)
C           if(torb(i).ne.torb(j))then
C                  nindx=nindx+1
C                  indx(1,nindx)=i
C                  indx(2,nindx)=j
C		  ilabel(i,j) = nindx
C                  kmax(i,j)=int(xj(i)+xj(j))
C                  kmin(i,j)=int(abs(xj(i)-xj(j)))
C             end if
C        end do
C      end do
C      write(6,*)'  nindx ',nindx
C
C---- compute phased Clebsch-Gordon array
      allocate(clb(0:j2max,1:nsps(1),1:nsps(1)))
      do i = 1,nsps(1)
         ia = spsqn(1,i,1)
         xji = float(spsqn(1,i,4))/2.
         xmi = float(spsqn(1,i,5))/2.
         do j = 1,nsps(1)
            ja = spsqn(1,j,1)
            xjj=float(spsqn(1,j,4))/2.
            xmj=float(spsqn(1,j,5))/2.
            iphase =  (spsqn(1,j,4)-spsqn(1,j,5))/2
            do k = kmin(ia,ja),kmax(ia,ja)
               xk = float(k)
               clb(k,i,j)=(-1)**iphase*
     &               clebr(xji,xmi,xjj,-xmj,xK,xmi-xmj)
            enddo
C            write(6,*)i,j,ia,ja,xji, xjj,kmin(ia,ja),kmax(ia,ja)        
         enddo
      enddo

      return
 3033 continue
      write(6,*)' I cannot understand this file. ',
     & 'Please choose another '
      goto 1

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine FindMaxJzPN(np,nn,Max2Jz)

C..............................................
C for a given space, finds the maximal 2*Jz that can be constructed
C (for both protons and neutrons) 
C..............................................
C
C   INPUT:	norb = # of j-orbits  norb(1) for p  norb(2) for n
C		jorb(i) = 2*j of ith orbit
C    		np = # of protons
C		nn = # of neutrons
C
C   OUTPUT: 	Max2Jz = 2 x maximal Jz for this space
C
C   CALLED BY 
C	main 
C   
C   SUBROUTINES CALLED
C	shell:  num. rec. routine for sorting
C
C   DESCRIPTION OF ALGORITHM:
C	for both protons and neutrons create an array MORB of m-values
C	order MORB and then simply fill
C

      use spspace
      implicit none
C      include 'gcb_dim.inc'

      integer,intent(IN)   :: np,nn		! # of protons, neutrons

      integer,intent(OUT)  :: Max2Jz


      integer,allocatable      ::  morb(:)
      integer i,m,ns
      
      
      Max2Jz = 0
      
      allocate(morb(nsps(1)))

      ns = 0
      do i = 1,norb(1)
	do m = -jorb(i),jorb(i),2
	    ns = ns+1
	    morb(ns) = m
	enddo
      enddo
      call shell(ns,nsps(1),morb)

      do i = 1,np
      	max2jz = Max2Jz + morb(i)
      enddo

      ns = 0
      do i = 1+norb(1),norb(2)
	do m = -jorb(i),jorb(i),2
	    ns = ns+1
	    morb(ns) = m
	enddo
      enddo
      call shell(ns,ns,morb)

      do i = 1,nn
      	max2jz = Max2Jz + morb(i)
      enddo
      
      max2jz = abs(max2jz)
      deallocate(morb)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine shell(n,np,a)
C
C  A is (integer) array of length n (dimensioned np)
C  sorted by Shell's method, found in Num Rec 8.1
C  into ascending numerical order
C

      implicit none
      
      integer n,np
      
      integer a(np)		! array to be sorted
      integer v			! dummy
      
      integer i,j,inc
C

      inc = 1
    1 inc =3*inc+1
      if(inc.le.n)goto 1
    2 continue
      inc = inc/3
      do i =inc+1,n
      	v = a(i)
      	j = i
    3 	if(a(j-inc).gt.v)then
      	    a(j) = a(j-inc)
      	    j=j-inc
      	    if(j.le.inc)goto 4
      	    goto 3
      	endif
    4 	a(j) = v
      enddo
      if(inc.gt.1)goto 2
      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
