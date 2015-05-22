       program SHERPA
C----------------------------------------------
C  SHEll-model RPA 
C
C  This program performs HF+RPA calculations 
C  in a shell-model (occupation number) framework, 
C  solving the matrix RPA equations.
C
C  version 1.0:  2001-2003 Ionel Stetcu, 
C				Louisiana State University
C			   Calvin Johnson, 
C				LSU & San Diego State University
C version 1.1:  April 2004  IS + CWJ
C              includes pnRPA, improved HF
C version 2.0:  May 2004 IS and CWJ
C               redesign for f90
C version 2.1   Feb 2005 IS
C version 2.2   Mar 2006 CWJ
C-----------------------------------------------
C
C SUBROUTINES CALLED:
C    HF_menu
C    Stablemenu
C    RPAmenu
C    TDAmenu
C
C..............MENU CHOICES......................

      use flags
      use observables
      implicit none
      character*1 menu_choice

      HFflag = .false. 
      hamflag=.false.
      nukeflag=.false.
      obsflag = .false.  ! redo observable each time
      jxyz=.false.       ! the matrices Jx,Jy,Jz have to be computed

      write(6,*)' '
      write(6,*)' '

      write(6,*)' ********************************** '
      write(6,*)' '
      write(6,*)' * Welcome to SHERPA:             * '
      write(6,*)' * the SHEll-model RPA code       * '
      write(6,*)' '
      write(6,*)' * 2014 Edition                   * '
      write(6,*)' * main author: Ionel Stetcu      * '
      write(6,*)' * (2nd author: Calvin Johnson)   * '
      write(6,*)' '
      write(6,*)' * Please reference:              * '
      write(6,*)' * I. Stetcu and C. W. Johnson    * '
      write(6,*)' * PRC 66, 034301(2002); 67, 043315 (2013) '
      write(6,*)' * 69, 024311 (2004)              *'

      write(6,*)' ********************************** '

      write(6,*)' '

1     continue
      RPAflag = .false.  !for now, redo rpa each time

      if(HFflag.and.nukeflag)then	! choose between menus

      write(6,*)' '
      write(6,*)'           MAIN MENU               '
      write(6,*)' '
      write(6,*)' Choose one of the following menus:'
      write(6,*)'(H) Hartree-Fock                   '
      write(6,*)'(R) Random phase approximation     '
      write(6,*)'(S) Stability matrix analysis      '
      write(6,*)'(P) pn random phase approximation  '
      write(6,*)'(X) Exit program                   '
      read(5,'(a)')menu_choice

	else

      menu_choice = 'HF'   ! must get HF state

      endif

C............. HARTREE-FOCK..............................

      if(menu_choice.eq.'H' .or. menu_choice.eq.'h')then
	call HF_menu
        goto 1
      endif

C..............RPA......................................

      if(menu_choice.eq.'R' .or. menu_choice.eq.'r')then

        call rpamenu
c        print*,'Sorry, this option not available yet in version 2.0'
c        print*,'Use a previous version.'
        goto 1
      endif

C..............STABILITY MATRIX........................

      if(menu_choice.eq.'S' .or. menu_choice.eq.'s')then

      	write(6,*)' '
 	write(6,*)' 	STABILITY-MATRIX MENU '
        call stablemenu
C        print*,'Sorry, this option not available yet in version 2.0'
C        print*,'Use a previous version.'
        goto 1
      endif

      if(menu_choice.eq.'P' .or. menu_choice.eq.'p')then
 	call pnRPAmenu
	goto 1
      endif

      if(menu_choice.eq.'X' .or. menu_choice.eq.'x')then
 	stop
      endif


      write(6,*)' That choice not valid '
      goto 1

      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine Stablemenu
C
C  
C
      use flags
      implicit none
      character*1 menu_choice

C..............FLAGS.............................




1     write(6,*)' '
      write(6,*)' 	STABILITY MENU '
      write(6,*)' (mostly for pedantic purposes)'
      write(6,*)' '
      write(6,*)' Choose one of the following: '
      write(6,*)' (S) Check stability, zeros           '
      write(6,*)' (W) Write out stability frequencies  '
      write(6,*)' (X) Exit this menu '
      
      read(5,'(a)')menu_choice


C------
      if(menu_choice.eq.'s' .or. menu_choice.eq.'W')then
 	if(.not.rpaflag)then
c	    call setupAB
	    rpaflag=.true.
        endif
        call checkstability
        goto 1
      endif


C------
      if(menu_choice.eq.'w' .or. menu_choice.eq.'W')then
 	if(.not.rpaflag)then
c	    call setupAB
	    rpaflag=.true.
        endif
        call writeoutstability
        goto 1
      endif
      
C-----------exit----------------------

      if(menu_choice.eq.'X' .or.menu_choice.eq.'x')then
	return
      endif

      write(6,*)' That choice not valid or not implemented'
      goto 1
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine writeoutstability
C
C  writes out eigenvalues of stability matrix
C
C  INPUT:
C  OUTPUT:
C
C  SUBROUTINES CALLED:
C         solve_stability
C
C
      use RPAstab
      implicit none
C      include 'gcb_dim.inc'

C..............FILE HANDLING..........................

      character*1 ychar
      character filename*15  		! 
      integer ilast
 

C.........
      logical errRPA
      integer i
      integer nzero
      character*60 title

      call solve_stability(neigen,errRPA)

C---------------CHECK STABILITY, ZEROES----------------

      if(errRPA) then 
	do i = 1,neigen/2
	  if(stabe(i).lt.-cutoff)then
	     write(6,*)' Instability in real p-h amplitudes: ',
     & stabe(i)
          endif
        enddo
	do i = neigen/2+1,neigen
	  if(stabe(i).lt.-cutoff)then
	     write(6,*)' Instability in imaginary p-h amplitudes: ',
     & stabe(i)
          endif
        enddo
      endif

      nzero=0 
      do i =1,neigen/2
	if(abs(stabe(i)).lt.cutoff)nzero=nzero+1
      enddo
      if(nzero.gt.0)then
	write(6,*)nzero,' zeros in real p-h amplitudes '
      endif
      nzero=0 
      do i =1+neigen/2,neigen
	if(abs(stabe(i)).lt.cutoff)nzero=nzero+1
      enddo
      if(nzero.gt.0)then
	write(6,*)nzero,' zeros in imaginary p-h amplitudes '
      endif


C-------------NOW WRITE OUT ---------------------------

311    continue
      write(6,*)' Enter output filename (.sta)'
      read(5,'(a)')filename
      ilast=index(filename,' ')
      if(ilast.ne.0)then
         ilast=ilast-1
      else
         ilast=15
      endif

      open(unit=2,file=filename(1:ilast)//'.sta',status='new',err=33)
      goto 44
33    continue
      write(6,*)' That file already exists;',
     &' do you want to overwrite (y/n/x=exit)? '
      read(5,'(a)')ychar
      if(ychar.eq.'x' .or. ychar.eq.'X')return
      
      if(ychar.eq.'y' .or.ychar.eq.'Y')then
      open(unit=2,file=filename(1:ilast)//'.sta',status='old',err=33)
      else
 	goto 311
      endif
44    continue      
	
      write(6,*)' Enter title line '
      read(5,'(60a)')title
      write(2,'(60a)')title
      do i = 1,neigen
	write(2,120)i,stabe(i)
      enddo

 120  format(I3,2x,F12.6)

      close(unit=2)
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine checkstability
C
C  examines eigenvalues of stability matrix 
C
C  INPUT:
C  OUTPUT:
C
C  SUBROUTINES CALLED:
C          solve_stability
C
C
      use RPAstab
      implicit none
C      include 'gcb_dim.inc'

C..............FILE HANDLING..........................

      character*1 ychar
      character filename*15  		!
      integer ilast
 

C.........
      logical errRPA
      integer i
      integer nzero

      call solve_stability(neigen,errRPA)

C---------------CHECK STABILITY, ZEROES----------------

      if(errRPA) then 
	do i = 1,neigen/2
	  if(stabe(i).lt.-cutoff)then
	     write(6,*)' Instability in real p-h amplitudes: ',
     & stabe(i)
          endif
        enddo
	do i = neigen/2+1,neigen
	  if(stabe(i).lt.-cutoff)then
	     write(6,*)' Instability in imaginary p-h amplitudes: ',
     & stabe(i)
          endif
        enddo
      endif

      nzero=0 
      do i =1,neigen/2
	if(abs(stabe(i)).lt.cutoff)nzero=nzero+1
      enddo
      if(nzero.gt.0)then
	write(6,*)nzero,' zeros in real p-h amplitudes '
      endif
      nzero=0 
      do i =1+neigen/2,neigen
	if(abs(stabe(i)).lt.cutoff)nzero=nzero+1
      enddo
      if(nzero.gt.0)then
	write(6,*)nzero,' zeros in imaginary p-h amplitudes '
      endif

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC













