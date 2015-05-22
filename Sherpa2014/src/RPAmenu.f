CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine RPAmenu
C
C
C

      use RPAstab
      use flags
      implicit none
      character*1 menu_choice

      logical errRPA
      integer n
      character*1 ychar
C..............FLAGS.............................
C      logical pnRPAflag

      logical rpasolved

      rpasolved=.false.

1     continue
      write(6,*)' '
      write(6,*)' 	RPA MENU '
      write(6,*)' '
      write(6,*)' Choose one of the following: '
      write(6,*)' (S) RPA excitation spectrum  '
      write(6,*)' (E) RPA correlation energy '
      write(6,*)' (O) RPA corrections to observable'
      write(6,*)' (T) RPA transition strengths '
      write(6,*)' (Y) Magnitude of RPA Y matrix '
      write(6,*)' (P) Parity of RPA modes '
      write(6,*)' (X) Exit this menu '

      read(5,'(a)')menu_choice
      noRPA=.false.


C---------------------------------------------

      if(menu_choice.eq.'s' .or. menu_choice.eq.'S')then
 	if(.not.rpaflag)then
	    call setupAB
            if(noRPA)then
              print*,' Please choose another SD. One cannot set RPA.'
              return
            endif
	    rpaflag=.true.
        endif
        if(.not.rpasolved)then
            call solve_stability(n,errRPA)
	    if(errRPA)then
		write(6,*)' RPA unstable; continue (y/n)?'
		read(5,'(a)')ychar
		if(ychar.eq.'n' .or. ychar.eq.'N')goto 1
            endif
            call solve_rpa
	    rpasolved=.true.
        endif
        call spectrum_rpa
        goto 1
      endif


C---------------------------------------------

      if(menu_choice.eq.'e' .or. menu_choice.eq.'E')then
 	if(.not.rpaflag)then
	    call setupAB
            if(noRPA)then
              print*,' Please choose another SD. One cannot set RPA.'
              return
            endif
	    rpaflag=.true.
        endif
        if(.not.rpasolved)then
            call solve_stability(n,errRPA)
	    if(errRPA)then
		write(6,*)' RPA unstable; continue (y/n)?'
		read(5,'(a)')ychar
		if(ychar.eq.'n' .or. ychar.eq.'N')goto 1
            endif
            call solve_rpa
	    rpasolved=.true.
        endif
        call rpa_corr_energy
        goto 1
      endif

C---------------------------------------------

      if(menu_choice.eq.'o' .or. menu_choice.eq.'O')then
	if(.not.obsflag)then
	   write(6,*)' Must first read in observable '
	   write(6,*)' (Go to Hartree-Fock menu please)'
           return
        endif
 	if(.not.rpaflag)then
	    call setupAB
            if(noRPA)then
              print*,' Please choose another SD. One cannot set RPA.'
              return
            endif
	    rpaflag=.true.
        endif
        if(.not.rpasolved)then
            call solve_stability(n,errRPA)
	    if(errRPA)then
		write(6,*)' RPA unstable; continue (y/n)?'
		read(5,'(a)')ychar
		if(ychar.eq.'n' .or. ychar.eq.'N')goto 1
            endif

            call solve_rpa
	    rpasolved=.true.
        endif
        call RPA_obs
        goto 1
      endif

C--------------------TRANSITIONS-------------------

      if(menu_choice.eq.'t' .or. menu_choice.eq.'T')then
 	if(.not.rpaflag)then
	    call setupAB
            if(noRPA)then
              print*,' Please choose another SD. One cannot set RPA.'
              return
            endif
	    rpaflag=.true.
        endif
        if(.not.rpasolved)then
            call solve_stability(n,errRPA)
	    if(errRPA)then
		write(6,*)' RPA unstable; continue (y/n)?'
		read(5,'(a)')ychar
		if(ychar.eq.'n' .or. ychar.eq.'N')goto 1
            endif

            call solve_rpa
	    rpasolved=.true.
        endif
        call rpa_transition
        goto 1
      endif

C---------------------------------------------

      if(menu_choice.eq.'y' .or. menu_choice.eq.'Y')then
 	if(.not.rpaflag)then
	    call setupAB
            if(noRPA)then
              print*,' Please choose another SD. One cannot set RPA.'
              return
            endif
	    rpaflag=.true.
        endif
        if(.not.rpasolved)then
            call solve_stability(n,errRPA)
	    if(errRPA)then
		write(6,*)' RPA unstable; continue (y/n)?'
		read(5,'(a)')ychar
		if(ychar.eq.'n' .or. ychar.eq.'N')goto 1
            endif
            call solve_rpa
	    rpasolved=.true.
        endif
        write(6,*)' '
        write(6,*)' This option prints out lowest ',
     & ' (Y_mi(lambda))^2,'
        write(6,*)' which measures the 2p-2h strength in gs.'
        call RPA_Ystrength
        goto 1
      endif

C----------------------------PARITY-----------------

      if(menu_choice.eq.'p' .or. menu_choice.eq.'P')then
 	if(.not.rpaflag)then
	    call setupAB
            if(noRPA)then
              print*,' Please choose another SD. One cannot set RPA.'
              return
            endif
	    rpaflag=.true.
        endif
        if(.not.rpasolved)then
            call solve_stability(n,errRPA)
	    if(errRPA)then
		write(6,*)' RPA unstable; continue (y/n)?'
		read(5,'(a)')ychar
		if(ychar.eq.'n' .or. ychar.eq.'N')goto 1
            endif
            call solve_rpa
	    rpasolved=.true.
        endif
        write(6,*)' '
        write(6,*)' This option prints out parity of RPA modes '
        call RPA_parity
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

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine spectrum_rpa
C
C  writes out RPA excitation spectrum
C
C  
C  SUBROUTINES CALLED
C
C
      use RPAsolutions
      use RPAmatrices
      implicit none
C      include 'gcb_dim.inc'

C..........RPA SOLUTIONS........................
      real wtemp(nph) 	! for sorting
      real,pointer      :: w(:)

C..............FILE HANDLING..........................

      character*1 ychar
      character filename*15  		! 
      integer ilast
      character title*60

C...............CHECK ZEROES 
      integer i,nzero
      real cutoff
      data cutoff/0.01/

      nzero = 0 
      w=>w0
      do i = 1,nph!-n0

	if(w(i).lt.-cutoff)then
	   write(6,*)' appears to have negative frequencies ',
     & i,w(i)
        endif
 	if(abs(w(i)).lt.cutoff)nzero=nzero+1
      enddo
      if(nzero.gt.0)then
	write(6,*)nzero,' zero modes '
      endif

C..........check agreement
      if(n0.ne.nzero)then
	write(6,*)' disagree # of zeros ',nzero,n0
C	stop
      endif

C-------------NOW WRITE OUT ---------------------------

311    continue
      write(6,*)' Enter output filename (.ex)'
      read(5,'(a)')filename
      ilast=index(filename,' ')
      if(ilast.ne.0)then
         ilast=ilast-1
      else
         ilast=15
      endif

      open(unit=2,file=filename(1:ilast)//'.ex',status='new',err=33)
      goto 44
33    continue
      write(6,*)' That file already exists;',
     &' do you want to overwrite (y/n/x=exit)? '
      read(5,'(a)')ychar
      if(ychar.eq.'x' .or. ychar.eq.'X')return

      if(ychar.eq.'y' .or.ychar.eq.'Y')then
      open(unit=2,file=filename(1:ilast)//'.ex',status='old',err=33)
      else
 	goto 311
      endif
44    continue
	
      write(6,*)' Enter title line '
      read(5,'(60a)')title
      write(2,'(60a)')title
      wtemp=0.0
      do i = 1,nph-n0
 	wtemp(i)=w(i)
      enddo
      call srt(wtemp,nph,nph)
      do i = 1,nph
	write(2,100)i,wtemp(i)
      enddo

 100  format(I4,2x,F12.6)

      close(unit=2)

      return
      end 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine rpa_corr_energy
      use obop

      real ecor


      call RPA_ecorr(ecor)
      write(6,*)' HF energy      RPA corr        total '
      write(6,99)ehf,ecor,ehf+ecor
   99 format(3x,f10.4,3x,f10.4,3x,f10.4)

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine RPA_ecorr(ecor)
C
C  computes RPA correlation energy
C
C
      use RPAmatrices
      use RPAsolutions
      implicit none
C      include 'gcb_dim.inc'


C--------OUTPUT

      real ecor       ! correlation energy

C-----Temporary variables
      real,pointer        :: w(:)
      real*8 trace1,trace2
      integer i

      w=>w0
C.......Calculate Tr(A)........
      trace1=0.d0
      do i=1,nph
         trace1=trace1+dble(a(i,i))
      enddo

      trace2=0.d0
      do i=1,nph-n0
         trace2=trace2+w(i)
      enddo
      ecor=sngl(trace2-trace1)/2.0

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine RPA_obs
C
C  computes RPA corrections to a scalar observable already read in
C

      use spspace
      use nuclide
      use phBASIS
      use nuclide
      use observables
      use RPAsolutions
      use RPAmatrices
      implicit none
C      include 'gcb_dim.inc' 



C.............EXPECTATION VALUES

      real cutoff
      data cutoff/0.01/
C..................RPA SOLUTIONS....................
      integer,pointer,dimension(:,:)  :: mapp,mapn

C...........OUTPUT
      real O_HF
      real*8 O_corr

      mapp=>map_prot
      mapn=>map_neutr

      call OpCorr(X,Y,nphp,nphn,mapp,mapn,
     &      vecp,vecn,nsps,n0,np,nn,O_HF,O_corr)

C...........CHECK THAT HF expectation values agree.........

      if(abs(O_HF).lt.cutoff)then
	if(abs(O_hf-obsval).gt.cutoff)then
	 	write(6,*)' observables do not match ',
     & obsval,o_Hf
c		return
        endif

      else
 	if( abs(O_hf-obsval)/abs(o_hf).gt.cutoff)then
	 	write(6,*)' observables do not match ',
     & obsval,o_Hf
		 	write(6,*)' observables do not match ',
     & obsval,o_Hf
c		return
        endif
	
        
      endif
C...............write out

      write(6,*)' HF value    RPA corr      HF+RPA  '
      write(6,99)O_hf,O_corr,o_hf+o_corr
99    format(2x,f9.4,2x,f9.4,2x,f9.4)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine rpa_transition

      use spspace
      use nuclide
      use phBASIS
      use nuclide
      use RPAsolutions
      use RPAmatrices
      use obop
      implicit none
C...................FILE OPERATIONS

      character ychar*1
      character*15 filename
      integer ilast


      integer,pointer,dimension(:,:)  :: mapp,mapn


C...Transitions strength calculations.........

C Transition operator
      ITROUT=22
      call readTransOp

C...........open a file.................

311    continue
      write(6,*)' Enter output filename (.str)'
      read(5,'(a)')filename
      ilast=index(filename,' ')
      if(ilast.ne.0)then
         ilast=ilast-1
      else
         ilast=15
      endif

      open(unit=ITROUT,file=filename(1:ilast)//'.str',status='new',
     & err=33)
      goto 44
33    continue
      write(6,*)' That file already exists;',
     &' do you want to overwrite (y/n/x=exit)? '
      read(5,'(a)')ychar
      if(ychar.eq.'x' .or. ychar.eq.'X')return
      
      if(ychar.eq.'y' .or.ychar.eq.'Y')then
      open(unit=ITROUT,file=filename(1:ilast)//'.str',status='old',
     & err=33)
      else
 	goto 311
      endif
44    continue      

      write(6,*)' '

      call main_trans

C      call hfAVcom(vecp,vecn,nphp,nphn,nsps,map_prot,map_neutr,app,ann,
C     &                 apn,Bpp,Bnn,Bpn)


      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine RPA_YstrengthOLD
C
C  computes total 2p-2h amplitudes
C  OLD VERSION
C  
C  SUBROUTINES CALLED
C
C
      use RPAsolutions
      use RPAmatrices
      implicit none

      integer i,lambda

      real ystrength


      ystrength = 0.0
 
      do lambda = 1,nph-n0
 	do i =1,nph
	   ystrength=ystrength+y(i,lambda)*y(i,lambda)
        enddo

      enddo

      write(6,*)' '
      write(6,*)'  Sum |Y|^2 = ',ystrength
    
      return
      end 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine RPA_Ystrength
C
C  writes out RPA excitation spectrum
C
C  
C  SUBROUTINES CALLED
C
C
      use RPAsolutions
      use RPAmatrices
      implicit none
C      include 'gcb_dim.inc'

C..........RPA SOLUTIONS..............
      integer i,lambda,j


      real wtemp(nph) 	! for sorting
      real ytemp(nph)
      real,pointer      :: w(:)

C..............FILE HANDLING..........................

      character*1 ychar
      character filename*15  		! 
      integer ilast
      character title*60


      real ystrength
C...............CHECK ZEROES 
      integer nzero
      real cutoff
      data cutoff/0.01/

      nzero = 0 
      w=>w0

      write(63,*)nph
      write(63,*)w
      write(63,*)w0
      do i = 1,nph!-n0
	if(w(i).lt.-cutoff)then
	   write(6,*)' appears to have negative frequencies ',
     & i,w(i)
        endif
 	if(abs(w(i)).lt.cutoff)nzero=nzero+1
      enddo
      if(nzero.gt.0)then
	write(6,*)nzero,' zero modes '
      endif

C..........check agreement
      if(n0.ne.nzero)then
	write(6,*)' disagree # of zeros ',nzero,n0
C	stop
      endif

C-------------NOW WRITE OUT ---------------------------

311    continue
      write(6,*)' Enter output filename (.y)'
      read(5,'(a)')filename
      ilast=index(filename,' ')
      if(ilast.ne.0)then
         ilast=ilast-1
      else
         ilast=15
      endif

      open(unit=2,file=filename(1:ilast)//'.y',status='new',err=33)
      goto 44
33    continue
      write(6,*)' That file already exists;',
     &' do you want to overwrite (y/n/x=exit)? '
      read(5,'(a)')ychar
      if(ychar.eq.'x' .or. ychar.eq.'X')return

      if(ychar.eq.'y' .or.ychar.eq.'Y')then
      open(unit=2,file=filename(1:ilast)//'.y',status='old')
      else
 	goto 311
      endif
44    continue
	
      ytemp(:) = 0.
      wtemp=0.0
      do i = 1,nph!-n0
 	wtemp(i)=w(i)
        if(wtemp(i) > cutoff)then
        ystrength = 0.0
        do j = 1,nph
	   ystrength=ystrength+y(j,i)*y(j,i)
        enddo
        ytemp(i)=ystrength
        endif
        
      enddo
      call srty(wtemp,ytemp,nph,nph)
C      call srt(wtemp,nph,nph)
      do i = 1,nph
	write(2,100)i,wtemp(i),ytemp(i)
      enddo

 100  format(I4,2x,2F10.5)

      close(unit=2)

      return
      end 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine RPA_parity
C
C  writes out parity of RPA modes
C
C  
C  SUBROUTINES CALLED
C
C
      use RPAsolutions
      use RPAmatrices
      use phBASIS
      use spspace

      implicit none
C      include 'gcb_dim.inc'

C..........RPA SOLUTIONS..............
      integer i,lambda,j


      real wtemp(nph) 	! for sorting
      real partemp(nph)
      integer parvalue(nph)
      real,pointer      :: w(:)
      integer i1,i2
      real tmp

C..............FILE HANDLING..........................

      character*1 ychar
      character filename*15  		! 
      integer ilast
      character title*60


      real parstrength
C...............CHECK ZEROES 
      integer nzero
      real cutoff
      data cutoff/0.01/

      nzero = 0 
      w=>w0

      write(63,*)nph
      write(63,*)w
      write(63,*)w0
      do i = 1,nph!-n0
	if(w(i).lt.-cutoff)then
	   write(6,*)' appears to have negative frequencies ',
     & i,w(i)
        endif
 	if(abs(w(i)).lt.cutoff)nzero=nzero+1
      enddo
      if(nzero.gt.0)then
	write(6,*)nzero,' zero modes '
      endif

C..........check agreement
      if(n0.ne.nzero)then
	write(6,*)' disagree # of zeros ',nzero,n0
C	stop
      endif

!------------------ COMPUTE PARITY VECTOR ------------

      do i = 1,nphp
        i1 = map_prot(i,1)
        i2 = map_prot(i,2)
        parvalue(i) = (-1)**( spsqn(1,i1,3) + spsqn(1,i2, 3))

      enddo

      do i = 1,nphn
        i1 = map_neutr(i,1)
        i2 = map_neutr(i,2)
        parvalue(i+nphp) = (-1)**( spsqn(2,i1,3) + spsqn(2,i2, 3))

      enddo

C-------------NOW WRITE OUT ---------------------------

311    continue
      write(6,*)' Enter output filename (.parity)'
      read(5,'(a)')filename
      ilast=index(filename,' ')
      if(ilast.ne.0)then
         ilast=ilast-1
      else
         ilast=15
      endif

      open(unit=2,file=filename(1:ilast)//'.parity',status='new',err=33)
      goto 44
33    continue
      write(6,*)' That file already exists;',
     &' do you want to overwrite (y/n/x=exit)? '
      read(5,'(a)')ychar
      if(ychar.eq.'x' .or. ychar.eq.'X')return

      if(ychar.eq.'y' .or.ychar.eq.'Y')then
      open(unit=2,file=filename(1:ilast)//'.parity',status='old')
      else
 	goto 311
      endif
44    continue
	
      partemp(:) = 0.
      wtemp=0.0
      do i = 1,nph!-n0
 	wtemp(i)=w(i)
        if(wtemp(i) > cutoff)then
        parstrength = 0.0
        tmp = 0.0
        do j = 1,nph
	   parstrength=parstrength+(x(j,i)*x(j,i) +y(j,i)*y(j,i)) 
     &   *parvalue(j)
	   tmp=tmp+x(j,i)*x(j,i) +y(j,i)*y(j,i)

        enddo
        partemp(i)=parstrength/tmp
        endif
        
      enddo
      call srty(wtemp,partemp,nph,nph)
C      call srt(wtemp,nph,nph)
      do i = 1,nph
	write(2,100)i,wtemp(i),partemp(i)
      enddo

 100  format(I4,2x,2F10.5)

      close(unit=2)

      return
      end 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE srt(d,n,np)
      INTEGER n,np
      REAL d(np)
      INTEGER i,j,k
      REAL p
      do 13 i=1,n-1
        k=i
        p=d(i)
        do 11 j=i+1,n
          if(d(j).le.p)then
            k=j
            p=d(j)
          endif
11      continue
        if(k.ne.i)then
          d(k)=d(i)
          d(i)=p
        endif
13    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software #<6V*1(.31..

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE srty(d,y,n,np)
      INTEGER n,np
      REAL d(np),y(np)
      INTEGER i,j,k
      REAL p
      do 13 i=1,n-1
        k=i
        p=d(i)
        do 11 j=i+1,n
          if(d(j).le.p)then
            k=j
            p=d(j)
          endif
11      continue
        if(k.ne.i)then
          d(k)=d(i)
          d(i)=p
          p = y(k)
          y(k) = y(i)
          y(i) = p
        endif
13    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software #<6V*1(.31..     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine pnRPA_ecorr(ecor)
C
C  computes pnRPA correlation energy
C
C
      use RPAmatrices
      use RPAsolutions
      implicit none
C      include 'gcb_dim.inc'


C--------OUTPUT

      real ecor       ! correlation energy

      real*8 trace1
      integer i

c      print*,nphpn,nphnp
C.......Calculate the correlation energy........
      trace1=0.0
      do i=1,nphnp
        trace1=trace1-anppn(i,i)
        trace1=trace1+wnp(i)
      enddo
      do i=1,nphpn
        trace1=trace1-apnnp(i,i)
        trace1=trace1+wpn(i)
      enddo
      ecor=sngl(trace1/2.d0)
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine pnrpa_corr_energy

      use obop

      real ecor

       ecor=0.0

      call pnRPA_ecorr(ecor)
      write(6,*)' HF energy      RPA corr        total '
      write(6,99)ehf,ecor,ehf+ecor
   99 format(3x,f10.4,3x,f10.4,3x,f10.4)

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine pnRPAmenu
C
C
C
      use RPAstab
      use flags
      implicit none
      character*1 menu_choice

      logical errRPA

      integer n
      character*1 ychar


      logical rpasolved

      rpasolved=.false.
      noRPA=.false.

1     continue
      write(6,*)' '
      write(6,*)'      pnRPA MENU '
      write(6,*)' '
      write(6,*)' Choose one of the following: '
      write(6,*)' (S) pnRPA excitation spectrum (neighbouring nuclei) '
      write(6,*)' (E) pnRPA correlation energy '
C      write(6,*)' (O) pnRPA corrections to observable'
      write(6,*)' (T) pnRPA transition strengths '
C      write(6,*)' (Y) Magnitude of RPA Y matrix '
      write(6,*)' (X) Exit this menu '

      read(5,'(a)')menu_choice


      if(menu_choice.eq.'s' .or. menu_choice.eq.'S')then
 	if(.not.pnrpaflag)then
	    call setupABpn
            if(noRPA)then
              print*,' Please choose another SD. One cannot set RPA.'
              return
            endif
	    pnrpaflag=.true.
        endif
        if(.not.rpasolved)then
            call solve_pnstability(errRPA)
	    if(errRPA)then
		write(6,*)' pnRPA possibly unstable; continue (y/n)?'
		read(5,'(a)')ychar
		if(ychar.eq.'n' .or. ychar.eq.'N')goto 1
            endif
            call solve_pnrpa
	    rpasolved=.true.
        endif
        call spectrum_pnrpa
        goto 1
      endif


C---------------------------------------------

      if(menu_choice.eq.'e' .or. menu_choice.eq.'E')then
 	if(.not.pnrpaflag)then
	    call setupABpn
            if(noRPA)then
              print*,' Please choose another SD. One cannot set RPA.'
              return
            endif
	    pnrpaflag=.true.
        endif
        if(.not.rpasolved)then
            call solve_pnstability(errRPA)
	    if(errRPA)then
		write(6,*)' pnRPA possibly unstable; continue (y/n)?'
		read(5,'(a)')ychar
		if(ychar.eq.'n' .or. ychar.eq.'N')goto 1
            endif

            call solve_pnrpa
	    rpasolved=.true.
        endif
        call pnrpa_corr_energy
        goto 1
      endif

C--------------------TRANSITIONS-------------------

      if(menu_choice.eq.'t' .or. menu_choice.eq.'T')then
 	if(.not.pnrpaflag)then
	    call setupABpn
            if(noRPA)then
              print*,' Please choose another SD. One cannot set RPA.'
              return
            endif
	    pnrpaflag=.true.
        endif
        if(.not.rpasolved)then
            call solve_pnstability(errRPA)
	    if(errRPA)then
		write(6,*)' pnRPA possibly unstable; continue (y/n)?'
		read(5,'(a)')ychar
		if(ychar.eq.'n' .or. ychar.eq.'N')goto 1
            endif

            call solve_pnrpa
	    rpasolved=.true.
        endif
C        print*,' This option not implemented yet.'
        call rpa_transitionpn
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




      subroutine spectrum_pnrpa
C
C  writes out pnRPA excitation spectrum
C
C  
C  SUBROUTINES CALLED
C
C
      use RPAsolutions
      use RPAmatrices
      implicit none
C      include 'gcb_dim.inc'

C..........RPA SOLUTIONS........................
      real,allocatable    :: wtemp(:) 	! for sorting
      real,pointer      :: w(:)

C..............FILE HANDLING..........................

      character*1 ychar
      character filename*15  		! 
      integer ilast
      character title*60

C...............CHECK ZEROES 
      integer i,nzero
      real cutoff
      data cutoff/0.01/

C-------------NOW WRITE OUT ---------------------------

311    continue
      write(6,*)' Enter output filename (.ex)'
      read(5,'(a)')filename
      ilast=index(filename,' ')
      if(ilast.ne.0)then
         ilast=ilast-1
      else
         ilast=15
      endif

      open(unit=2,file=filename(1:ilast)//'.ex',status='new',err=33)
      goto 44
33    continue
      write(6,*)' That file already exists;',
     &' do you want to overwrite (y/n/x=exit)? '
      read(5,'(a)')ychar
      if(ychar.eq.'x' .or. ychar.eq.'X')return

      if(ychar.eq.'y' .or.ychar.eq.'Y')then
      open(unit=2,file=filename(1:ilast)//'.ex',status='old',err=33)
      else
 	goto 311
      endif
44    continue
	
      write(6,*)' Enter title line '
      read(5,'(60a)')title
      write(2,'(60a)')title
      allocate(wtemp(nphpn))
      write(2,*)' pn spectrum'
      wtemp=wpn
      call srt(wtemp,nphpn,nphpn)
      do i = 1,nphpn
	write(2,100)i,wtemp(i)
      enddo
      deallocate(wtemp)
      allocate(wtemp(nphnp))
      write(2,*)' pn spectrum'
      wtemp=wnp
      call srt(wtemp,nphnp,nphnp)
      do i = 1,nphnp
	write(2,100)i,wtemp(i)
      enddo
      deallocate(wtemp)

 100  format(I4,2x,F12.6)

      close(unit=2)

      return
      end 
