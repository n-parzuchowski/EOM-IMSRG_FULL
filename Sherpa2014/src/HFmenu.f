CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine HF_menu
C
C  controls all Hartree-Fock procedures
C
C  INPUT:
C
C  SUBROUTINES CALLED:
C  	ChooseNuclide: 
C 	HFrandomstart   grad descent with random starting vectors
C	read_sd
C	save_sd
      use flags
      implicit none


C.............MENU CONTOL...................

      character*1 menu_choice
C.	

C      include 'gcb_dim.inc'
      

C..............FILE HANDLING..........................

      character*1 ychar
     
C===============================================================

1     continue

      if(nukeflag)then		! read in options
      write(6,*)' '
      write(6,*)' 	HARTREE-FOCK MENU '
      write(6,*)' '
      write(6,*)' Choose one of the following: '
      write(6,*)' (C) Choose new nucleus/s.p. space(C)'
      write(6,*)' (N) Generate new HF state(s)     (N)'
      write(6,*)' (R) Read in HF state from file   (R)'
      if(HFflag)then
	write(6,*)' (E) Recompute energy             (E)' 
	write(6,*)' (W) Write HF state to file       (W)'
        write(6,*)' (D) Further gradient Descent     (D)'
        write(6,*)' (O) Compute scalar observable    (O)'
        write(6,*)' (X) Exit this menu               (X)'
 
      endif
      read(5,'(a)')menu_choice
      else	! FORCE GOING TO 'CHOOSE NUCLEUS '
	menu_choice='C'
      endif
C.........................CHOOSE NUCLEUS...............

      if(menu_choice.eq.'C'.or.menu_choice.eq.'c')then
      	call ChooseNuclide
        obsflag=.false.
 	nukeflag=.true.
        hfflag=.false.
      	goto 1
      endif

C------------------------ GENERATE NEW STATES ------------
      if(menu_choice.eq.'N'.or.menu_choice.eq.'n')then
	call HF(hamflag,.true.)
        HFflag=.true.        
     	goto 1
      endif

C-----------------NEW GRADIENT DESCENT ----------------------------

      if(menu_choice.eq.'D' .or.menu_choice.eq.'d')then
	if(.not.HFflag)then
   	  write(6,*)' Must create HF state first '
          goto 1
        endif
        call HF(hamflag,.false.)
        HFflag = .true.
        obsflag=.false.
      endif

C---------------RECOMPUTE ENERGY-----------------------------------

      if(menu_choice.eq.'E' .or. menu_choice.eq.'e')then
	if(.not.HFflag)then
   	  write(6,*)' Must create HF state first '
          goto 1
        endif
        call check_energy  !(hamflag)
	goto 1    
     	
      endif

C------------------READ IN SLATER DETERMINANT ---------------------

      if(menu_choice.eq.'R' .or.menu_choice.eq.'r')then
	if(.not.nukeflag)then
   	  write(6,*)' Must choose nucleus first '
          goto 1
        endif
        call read_sd
    
  	HFflag = .true.
        goto 1
      endif


C-----------------WRITE OUT SLATER DETERMINANT -------------------

      if(menu_choice.eq.'W' .or.menu_choice.eq.'w')then
	if(.not.HFflag)then
   	  write(6,*)' Must create HF state first '
          goto 1
        endif
        call save_sd(0)
	goto 1    
      endif

C.....................SCALAR OBSERVABLE..........................

      if(menu_choice.eq.'o' .or. menu_choice.eq.'O')then
	if(.not.HFflag)then
   	  write(6,*)' Must create/read in HF state first '
          goto 1
        endif

	call HFobs
        obsflag=.true.
	goto 1
      endif

C---------------------------------

      if(menu_choice.eq.'X' .or.menu_choice.eq.'x')then
	return
      endif

      write(6,*)' That choice not valid or not implemented '
      goto 1

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine ChooseNuclide
C
C  reads in s.p. space, Z,N of nuclide
C
C  SUBROUTINES CALLED
C	GetSPS		reads in s.p. space
C
      use spspace
      use nuclide
      use wfn
      implicit none

      write(6,*)' '
      write(6,*)' 	Specify nucleus '
      write(6,*)' '
C...........READ IN SINGLE-PARTICLE SPACE.........................

      call GetSPS

      write(6,*)nsps(1),nsps(2)
2     continue
      write(6,*)' Enter number of valence protons '
      read(5,*)np
      if(np.gt.nsps(1))then
      	write(6,*)' too many protons ',np,nsps(1)
      	goto 2
      endif
22     continue
      write(6,*)' Enter number of valence neutrons '
      read(5,*)nn
      if(nn.gt.nsps(2))then
      	write(6,*)' too many neutrons ',nn,nsps(2)
      	goto 22
      endif


      if(allocated(psd))deallocate(psd)
      if(allocated(nsd))deallocate(nsd)
      if(allocated(rhop))deallocate(rhop)
      if(allocated(rhon))deallocate(rhon)
      if(allocated(qprot))deallocate(qprot)
      if(allocated(qneutr))deallocate(qneutr)

      allocate(psd(nsps(1),np),nsd(nsps(2),nn))
      allocate(rhop(nsps(1),nsps(1)),rhon(nsps(2),nsps(2)))
      allocate(qprot(nsps(1),nsps(1)),qneutr(nsps(2),nsps(2)))

c      write(6,*)' '
c      write(6,*)'Enter r2 for core nucleons as follows:'
c      write(6,*)'   36 for sd shell'
c      write(6,*)'  120 for pf shell'
c      write(6,*)'   all the other shells you have to calculate yourself'
c      read(5,*)r20
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine HF(hamflag,randomflag)
C
C  creates HF states starting either from random Slater Determinants
C   or from previously existing SD
C
C  SUBROUTINES CALLED:
C   	random_sd:      generates a random slater determinant to start with
C       make_hvec_iso:      creates vectors needed for computing
C                            Hamiltonian matrix elements 
C	save_sd		: saves slater determinants
C

      use spspace
      use HAMILTONIAN
      use wfn
      use nuclide
      use observables
      use obop
      use hfdiagold
      implicit none

      logical randomflag

      integer ndiag, idiag
C------------J-PROJECTION---------------------------------------
      real Pi
      integer i,j,n,m  !,m,mp		! dummies


      real hme
      real,allocatable          :: en(:)
      real                      :: psdmin(nsps(1),np),nsdmin(nsps(2),nn)
      real emin
      integer imin

      real over,protons


C--------Random generation of SD----------------------------------
      integer iseed,iseed_keep
      integer maxseed
      integer i_seed

C-------- Parameters----------------------------------------
      real gtol    ! the smallness parameter for gradient
      real test
      integer iter,err,irot,iwrite
      logical hamflag
      integer pade
C--------J^2------------------------------------------------
      real jzp_av,jzn_av
      real jx2,jy2,jz2
      real J2

C--------Deformation---------------------------------------------
      real beta,gamma
c      real r20

C..............FILE HANDLING..........................

      character*1 ychar
      logical writeout,writesd
      character filename*15  		! 
      integer ilast
      integer iout 			! choice of SD to keep
      integer tempfile			! location of temporary file
      logical saveflag
      logical errflag
      data tempfile/99/
     
C===============================================================

      rewind(tempfile)
      writeout=.false.

C...........READ IN HAMILTONIAN...................................
      if(.not.hamflag)then
         write(6,*)' Hamiltonian for minimization: '
         call  make_hvec_iso(.true.)
         hamflag=.true.
      endif
C.......Calculate the matrix corresponding to Jz......................
      if(.not.jxyz)then
         allocate(jx_matrix(nsps(1),nsps(1)),jy_matrix(nsps(1),nsps(1)))
         allocate(jz_matrix(nsps(1),nsps(1)))
         call make_Jz(jz_matrix,spsqn,nsps(1))
C.......Also, calculate Jx and Jy.....................................
         call make_Jxy(jx_matrix,jy_matrix,spsqn,nsps(1))
         jxyz=.true.
      endif

C......Read in a few options..........................................
      if(randomflag)then
      write(*,*)'Insert iseed'
      read(5,*)iseed_keep
      endif

      write(6,*)' Enter number of diagonalization steps '
      print*,'(If negative, use gradient descent)'
      read(5,*)ndiag

      if(ndiag<0)then
         write(6,*)'Insert Choice: 1-Pade Approx, 0-Exact'
         read(5,*)pade
      endif

      if(randomflag)then
      write(6,*)'Insert the number of different random SD'
      read(5,*)maxseed
      else
         maxseed=1
      endif
C      Jz = 0.0

3      continue

      write(6,108)
  108 format('  #    Ehf     J^2 (  Jx^2  Jy^2  Jz^2 ) ',
     &' beta  gamma (#iter) ')
      if(writeout)write(71,108)

C      do i=1,nsps(2)
C         do j=1,nsps(2)
C	    irhon(i,j)=0.0
C	 enddo
C      enddo
   
      i_seed=0
      iwrite=0
      allocate(en(maxseed))

      allocate(gammap_old(nsps(1),nsps(1)),gamman_old(nsps(2),nsps(2)))

      do while(i_seed.lt.maxseed)  ! start loop for generating different SDs
      iseed_keep=iseed_keep+12
      iseed=iseed_keep


C              INTIAL SLATER DETERMINANT
C................Generate a random Slater Determinant.................
      if(randomflag)then
      call random_sd(iseed,psd,nsd)
      call orth_SVD(psd,nsps(1),np)
      call orth_SVD(nsd,nsps(2),nn)
      endif

C......................INITIAL INFORMATION...................
      call make_projectors(psd,rhop,qprot,np,nsps(1),overp)
      call make_projectors(nsd,rhon,qneutr,nn,nsps(2),overn)

      over=overp*overn

      call hmult(nsps,nmatpp,nmatnn,nmatpn,
     &  map_combpp,map_combnn,map_combpn,
     &  hvecpp,hvecnn,hvecpn,e_spepp,e_spenn,
     &  overp,overn,rhop,rhon,over,hme)

C      print*,hme

C...J^2, another method
      call calculateJ2(Jx_matrix,Jy_matrix,Jz_matrix,nsps,
     &                        rhop,rhon,Jx2,jy2,jz2)
      j2 = jx2+jy2+jz2
c      print*,j2

C..........PRINT OUT THE RESULTS.......................................

      call deform(beta,gamma)
      if(maxseed.eq.1)then
      write(6,*)' Initial state: '
      write(6,299)0,hme,j2,jx2,jy2,jz2,beta,gamma,
     &                                        0
      endif
C.................MINIMIZATION BY DIAGONALIZATION ......................

      diag1=.false.
      if(ndiag.le.0)then
         do idiag = 1,20
	   call HFdiag(hme)
C            write(60,*)idiag,hme
         enddo
C         print*,hme
C        call HFgrdesc(hme,pade)
        call prjgrd(pade,hme,iter,err)
        if(err<0)goto 10
      else
         iter=ndiag
         do idiag = 1,ndiag
	   call HFdiag(hme)
C            write(60,*)idiag,hme
         enddo
      endif


33    continue
      i_seed=i_seed+1
      call make_projectors(psd,rhop,qprot,np,nsps(1),overp)
      call make_projectors(nsd,rhon,qneutr,nn,nsps(2),overn)
C...J^2
      call calculateJ2(Jx_matrix,Jy_matrix,Jz_matrix,nsps,
     &                        rhop,rhon,Jx2,jy2,jz2)
      j2 = jx2+jy2+jz2

C..........PRINT OUT THE RESULTS.......................................

      en(i_seed)=hme
      call deform(beta,gamma)
      write(6,299)i_seed,hme,j2,jx2,jy2,jz2,beta,gamma,
     &                                        iter
      if(writeout)then
           write(71,299)i_seed,hme,j2,jx2,jy2,jz2,beta,gamma,iter
      endif
C...............WRITE TO TEMPORARY FILE.................

	call writeoutsd(nsps(1),np,tempfile,psd)
	call writeoutsd(nsps(2),nn,tempfile,nsd)

10    continue
      enddo  ! loop over i_seed
      deallocate(gammap_old,gamman_old)

98    format(2(G12.6,1X))
99    format('i_seed=',I3,1X,'iter=',I5,1X,'Energy=',G12.6,1X,'J=',F7.4,
     &       1X,'Jz_p=',F7.4,1X,'Jz_n=',F7.4)
  199 format(' E = ',g12.6,1x,' J = ',f7.4,2x,
     7 ' beta= ',F7.4,1X,'  gamma= ',F7.4,1x,' (#iter = ',i5,' )')
  299 format(i3,1x,f8.3,f7.3,'(',3(1x,f5.2),' )',2(f6.3,1X),'(',i5,')')
  119 format(10(2x,F10.5))
      close(unit=71)
      close(unit=77)
      close(unit=78)
C..................................................choose a state

333   continue
      saveflag=.false.

      if(maxseed.eq.1)then
	iout = 1
      else
      write(6,*)' ' 
      write(6,*)' Choose one of these for further calculation '
      if(.not.saveflag)then
      write(6,*)' (enter [-1] to choose lowest, or [0] to save ',
     & 'all to file) '
      endif
      read(5,*)iout 
      endif
      if(iout.eq.-1)then
C           print*,'Choosing minimum energy...'
           emin=en(1)
           iout=1
          do i = 2,maxseed
               if(en(i)<emin)then
                  iout=i
                  emin=en(i)
               endif
          enddo
      endif
C      print*,en

      if(iout.eq.0)then
        if(saveflag)then
	   write(6,*)' Already saved! '
	   goto 333
	endif

	call save_sd(tempfile)
	saveflag=.true.
	write(6,*)' Must still choose one for further calculations'
	goto 333
      else
	if(iout.lt.0 .or.iout.gt. maxseed )then
	    write(6,*)' Must be less than or equal to ',maxseed
	    goto 333
        endif
        rewind(tempfile)

C        print*,'IOUT=',iout
        do i = 1,iout
	  call readinsd(nsps(1),np,tempfile,psd,errflag)
	  if(errflag)stop
	  call readinsd(nsps(2),nn,tempfile,nsd,errflag)
	  if(errflag)stop

        enddo

      endif

      Ehf=En(iout)

      call make_projectors(psd,rhop,qprot,np,nsps(1),overp)
      call make_projectors(nsd,rhon,qneutr,nn,nsps(2),overn)
      deallocate(en)
C      return

        call hmult(nsps,nmatpp,nmatnn,nmatpn,
     &  map_combpp,map_combnn,map_combpn,
     &  hvecpp,hvecnn,hvecpn,e_spepp,e_spenn,
     &  overp,overn,rhop,rhon,over,hme)
C      print*,Ehf,hme


      return

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine save_sd(tempfile)
C
C  saves slater determinant wfn to a file
C  also saves s.p. information, plus a header
C  has option to write out multiple SDs from tempfile
C
C  INPUT
C     tempfile			! where temp SDs saved
C				! if = 0, no temp file
C  SUBROUTINES CALLED
C	writeoutsd
C	readinsd
C
      use spspace
      use wfn
      use nuclide

      implicit none

      real,allocatable           :: sdtemp(:,:)
C..............FILE HANDLING..........................

      character*1 ychar
      character filename*15  		! 
      integer ilast
      integer tempfile			! location of temporary file
      logical append			! controls appending
      logical errflag

C..............MISC........................

      integer i,ii,j,n,z
      character title*60

311    continue
      write(6,*)' Enter output filename (.sd)'
      read(5,'(a)')filename
      ilast=index(filename,' ')
      if(ilast.ne.0)then
         ilast=ilast-1
      else
         ilast=15
      endif

      append = .false.
      open(unit=2,file=filename(1:ilast)//'.sd',status='new',err=33,
     &      form='unformatted')
      goto 44
33    continue
      write(6,*)' That file already exists;',
     & ' overwrite (o) or append (a) or new (n) file?(o/n/a) '
C
C IF A FILE ALREADY EXISTS YOU CAN OVERWRITE OR APPEND

      read(5,'(a)')ychar
      if(ychar.eq.'O' .or. ychar.eq.'o')then
	open(unit=2,file=filename(1:ilast)//'.sd',status='unknown' ,
     &       form='unformatted')
      elseif(ychar.eq.'A' .or. ychar.eq.'a')then
        append = .true.
	open(unit=2,file=filename(1:ilast)//'.sd',status='old',
     &       form='unformatted')
      else
	goto 311
      endif
44    continue

C
C AS LONG AS NOT APPENDING, WRITE DOWN S.P. INFORMATION
C ASSUME PROTON, NEUTRON SPACES THE SAME
C
      if(.not.append)then
      	write(6,*)'(I assume the proton, neutron spaces the same)'
	write(2)norb(1)
	do i =1,norb(1)
	   write(2)i,nrorb(i),jorb(i)
        enddo
	write(2)np,nn
        write(6,*)' Enter a comment '
        read(5,'(60a)')title
C        ilast=index(title,' ')
        write(2)title
      else
C
C  IF APPEND, THEN CHECK TO MAKE SURE INFO MATCHES
C
 	read(2)n
	if(n.ne.norb(1))then
	    write(6,*)' inconsistent orbit info ',
     &  norb(1),n
	    close(unit=2)
	    goto 311	    
        endif
	do i =1,norb(1)
	   read(2)ii,n,j
	   if(n.ne.nrorb(i).or. j.ne.jorb(i))then
		write(6,*)' inconsistent orbits ',
     & i,nrorb(i),n,jorb(i),j
		close(unit=2)
		goto 311
	   endif
        enddo
	read(2)z,n
        if(z.ne.np .or. n.ne.nn)then
		write(6,*)' wrong nuclide ',
     & np,z,nn,n
		close(unit=2)
		goto 311
	endif
        read(2)title
        write(6,*)title
C............ now readout to dummy array.............
	do
          allocate(sdtemp(nsps(1),np))
	  call readinsd(nsps(1),np,2,sdtemp,errflag)
          deallocate(sdtemp)
	  if(errflag)goto 443
          allocate(sdtemp(nsps(1),nn))
	  call readinsd(nsps(2),nn,2,sdtemp,errflag)
          deallocate(sdtemp)
	  if(errflag)goto 443
        enddo
443     continue

      endif
      if(tempfile.eq.0)then
	call writeoutsd(nsps(1),np,2,psd)
	call writeoutsd(nsps(2),nn,2,nsd)
      else
        rewind(tempfile)
	do i = 1,1000
	  call readinsd(nsps(1),np,tempfile,psd,errflag)
	  if(errflag)goto 444
 	  call writeoutsd(nsps(1),np,2,psd)
	  call readinsd(nsps(2),nn,tempfile,nsd,errflag)
	  if(errflag)goto 444
 	  call writeoutsd(nsps(2),nn,2,nsd)

        enddo
444     continue
      endif

      close(unit=2)
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine read_sd
C
C  reads slater determinant wfn from file
C
C  INPUT
C     hamflag			! a Hamiltonian already read in
C
C  SUBROUTINES CALLED
C	readinsd
C
      use spspace
      use nuclide
      use wfn
      use Hamiltonian
      use observables
      use flags
      implicit none

C--------Jz averages------------------------------------------------
      real jzp_av,jzn_av

C--------Jx and Jy (J^2)-------------------------------------------
      real jx2,jy2,jz2
      real J2

C--------Deformation---------------------------------------------
      real                   :: beta,gamma

C..............FILE HANDLING..........................

      character*1 ychar
      character filename*15  		! 
      integer ilast
      integer tempfile			! location of temporary file
      data tempfile/99/
      logical errflag


C..............MISC........................

      integer i,ii,j,n,m,iout,z
      character title*60
      real,allocatable     :: psdtemp(:,:),nsdtemp(:,:)
      real                 :: over,hme

311    continue
      write(6,*)' Enter input filename (.sd)'
      read(5,'(a)')filename
      ilast=index(filename,' ')
      if(ilast.ne.0)then
         ilast=ilast-1
      else
         ilast=15
      endif

      open(unit=2,file=filename(1:ilast)//'.sd',status='old',err=33,
     &      form='unformatted')
      goto 44
33    continue
      write(6,*)' That file does not exist; ',
     &'do you wish to try another file (y/n)?'
      read(5,'(a)')ychar
      if(ychar.eq.'n' .or.ychar.eq.'N')then
	return
      else
      	goto 311
      endif
44    continue

C..............READ IN HEADER INFO...............

      read(2)n
      if(n.ne.norb(1))then
	write(6,*)' # of orbits mismatch ',n,norb(1)    
        goto 302
      endif
      do i = 1,norb(1)
	read(2)ii,n,j
	if(n.ne.nrorb(i) .or. j.ne.jorb(i))then
	  write(6,*)' mismatch n,j:',n,nrorb(i),j,jorb(i)
	  goto 302
        endif
      enddo
      goto 336

302     continue
        write(6,*)' The single particle space does not match ',
     &' with that previously chosen.'
322     continue
        write(6,*)' Exit (x) or choose another file (c)?'
        read(5,'(a)')ychar
        if(ychar.eq.'x' .or. ychar.eq.'X')then
	    close(unit=2)
	    return
        endif
        if(ychar.eq.'c' .or. ychar.eq.'C')then
	  close(unit=2)
	  goto 311
        else 		
 	  write(6,*)' That selection not valid '
 	  goto 322

        endif
336   continue

C...............CHECK IF N,Z match........

      read(2)z,n
      if(n.ne.nn .or. z.ne.np)then
	write(6,*)' mismatch Z,N. Old: ',np,nn,
     & ', new: ',z,n
	write(6,*)' Do you want to use new values?'
	read(5,'(a)')ychar
	if(ychar.eq.'y' .or.ychar.eq.'Y')then
	    np=z
	    nn=n
        else
	    goto 322
        endif      
      endif

      read(2)title

      ilast=index(title,' ')
      if(ilast>1)then
        ilast=ilast-1
        write(6,*)'Title card: "',title(1:ilast),'"'
      else
       print*, 'no title card'
      endif

C.............CHOOSING A NEW HAMILTONIAN...............

      ychar='n'
      if(.not.hamflag)then
	write(6,*)' You MUST choose a Hamiltonian '
      else
	write(6,*)' A Hamiltonian has been read in. ',
     & 'Do you want to use a new one? (y/n)?'
	read(5,'(a)')ychar
      endif

      if(.not.hamflag .or. ychar.eq.'Y' .or.ychar.eq.'y')then

      call  make_hvec_iso(.true.)
	hamflag=.true.
      endif 

C................NOW READ IN...FIRST PREPARE

C.......Calculate the matrix corresponding to Jz......................
      if(.not.jxyz)then
         allocate(jz_matrix(nsps(1),nsps(1)),jx_matrix(nsps(1),nsps(1)),
     &                    jy_matrix(nsps(1),nsps(1)))
         call make_Jz(jz_matrix,spsqn,nsps(1))

C.......Also, calculate Jx and Jy.....................................
         call make_Jxy(jx_matrix,jy_matrix,spsqn,nsps(1))
         jxyz=.true.
      endif

C.................FINALLY READ IN UNTIL DONE............

      iout = 0
      rewind(tempfile)

      write(6,108)
  108 format('  #     <H>    J^2 (  Jx^2  Jy^2  Jz^2 ) ',
     &' beta  gamma ')

      allocate(psdtemp(nsps(1),np),nsdtemp(nsps(2),nn))
      do
	call readinsd(nsps(1),np,2,psdtemp,errflag)
	if(errflag)goto 443
	call readinsd(nsps(2),nn,2,nsdtemp,errflag)
	if(errflag)goto 443
        iout = iout+1

      call make_projectors(psdtemp,rhop,qprot,np,nsps(1),overp)
      call make_projectors(nsdtemp,rhon,qneutr,nn,nsps(2),overn)

        call hmult(nsps,nmatpp,nmatnn,nmatpn,
     &  map_combpp,map_combnn,map_combpn,
     &  hvecpp,hvecnn,hvecpn,e_spepp,e_spenn,
     &  overp,overn,rhop,rhon,over,hme)


C...J^2...............
      call calculateJ2(Jx_matrix,Jy_matrix,Jz_matrix,nsps,
     &                        rhop,rhon,Jx2,jy2,jz2)
      j2 = jx2+jy2+jz2

C..........PRINT OUT THE RESULTS.......................................
      call deform(beta,gamma)
      write(6,299)iout,hme,j2,jx2,jy2,jz2,beta,gamma
  299 format(i3,1x,f8.3,f7.3,'(',3(1x,f5.2),' )',2(f6.3,1X))

C...............WRITE TO TEMPORARY FILE.................

	call writeoutsd(nsps(1),np,tempfile,psdtemp) 
	call writeoutsd(nsps(2),nn,tempfile,nsdtemp) 

      enddo
      write(6,*)' should not have gotten here '
      stop
443   continue
      if(iout.gt.1)then
333     continue
	write(6,*)' Which state do you wish? '
	read(5,*)j

	if(j.lt.1 .or.j.gt. iout )then
	    write(6,*)' Must be less than or equal to ',iout
	    goto 333
        endif
        rewind(tempfile)
        do i = 1,j
           call readinsd(nsps(1),np,tempfile,psd,errflag)
           call readinsd(nsps(2),nn,tempfile,nsd,errflag)

        enddo

      else
	   do n = 1,nsps(1)
		do m = 1,np
		  psd(n,m)=psdtemp(n,m)
		enddo
           enddo
	   do n = 1,nsps(2)
		do m = 1,nn
		   nsd(n,m)=nsdtemp(n,m)
 		enddo
	   enddo
	
      endif
      call make_projectors(psd,rhop,qprot,np,nsps(1),overp)
      call make_projectors(nsd,rhon,qneutr,nn,nsps(2),overn)

      deallocate(psdtemp,nsdtemp)

        call hmult(nsps,nmatpp,nmatnn,nmatpn,
     &  map_combpp,map_combnn,map_combpn,
     &  hvecpp,hvecnn,hvecpn,e_spepp,e_spenn,
     &  overp,overn,rhop,rhon,over,hme)
      print*,hme


      close(unit=2)
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine check_energy
C
C  reads energy, etc. of current wfn
C
C  INPUT
C     hamflag			! a Hamiltonian already read in
C
C  SUBROUTINES CALLED
C
      use spspace
      use Hamiltonian
      use wfn
      use observables
      use flags
      use nuclide
      implicit none
C      include 'gcb_dim.inc'
      
      real hme,over

C-------------SLATER DETERMINANTS------------------------------
   
C      integer,intent(IN)              :: np,nn		! # of protons, neutrons

C--------Jz averages------------------------------------------------
      real jzp_av,jzn_av

C--------Jx and Jy (J^2)--------------------------------------------
      real jx2,jy2,jz2
      real J2

C--------Deformation---------------------------------------------
      real beta,gamma

C..............FILE HANDLING..........................

      character*1 ychar
      character filename*15  		! 
      integer ilast
      integer tempfile			! location of temporary file
      data tempfile/99/
      logical errflag

C..............MISC........................

      integer i,ii,j,n,m,iout,z
      character title*60
     

C + + + + + + + + + + + + + + + + + + + + + + +



C.............CHOOSING A NEW HAMILTONIAN...............

      ychar='n'
      if(.not.hamflag)then
	write(6,*)' You MUST choose a Hamiltonian '
      else
	write(6,*)' A Hamiltonian has been read in. ',
     & 'Do you want to use a new one? (y/n)?'
	read(5,'(a)')ychar
      endif

      if(.not.hamflag .or. ychar.eq.'Y' .or.ychar.eq.'y')then
        call  make_hvec_iso(.true.)
	hamflag=.true.
      endif 

c      call make_projectors(psd,rhop,qprot,np,nsps(1),overp)
c      call make_projectors(nsd,rhon,qneutr,nn,nsps(2),overn)

        call hmult(nsps,nmatpp,nmatnn,nmatpn,
     &  map_combpp,map_combnn,map_combpn,
     &  hvecpp,hvecnn,hvecpn,e_spepp,e_spenn,
     &  overp,overn,rhop,rhon,over,hme)


C...J^2, another method
      call calculateJ2(Jx_matrix,Jy_matrix,Jz_matrix,nsps,
     &                        rhop,rhon,Jx2,jy2,jz2)
      j2 = jx2+jy2+jz2

C..........PRINT OUT THE RESULTS.......................................
      call deform(beta,gamma)

      write(6,108)
  108 format('        <H>    J^2 (  Jx^2  Jy^2  Jz^2 ) ',
     &' beta  gamma ')

      write(6,299)hme,j2,jx2,jy2,jz2,beta,gamma
  299 format(4x,f8.3,f7.3,'(',3(1x,f5.2),' )',2(f6.3,1X))


      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine HFobs
C
C  computes HF expectation value of a scalar 1+2 body observable
C
      use observables
      use wfn
      use spspace
      use nuclide
      use flags
C      use HAMILTONIAN
      implicit none

      real hme,over
      character*1 obschar

c      print*,orb_qn

      print*,nmatpp_O,nmatnn_O,nmatpn_O

      if(.not.obsflag)then
         write(6,*)' Enter file for scalar observable '
         call  make_hvec_iso(.false.)
         obsflag=.true.
      else
         print*,' An observable has been alsready read. Do you',
     &          ' want to read another one? (y/n)'
         read(5,*)obschar
         if(obschar=='y' .or. obschar=='Y')then
             call  make_hvec_iso(.false.)
             obsflag=.true.
         endif
      endif
      
c      print*,nmatpp_O,nmatnn_O,nmatpn_O

c      call make_projectors(psd,rhop,qprot,np,nsps(1),overp)
c      call make_projectors(nsd,rhon,qneutr,nn,nsps(2),overn)


c      print200,hvecpp(1:10)
c      print200,hvecpp_O(1:10)
c      print200,hvecnn(1:10)
c      print200,hvecnn_O(1:10)
C      print200,hvecpn(2000:2010)
c      print200,hvecpn_O(2000:2010)

c       print*,nsps

 200  format(10F12.4)

c      write(21,*)rhop

        call hmult(nsps,nmatpp_O,nmatnn_O,nmatpn_O,
     &  map_combpp_O,map_combnn_O,map_combpn_O,
     &  hvecpp_O,hvecnn_O,hvecpn_O,e_spepp_O,e_spenn_O,
     &  overp,overn,rhop,rhon,over,obsval)

c        call hmult(nsps,nmatpp,nmatnn,nmatpn,
c     &  map_combpp,map_combnn,map_combpn,
c     &  hvecpp,hvecnn,hvecpn,e_spepp,e_spenn,
c     &  overp,overn,rhop,rhon,over,obsval)

      write(6,122)obsval
 122  format(' Expectation Value=',F9.3)
C      print*,orb_qn

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC









