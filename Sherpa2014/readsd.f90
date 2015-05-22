!
!  this code translates a .sd file between formatted and 
!  unformatted
!

module sddata

implicit none

character title*60
integer norb, np,nn,nsps
integer, allocatable :: nrorb(:), jorb(:)
real, allocatable :: psd(:,:),nsd(:,:)

end module sddata

!=================================================

      implicit none

      real,allocatable           :: sdtemp(:,:)
!..............FILE HANDLING..........................

      character*1 ychar
      character filename*15             !
      integer ilast
      integer tempfile                  ! location of temporary file
      logical append                    ! controls appending
      logical errflag
      logical formatit
!..............MISC........................

      integer i,ii,j,n,z
      character title*60

      write(6,*)' Enter filename (.sd /.fsd)'
      read(5,'(a)')filename
      ilast=index(filename,' ')
      if(ilast.ne.0)then
         ilast=ilast-1
      else
         ilast=15
      endif

      print*,' Do you want to convert unformatted to formatted? (y/n) '
      print*,' (If no (n) then convert formatted to unformatted '
      read(5,'(a)')ychar
      if(ychar == 'y' .or. ychar == 'Y')then
         formatit = .true.
         print*,' Reading from .sd file and formatting as .fsd file '
      else
         formatit = .false.
         print*,' Reading from formatted .fsd file and writing as  unformatted.sd file '

      endif

      if(formatit)then

      open(unit=2,file=filename(1:ilast)//'.sd',status='old',err=33,form='unformatted')
      open(unit=22,file=filename(1:ilast)//'.fsd',status='unknown')
      call readinsd(2)

       else

      open(unit=2,file=filename(1:ilast)//'.sd',status='unknown',form='unformatted')
      open(unit=22,file=filename(1:ilast)//'.fsd',status='old',err=33)
      call writeoutsd(2)

      endif

      stop
33    continue  
      print*,' File does not exist '
      end
!===========================================
      subroutine readinsd(ifile)

      use sddata
      implicit none
!      include 'gcb_dim.inc'


      integer i,j
      integer ifile
      logical errflag

      errflag=.false.

	read(ifile)norb
        write(22,*)norb
        allocate(nrorb(norb), jorb(norb))
        nsps = 0
	do i =1,norb
	   read(ifile)j,nrorb(i),jorb(i)
           write(22,*)i,nrorb(i),jorb(i)
           nsps = nsps +jorb(i)+1

        enddo
	read(ifile)np,nn
        read(ifile)title
        write(22,'(60a)')title
        print*,title


        if(np > 0)allocate(psd(nsps,np))
        if(nn > 0)allocate(nsd(nsps,nn))

      do i = 1,np
        read(ifile,err=103,end=103)(psd(j,i),j=1,nsps)
        write(22,100)(psd(j,i),j=1,nsps)
100   format(20(F12.6,1x))

      enddo
      do i = 1,nn
        read(ifile,err=103,end=103)(nsd(j,i),j=1,nsps)
        write(22,100)(nsd(j,i),j=1,nsps)
      enddo
      return
  103 continue
      print*,' oops '
      stop

      end

!----------------------------------------------------

      subroutine writeoutsd(ifile) 

      use sddata
      implicit none
!      include 'gcb_dim.inc'

      

      integer i,j
      integer ifile
      
      read(22,*)norb
	write(ifile)norb
      allocate(nrorb(norb),jorb(norb))
        nsps = 0
	do i =1,norb
           read(22,*)j,nrorb(i),jorb(i)
            nsps = nsps + jorb(i)  + 1
	   write(ifile)i,nrorb(i),jorb(i)
        enddo
        read(22,*)np,nn
	write(ifile)np,nn
        read(22,'(60a)')title
        print*,title
        write(ifile)title
        if(np > 0)allocate(psd(nsps,np))
        if(nn > 0)allocate(nsd(nsps,nn))
      do i = 1,np
        read(22,*)(psd(j,i),j=1,nsps)
        write(ifile)(psd(j,i),j=1,nsps)
      enddo
      do i = 1,nn
        read(22,*)(nsd(j,i),j=1,nsps)

        write(ifile)(nsd(j,i),j=1,nsps)
      enddo
100   format(20(F12.6,1x))
      close(22)
      close(2)
      return
      end
      

 