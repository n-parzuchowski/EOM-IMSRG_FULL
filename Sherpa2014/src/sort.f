      program sortStr

      implicit none
      integer nmax
      parameter(nmax=5000)
      real en(nmax),str(nmax)
      integer ilast
      character*20 filename
      integer n
      character smRPA*1
      integer nfiles
      real*8 sumStr
      real*8 S1,S2
      real*8 CTR,WDT

      integer i,nfc,j
      real e,s,egs
 
      nfiles=1
      nfc=1
      n=0

    2 continue
      write(6,*)'SM or RPA distribution strengths? (s/r)'
      read(5,*)smRPA

    3 continue
      if(smRPA.eq.'s' .or. smRPA.eq.'S')then

        write(6,*)'Introduce 1st filename with distribution',
     1                                       ' from SM(.str)'
        if(nfc.eq.1)write(6,*)'WARNING: this file should contain',
     1                ' the ground state energy on the last line'
      elseif(smRPA.eq.'r' .or. smRPA.eq.'R')then
        write(6,*)'Introduce filename with RPA distribution (.str).'
      else
        write(6,*)'Invalid choice. Please try again.'
        goto 2
      endif
1     continue
      read(5,'(a)')filename
      ilast=index(filename,' ')
      if(ilast.ne.0)then
         ilast=ilast-1
      else
         ilast=15
      endif
      open(unit=61,file=filename(1:ilast)//'.str',status='old',
     &err=131)
      if(nfc.ne.1)goto 132
      open(unit=62,file=filename(1:ilast)//'.smo.dat',status='unknown',
     &err=131)
        goto 132

  131 continue
        write(6,*)' That file does not exist '

          goto 1
  132 continue

      do i=1,1000
        read(61,*,ERR=10,END=20)e,s
        n=n+1
        if(n.gt.nmax)stop 'Increase nmax.'
        en(n) =e
        str(n)=s
C        if(smRPA.eq.'s' .or. smRPA.eq.'S')str(n)=str(n)/2.0
   10   continue
      enddo
   20 continue
      close(61)
      if(nfiles.eq.1)then
      if(smRPA.eq.'S' .or. smRPA.eq.'s')then
         Egs=en(n)
      else
         Egs=0.0
      endif

      write(6,*)'First file processed. How many others do you want?'
      read(5,*)nfiles
      endif
      nfc=nfc+1
      if(nfc.le.nfiles)goto 3

      call bubbleMethod(en,str,nmax,n)

      sumStr=0.0
      S1=0.0
      S2=0.0
      write(62,100)en(1)-Egs,sumStr,sumStr
      do i=1,n
        write(62,100)en(i)-Egs,str(i),sumStr
        sumStr=sumStr+str(i) !/2.
        S1=S1+(en(i)-Egs)*str(i) !/2.
        S2=S2+((en(i)-Egs)**2)*str(i) ! /2.
        write(62,100)en(i)-Egs,str(i),sumStr
      enddo
100   format(3(2x,F10.3))
      write(6,101)sumStr,S1,S2
101   format('B(GT)=',F10.2,', S1=',F10.2,', S2=',F10.2)
      CTR=S1/sumStr
      WDT=S2/sumStr-CTR**2
      WDT=dsqrt(WDT)
      write(6,102)CTR,WDT
102   format('CENTROID=',F10.2,', WIDTH=',F10.2)

      close(62)

      end



      subroutine bubbleMethod(array,val,nmax,ndim)

      implicit none
      integer ndim,nmax,i,st
      real array(nmax),val(nmax),tmp,tmpval

      st=1
      do while(st.gt.0)
         st=0
         do i=1,ndim-1
            if(array(i).gt.array(i+1))then
               tmp = array(i)
               tmpval=val(i)
               array(i) = array(i+1)
               array(i+1) = tmp
               val(i)=val(i+1)
               val(i+1)=tmpval
               st=st+1
	    endif
         enddo
      enddo
      
      return
      end
      
