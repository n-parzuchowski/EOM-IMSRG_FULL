CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  original routines written by W.E. Ormand
C  modified by CWJ   March/April 1998  LSU
C
C  revisions for isospin  CWJ April 2001 LSU
C
C  revision for f90  IS May 2004 UA
C
C  2 small bugs fixed in make_hvec_iso Feb 2008 by CWJ
C

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine make_hvec_iso(ham)
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C          read in hamiltonian file IN ISOSPIN FORMALISM
C
c          Subroutine to put Hamiltonian into the m-scheme
c          returns only non-zero combinations in the m-scheme
c          hvec returns the m-schem matrix element <i1,i2|V|i3,i4>
c          nmat is the number of nonzero elements
c          map_comb=ipack(i,j,k,l) is a packed integer labeling the M.E.
c          Matrix element of two-body part of Ham is 
c          <SD_1|H|SD_2>=sum[(i<j),(k<l)]<i,j|V|k,l>*
c              <SD_1|a+(i)a+(j)a(k)a(l)|SD_2>
c          one-body part of Ham is stored in the vector e_spe
C
C  INPUT:
C	norb = # of j-shells
C	nsps = # of m-states
C	orbqn = quantum #'s of j-shells
C	spsqn  = quantum #'s of m-states
C
C  OUTPUT:	
C	nmatpp,nmatnn,nmatpn:
C   -- the routine decodes the hamiltonian into all m-scheme TBME's
C   nmatpp = # of nonzero pp matrix elements; similar for nn, pn
C      
C	map_combpp(i): list of (packed) nonzero indices, i = 1 to nmatpp
C   to get indices IJKL of i'th matrix element do UNPACK(map_combpp)
C	hvecpp(i) = ith TBME associated with indices unpacked from map_combpp(i)
C	similar for map_combnn, map_combpn, hvecnn, hvecpn
C	e_spepp(i,j) = one-body matrix element between i,j m-states
C	e_spenn(i,j)
 
C  SUBROUTINES CALLED:
C	ipack: packs 4 small integers ijkl into one integer
C 	unpack: turns integer X into 4 integers ijkl 
C

      use SPspace
      use chf_dim
      use HAMILTONIAN
      use observables
      implicit none

C OUTPUT

      logical,intent(IN)    :: ham   ! = true, computes the Hamiltonan
                        ! = false computes observables
C-----------HAMILTONIAN------------------------------------------
      integer    :: nmatpp_p		! # of pp matrix elements 
      integer    :: nmatnn_p		! # of nn matrix elements
      integer    :: nmatpn_p		! # of np matrix elements

C-------- list of packed 4-indices with nonzero TBME     
c      integer,pointer,dimension(:)  :: map_combpp_p,map_combnn_p
c      integer,pointer,dimension(:)  :: map_combpn_p

C---------- nonzero TBMEs associated with above list   
c      real,pointer,dimension(:)     :: hvecpp_p,hvecnn_p,hvecpn_p
      
C---------one-body matrix elements      
c      real,pointer,dimension(:,:)  :: e_spepp_p,e_spenn_p

c==============  data used internally
c==============  interaction file name
      character input*15
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c============== define dimensions of ej
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real,allocatable,dimension(:,:,:,:,:,:)  :: ej   ! in pn formalism
      real,allocatable,dimension(:,:,:,:,:,:)  :: vjt  ! in isospin
      
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c============== single-particle enrgies
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real,allocatable  :: spe(:)

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C                 external field (NEW 2007)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      logical success
      character ychar*1,title*70
      real, allocatable :: vextp(:,:),vextn(:,:)
      real :: vstrength
      integer nspstest,jtest,mtest
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c============== phases and matrix element read in
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real phab,phcd,vvv
      integer iaabb,iccdd
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c============== integers 
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer nint,nspe
      integer i1,i2,i3,i4,ia,ib,ic,id,jj,iso
      integer i,j,k,l,ii,ilast,index
      integer jabmax,jabmin,jcdmax,jcdmin,jmax,jmin,jr,is,ismin
      integer mij,mkl,iij,ikl,ma,mb,mc,md,iza,izb,izc,izd
      integer ja,jb,jc,jd
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c==============  Integers to check if orbits are the same for pn part
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer ichkab,ichkcd
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c============== scaling parameters for Ham. (A/B)**X 
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real A,B,X,scale,xspe,sum,sum2
      real,allocatable    :: temp(:)

      integer n_ham	! number of interactions to be read in
      integer i_ham
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c============== External functions
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer ipack
      real cleb                           ! clebsch-gordon coeff.

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c==============  isospin to pn conversion
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C---- SINGLE-PARTICLE LIST------

      integer,allocatable    :: jt(:),tt(:)
      integer,allocatable    :: lt(:),nt(:)

      integer,allocatable    :: jtot(:)		! spin
      integer,allocatable    :: mout(:,:)	! set of 4 states
      integer                :: nmepn

      real xme
      integer ttot
      integer,allocatable    :: map(:)
      
      integer ta,tb,tc,td,tkl,tij   ! ******see below***
      real factkeep

      logical          :: countonly


c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c==============  Standard useful constants
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real zero,one,two,fact,factor,small
C      integer j99,i99,l99,k99

      zero=0.0e0
      one=1.0e0
      two=2.0e0
      small=1.0e-6

C-------------- ERROR TRAP -------
C -------------- # of proton and neutron single-particle energies must be =

      if(nsps(1).ne.nsps(2))then
        write(6,*)' proton states neq neutron states ',nsps(1),nsps(2)
        stop
      endif
C---------------END OF ERROR TRAP --------------------------
      nspe=norb(1)
C      print*,'nspe=',nspe
      allocate(vjt(0:j2max,0:1,1:nspe,1:nspe,1:nspe,1:nspe))
      allocate(ej(0:j2max,0:1,1:2*nspe,1:2*nspe,1:2*nspe,1:2*nspe))
      allocate(spe(norb(1)+norb(2)))

C      deassociate(hvecpp_p,hvecnn_p,hvecpn_p)
C      deassociate(map_combpp_p,map_combnn_p,map_combpn_p)

c      deallocate(e_spepp_p,e_spenn_p,hvecpp_p,hvecnn_p,hvecpn_p)
c      deallocate(map_combpp_p,map_combnn_p,map_combpn_p)


C....Make sure the matrices to store Hamiltonian/Observable are not already allocated.....
      if(ham)then
         if(allocated(e_spepp))deallocate(e_spepp)
         if(allocated(e_spenn))deallocate(e_spenn)
         if(allocated(hvecpp))deallocate(hvecpp)
         if(allocated(hvecnn))deallocate(hvecnn)
         if(allocated(hvecpn))deallocate(hvecpn)
         if(allocated(map_combpp))deallocate(map_combpp)
         if(allocated(map_combnn))deallocate(map_combnn)
         if(allocated(map_combpn))deallocate(map_combpn)
      else
         if(allocated(e_spepp_O))deallocate(e_spepp_O)
         if(allocated(e_spenn_O))deallocate(e_spenn_O)
         if(allocated(hvecpp_O))deallocate(hvecpp_O)
         if(allocated(hvecnn_O))deallocate(hvecnn_O)
         if(allocated(hvecpn_O))deallocate(hvecpn_O)
         if(allocated(map_combpp_O))deallocate(map_combpp_O)
         if(allocated(map_combnn_O))deallocate(map_combnn_O)
         if(allocated(map_combpn_O))deallocate(map_combpn_O)
      endif
C      nullify(e_spepp_p,e_spenn_p,hvecpp_p,hvecnn_p,hvecpn_p)
C      nullify(map_combpp_p,map_combnn_p,map_combpn_p)

      spe=0.0
      vjt=0.0

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c==============  get the two-body matrix elements
c                fill in the array ej(J,T,ia,ib,ic,id)
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      nspe=norb(1)
      allocate(temp(nspe))

      write(6,*)' Enter number of interactions to read in '
      read(5,*)n_ham
      if(n_ham .le. 0)n_ham = 1

      do i_ham = 1,n_ham

1     continue
      write(6,*)'Enter interaction file name (.int)'
      write(6,*)' ATTN: I ASSUME INTERACTION IS IN ISOSPIN FORM '

      if(i_ham.gt.1)then
           write(6,*)' Enter "null" if you do not want another '
      endif

      read(5,'(a)')input
      if(input.eq.'null')goto 22
      ilast=index(input,' ')-1
      open(unit=1,file=input(1:ilast)//'.int',status='old',err=2)
      goto 3
2     continue
      write(6,*)' That file does not exist '
      goto 1
3     continue

      write(6,*)'Enter parameters for scaling (A/B)**X,'
     &, 'single-particle energies '
      write(6,*)' (if B=0 then scale = 1) '
      read(5,*)A,B,X,xspe
      if(b.eq.0.)then
	scale = 1.0
      else
       scale=(A/B)**X       
      endif

c      scale = 1.0
c      if(scale.ne.1.0)write(6,*)' scale = ',scale,xspe
c      if(xspe.eq.0.)then
c	write(6,*)' single-particle set = 0'
c      endif

      read(1,*)nint,(temp(i),i=1,nspe)
      write(6,*)nspe,(temp(i)*xspe,i=1,nspe)
      do i =1,nspe
	spe(i)=spe(i)+temp(i)*xspe
      enddo



      do i = 1,norb(2)
        spe(i+norb(1))=spe(i+norb(1))+temp(i)*xspe   ! duplicate single particle energies
      enddo
C++++++++++++++++++++++++++++END OF ERROR TRAP++++++++++++++++++++++++

      do i=1,nint
         read(1,*)ia,ib,ic,id,jj,iso,vvv
         if(jj.gt.j2max)then
            write(6,*)'You shoud not be here: ',jj,'>',j2max
            write(6,*)'Check dimensions or interaction file.'
            stop
         end if
         vvv=vvv*scale
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C =========== put into "standard order":   a < = b, c < = d,  a <= c ====
C====================================================
         iaabb=int(orb_qn(ia,3)+orb_qn(ib,3)-dfloat(jj+iso))
         iccdd=int(orb_qn(ic,3)+orb_qn(id,3)-dfloat(jj+iso))
         phab=(-1.)**iaabb
         phcd=(-1.)**iccdd

      if(ia.gt.ib)then 
              vvv=vvv*phab
              i1 = ia
              ia=ib
              ib=i1
      endif
      if(ic.gt.id)then
	vvv=vvv*phcd
	i1 =ic
	ic=id
	id = i1
      endif
      if(ia.gt.ic)then
	i1=ia
	i2 =ib
	ia=ic
	ib =id
	ic = i1
	id =i2
      endif
      if(ia == ic .and. ib < id)then
        i2 = ib
        ib = id
        id = i2
      endif

C      print*,ia,ib,ic,id
      vjt(jj,iso,ia,ib,ic,id)=vjt(jj,iso,ia,ib,ic,id)+vvv
      enddo		! end loop over nint
      close(1)
      enddo		! end loop over n_ham
22     continue      

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C                 option to read in an external field   NEW
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      allocate(vextp(nsps(1),nsps(1)),vextn(nsps(1),nsps(1)))
      vextp = 0.
      vextn = 0.

      print*,' Do you want to read in an external field (y/n)?'
      read(5,'(a)')ychar
      if(ychar == 'y' .or. ychar =='Y')then


C-------------- OPEN A FILE ----------------------------------
      success = .false.

      do while(.not.success)
          print*,' Enter file (.exf)'
          read(5,'(a)')input
          ilast = index(input,' ')-1
C........... ATTEMPT TO OPEN .exf FILE.................
          open(unit=3,file=input(1:ilast)//'.exf',status='old', 
     &     err=101)
          success = .true.
          cycle
101       continue

          print*,input(1:ilast),'.exf does not exist '
      enddo

C--------------- CHECK S.P. Space is the same ---------

      read(3,*)title
      write(6,*)title
      read(3,*)nspstest
      if(nspstest /= nsps(1))then
        print*,' Mismatch in size of space ',nspstest,nsps(1)
        stop
      endif
      do i = 1,nsps(1)
         read(3,*)jtest,mtest
         if(jtest /= spsqn(1,i,4))then
            print*,' mismatch in js ',jtest,spsqn(1,i,4)
            stop
         endif
         if(mtest /= spsqn(1,i,5))then
            print*,' mismatch in ms ',mtest,spsqn(1,i,5)
            stop
         endif
      enddo
C---------------ALLOCATE AND READ IN --------------


      do i = 1,nsps(1)
         read(3,*)(vextp(i,j),j=1,nsps(1))
      enddo
      do i = 1,nsps(1)
         read(3,*)(vextn(i,j),j=1,nsps(1))
      enddo

C---------------- ASK FOR STRENGTH --------------------

      print*,' Enter strength of field '
      read*,vstrength
      vextp(:,:) = vstrength*vextp(:,:)
      vextn(:,:) = vstrength*vextn(:,:)

      close(unit=3)
      endif

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C==============    Complete the whole matrix via phase arguments
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      do ia = 1,nspe
         do ib =ia,nspe
            do ic = ia,nspe
               do id = ic,nspe
                  do jj = 0,j2max
                    do iso =0,1
         vvv= vjt(jj,iso,ia,ib,ic,id)
         if(vvv.ne.0.0)then

 
         vjt(jj,iso,ic,id,ia,ib)=vvv
 
         iaabb=int(orb_qn(ia,3)+orb_qn(ib,3)-dfloat(jj+iso))
         iccdd=int(orb_qn(ic,3)+orb_qn(id,3)-dfloat(jj+iso))
         phab=(-1.)**iaabb
         phcd=(-1.)**iccdd
         if(ia.ne.ib)then
            vjt(jj,iso,ib,ia,ic,id)=phab*vvv
            
            vjt(jj,iso,ic,id,ib,ia)=phab*vvv
         end if
         if(ic.ne.id)then
            vjt(jj,iso,ia,ib,id,ic)=phcd*vvv
            vjt(jj,iso,id,ic,ia,ib)=phcd*vvv
         end if
         if(ia.ne.ib.and.ic.ne.id)then
           vjt(jj,iso,ib,ia,id,ic)=phcd*phab*vvv
           vjt(jj,iso,id,ic,ib,ia)=phcd*phab*vvv
         end if
      endif
      end do
      enddo
      enddo
      enddo
      enddo
      enddo
      


C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C                 convert from isospin to pn formalism 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    

C-------------NOW I HAVE TO CREATE ALL THE POSSIBLE p+n ORBITS----

C----- first, protons
      j = 0
      allocate(jt(norb(1)+norb(2)),lt(norb(1)+norb(2)),
     &         nt(norb(1)+norb(2)),tt(norb(1)+norb(2)))
      allocate(map(2*nspe))
      do i = 1,norb(1)
      	    j = j+1
      	    map(j) = i
      	    jt(j)=int(orb_qn(i,3)*2)
      	    lt(j)=orb_qn(i,2)
            nt(j)=orb_qn(i,1)
      	    tt(j)=1
C	    if(i.gt.norb(1)/2)tt(j)=-1
      enddo
C--------then neutrons
      do i = 1,norb(1)
      	    j = j+1
      	    map(j) = i
      	    jt(j)=int(2*orb_qn(i,3))
            nt(j)=orb_qn(i,1)
      	    lt(j)=orb_qn(i,2)
      	    tt(j) = -1
      enddo

      call TBMElist_count(norb(1)+norb(2),jt,lt,tt,nmepn,nt)
      allocate(mout(nmepn,4),jtot(nmepn))
      call TBMElist(norb(1)+norb(2),jt,lt,tt,nmepn,mout,jtot,nt)

C      stop

      ej=0.0
      do i = 1,nmepn
         ttot = tt(mout(i,1))+tt(mout(i,2))
         ia=mout(i,1)
         ib=mout(i,2)
         ic=mout(i,3)
         id=mout(i,4)
         jj=jtot(i)
      	 if(abs(ttot).gt.0)then 
      	    iso =1
      	    xme = vjt(jj,1,map(ia),map(ib),map(ic),map(id))
         else
            iso = 0
      	    xme = (vjt(jj,0,map(ia),map(ib),map(ic),map(id))
     &  +vjt(jj,1,map(ia),map(ib),map(ic),map(id)))       

         endif
C         print*,ia,ib,ic,id
         ej(jj,iso,ia,ib,ic,id)=xme
         ej(jj,iso,ic,id,ia,ib)=xme
         iaabb=int(orb_qn(ia,3)+orb_qn(ib,3)-dfloat(jj+iso))
         iccdd=int(orb_qn(ic,3)+orb_qn(id,3)-dfloat(jj+iso))
         phab=(-1.)**iaabb
         phcd=(-1.)**iccdd
         if(ia.ne.ib)then
            ej(jj,iso,ib,ia,ic,id)=phab*xme
            ej(jj,iso,ic,id,ib,ia)=phab*xme
         end if
         if(ic.ne.id)then
            ej(jj,iso,ia,ib,id,ic)=phcd*xme
            ej(jj,iso,id,ic,ia,ib)=phcd*xme
         end if
         if(ia.ne.ib.and.ic.ne.id)then
            ej(jj,iso,ib,ia,id,ic)=phcd*phab*xme
            ej(jj,iso,id,ic,ib,ia)=phcd*phab*xme
         end if         
      enddo


  
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                 First construct the pp part of the Hamiltonian
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c===============  fill array for single-particle energies
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C      if(allocated(e_spepp))deallocate(e_spepp)
C      if(allocated(e_spenn))deallocate(e_spenn)

c      allocate(e_spepp_p(nsps(1),nsps(1)))
      if(ham)then
         allocate(e_spepp(nsps(1),nsps(1)))
c         e_spepp_p=>e_spepp
         e_spepp=0.0
         do i=1,nsps(1)
            ia=spsqn(1,i,1)
            e_spepp(i,i)=spe(ia)
            do j = 1,nsps(1)
              e_spepp(i,j) = e_spepp(i,j)+vextp(i,j)
            enddo
         end do
      else
         allocate(e_spepp_O(nsps(1),nsps(1)))
         e_spepp_O=0.0
         do i=1,nsps(1)
            ia=spsqn(1,i,1)
            e_spepp_O(i,i)=spe(ia)
            do j = 1,nsps(1)
              e_spepp_O(i,j) = e_spepp_O(i,j)+vextp(i,j)
            enddo
         end do
c         e_spepp_p=>e_spepp_O
      endif


c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c===============  Now make the Hamiltonian vector
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      countonly=.true.
 1999 continue
      nmatpp_p=0
      do i=1,nsps(1)
         ta=spsqn(1,i,6)
         do j=i+1,nsps(1)
            tb=spsqn(1,j,6)
            do k=1,nsps(1)
               tc=spsqn(1,k,6)
               do l=k+1,nsps(1)
                  ma=spsqn(1,i,5)               ! 2*jz-value for state #1
                  mb=spsqn(1,j,5)               ! 2*jz-value for state #2
                  mc=spsqn(1,k,5)               ! 2*jz-value for state #3
                  md=spsqn(1,l,5)               ! 2*jz-value for state #4
                  mij=ma+mb                     ! 2*jz for creat. ops
                  mkl=mc+md                     ! 2*jz for annih. ops
                  td=spsqn(1,l,6)
		  tij=ta+tb
		  tkl=tc+td
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c===============  Check if z-proj. of J is the same. if so, compute
c===============  a matrix element
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                  if((mij.eq.mkl).and.(tij.eq.tkl))then
                     ia=spsqn(1,i,1)            ! orbit label for state #1
                     ib=spsqn(1,j,1)            ! orbit label for state #2
                     ic=spsqn(1,k,1)            ! orbit label for state #3
                     id=spsqn(1,l,1)            ! orbit label for state #4
                     ja=spsqn(1,i,4)            ! 2*j-value for state #1
                     jb=spsqn(1,j,4)            ! 2*j-value for state #2
                     jabmax=(ja+jb)             ! 2*Max J for states #1&#2
                     jabmin=iabs(ja-jb)         ! 2*Min J for states #1&#2
                     jc=spsqn(1,k,4)            ! 2*j-value for state #3
                     jd=spsqn(1,l,4)            ! 2*j-value for state #4
                     jcdmax=(jc+jd)             ! 2*Max J for states #3&#4
                     jcdmin=iabs(jc-jd)         ! 2*Min J for states #3&#4
                     jmin=max(jabmin,jcdmin)
                     jmax=min(jabmax,jcdmax)
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c=============== Sqrt. factor if(orbit#1=orbit#2) or if(orbit#3=orbit#4)
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                     fact=one
                     if(ia.eq.ib)fact=fact*sqrt(two)
                     if(ic.eq.id)fact=fact*sqrt(two)
		     if(ja.eq.jb)then
		      if(jc.eq.jd)then
		       if(ta*tb.lt.0)fact=two*fact
		      endif
		     endif
                     sum=zero
                     do jr=jmin,jmax,2       !  sum over J-values
                        jj=jr/2

                        sum=sum+cleb(ja,ma,jb,mb,jr,mij)*
     &                          cleb(jc,mc,jd,md,jr,mkl)*
     &                  (ej(jj,1,ia,ib,ic,id)*cleb(1,ta,1,tb,2,tij)
     &                                       *cleb(1,tc,1,td,2,tkl))
                     end do
                     sum=-fact*sum
                     if(abs(sum).gt.small)then
                        nmatpp_p=nmatpp_p+1
                        if(.not.countonly)then
                          if(ham)then
                           hvecpp(nmatpp_p)=sum
                           map_combpp(nmatpp_p)=ipack(i,j,k,l)
                          else
                           hvecpp_O(nmatpp_p)=sum
                           map_combpp_O(nmatpp_p)=ipack(i,j,k,l)
                          endif
                        endif
                    end if
                  end if
               end do
            end do
         end do
      end do
      if(countonly)then
         countonly=.false.
         print*,' # of pp TBMEs is:',nmatpp_p
c         allocate(hvecpp_p(nmatpp_p),map_combpp_p(nmatpp_p))
         if(ham)then
            nmatpp=nmatpp_p
            allocate(hvecpp(nmatpp),map_combpp(nmatpp))
c            hvecpp_p=>hvecpp
c            map_combpp_p=>map_combpp
         else
            nmatpp_O=nmatpp_p
            allocate(hvecpp_O(nmatpp_O),map_combpp_O(nmatpp_O))
c            hvecpp_p=>hvecpp_O
c            map_combpp_p=>map_combpp_O
         endif
         goto 1999
      endif


c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                 Now construct the nn part of the Hamiltonian
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c===============  fill array for single-particle energies
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c      allocate(e_spenn_p(nsps(1),nsps(1)))
      if(ham)then
         allocate(e_spenn(nsps(2),nsps(2))) 
         e_spenn=0.0
         do i=1,nsps(2)
            ia=spsqn(2,i,1)
            e_spenn(i,i)=spe(ia)
            do j = 1,nsps(2)
              e_spenn(i,j) = e_spenn(i,j)+vextn(i,j)
            enddo
         end do
c          e_spenn_p=>e_spenn
      else
         allocate(e_spenn_O(nsps(2),nsps(2)))
         e_spenn_O=0.0
         do i=1,nsps(2)
            ia=spsqn(2,i,1)
            e_spenn_O(i,i)=spe(ia)
            do j = 1,nsps(2)
              e_spenn_O(i,j) = e_spenn_O(i,j)+vextn(i,j)
            enddo
c         print200,e_spenn_p(i,i)
         end do
C          e_spenn_p=>e_spenn_O
      endif
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c===============  Now make the Hamiltonian vector
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      countonly=.true.
 1998 continue
      nmatnn_p=0
      do i=1,nsps(2)
         do j=i+1,nsps(2)
            do k=1,nsps(2)
               do l=k+1,nsps(2)
                  ma=spsqn(2,i,5)               ! 2*jz-value for state #1
                  mb=spsqn(2,j,5)               ! 2*jz-value for state #2
                  mc=spsqn(2,k,5)               ! 2*jz-value for state #3
                  md=spsqn(2,l,5)               ! 2*jz-value for state #4
                  mij=ma+mb                     ! 2*jz for creat. ops
                  mkl=mc+md                     ! 2*jz for annih. ops
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c===============  Check if z-proj. of J is the same. if so, compute
c===============  a matrix element
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                  if((mij.eq.mkl))then
                     ia=spsqn(2,i,1)            ! orbit label for state #1
                     ib=spsqn(2,j,1)            ! orbit label for state #2
                     ic=spsqn(2,k,1)            ! orbit label for state #3
                     id=spsqn(2,l,1)            ! orbit label for state #4
                     ja=spsqn(2,i,4)            ! 2*j-value for state #1
                     jb=spsqn(2,j,4)            ! 2*j-value for state #2
                     jabmax=(ja+jb)             ! 2*Max J for states #1&#2
                     jabmin=iabs(ja-jb)         ! 2*Min J for states #1&#2
                     jc=spsqn(2,k,4)            ! 2*j-value for state #3
                     jd=spsqn(2,l,4)            ! 2*j-value for state #4
                     jcdmax=(jc+jd)             ! 2*Max J for states #3&#4
                     jcdmin=iabs(jc-jd)         ! 2*Min J for states #3&#4
                     jmin=max(jabmin,jcdmin)
                     jmax=min(jabmax,jcdmax)
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c=============== Sqrt. factor if(orbit#1=orbit#2) or if(orbit#3=orbit#4)
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                     fact=one
                     if(ia.eq.ib)fact=fact*sqrt(two)
                     if(ic.eq.id)fact=fact*sqrt(two)
                     sum=zero
                     do jr=jmin,jmax,2       !  sum over J-values
                        jj=jr/2
                        sum=sum+cleb(ja,ma,jb,mb,jr,mij)*
     &                          cleb(jc,mc,jd,md,jr,mkl)*
     &                          ej(jj,1,ia,ib,ic,id)
                     end do
                     sum=-fact*sum
                     if(abs(sum).gt.small)then
                        nmatnn_p=nmatnn_p+1
                        if(.not.countonly)then
                         if(ham)then
                           hvecnn(nmatnn_p)=sum
                           map_combnn(nmatnn_p)=ipack(i,j,k,l)
                         else
                           hvecnn_O(nmatnn_p)=sum
                           map_combnn_O(nmatnn_p)=ipack(i,j,k,l)
                         endif
                        endif
                    end if
                  end if
               end do
            end do
         end do
      end do
      if(countonly)then
         print*,' # of nn TBMEs is:',nmatnn_p
         countonly=.false.
c         allocate(hvecnn_p(nmatnn_p),map_combnn_p(nmatnn_p))
         if(ham)then
            nmatnn=nmatnn_p
            allocate(hvecnn(nmatnn),map_combnn(nmatnn))
c            hvecnn_p=>hvecnn
c            map_combnn_p=>map_combnn
         else
            nmatnn_O=nmatnn_p
            allocate(hvecnn_O(nmatnn_O),map_combnn_O(nmatnn_O))
c            hvecnn_p=>hvecnn_O
c            map_combnn_p=>map_combnn_O
         endif
         goto 1998
      endif
     
      
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                 Now, construct the pn part of the Hamiltonian
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c===============  No single-particle part
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c===============  Now make the Hamiltonian vector
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      countonly=.true.
 1997 continue
      nmatpn_p=0
      do j=1,nsps(2)
         do l=1,nsps(2)
            do i=1,nsps(1)
               do k=1,nsps(1)
                  mb=spsqn(2,j,5)               ! 2*jz-value for state #1
                  md=spsqn(2,l,5)               ! 2*jz-value for state #2
                  ma=spsqn(1,i,5)               ! 2*jz-value for state #3
                  mc=spsqn(1,k,5)               ! 2*jz-value for state #4
                  mij=ma+mb                     ! 2*jz for creat. ops
                  mkl=mc+md                     ! 2*jz for annih. ops
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c===============  Check if z-proj. of J is the same. if so, compute
c===============  a matrix element
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                  if((mij.eq.mkl))then
                     ia=spsqn(1,i,1)            ! orbit label for state #1
                     ib=spsqn(2,j,1)            ! orbit label for state #2
                     ic=spsqn(1,k,1)            ! orbit label for state #3
                     id=spsqn(2,l,1)            ! orbit label for state #4
                     ja=spsqn(1,i,4)            ! 2*j-value for state #1
                     jb=spsqn(2,j,4)            ! 2*j-value for state #2
                     jabmax=(ja+jb)             ! 2*Max J for states #1&#2
                     jabmin=iabs(ja-jb)         ! 2*Min J for states #1&#2
                     jc=spsqn(1,k,4)            ! 2*j-value for state #3
                     jd=spsqn(2,l,4)            ! 2*j-value for state #4
                     jcdmax=(jc+jd)             ! 2*Max J for states #3&#4
                     jcdmin=iabs(jc-jd)         ! 2*Min J for states #3&#4
                     jmin=max(jabmin,jcdmin)
                     jmax=min(jabmax,jcdmax)
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c=============== Sqrt. factor if(orbit#1=orbit#2) or if(orbit#3=orbit#4)
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C corrected Feb 2008 by CWJ
C
                     ichkab=abs(spsqn(1,i,2)-spsqn(2,j,2))+
     &                      abs(spsqn(1,i,3)-spsqn(2,j,3))+
     &                      abs(spsqn(1,i,4)-spsqn(2,j,4))
                     ichkcd=abs(spsqn(1,k,2)-spsqn(2,l,2))+
     &                      abs(spsqn(1,k,3)-spsqn(2,l,3))+
     &                      abs(spsqn(1,k,4)-spsqn(2,l,4))
                     fact=one
                     if(ichkab.eq.0)fact=fact*sqrt(two)
                     if(ichkcd.eq.0)fact=fact*sqrt(two)
                     sum=zero
                     do jr=jmin,jmax,2       !  sum over J-values
                        jj=jr/2
                        sum=sum+cleb(ja,ma,jb,mb,jr,mij)*
     &                          cleb(jc,mc,jd,md,jr,mkl)*
     &                          (ej(jj,1,ia,ib,ic,id)+
     &                           ej(jj,0,ia,ib,ic,id))
                     end do
                     sum=0.5*fact*sum
                     if(abs(sum).gt.small)then
                        nmatpn_p=nmatpn_p+1
                        if(.not.countonly)then
                          if(ham)then
                           hvecpn(nmatpn_p)=sum
                           map_combpn(nmatpn_p)=ipack(j,l,i,k)
                          else
                           hvecpn_O(nmatpn_p)=sum
                           map_combpn_O(nmatpn_p)=ipack(j,l,i,k)
                          endif
                        endif
                    end if
                  end if
               end do
            end do
         end do
      end do
      if(countonly)then
         countonly=.false.
         print*,' # of pn TBMEs is:',nmatpn_p
c         allocate(hvecpn_p(nmatpn_p),map_combpn_p(nmatpn_p))
         if(ham)then
            nmatpn=nmatpn_p
            allocate(hvecpn(nmatpn),map_combpn(nmatpn))
c            hvecpn_p=>hvecpn
c            map_combpn_p=>map_combpn
         else
            nmatpn_O=nmatpn_p
            allocate(hvecpn_O(nmatpn_O),map_combpn_O(nmatpn_O))
c            hvecpn_p=>hvecpn_O
c            map_combpn_p=>map_combpn_O
         endif
         goto 1997
      endif
c      print200,hvecpn_p(1:10)
 200  format(10F12.4)
      deallocate(vjt,spe,ej,temp,mout,map,jtot,jt,lt,nt,tt)
C      nullify(e_spepp_p,e_spenn_p,hvecpp_p,hvecnn_p,hvecpn_p)
C      nullify(map_combpp_p,map_combnn_p,map_combpn_p)
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

