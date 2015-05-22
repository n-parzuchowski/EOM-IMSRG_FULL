CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  original routines written by W.E. Ormand
C  modified by CWJ   March/April 1998  LSU
C
C  revisions for isospin  CWJ April 2001 LSU
C
C  revision for f90  IS May 2004 UA


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine make_hvec_iso(nmatpp,nmatnn,nmatpn,
     &  map_combpp,map_combnn,map_combpn,
     &  hvecpp,hvecnn,hvecpn,e_spepp,e_spenn)
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
c      use HAMILTONIAN
c      use observables
      implicit none

C OUTPUT

c      logical     ham   ! = true, computes the Hamiltonan
                        ! = false computes observables
C-----------HAMILTONIAN------------------------------------------
      integer    :: nmatpp		! # of pp matrix elements 
      integer    :: nmatnn		! # of nn matrix elements
      integer    :: nmatpn		! # of np matrix elements

C-------- list of packed 4-indices with nonzero TBME     
      integer,pointer,dimension(:)  :: map_combpp,map_combnn
      integer,pointer,dimension(:)  :: map_combpn

C---------- nonzero TBMEs associated with above list   
      real,pointer,dimension(:)     :: hvecpp,hvecnn,hvecpn
      
C---------one-body matrix elements      
      real,pointer,dimension(:,:)  :: e_spepp,e_spenn

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
      allocate(vjt(0:j2max,0:1,1:nspe,1:nspe,1:nspe,1:nspe))
      allocate(ej(0:j2max,0:1,1:2*nspe,1:2*nspe,1:2*nspe,1:2*nspe))
      allocate(spe(norb(1)+norb(2)))

C      deassociate(hvecpp_p,hvecnn_p,hvecpn_p)
C      deassociate(map_combpp_p,map_combnn_p,map_combpn_p)

C....Make sure the matrices to store Hamiltonian/Observable are not already allocated.....
c      if(allocated(e_spepp))deallocate(e_spepp)
c      if(allocated(e_spenn))deallocate(e_spenn)
c      if(allocated(hvecpp))deallocate(hvecpp)
c      if(allocated(hvecnn))deallocate(hvecnn)
c      if(allocated(hvecpn))deallocate(hvecpn)
c      if(allocated(map_combpp))deallocate(map_combpp)
c      if(allocated(map_combnn))deallocate(map_combnn)
c      if(allocated(map_combpn))deallocate(map_combpn)

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
      vjt(jj,iso,ia,ib,ic,id)=vjt(jj,iso,ia,ib,ic,id)+vvv
      enddo		! end loop over nint
      close(1)
      enddo		! end loop over n_ham
22     continue      

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

      allocate(e_spepp(nsps(1),nsps(1)))
      do i=1,nsps(1)
         ia=spsqn(1,i,1)
         e_spepp(i,i)=spe(ia)
      end do
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c===============  Now make the Hamiltonian vector
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      countonly=.true.
 1999 continue
      nmatpp=0
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
                        nmatpp=nmatpp+1
                        if(.not.countonly)then
                           hvecpp(nmatpp)=sum
                           map_combpp(nmatpp)=ipack(i,j,k,l)
                        endif
                    end if
                  end if
               end do
            end do
         end do
      end do
      if(countonly)then
         countonly=.false.
         print*,' # of pp TBMEs is:',nmatpp
         allocate(hvecpp(nmatpp),map_combpp(nmatpp))
         goto 1999
      endif


c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                 Now construct the nn part of the Hamiltonian
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c===============  fill array for single-particle energies
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      allocate(e_spenn(nsps(2),nsps(2)))
      do i=1,nsps(2)
         ia=spsqn(2,i,1)
         e_spenn(i,i)=spe(ia)
c         print200,e_spenn_p(i,i)
      end do
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c===============  Now make the Hamiltonian vector
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      countonly=.true.
 1998 continue
      nmatnn=0
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
                        nmatnn=nmatnn+1
                        if(.not.countonly)then
                           hvecnn(nmatnn)=sum
                           map_combnn(nmatnn)=ipack(i,j,k,l)
                        endif
                    end if
                  end if
               end do
            end do
         end do
      end do
      if(countonly)then
         print*,' # of nn TBMEs is:',nmatnn
         countonly=.false.
         allocate(hvecnn(nmatnn),map_combnn(nmatnn))
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
      nmatpn=0
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
                     ichkab=spsqn(1,i,2)-spsqn(2,j,2)+
     &                      spsqn(1,i,3)-spsqn(2,j,3)+
     &                      spsqn(1,i,4)-spsqn(2,j,4)
                     ichkcd=spsqn(1,k,2)-spsqn(2,l,2)+
     &                      spsqn(1,k,3)-spsqn(2,l,3)+
     &                      spsqn(1,k,4)-spsqn(2,l,4)
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
                        nmatpn=nmatpn+1
                        if(.not.countonly)then
                           hvecpn(nmatpn)=sum
                           map_combpn(nmatpn)=ipack(j,l,i,k)
                        endif
                    end if
                  end if
               end do
            end do
         end do
      end do
      if(countonly)then
         countonly=.false.
         print*,' # of pn TBMEs is:',nmatpn
         allocate(hvecpn(nmatpn),map_combpn(nmatpn))
         goto 1997
      endif
c      print200,hvecpn_p(1:10)
 200  format(10F12.4)
      deallocate(vjt,spe,ej,temp,mout,map,jtot,jt,lt,nt,tt)
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine TBMElist(nsp,jt,lt,tt,nme,m,jtot,nt)

C  creates list of all possible TBME matrix elements
C
C  INPUT:
C       NSP:   # of s.p. states
C	JT(i): 2 x j of i'th s.p. state
C       LT(i)  = orb ang moment L of s.p. state (necessary for parity)
C       TT(i): 2 x t of i'th s.p. states
C	NT(i): radial q# N of i'th s.p. state (principle q#)
C
C  OUTPUT:
C	NME:  # of TBMEs
C	M(i,x): for i'th TBME, for x = 1,4 what is s.p. state
C     		< 1 2 | V | 3 4 >
C       JTOT(i): total J for i'th TBME
C


      implicit none
C      include 'gcb_dim.inc'

      integer,intent(IN)    :: nsp
C---- SINGLE-PARTICLE LIST------

      integer jt(nsp),tt(nsp),nt(nsp)
      integer lt(nsp),ltot
      integer,intent(INOUT)            :: nme		! # of matrix elements
      integer,intent(OUT)              :: m(nme,4)	! set of 4 states
      integer,intent(OUT)              :: jtot(nme)     ! spin
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      integer i,l			! count over orbits

      integer,allocatable              :: jpair(:),tpair(:),pair(:,:)
      integer npairs,jmin,jmax,nme0



      integer k
C--------------------- sort through and create all possible pairs

C... first just count the possible pairs

      nme0=nme
      npairs = 0
      do i = 1,nsp
	do l = i,nsp
	    jmin = (abs(jt(i)-jt(l)))/2
	    jmax = (jt(i)+jt(l))/2
            do k = jmin,jmax
		if(i.eq.l )then
		  if( mod(k+2,2).eq.0)then

		   npairs = npairs +1
		   endif

		else
		   npairs = npairs + 1
		endif
	    enddo
	enddo
      enddo


      allocate(pair(npairs,2),jpair(npairs),tpair(npairs))
      npairs = 0
      do i = 1,nsp
	do l = i,nsp
	    jmin = (abs(jt(i)-jt(l)))/2
	    jmax = (jt(i)+jt(l))/2
            do k = jmin,jmax
		if(i.eq.l )then
		  if( mod(k+2,2).eq.0)then

		   npairs = npairs +1
		   pair(npairs,1)=i
		   pair(npairs,2)=l
		   jpair(npairs)=k
		   tpair(npairs)=tt(i)+tt(l)
		   endif

		else
		   npairs = npairs + 1
		   pair(npairs,1)=i
		   pair(npairs,2)=l

		   jpair(npairs) = k
		   tpair(npairs)=tt(i)+tt(l)

		endif
	    enddo
	enddo
      enddo

C--------------------- now list all possible matrix elements

      nme = 0
      do i = 1,npairs
	do k = i,npairs
	   if(jpair(i).eq.jpair(k).and.tpair(i).eq.tpair(k))then

C--- check parity
             ltot = lt(pair(i,1))+lt(pair(i,2))+
     &             lt(pair(k,1))+lt(pair(k,2))
             if(ltot.eq.2*(ltot/2))then
	  	 nme= nme+1
                 if(nme>nme0)stop 'Something wrong...'
		 jtot(nme)=jpair(i)
  		 m(nme,1)=pair(i,1)
		 m(nme,2)=pair(i,2)
		 m(nme,3)=pair(k,1)
		 m(nme,4)=pair(k,2)
             endif
	    endif
	enddo
      enddo

      deallocate(pair,jpair,tpair)
C--------------------- now sort into proper order

      if(nme0/=nme)stop 'Problem with the # of TBME'
C      write(6,*)'There are ', nme,' two-body matrix elements.'

      return

      end

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine hmult(nsps,nmatpp,nmatnn,nmatpn,
     &  map_combpp,map_combnn,map_combpn,
     &  hvecpp,hvecnn,hvecpn,e_spepp,e_spenn,
     &  overp,overn,rhopij,rhonij,over,hme)
      
C
C IS: corrected to account for non-diagonal elements in the
C     one-body term
C  routine to compute hamiltonian matrix elements
C  < psil | H | psir >
C
C  INPUT:
C	nsps:	# of m-states
C	nmatpp,nmatnn,nmatnp: # of 2-body matrix elements
C	map_combpp,map_combnn,map_combpn: mapping matrix elements
C	hvecpp(i) = ith TBME associated with indices unpacked from map_combpp(i)
C	e_spepp(i,j) = one-body matrix element between i,j m-states
C	overp,overn = <psil | psir> for protons, neutrons
C	rhopij,rhonij = one-body density matrices for protons, neutrons
C
C  OUTPUT:
C	over = total overlap
C	hme = total matrix element
C   
      implicit none


C-----------HAMILTONIAN------------------------------------------
      integer nsps(2)		! # of single-particle m-states

      integer nmatpp		! # of pp matrix elements 
      integer nmatnn		! # of nn matrix elements
      integer nmatpn		! # of np matrix elements

C-------- list of packed 4-indices with nonzero TBME     
      integer map_combpp(nmatpp)
      integer map_combnn(nmatnn)
      integer map_combpn(nmatpn)
      
C---------- nonzero TBMEs associated with above list   
      real hvecpp(nmatpp)
      real hvecnn(nmatnn)
      real hvecpn(nmatpn)
      
C---------one-body matrix elements
      real e_spepp(nsps(1),nsps(1))
      real e_spenn(nsps(2),nsps(2))
      
     
      integer n_one
      real,intent(OUT)      :: hme
      
C------One-body densities----------------------------------------
      real*8     :: rhopij(nsps(1),nsps(1)),rhonij(nsps(2),nsps(2))
      real*8     :: overp,overn
      real       :: over

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c============  Often used integers
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer i,j,k,n,l,m
      integer ilast
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C============  Declare the type for some useful constants
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real zero,one,two,pi,temp,cb,small
      real*8 den_test
      real*8 obp,tbpp,obn,tbnn,tbpn
      
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c============  Initialize the useful constants
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      zero=0.0e0
      one=1.0e0
      two=2.0e0
      small=1.0e-7
      pi=two*asin(one)      

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c============   Compute the matrix element of the Hamiltonian
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c============   First the proton-proton part
c      obp=zero       ! Initially was assumed diagonal interaction
c      do n=1,nsps(1) 
c        obp=obp+e_spepp(n,n)*rhopij(n,n)*overp
c      end do
      obp=zero
      do n=1,nsps(1)
       do j=1,nsps(1)
        obp=obp+e_spepp(n,j)*rhopij(j,n)*overp
       enddo
      end do
      tbpp=zero
      do n=1,nmatpp
         call unpack0(map_combpp(n),i,j,k,l)
         den_test=rhopij(i,l)*rhopij(j,k)-rhopij(i,k)*rhopij(j,l)
         tbpp=tbpp+den_test*hvecpp(n)*overp
      end do
c============   For proton part, multiply by neutron overlap
      obp=obp*overn
      tbpp=tbpp*overn
c============   Next the neutron-neutron part
C      obn=zero
c      do n=1,nsps(2)
c         obn=obn+e_spenn(n,n)*rhonij(n,n)*overn
c      end do
      obn=zero
      do n=1,nsps(2)
       do j=1,nsps(2)
         obn=obn+e_spenn(n,j)*rhonij(j,n)*overn
       enddo
      end do
      tbnn=zero
      do n=1,nmatnn
         call unpack0(map_combnn(n),i,j,k,l)
         den_test=rhonij(i,l)*rhonij(j,k)-rhonij(i,k)*rhonij(j,l)
         tbnn=tbnn+den_test*hvecnn(n)*overn
      end do
c============   For neutron part, multiply by proton overlap
      obn=obn*overp
      tbnn=tbnn*overp
c============   Finally, the proton neutron part
      tbpn=zero
      do n=1,nmatpn
          call unpack0(map_combpn(n),j,l,i,k)
          den_test=rhonij(j,l)*rhopij(i,k)
          tbpn=tbpn+den_test*hvecpn(n)*overn*overp
      end do
            
      hme=real(obp+tbpp+obn+tbnn+tbpn)
      over=real(overp*overn)      
      
      return
      end subroutine hmult
      



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      integer function ipack(i1,i2,i3,i4)
c==============   Function to pack four small integers into one
      
      implicit none
      integer i1,i2,i3,i4
      ipack=ior(I4,ishft(ior(I3,ishft(ior(I2,ishft(I1,8)),8)),8))
      return
      end function ipack
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine unpack0(packed,i,j,k,l)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c      Subroutine to unpack the four integers i,j,k,l stored in ipack using
c      ipack(I1,I2,I3,I4)=
c     & ior(I4,ishft(ior(I3,ishft(ior(I2,ishft(I1,8)),8)),8))
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit none
      integer i,j,k,l
      integer packed,i1,itemp
      i1=255               ! 1111 1111 0000 0000 0000 0000 0000 0000
      itemp=packed
      l=iand(itemp,i1)
      itemp=ishft(itemp,-8)
      k=iand(itemp,i1)
      itemp=ishft(itemp,-8)
      j=iand(itemp,i1)
      itemp=ishft(itemp,-8)
      i=iand(itemp,i1)
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



      subroutine TBMElist_count(nsp,jt,lt,tt,nme,nt)

C  creates list of all possible TBME matrix elements
C
C  INPUT:
C       NSP:   # of s.p. states
C	JT(i): 2 x j of i'th s.p. state
C       LT(i)  = orb ang moment L of s.p. state (necessary for parity)
C       TT(i): 2 x t of i'th s.p. states
C	NT(i): radial q# N of i'th s.p. state (principle q#)
C
C  OUTPUT:
C	NME:  # of TBMEs
C	M(i,x): for i'th TBME, for x = 1,4 what is s.p. state
C     		< 1 2 | V | 3 4 >
C       JTOT(i): total J for i'th TBME
C


      implicit none
C      include 'gcb_dim.inc'

      integer,intent(IN)    :: nsp
C---- SINGLE-PARTICLE LIST------

      integer jt(nsp),tt(nsp),nt(nsp)
      integer lt(nsp),ltot
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      integer i,l			! count over orbits

      integer,allocatable              :: jtot(:)		! spin
      integer,allocatable              :: jpair(:),tpair(:),pair(:,:)
      integer npairs,jmin,jmax

      integer,intent(OUT)              :: nme		! # of matrix elements

      integer k
      logical countonly

C--------------------- sort through and create all possible pairs

C... first just count the possible pairs
      npairs = 0
      do i = 1,nsp
	do l = i,nsp
	    jmin = (abs(jt(i)-jt(l)))/2
	    jmax = (jt(i)+jt(l))/2
            do k = jmin,jmax
		if(i.eq.l )then
		  if( mod(k+2,2).eq.0)then

		   npairs = npairs +1
		   endif

		else
		   npairs = npairs + 1
		endif
	    enddo
	enddo
      enddo


      allocate(pair(npairs,2),jpair(npairs),tpair(npairs))
      npairs = 0
      do i = 1,nsp
	do l = i,nsp
	    jmin = (abs(jt(i)-jt(l)))/2
	    jmax = (jt(i)+jt(l))/2
            do k = jmin,jmax
		if(i.eq.l )then
		  if( mod(k+2,2).eq.0)then

		   npairs = npairs +1
		   pair(npairs,1)=i
		   pair(npairs,2)=l
		   jpair(npairs)=k
		   tpair(npairs)=tt(i)+tt(l)
		   endif

		else
		   npairs = npairs + 1
		   pair(npairs,1)=i
		   pair(npairs,2)=l

		   jpair(npairs) = k
		   tpair(npairs)=tt(i)+tt(l)

		endif
	    enddo
	enddo
      enddo

      write(6,*)'There are ',npairs,' two-particle pairs possible '

C--------------------- now list all possible matrix elements
      countonly=.true.
C      countonly=.false.
 2004 continue
      nme = 0
      do i = 1,npairs
	do k = i,npairs
	   if(jpair(i).eq.jpair(k).and.tpair(i).eq.tpair(k))then

C--- check parity
             ltot = lt(pair(i,1))+lt(pair(i,2))+
     &             lt(pair(k,1))+lt(pair(k,2))
             if(ltot.eq.2*(ltot/2))then
		nme= nme+1
                endif
             endif
	enddo
      enddo
      deallocate(pair,jpair,tpair)
C--------------------- now sort into proper order

      write(6,*)'There are ', nme,' two-body matrix elements.'

      return

      end

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++







