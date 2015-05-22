
      subroutine HFdiag(energy)
C
C  INPUT: np/nn number of protons/neutrons
C
C  OUTPUT:
C

C      include 'gcb_dim.inc'
      use spspace
      use HAMILTONIAN
      use wfn
      use nuclide
      use HFdiagold
      use old_rho

      implicit none

C      real*8,allocatable    :: rhopij_comp(:,:),rhonij_comp(:,:)
C      real,allocatable      :: psdMin(:,:),nsdMin(:,:)
      real,dimension(nsps(1),nsps(1)) :: hpp,gammap
      real,dimension(nsps(2),nsps(2)) :: hnn,gamman
      real                            :: energp(nsps(1)),energn(nsps(2))
      real,allocatable      :: vec(:,:),work(:)
      real*8,allocatable    :: psil(:,:),psir(:,:)
      real                  :: gammap1(nsps(1),nsps(1)),
     &                         gamman1(nsps(2),nsps(2))
c      real*8 overp,overn
      real*8 denergy,sum
      real   over
      real,intent(OUT)      :: energy
      real hme

      logical replace
      
      integer i,m,j
      

      call make_projectors(psd,rhop,qprot,np,nsps(1),overp)
      call make_projectors(nsd,rhon,qneutr,nn,nsps(2),overn)
      call make_Gamma(gammap,gamman)

!      do i = 1,nsps(1)
!        write(18,103)(rhop(i,j),j=1,nsps(1))
!103     format(6f10.4)
!      enddo

      call make_hmatrix(gammap,e_spepp,hpp,nsps(1))      
      call make_hmatrix(gamman,e_spenn,hnn,nsps(2))      

      if(.not.allocated(oldrhop)) allocate(oldrhop(nsps(1),nsps(1)))
      if(.not.allocated(oldrhon)) allocate(oldrhon(nsps(2),nsps(2)))
!      call commuter(nsps(1),oldrhop,hpp)

      do i = 1,nsps(1)
         do j = 1,nsps(1)
            oldrhop(i,j) = rhop(i,j)
         enddo
      enddo
!      call commuter(nsps(1),oldrhop,hpp)

      if(diag1)then
         gammap=0.9*gammap+0.1*gammap_old
         gamman=0.9*gamman+0.1*gamman_old
         diag1=.true.
      endif
      gammap_old=gammap
      gamman_old=gamman

      allocate(vec(nsps(1),nsps(1)),work(nsps(1)))
      call eig(hpp,nsps(1),nsps(1),energp,vec,work)
      call eigsrt(energp,vec,nsps(1),nsps(1))
C      rewind(15)
      do i=1,nsps(1),2
C      write(15,*)energp(i)-energp(1)
      enddo
C      rewind(14)
      do i=1,np
      sum=0.0
      do m=1,8
         sum=sum+vec(m,i)**2
      enddo
C      write(14,*)i,sum
C      write(14,'(40(F4.2,1x))')vec(:,i)
      enddo



      do i=1,np                                ! Fill SD ---- protons
         do m=1,nsps(1)
            psd(m,i)=vec(m,i)
	 enddo
      enddo  
C      rewind(10)
C      write(10,1000)energp    
      deallocate(vec,work)
        
      allocate(vec(nsps(2),nsps(2)),work(nsps(2)))
      call eig(hnn,nsps(2),nsps(2),energn,vec,work)
      call eigsrt(energn,vec,nsps(2),nsps(2))
     
      do i=1,nn                                ! Fill SD ---- neutrons
         do m=1,nsps(2)
	    nsd(m,i)=vec(m,i)
	 enddo
      enddo

      deallocate(vec,work)


C      call make_projectors(psd,rhop,qprot,np,nsps(1),overp)
C      call make_projectors(nsd,rhon,qneutr,nn,nsps(2),overn)
C      call make_Gamma(gammap,gamman)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C E_HF=Tr[rho*(T+Gamma/2)]
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      denergy = 0.d0
      do i=1,nsps(1)
         do j=1,nsps(1)
	    denergy=denergy+rhop(i,j)*dble(e_spepp(j,i))
	    denergy=denergy+0.5*rhop(i,j)*dble(gammap(j,i))
	 enddo
      enddo
      do i=1,nsps(2)
         do j=1,nsps(2)
	    denergy=denergy+rhon(i,j)*dble(e_spenn(j,i))
	    denergy=denergy+0.5*rhon(i,j)*dble(gamman(j,i))
	 enddo
      enddo
      energy=sngl(denergy)

      
!      return

CIS I don't believe the next lines are necessary
            
      call make_projectors(psd,rhop,qprot,np,nsps(1),overp)
      call make_projectors(nsd,rhon,qneutr,nn,nsps(2),overn)
      call make_Gamma(gammap,gamman)

!      do i = 1,nsps(1)
!        write(18,103)(rhop(i,j),j=1,nsps(1))
!103     format(6f10.4)
!      enddo

      call make_hmatrix(gammap,e_spepp,hpp,nsps(1))      
      call make_hmatrix(gamman,e_spenn,hnn,nsps(2))      
      overp=1.0
      overn=1.0
      
!      call commuter(nsps(1),oldrhop,hpp)
!      print*,' '
!      call hmult(nsps,nmatpp,nmatnn,nmatpn,
!     &  map_combpp,map_combnn,map_combpn,
!     &  hvecpp,hvecnn,hvecpn,e_spepp,e_spenn,
!     &  overp,overn,rhop,rhon,over,hme)


  
c      write(6,*)energy,hme
!      if(abs(energy-hme)>1.E-3)then
!          print*,'Warning:',hme,energy
c          stop 'Energies do not match'
!      endif
! 1000 format(12F12.6)
!      energy=hme

      return
      end

!=================================================================

      subroutine commuter(n,array1,array2)

      implicit none
      integer n
      real :: array1(n,n)
      real ::   array2(n,n)

      integer i,j,k
      real :: tmp,sum2

      sum2 = 0.0
      do i = 1,n
         do j = 1,n
            tmp = 0.0
            do k = 1,n
               tmp = tmp + array1(i,k)*array2(k,j) 
     &                  - array2(i,k)*array1(k,j)
            enddo
            sum2 = sum2 + tmp*tmp
         enddo ! j

      enddo !i

      sum2 = sum2/n**2
      sum2 = sqrt(sum2)
      print*,' rms of [rho,h] = ',sum2
      return

      end
!==============================================================


