      module flags
         implicit none
         logical HFFLAG	    ! = true if hf state exists
			    ! = false otherwise
         logical nukeflag   ! = true if nucleus chosen
			    ! (s.p. space, valence N,Z specified)
			    ! = false otherwise
         logical hamflag    ! if hamiltonian read in
         logical RPAflag    ! if RPA matrices computed
         logical obsflag
      end module flags

      module RPAstab
         implicit none
         real,allocatable    :: stabe(:) ! stability matrix eigenvalues
         integer             :: neigen   ! # of eigenvalues
         real,parameter      :: cutoff=-0.01
      end module RPAstab

      module spspace
         implicit none      
C-------------SINGLE-PARTICLE STATES---------------------------   

         integer norb(2)		! # of single-particle j-orbits
         				! 1 = p,  2 = n
         integer nsps(2)		! # of single-particle m-states
      		   		        ! 1 = p, 2 = n
      
         integer,allocatable  :: jorb(:)	! 2 x J of orbit (note: p,n different)
         integer,allocatable  :: nrorb(:)	! radial quantum number N of j-orbit
         integer,allocatable  :: torb(:)	! 2 x Tz of orbit ( p/n = +/-1 )
      
         real,allocatable     :: orb_qn(:,:)	! quantum #'s of j-orbits  -- from WE0
      				! note: p,n have different orbits
      				! orb_qn(i,1) = N of ith orbit
      				! orb_qn(i,2) = L of ith orbit
      				! orb_qn(i,3) = J of ith orbit

         integer,allocatable  :: spsqn(:,:,:)
      				! quantum #'s of s.p. m-states - from WEO
      				! spsqnp(i,j)  
      				! i = label of s.p. m-state
      				! j = 1 -> label of j-orbit (note p,n different)
      				! j = 2 -> radial quantum number N
      				! j = 3 -> L 
      				! j = 4 -> 2 x J of state
      				! j = 5 -> 2 x M of state
      				! j = 6 -> 2 x Tz of state ( p/n=+/-1 )
         logical spinless
      end module spspace

      module chf_dim
         implicit none
         integer              :: nindx
         integer              :: j2max     ! maximum 2 X J for s.p.'s
         integer,allocatable  :: kmax(:,:) ! given jl1 x jl2 = ij, ang mom k of ij 
         integer,allocatable  :: kmin(:,:) ! must be between kmin(1,2) and kmax(1,2)
         integer,allocatable  :: indx(:,:)
         real*8,allocatable   :: clb(:,:,:)    ! Clebsch-Gordon
         real*8,allocatable   :: jphase(:)     ! useful phase
      end module chf_dim

      module nuclide
         implicit none
         integer                          :: np,nn   ! # of protons/neutrons
      end module nuclide


      module HAMILTONIAN
      interface
      subroutine make_hvec_iso(nmatpp,nmatnn,nmatpn,
     &  map_combpp,map_combnn,map_combpn,
     &  hvecpp,hvecnn,hvecpn,e_spepp,e_spenn)
C-----------HAMILTONIAN------------------------------------------
         integer                  :: nmatpp		! # of pp matrix elements 
         integer                  :: nmatnn		! # of nn matrix elements
         integer                  :: nmatpn		! # of np matrix elements

C-------- list of packed 4-indices with nonzero TBME     
         integer,pointer,dimension(:):: map_combpp,map_combnn
         integer,pointer,dimension(:):: map_combpn

C---------- nonzero TBMEs associated with above list   
         real,pointer,dimension(:)  :: hvecpp,hvecnn,hvecpn
      
C---------one-body matrix elements      
         real,pointer,dimension(:,:)   :: e_spepp,e_spenn
      end subroutine make_hvec_iso
      end interface
C-----------HAMILTONIAN------------------------------------------
      integer                  :: nmatpp		! # of pp matrix elements 
      integer                  :: nmatnn		! # of nn matrix elements
      integer                  :: nmatpn		! # of np matrix elements

C-------- list of packed 4-indices with nonzero TBME     
      integer,pointer,dimension(:):: map_combpp,map_combnn
      integer,pointer,dimension(:):: map_combpn

C---------- nonzero TBMEs associated with above list   
      real,pointer,dimension(:)  :: hvecpp,hvecnn,hvecpn
      
C---------one-body matrix elements      
      real,pointer,dimension(:,:)   :: e_spepp,e_spenn
      end module HAMILTONIAN

      module wfn
         implicit none
         real,allocatable,dimension(:,:)  :: pSD,nSD ! proton/neutron SD's      
          				         ! columns = particles, rows = m-states
         real*8,allocatable,dimension(:,:):: rhop, rhon
                                                 ! density matrix for protons/neutrons
         real*8 overp,overn
         real*8,allocatable,dimension(:,:):: qprot,qneutr ! Q operator
                                                          ! Q=1-rho 
      end module wfn

      module observables
      implicit none

C-----------J operator--------------------------------------------------
         real,allocatable,dimension(:,:) ::jx_matrix,jy_matrix,jz_matrix
         logical                         :: jxyz ! true if the Jx,Jy,Jz matrices have
                                                 ! been already computed
      interface
      subroutine make_hvec_iso(nmatpp_O,nmatnn_O,nmatpn_O,
     &  map_combpp_O,map_combnn_O,map_combpn_O,
     &  hvecpp_O,hvecnn_O,hvecpn_O,e_spepp_O,e_spenn_O)
C-----------OBSERVABLE--------------------------------------------------
         integer                      :: nmatpp_O		! # of pp matrix elements 
         integer                      :: nmatnn_O		! # of nn matrix elements
         integer                      :: nmatpn_O		! # of np matrix elements

C-------- list of packed 4-indices with nonzero TBME     
         integer,pointer   :: map_combpp_O(:),
     &                                   map_combnn_O(:),
     &                                   map_combpn_O(:)
      
C---------- nonzero TBMEs associated with above list   
         real,pointer        :: hvecpp_O(:),hvecnn_O(:),
     &                              hvecpn_O(:)
 
C---------one-body matrix elements      
         real,pointer        :: e_spepp_O(:,:),e_spenn_O(:,:)
      end subroutine make_hvec_iso
      end interface
C-----------OBSERVABLE--------------------------------------------------
         integer                      :: nmatpp_O		! # of pp matrix elements 
         integer                      :: nmatnn_O		! # of nn matrix elements
         integer                      :: nmatpn_O		! # of np matrix elements

C-------- list of packed 4-indices with nonzero TBME     
         integer,pointer   :: map_combpp_O(:),
     &                                   map_combnn_O(:),
     &                                   map_combpn_O(:)
      
C---------- nonzero TBMEs associated with above list   
         real,pointer        :: hvecpp_O(:),hvecnn_O(:),
     &                              hvecpn_O(:)
 
C---------one-body matrix elements      
         real,pointer        :: e_spepp_O(:,:),e_spenn_O(:,:)
      end module observables












