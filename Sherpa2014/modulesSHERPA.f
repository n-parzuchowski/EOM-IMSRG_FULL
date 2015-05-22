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
         logical pnRPAflag
         logical noRPA      ! true if RPA/pnRPA cannot be performed
                            !      because one cannot set the RPA matrices
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
      				! spsqnp(it,i,j)  it = 1,2 proton/neutron  
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
         implicit none
C-----------HAMILTONIAN------------------------------------------
         integer                  :: nmatpp		! # of pp matrix elements 
         integer                  :: nmatnn		! # of nn matrix elements
         integer                  :: nmatpn		! # of np matrix elements

C-------- list of packed 4-indices with nonzero TBME     
         integer,allocatable,dimension(:):: map_combpp,map_combnn
         integer,allocatable,dimension(:):: map_combpn

C---------- nonzero TBMEs associated with above list   
         real,allocatable,dimension(:)  :: hvecpp,hvecnn,hvecpn
      
C---------one-body matrix elements      
         real,allocatable,target,dimension(:,:)   :: e_spepp,e_spenn
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
C-----------OBSERVABLE--------------------------------------------------
         integer                      :: nmatpp_O		! # of pp matrix elements 
         integer                      :: nmatnn_O		! # of nn matrix elements
         integer                      :: nmatpn_O		! # of np matrix elements

C-------- list of packed 4-indices with nonzero TBME     
         integer,allocatable          :: map_combpp_O(:),
     &                                   map_combnn_O(:),
     &                                   map_combpn_O(:)
      
C---------- nonzero TBMEs associated with above list   
         real,allocatable        :: hvecpp_O(:),hvecnn_O(:),
     &                              hvecpn_O(:)
 
C---------one-body matrix elements      
         real,allocatable        :: e_spepp_O(:,:),e_spenn_O(:,:)
         real                    :: obsval
      end module observables

      module RPAmatrices
         integer                         :: nphp,nphn,nph ! number of proton/neutron matrix elements
                                                          ! in ph basis, as well as their sum
         real,allocatable,dimension(:,:) :: App,Ann,Apn
         real,allocatable,dimension(:,:) :: Bpp,Bnn,Bpn
         real,allocatable,dimension(:,:) :: A,B
         integer                         :: nphpn,nphnp,nph_pn  ! number of proton/neutron matrix elements
                                                          ! in ph basis, as well as their sum
         real,allocatable,dimension(:,:) :: Apnnp,Anppn
         real,allocatable,dimension(:,:) :: Bpnnp,Bnppn
      end module RPAmatrices

      module RPAsolutions
         real,allocatable,dimension(:),target   :: W0
         real,pointer,dimension(:)              :: W1
         real,allocatable,dimension(:,:)        :: X,Y
         integer                                :: n0
         real,allocatable,dimension(:,:),target :: xpn,xnp,ypn,ynp
         real,allocatable,dimension(:)          :: wpn,wnp
         integer                                :: dim_pn,dim_np
      end module RPAsolutions

      module phBASIS
      real,allocatable    :: vecp(:,:),vecn(:,:)
                                                             ! hpp/hnn eigenvectors
      real,allocatable    :: energn(:),energp(:) ! s.p. energies in ph basis
      integer,allocatable,target :: map_prot(:,:),map_neutr(:,:)
      integer,allocatable        :: map_pn(:,:),map_np(:,:)
                            ! map of ph states for protons/neutrons
      end module phBASIS


      module obop        ! one body operator for transitions
         integer ITROUT             ! the file to be used to write out the
                                    ! sum rule
         real cp,cn                 ! effective charges for proton, neutron
         real EHF                   ! the HF energy
         integer JJ                 ! L of Multipole operator
         integer TT                 ! T of the operator
         real,allocatable             :: opp(:,:),onn(:,:)
         real,allocatable             :: opn(:,:)   ! for beta decays only
      end module obop

      module hfdiagold
         real,allocatable             :: gammap_old(:,:),gamman_old(:,:)
         logical                      :: diag1
      end module hfdiagold

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      module old_rho
      implicit none

      real, allocatable :: oldrhop(:,:),oldrhon(:,:)

      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC










