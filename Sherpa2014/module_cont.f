      module RPAmatrices
         integer                         :: nphp,nphn,nph
         real,allocatable,dimension(:,:) :: App,Ann,Apn
         real,allocatable,dimension(:,:) :: Bpp,Bnn,Bpn
         real,allocatable,dimension(:,:) :: A,B
      end module RPAmatrices

      module RPAsolutions
         real,allocatable,dimension(:)   :: X,Y,W0,w1
         integer                         :: n0
      end module RPAsolutions


      module RPAstab
         implicit none
         real,allocatable    :: stabe(:) ! stability matrix eigenvalues
         integer             :: neigen   ! # of eigenvalues
         real,parameter      :: cutoff=-0.01
      end module RPAstab
