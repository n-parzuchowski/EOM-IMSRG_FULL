module deuteron
  use basic_IMSRG
  implicit none 

contains
  subroutine compute_deuteron_ground_state(HS,jbas,rho21)    
!!! diagonalization of the 1+ Tz=0 block
    implicit none

    type(sq_op) :: HS
    type(sq_op),optional :: rho21
    type(spd) :: jbas
    integer ::dm,q,np,nb,nh,a,b,c,d,II,JJ
    real(8),allocatable,dimension(:) :: PSI
    real(8) :: Egs,dstate
    real(8),allocatable,dimension(:,:) :: H,rho
    integer,allocatable,dimension(:,:) :: qnums 
    real(8),allocatable,dimension(:) :: workl,DX,QX,eigs,resid,work,workD
    real(8),allocatable,dimension(:,:) :: V,Z
    integer :: lwork,info,ido,ncv,ldv,iparam(11),ipntr(11)
    integer :: ishift,mxiter,nconv,mode,lworkl,ldz,nev,inc
    real(8) :: tol,sigma,sm
    character(1) :: BMAT,HOWMNY 
    character(2) :: which
    logical :: rvec
    logical,allocatable,dimension(:) :: selct


    print*, "SETTING UP DEUTERON MATRIX"
!!! block for the deuteron ground state
    q = block_index(2,0,0)

    np = HS%mat(q)%npp
    nb = HS%mat(q)%nph
    nh = HS%mat(q)%nhh
    dm = np+nb+nh
    allocate(PSI(dm),H(dm,dm),qnums(dm,2))

!!! keep track of the two-body basis 
    qnums(1:nh,:) = HS%mat(q)%qn(3)%Y
    qnums(nh+1:nh+nb,:) = HS%mat(q)%qn(2)%Y
    qnums(nh+nb+1:dm,:) = HS%mat(q)%qn(1)%Y    

!!! construct matrix from two-body matrix elements
    H(1:nh,1:nh) = HS%mat(q)%gam(5)%X
    H(1:nh,nh+1:nh+nb) = transpose(HS%mat(q)%gam(6)%X)
    H(nh+1:nh+nb,1:nh) = HS%mat(q)%gam(6)%X    
    H(nh+1:nh+nb,nh+1:nh+nb) = HS%mat(q)%gam(4)%X
    H(nh+nb+1:dm,1:nh) = HS%mat(q)%gam(3)%X
    H(1:nh,nh+nb+1:dm) = transpose(HS%mat(q)%gam(3)%X)
    H(nh+nb+1:dm,nh+1:nb) = HS%mat(q)%gam(2)%X
    H(nh+1:nb,nh+nb+1:dm) = transpose(HS%mat(q)%gam(2)%X)
    H(nh+nb+1:dm,nh+nb+1:dm) = HS%mat(q)%gam(1)%X

!!! add in the one-body operators

    do II = 1, dm
       a = qnums(II,1)
       b = qnums(II,2)
       do JJ = II,dm
          c = qnums(JJ,1)
          d = qnums(JJ,2)

          H(II,JJ) = H(II,JJ) + T_twobody(a,b,c,d,2,HS,jbas)
          H(JJ,II) = H(II,JJ) 

       end do
    end do
    
    print *, "DIAGONALIZING MATRIX OF SIZE: ",dm

    if (dm > 2000) then
!!! Diagonalize?
       nev = 5 ! I only care about the ground state right now. 
       ido = 0  ! status integer is 0 at start
       BMAT = 'I' ! standard eigenvalue problem (N for generalized) 
       which = 'SM' ! compute smallest eigenvalues in magnitude ('SA') is algebraic. 
       tol = 0.0 ! error tolerance? (wtf zero?) 
       info = 0
       ncv = 5*nev ! number of lanczos vectors I guess
       lworkl = ncv*(ncv+8) 
       allocate(V(dm,NCV),workl(lworkl))
       LDV = dm  
       ishift = 1
       mxiter = 500 
       mode = 1

       allocate(eigs(dm),resid(dm),work(10*dm),workD(3*dm)) 

       iparam(1) = ishift
       iparam(3) = mxiter
       iparam(7) = mode
       ii = 0

       inc = 1
       do 
          ! so V is the krylov subspace matrix that is being diagonalized
          ! it does not need to be initialized, so long as you have the other 
          ! stuff declared right, the code should know this. 
          call dsaupd ( ido, bmat, dm, which, nev, tol, resid, &
               ncv, v, ldv, iparam, ipntr, workd, workl, &
               lworkl, info )
          ! The actual matrix only gets multiplied with the "guess" vector in "matvec_prod" 


          ii=ii+1 

          if ( ido /= -1 .and. ido /= 1 ) then
             exit
          end if

          call dgemv('N',dm,dm,al,H,dm,workd(ipntr(1)),inc,bet,workd(ipntr(2)),inc) 


       end do
       print *, "converged after", ii, "iterations"
       write(6,*) 
       ! the ritz values are out of order right now. Need to do post
       ! processing to fix this, and get the eigenvectors
       rvec= .true. 
       howmny = 'A'

       allocate(selct(NCV)) 
       allocate(DX(NEV)) 
       allocate(Z(dm,NEV)) 
       ldz = dm  
       call dseupd( rvec, howmny, selct, DX, Z, ldv, sigma, &
            bmat, dm, which, nev, tol, resid, ncv, v, ldv, &
            iparam, ipntr, workd, workl, lworkl, info )


       Egs = DX(1)
       PSI = Z(:,1)
       PRINT*, Egs
       
       deallocate(DX)

    else
       allocate(DX(dm))
       allocate(qx(10*dm))
       info = 0
       DX = 0.0
       QX = 0.0
       call dsyev('V','U',dm,H,dm,DX,QX,10*dm,info)
       
       open(unit=36,file = 'wavefunction.dat')
       do ii = 1, dm 
          write(36,*) H(ii,1)
       end do
       close(36)
       PSI = H(:,1) 
    end if

    open(unit=29,file=trim(OUTPUT_DIR)//"deuteron_energy.dat",position="append")
    write(29,*)  HS%hospace, HS%eMax , DX(1)
    close(29)

    print*, 'Energy: ', DX(1)
    if (present(rho21)) then
       ! calculate short range pair density
       
       allocate(rho(dm,dm))
       
!!! construct matrix from two-body matrix elements
       rho(1:nh,1:nh) = rho21%mat(q)%gam(5)%X
       rho(1:nh,nh+1:nh+nb) = transpose(rho21%mat(q)%gam(6)%X)
       rho(nh+1:nh+nb,1:nh) = rho21%mat(q)%gam(6)%X    
       rho(nh+1:nh+nb,nh+1:nh+nb) = rho21%mat(q)%gam(4)%X
       rho(nh+nb+1:dm,1:nh) = rho21%mat(q)%gam(3)%X
       rho(1:nh,nh+nb+1:dm) = transpose(rho21%mat(q)%gam(3)%X)
       rho(nh+nb+1:dm,nh+1:nb) = rho21%mat(q)%gam(2)%X
       rho(nh+1:nb,nh+nb+1:dm) = transpose(rho21%mat(q)%gam(2)%X)
       rho(nh+nb+1:dm,nh+nb+1:dm) = rho21%mat(q)%gam(1)%X
       
       sm = 0.d0 
       do II = 1, dm
          do JJ = 1,dm

             sm = sm + PSI(II)*PSI(JJ) * rho(II,JJ)
          
          end do
       end do

       ! call print_matrix(rho(42:54,120:132)) 

       ! print*, qnums(42:54,1)
       ! print*, qnums(42:54,2)
       ! print*
       ! print*, qnums(120:132,1)
       ! print*, qnums(120:132,2)


       print*
       print*, 'Density: ',sm

    open(unit=29,file=trim(OUTPUT_DIR)//"deuteron_density.dat",position="append")
    write(29,*)  HS%hospace, HS%eMax , sm
    close(29)



    end if
       

  end subroutine compute_deuteron_ground_state


          
          
          
          
          
    

    
end module deuteron
