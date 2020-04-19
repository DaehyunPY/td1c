!######################################################################
subroutine hprod_getprojhigh()

  use, intrinsic :: iso_c_binding
  use mod_bas, only : tmat
  use mod_rad, only : nrad
  use mod_sph, only : lmax1
  use mod_hprod, only : projhigh_orbs, projhigh_eigs, projhigh_cutoff, projhigh_nfun, projhigh_ncut

  implicit none
  integer(c_int) :: l, dim,idim,ifun
  complex(c_double_complex), allocatable :: fock(:,:,:)
  complex(c_double_complex), allocatable :: orbs(:,:,:)
  real(c_double_complex), allocatable :: eigs(:,:)

  dim = nrad - 1
  allocate(fock(1:dim, 1:dim, 0:lmax1))
  allocate(orbs(1:dim, 1:dim, 0:lmax1))
  allocate(eigs(1:dim, 0:lmax1))
  fock(1:dim, 1:dim, 0:lmax1) = 0d0
  orbs(1:dim, 1:dim, 0:lmax1) = 0d0
  eigs(1:dim, 0:lmax1) = 0d0
  call hprod_fulldiag_fock1(fock)
  call hprod_getprojhigh_diag(fock, orbs, eigs)

  allocate(projhigh_ncut(0:lmax1))
  do l = 0, lmax1
     projhigh_ncut(l) = 1
     do idim = 1, dim
        if (eigs(idim,l) > projhigh_cutoff) then
           projhigh_ncut(l) = idim
           exit
        end if
     end do
  end do

  projhigh_nfun = dim - minval(projhigh_ncut) + 1
  allocate(projhigh_orbs(1:dim,1:projhigh_nfun,0:lmax1))
  allocate(projhigh_eigs(      1:projhigh_nfun,0:lmax1))

  projhigh_orbs(1:dim, 1:projhigh_nfun, 0:lmax1) = 0d0
  projhigh_eigs(       1:projhigh_nfun, 0:lmax1) = 0d0
  do l = 0, lmax1
     ifun = 0
     do idim = projhigh_ncut(l), dim
        ifun = ifun + 1
        projhigh_orbs(:,ifun,l) = orbs(:,idim,l)
        projhigh_eigs(  ifun,l) = eigs(  idim,l)
     end do
     projhigh_ncut(l) = dim - projhigh_ncut(l) + 1
  end do

!debug
  write(6, "('hprod_getprojhigh_diag: removed high-energy states')")
  do l = 0, lmax1
     write(6,"('# l = ',i5,', ncut = ',i5)") l,projhigh_ncut(l)
     do ifun = 1, projhigh_ncut(l)
        write(6,"(i5,e15.7)") ifun,projhigh_eigs(ifun,l)
     end do
     write(6,*)
     write(6,*)
  end do
!debug
  
  deallocate(fock)
  deallocate(orbs)
  deallocate(eigs)

end subroutine hprod_getprojhigh
!######################################################################
subroutine hprod_getprojhigh_diag(fock, orbs, eigs)
!
  use, intrinsic :: iso_c_binding
  use mod_const, only : czero, runit, iunit
  use mod_rad, only : nrad
  use mod_sph, only : lmax1

  implicit none
  ! ### args ###
  complex(c_double_complex), intent(in) :: fock(1:(nrad-1), 1:(nrad-1), 0:lmax1)
  complex(c_double_complex), intent(out) :: orbs(1:(nrad-1), 1:(nrad-1), 0:lmax1)
  real(c_double), intent(out) :: eigs(1:(nrad-1), 0:lmax1)

  ! ### local ###
  character(len = 9), parameter :: fname = "HFOrb.dat"
  integer(c_int), parameter :: iow = 10
  integer(c_int) :: dim, l, ifun, irad
  complex(c_double_complex), allocatable :: ftmp(:,:), utmp(:,:)

  dim = nrad - 1
  allocate(ftmp(1:dim, 1:dim))
  allocate(utmp(1:dim, 1:dim))

  open(unit = iow, file = trim(fname), status = 'unknown', form = 'formatted')

  do l = 0, lmax1
     utmp(1:dim, 1:dim) = czero
     ftmp(1:dim, 1:dim) = fock(1:dim, 1:dim, l)

     call lapack_zheev(dim, ftmp, utmp)
     do ifun = 1, dim
        eigs(ifun,l) = dble(ftmp(ifun,ifun))
     end do
     orbs(:,:,l) = utmp(:,:)
  end do

  !debug
  write(6, "('hprod_getprojhigh_diag: Eigenvalues')")
  do ifun = 1, dim
     write(6, "(i5)", advance ='no') ifun
     do l = 0, min(10,lmax1)
        write(6, "(e15.7)", advance='no') eigs(ifun,l)
     end do
     write(6, *)
  end do
  !debug

  deallocate(utmp)
  deallocate(ftmp)

end subroutine hprod_getprojhigh_diag
!######################################################################
