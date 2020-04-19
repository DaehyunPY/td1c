!######################################################################
subroutine h1exp_prop(icomp, isplit, igauge, lfield, dtime, maxcyc, thresh, wfn)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : mmax1
  use mod_bas, only : nbas
  use mod_ormas, only : nfun, neltot

  implicit none
  integer(c_int), intent(in) :: icomp, isplit, igauge, maxcyc
  real(c_double), intent(in) :: lfield(3), dtime, thresh
  complex(c_double_complex), intent(inout) :: wfn(1:nbas, 1:nfun)

  integer(c_int), allocatable :: norb(:)
  integer(c_int), allocatable :: mmap(:,:)
  real(c_double), allocatable :: alph(:), beta(:)
  complex(c_double_complex), allocatable :: torb(:,:)
  integer(c_int) :: m, maxorb, maxvec, ntot

  ! #### for hydrogen only ####
  if (neltot(3) > 1) stop "h1exp_prop: only for one electron."

  allocate(norb(-mmax1:mmax1))
  allocate(mmap(1:nfun, -mmax1:mmax1))
  call h1exp_mmap(norb, mmap, maxorb)

  maxvec = maxcyc * maxorb
  allocate(alph(1:maxvec))
  allocate(beta(1:maxvec))
  allocate(torb(1:nbas, 1:maxvec))

  do m = -mmax1, mmax1
     call h1exp_krylov(isplit, igauge, lfield, dtime, maxcyc, thresh, wfn, ntot, alph, beta, torb)
     call h1exp_exporb(icomp, dtime, maxcyc, ntot, alph, beta, torb, wfn)
  end do

  deallocate(torb)
  deallocate(beta)
  deallocate(alph)
  deallocate(mmap)
  deallocate(norb)

!  stop "@ h1exp_prop..."

end subroutine h1exp_prop
!######################################################################
subroutine h1exp_prop1(icomp, isplit, igauge, lfield, dtime, maxcyc, wfn)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : mmax1
  use mod_bas, only : nbas
  use mod_ormas, only : nfun, neltot

  implicit none
  integer(c_int), intent(in) :: icomp, isplit, igauge, maxcyc
  real(c_double), intent(in) :: lfield(3), dtime
  complex(c_double_complex), intent(inout) :: wfn(1:nbas, 1:nfun)

  integer(c_int), allocatable :: norb(:)
  integer(c_int), allocatable :: mmap(:,:)
  real(c_double), allocatable :: norm(:)
  complex(c_double_complex), allocatable :: hmat(:,:)
  complex(c_double_complex), allocatable :: torb(:,:)
  integer(c_int) :: m, maxorb, maxvec, ntot

  ! #### for hydrogen only ####
  if (neltot(3) > 1) stop "h1exp_prop: only for one electron."

  allocate(norb(-mmax1:mmax1))
  allocate(mmap(1:nfun, -mmax1:mmax1))
  call h1exp_mmap(norb, mmap, maxorb)

  maxvec = maxcyc * maxorb
  allocate(norm(1:maxvec))
  allocate(hmat(1:maxvec, 1:maxvec))
  allocate(torb(1:nbas, 1:maxvec))

  do m = -mmax1, mmax1
     ! ##### subspace construction #####
     call h1exp_krylov1(isplit, igauge, lfield, maxcyc, maxvec, norb(m), mmap(1,m), wfn, &
          ntot, norm, torb, hmat)
     ! ##### exponential propagation within the subspace #####
     call h1exp_exporb1(icomp, dtime, maxvec, ntot, norm, torb, hmat, norb(m), mmap(1,m), wfn)
  end do

  deallocate(torb)
  deallocate(hmat)
  deallocate(norm)
  deallocate(mmap)
  deallocate(norb)

!  stop "@ h1exp_prop..."

end subroutine h1exp_prop1
!######################################################################
subroutine h1exp_mmap(norb, mmap, maxorb)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mval
  use mod_sph, only : mmax1
  use mod_ormas, only : nfcore, nfun

  implicit none
  integer(c_int), intent(out) :: norb(-mmax1:mmax1)
  integer(c_int), intent(out) :: mmap(1:nfun, -mmax1:mmax1)
  integer(c_int), intent(out) :: maxorb

  integer(c_int) :: m, ifun

  do m = -mmax1, mmax1
     norb(m) = 0
     do ifun = nfcore + 1, nfun
        if (mval(ifun) == m) then
           norb(m) = norb(m) + 1
           mmap(norb(m), m) = ifun
!          write(6, "('h1exp_mmap: ', 3i5)") m, norb(m), mmap(norb(m), m)
        end if
     end do
  end do
  maxorb = maxval(norb(-mmax1:mmax1))

end subroutine h1exp_mmap
!######################################################################
