!################################################################################
subroutine lapack_dsyev(n, a, u)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_int), intent(in) :: n
  real(c_double), intent(inout) :: a(1:n, 1:n)
  real(c_double), intent(out) :: u(1:n, 1:n)

  integer(c_int) :: nf90, info, i, len, lwork
  real(c_double) :: rlwork

  logical, parameter :: dsc = .false.
  real(c_double), allocatable :: tmpa(:,:)
  real(c_double), allocatable :: work(:)
  real(c_double), allocatable :: eig(:)

  nf90 = n
!  write(6, "('n =         ', i5)") n
!  write(6, "('nf90 =      ', i5)") nf90
!  write(6, "('KIND(0) =   ', i5)") c_long
!  write(6, "('KIND(10) =  ', i5)") kind(10)
!  write(6, "('KIND(0E0) = ', i5)") kind(0e0)
!  write(6, "('KIND(0D0) = ', i5)") kind(0d0)
!  write(6, "('KIND(1.0) = ', i5)") kind(1.0)

  len = max(nf90*nf90, 3*nf90-2)
  allocate(tmpa(1:nf90, 1:nf90))
  allocate(eig(1:nf90))
  tmpa(1:nf90, 1:nf90) = a(1:nf90, 1:nf90)
!debug  tmpa(1:n, 1:n) = 0.d+0
!debug  do i = 1, n
!debug     tmpa(i, i) = 1.d+0
!debug  end do
!debug  tmpa(1, 2) = 1.d-1
!debug  tmpa(2, 1) = 1.d-1

  call dsyev("v", "u", nf90, tmpa, nf90, eig, rlwork, -1, info)
  lwork = int(rlwork)
!debug  write(6, "('info:          ', i20)") info
!debug  write(6, "('rlwork, lwork: ', f20.10, i20)") rlwork, lwork

  allocate(work(lwork))
  call dsyev("v", "u", nf90, tmpa, nf90, eig, work, lwork, info)
!debug  write(6, "('info:          ', i20)") info

  if (info /= 0) then
     write(6, "('util_diag: info = ', i20)") info
     stop 'error in util_diag.'
  else
     a(1:nf90, 1:nf90) = 0.d+0
     u(1:nf90, 1:nf90) = 0.d+0
     if (dsc) then
        do i = 1, nf90
           a(i, i) = eig(nf90 - i + 1)
           u(1:nf90, i) = tmpa(1:nf90, nf90 - i + 1)
        end do
     else
        do i = 1, nf90
           a(i, i) = eig(i)
           u(1:nf90, i) = tmpa(1:nf90, i)
        end do
     end if
  end if

!debug
!write(6,"('diag: eigenvalues:')")
!write(6,"(F20.10)") eig(1:nf90)
!debug

  deallocate(work)
  deallocate(eig)
  deallocate(tmpa)

!stop 'for debug in lapack_ev.'
end subroutine lapack_dsyev
!################################################################################
!################################################################################
subroutine lapack_zheev(n, a, u)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_int), intent(in) :: n
  complex(c_double_complex), intent(inout) :: a(1:n, 1:n)
  complex(c_double_complex), intent(out) :: u(1:n, 1:n)

  integer(c_int) :: nf90, info, lwork, i, len
  complex(c_double_complex) :: clwork

  logical, parameter :: dsc = .false.
  complex(c_double_complex), allocatable :: atmp(:, :)
  complex(c_double_complex), allocatable :: work(:)
  real(c_double), allocatable :: eig(:)
  real(c_double), allocatable :: rwork(:)
!debug
!debug  integer(c_int) :: j, k
!debug

  nf90 = n
  len = max(nf90*nf90, 3*nf90-2)
  allocate(rwork(len))
  allocate(atmp(nf90, nf90))
  allocate(eig(nf90))
  atmp(1:nf90, 1:nf90) = a(1:nf90, 1:nf90)

!debug
!debugwrite(6, "('A: A')")
!debugwrite(6, "(2f20.10)") atmp(1:nf90, 1:nf90)
!debug
  call zheev("v", "u", nf90, atmp, nf90, eig, clwork, -1, rwork, info)
  lwork = int(clwork)

  allocate(work(lwork))
  call zheev("v", "u", nf90, atmp, nf90, eig, work, lwork, rwork, info)
!debug
!debugwrite(6, "('A: A = U * a * U+')")
!debugdo i = 1, nf90
!debug   do j = 1, nf90
!debug      a(i, j) = (0.d+0, 0.d+0)
!debug      do k = 1, nf90
!debug         a(i, j) = a(i, j) + atmp(i, k) * eig(k) * conjg(atmp(j, k))
!debug      end do
!debug   end do
!debugend do
!debugwrite(6, "(2f20.10)") a(1:nf90, 1:nf90)
!debug

  if (info /= 0) then
     write(6, "('util_diag: info = ', i20)") info
     stop 'error in util_diag_comp.'
  else
     u(1:nf90, 1:nf90) = (0.d+0, 0.d+0)
     a(1:nf90, 1:nf90) = (0.d+0, 0.d+0)
     if (dsc) then
        do i = 1, nf90
           a(i, i) = eig(nf90 - i + 1)
           u(1:nf90, i) = atmp(1:nf90, nf90 - i + 1)
        end do
     else
        do i = 1, nf90
           a(i, i) = eig(i)
           u(1:nf90, i) = atmp(1:nf90, i)
        end do
     end if
  end if

!debug
!write(6,"('diag: eigenvalues:')")
!write(6,"(F20.10)") eig(1:nf90)
!debug

  deallocate(work)
  deallocate(eig)
  deallocate(atmp)
  deallocate(rwork)

end subroutine lapack_zheev
!################################################################################
subroutine lapack_zheev_dsc(n, a, u)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_int), intent(in) :: n
  complex(c_double_complex), intent(inout) :: a(1:n, 1:n)
  complex(c_double_complex), intent(out) :: u(1:n, 1:n)

  integer(c_int) :: nf90, info, lwork, i, len
  complex(c_double_complex) :: clwork

  logical, parameter :: dsc = .true.
  complex(c_double_complex), allocatable :: atmp(:, :)
  complex(c_double_complex), allocatable :: work(:)
  real(c_double), allocatable :: eig(:)
  real(c_double), allocatable :: rwork(:)
!debug
!debug  integer(c_int) :: j, k
!debug

  nf90 = n
  len = max(nf90*nf90, 3*nf90-2)
  allocate(rwork(len))
  allocate(atmp(nf90, nf90))
  allocate(eig(nf90))
  atmp(1:nf90, 1:nf90) = a(1:nf90, 1:nf90)

!debug
!debugwrite(6, "('A: A')")
!debugwrite(6, "(2f20.10)") atmp(1:nf90, 1:nf90)
!debug
  call zheev("v", "u", nf90, atmp, nf90, eig, clwork, -1, rwork, info)
  lwork = int(clwork)

  allocate(work(lwork))
  call zheev("v", "u", nf90, atmp, nf90, eig, work, lwork, rwork, info)
!debug
!debugwrite(6, "('A: A = U * a * U+')")
!debugdo i = 1, nf90
!debug   do j = 1, nf90
!debug      a(i, j) = (0.d+0, 0.d+0)
!debug      do k = 1, nf90
!debug         a(i, j) = a(i, j) + atmp(i, k) * eig(k) * conjg(atmp(j, k))
!debug      end do
!debug   end do
!debugend do
!debugwrite(6, "(2f20.10)") a(1:nf90, 1:nf90)
!debug

  if (info /= 0) then
     write(6, "('util_diag: info = ', i20)") info
     stop 'error in util_diag_comp.'
  else
     u(1:nf90, 1:nf90) = (0.d+0, 0.d+0)
     a(1:nf90, 1:nf90) = (0.d+0, 0.d+0)
     if (dsc) then
        do i = 1, nf90
           a(i, i) = eig(nf90 - i + 1)
           u(1:nf90, i) = atmp(1:nf90, nf90 - i + 1)
        end do
     else
        do i = 1, nf90
           a(i, i) = eig(i)
           u(1:nf90, i) = atmp(1:nf90, i)
        end do
     end if
  end if

!debug
!write(6,"('diag: eigenvalues:')")
!write(6,"(F20.10)") eig(1:nf90)
!debug

  deallocate(work)
  deallocate(eig)
  deallocate(atmp)
  deallocate(rwork)

end subroutine lapack_zheev_dsc
!################################################################################
