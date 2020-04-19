!################################################################################
subroutine util_diag_comp(dsc, n, a, u)

  use, intrinsic :: iso_c_binding

  implicit none
  logical(c_bool), intent(in) :: dsc
  integer(c_int), intent(in) :: n
  complex(c_double_complex), intent(inout) :: a(1:n, 1:n)
  complex(c_double_complex), intent(out) :: u(1:n, 1:n)

  integer(c_int) :: info, lwork, i, len
  complex(c_double_complex) :: clwork

  complex(c_double_complex), allocatable :: tmpa(:, :)
  complex(c_double_complex), allocatable :: work(:)
  real(c_double), allocatable :: eig(:), rwork(:)
!debug
!debug  integer(c_int) :: j, k
!debug

  len = max(n*n, 3*n-2)
  allocate(rwork(len))
  allocate(tmpa(n, n))
  allocate(eig(n))
  tmpa(1:n, 1:n) = a(1:n, 1:n)

!debug
!debugwrite(6, "('A: A')")
!debugwrite(6, "(2f20.10)") tmpa(1:n, 1:n)
!debug
  call zheev("v", "u", n, tmpa, n, eig, clwork, -1, rwork, info)
  lwork = int(clwork)

  allocate(work(lwork))
  call zheev("v", "u", n, tmpa, n, eig, work, lwork, rwork, info)
!debug
!debugwrite(6, "('A: A = U * a * U+')")
!debugdo i = 1, n
!debug   do j = 1, n
!debug      a(i, j) = (0.d+0, 0.d+0)
!debug      do k = 1, n
!debug         a(i, j) = a(i, j) + tmpa(i, k) * eig(k) * conjg(tmpa(j, k))
!debug      end do
!debug   end do
!debugend do
!debugwrite(6, "(2f20.10)") a(1:n, 1:n)
!debug

  if (info /= 0) then
     write(6, "('util_diag: info = ', i20)") info
     stop 'error in util_diag_comp.'
  else
     u(1:n, 1:n) = (0.d+0, 0.d+0)
     a(1:n, 1:n) = (0.d+0, 0.d+0)
     if (dsc) then
        do i = 1, n
           a(i, i) = eig(n - i + 1)
           u(1:n, i) = tmpa(1:n, n - i + 1)
        end do
     else
        do i = 1, n
           a(i, i) = eig(i)
           u(1:n, i) = tmpa(1:n, i)
        end do
     end if
  end if

!debug
!write(6,"('diag: eigenvalues:')")
!write(6,"(F20.10)") eig(1:n)
!debug

  deallocate(work)
  deallocate(eig)
  deallocate(tmpa)
  deallocate(rwork)

end subroutine util_diag_comp
!################################################################################
!################################################################################
subroutine util_diag_comp2(dsc, n, s, a, u)

  use, intrinsic :: iso_c_binding

  implicit none
  logical(c_bool), intent(in) :: dsc
  integer(c_int), intent(in) :: n
  complex(c_double_complex), intent(in) :: s(1:n, 1:n)
  complex(c_double_complex), intent(inout) :: a(1:n, 1:n)
  complex(c_double_complex), intent(out) :: u(1:n, 1:n)

  integer(c_int) :: info, lwork, i, len
  complex(c_double_complex) :: clwork

  complex(c_double_complex), allocatable :: tmpa(:, :)
  complex(c_double_complex), allocatable :: work(:)
  real(c_double), allocatable :: eig(:), rwork(:)
!debug
!debug  integer(c_int) :: j, k
!debug

  len = max(n*n, 3*n-2)
  allocate(rwork(len))
  allocate(tmpa(n, n))
  allocate(eig(n))
  tmpa(1:n, 1:n) = a(1:n, 1:n)

!debug
!debugwrite(6, "('A: A')")
!debugwrite(6, "(2f20.10)") tmpa(1:n, 1:n)
!debug
  call zhegv(1, "v", "u", n, tmpa, n, s, n, eig, &
       & clwork, -1, rwork, info)
  lwork = int(clwork)

  allocate(work(lwork))
  call zhegv(1, "v", "u", n, tmpa, n, s, n, eig, &
       & work, lwork, rwork, info)
!debug
!debugwrite(6, "('A: A = U * a * U+')")
!debugdo i = 1, n
!debug   do j = 1, n
!debug      a(i, j) = (0.d+0, 0.d+0)
!debug      do k = 1, n
!debug         a(i, j) = a(i, j) + tmpa(i, k) * eig(k) * conjg(tmpa(j, k))
!debug      end do
!debug   end do
!debugend do
!debugwrite(6, "(2f20.10)") a(1:n, 1:n)
!debug

  if (info /= 0) then
     write(6, "('util_diag: info = ', i20)") info
     stop 'error in util_diag_comp.'
  else
     u(1:n, 1:n) = (0.d+0, 0.d+0)
     a(1:n, 1:n) = (0.d+0, 0.d+0)
     if (dsc) then
        do i = 1, n
           a(i, i) = eig(n - i + 1)
           u(1:n, i) = tmpa(1:n, n - i + 1)
        end do
     else
        do i = 1, n
           a(i, i) = eig(i)
           u(1:n, i) = tmpa(1:n, i)
        end do
     end if
  end if

!debug
!write(6,"('diag: eigenvalues:')")
!write(6,"(F20.10)") eig(1:n)
!debug

  deallocate(work)
  deallocate(eig)
  deallocate(tmpa)
  deallocate(rwork)

end subroutine util_diag_comp2
!################################################################################
subroutine util_diag_real(dsc, n, a, u)

  use, intrinsic :: iso_c_binding

  implicit none
  logical(c_bool), intent(in) :: dsc
  integer(c_int), intent(in) :: n
  real(c_double), intent(inout) :: a(1:n, 1:n)
  real(c_double), intent(out) :: u(1:n, 1:n)

  integer(c_int) :: info, lwork, i, len
  real(c_double) :: rlwork

  real(c_double), allocatable :: tmpa(:, :)
  real(c_double), allocatable :: work(:)
  real(c_double), allocatable :: eig(:)

  len = max(n*n, 3*n-2)
  allocate(tmpa(n, n))
  allocate(eig(n))
  tmpa(1:n, 1:n) = a(1:n, 1:n)

  call dsyev("v", "u", n, tmpa, n, eig, rlwork, -1, info)
  lwork = int(rlwork)

  allocate(work(lwork))
  call dsyev("v", "u", n, tmpa, n, eig, work, lwork, info)

  if (info /= 0) then
     write(6, "('util_diag: info = ', i20)") info
     stop 'error in util_diag.'
  else
     a(1:n, 1:n) = 0.d+0
     u(1:n, 1:n) = 0.d+0
     if (dsc) then
        do i = 1, n
           a(i, i) = eig(n - i + 1)
           u(1:n, i) = tmpa(1:n, n - i + 1)
        end do
     else
        do i = 1, n
           a(i, i) = eig(i)
           u(1:n, i) = tmpa(1:n, i)
        end do
     end if
  end if
  deallocate(work)

  deallocate(eig)
  deallocate(tmpa)

end subroutine util_diag_real
!################################################################################
!################################################################################
subroutine util_diag_real2(dsc, n, s, a, u)

  use, intrinsic :: iso_c_binding

  implicit none
  logical(c_bool), intent(in) :: dsc
  integer(c_int), intent(in) :: n
  real(c_double), intent(in) :: s(1:n, 1:n)
  real(c_double), intent(inout) :: a(1:n, 1:n)
  real(c_double), intent(out) :: u(1:n, 1:n)

  integer(c_int) :: info, lwork, i, len
  real(c_double) :: rlwork

  real(c_double), allocatable :: tmpa(:, :)
  real(c_double), allocatable :: work(:)
  real(c_double), allocatable :: eig(:)

  len = max(n*n, 3*n-2)
  allocate(tmpa(n, n))
  allocate(eig(n))
  tmpa(1:n, 1:n) = a(1:n, 1:n)

  call dsygv(1, "v", "u", n, tmpa, n, s, n, eig, rlwork, -1, info)
  lwork = int(rlwork)

  allocate(work(lwork))
  call dsygv(1, "v", "u", n, tmpa, n, s, n, eig, work, lwork, info)

  if (info /= 0) then
     write(6, "('util_diag_real2: info = ', i20)") info
     stop 'error in util_diag_real2'
  else
     u(1:n, 1:n) = 0.d+0
     a(1:n, 1:n) = 0.d+0
     if (dsc) then
        do i = 1, n
           a(i, i) = eig(n - i + 1)
           u(1:n, i) = tmpa(1:n, n - i + 1)
        end do
     else
        do i = 1, n
           a(i, i) = eig(i)
           u(1:n, i) = tmpa(1:n, i)
        end do
     end if
  end if

!debug
!write(6,"('diag: eigenvalues:')")
!write(6,"(F20.10)") eig(1:n)
!debug

  deallocate(work)
  deallocate(eig)
  deallocate(tmpa)

end subroutine util_diag_real2
!################################################################################
