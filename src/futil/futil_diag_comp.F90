!################################################################################
subroutine futil_diag_comp(dsc, n, a, u)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero

  implicit none
  logical, intent(in) :: dsc
  integer(c_long), intent(in) :: n
  complex(c_double_complex), intent(inout) :: a(1:n, 1:n)
  complex(c_double_complex), intent(out) :: u(1:n, 1:n)

  integer(c_long) :: info, lwork, i, len
  complex(c_double_complex) :: clwork

  complex(c_double_complex), allocatable :: tmpa(:, :)
  complex(c_double_complex), allocatable :: work(:)
  real(c_double), allocatable :: eig(:), rwork(:)

  len = max(n*n, 3*n-2)
  allocate(rwork(len))
  allocate(tmpa(n, n))
  allocate(eig(n))
  tmpa(1:n, 1:n) = a(1:n, 1:n)

  call zheev("v", "u", n, tmpa, n, eig, clwork, -1, rwork, info)
  lwork = int(clwork)

  allocate(work(lwork))
  call zheev("v", "u", n, tmpa, n, eig, work, lwork, rwork, info)
!debug
!debugwrite(6, "('A: A = U * a * U+')")
!debugdo i = 1, n
!debug   do j = 1, n
!debug      a(i, j) = czero
!debug      do k = 1, n
!debug         a(i, j) = a(i, j) + tmpa(i, k) * eig(k) * conjg(tmpa(j, k))
!debug      end do
!debug   end do
!debugend do
!debugwrite(6, "(2f20.10)") a(1:n, 1:n)
!debug

  if (info /= 0) then
     write(6, "('futil_diag: info = ', i20)") info
     stop 'error in futil_diag_comp.'
  else
     u(1:n, 1:n) = czero
     a(1:n, 1:n) = czero
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

end subroutine futil_diag_comp
!################################################################################
