!################################################################################
subroutine util_zexpax(n, a, x, expax)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_long), intent(in) :: n
  complex(c_double_complex), intent(in) :: a
  complex(c_double_complex), intent(inout) :: x(1:n, 1:n)
  complex(c_double_complex), intent(out) :: expax(1:n, 1:n)

  integer(c_long) :: i, j, k
  complex(c_double_complex) :: fac
  complex(c_double_complex), allocatable :: u(:,:)
  complex(c_double_complex), parameter :: czero = (0.d+0, 0.d+0)

  allocate(u(1:n, 1:n))

  call lapack_zheev(n, x, u)

  expax(1:n, 1:n) = czero
  do k = 1, n
     fac = exp(a * x(k, k))
     do i = 1, n
        do j = 1, n
           expax(j, i) = expax(j, i) + u(j, k) * fac * conjg(u(i, k))
        end do
     end do
  end do

  deallocate(u)

end subroutine util_zexpax
!################################################################################
