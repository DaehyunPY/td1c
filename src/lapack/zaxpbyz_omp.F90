!######################################################################
subroutine zaxpbyz_omp(n, a, x, b, y, z)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_long), intent(in) :: n
  complex(c_double_complex), intent(in) :: a, b
  complex(c_double_complex), intent(in) :: x(1:n)
  complex(c_double_complex), intent(in) :: y(1:n)
  complex(c_double_complex), intent(out) :: z(1:n)

  integer(c_long) :: i

  !$omp parallel default(shared) private(i)
  !$omp do
  do i = 1, n
     z(i) = a * x(i) + b * y(i)
  end do
  !$omp end do
  !$omp end parallel

end subroutine zaxpbyz_omp
!######################################################################
