!######################################################################
subroutine zadd_omp(n, x, y)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_long), intent(in) :: n
  complex(c_double_complex), intent(in) :: x(1:n)
  complex(c_double_complex), intent(inout) :: y(1:n)

  integer(c_long) :: i

  !$omp parallel default(shared) private(i)
  !$omp do
  do i = 1, n
     y(i) = y(i) + x(i)
  end do
  !$omp end do
  !$omp end parallel

end subroutine zadd_omp
!######################################################################
