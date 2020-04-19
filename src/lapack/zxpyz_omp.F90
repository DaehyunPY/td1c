!######################################################################
subroutine zxpyz_omp(n, x, y, z)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_int), intent(in) :: n
  complex(c_double_complex), intent(in) :: x(1:n)
  complex(c_double_complex), intent(in) :: y(1:n)
  complex(c_double_complex), intent(out) :: z(1:n)

  integer(c_int) :: i

  !$omp parallel default(shared) private(i)
  !$omp do
  do i = 1, n
     z(i) = x(i) + y(i)
  end do
  !$omp end do
  !$omp end parallel

end subroutine zxpyz_omp
!######################################################################
