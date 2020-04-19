!######################################################################
subroutine zaxpyz_omp(n, a, x, y, z)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_int), intent(in) :: n
  complex(c_double_complex), intent(in) :: a
  complex(c_double_complex), intent(in) :: x(1:n)
  complex(c_double_complex), intent(in) :: y(1:n)
  complex(c_double_complex), intent(out) :: z(1:n)

  integer(c_int) :: i

! use blas
! call zaxpy(n, a, x, 1, y, 1)

! in case inc = 0
! do i = 1, n
!    y(i) = y(i) + a * x(1)
! end do

  !$omp parallel default(shared) private(i)
  !$omp do
  do i = 1, n
     z(i) = y(i) + a * x(i)
  end do
  !$omp end do
  !$omp end parallel

end subroutine zaxpyz_omp
!######################################################################
