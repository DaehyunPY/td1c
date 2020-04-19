!######################################################################
subroutine zaxpy_omp(n, a, x, y)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_long), intent(in) :: n
  complex(c_double_complex), intent(in) :: a
  complex(c_double_complex), intent(in) :: x(1:n)
  complex(c_double_complex), intent(inout) :: y(1:n)

  integer(c_long) :: i

! use blas
! call zaxpy(n, a, x, 1, y, 1)

! in case inc = 0
! do i = 1, n
!    y(i) = y(i) + a * x(1)
! end do

  !$omp parallel default(shared) private(i)
  !$omp do
  do i = 1, n
     y(i) = y(i) + a * x(i)
  end do
  !$omp end do
  !$omp end parallel

end subroutine zaxpy_omp
!######################################################################
