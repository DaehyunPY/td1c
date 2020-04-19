!######################################################################
subroutine zscal_omp(n, a, x)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_int), intent(in) :: n
  complex(c_double_complex), intent(in) :: a
  complex(c_double_complex), intent(inout) :: x(1:n)

  integer(c_int) :: i

! use blas
! call zscal(n, a, x, 1, y, 1)

  !$omp parallel default(shared) private(i)
  !$omp do
  do i = 1, n
     x(i) = a * x(i)
  end do
  !$omp end do
  !$omp end parallel

end subroutine zscal_omp
!######################################################################
subroutine zscal2_omp(n, a, x, y)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_int), intent(in) :: n
  complex(c_double_complex), intent(in) :: a
  complex(c_double_complex), intent(in) :: x(1:n)
  complex(c_double_complex), intent(out) :: y(1:n)

  integer(c_int) :: i

! use blas
! call zscal(n, a, x, 1, y, 1)

  !$omp parallel default(shared) private(i)
  !$omp do
  do i = 1, n
     y(i) = a * x(i)
  end do
  !$omp end do
  !$omp end parallel

end subroutine zscal2_omp
!######################################################################
subroutine dscal_omp(n, a, x)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_int), intent(in) :: n
  real(c_double), intent(in) :: a
  real(c_double), intent(inout) :: x(1:n)

  integer(c_int) :: i

! use blas
! call dscal(n, a, x, 1, y, 1)

  !$omp parallel default(shared) private(i)
  !$omp do
  do i = 1, n
     x(i) = a * x(i)
  end do
  !$omp end do
  !$omp end parallel

end subroutine dscal_omp
!######################################################################
