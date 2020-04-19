!######################################################################
subroutine zclear_omp(n, x)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_int), intent(in) :: n
  complex(c_double_complex), intent(out) :: x(1:n)

  integer(c_int) :: i
  complex(c_double_complex) :: czero = (0.d+0, 0.d+0)

!!NOOMP  !$omp parallel default(shared) private(i)
!!NOOMP  !$omp do
  do i = 1, n
     x(i) = czero
  end do
!!NOOMP  !$omp end do
!!NOOMP  !$omp end parallel

end subroutine zclear_omp
!######################################################################
subroutine dclear_omp(n, x)

  use, intrinsic :: iso_c_binding
  use mod_const, only : zero

  implicit none
  integer(c_int), intent(in) :: n
  real(c_double), intent(out) :: x(1:n)

  integer(c_int) :: i

  !$omp parallel default(shared) private(i)
  !$omp do
  do i = 1, n
     x(i) = zero
  end do
  !$omp end do
  !$omp end parallel

end subroutine dclear_omp
!######################################################################
