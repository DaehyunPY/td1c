!######################################################################
subroutine zcopy_omp(n, x, y)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_long), intent(in) :: n
  complex(c_double_complex), intent(in) :: x(1:n)
  complex(c_double_complex), intent(out) :: y(1:n)

  integer(c_long) :: i

!!NOOMP  !$omp parallel default(shared) private(i)
!!NOOMP  !$omp do
  do i = 1, n
     y(i) = x(i)
  end do
!NOOMP  !$omp end do
!NOOMP  !$omp end parallel

end subroutine zcopy_omp
!######################################################################
