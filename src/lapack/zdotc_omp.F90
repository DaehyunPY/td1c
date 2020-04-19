!######################################################################
subroutine zdotc_omp(n, x, y, val)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_long), intent(in) :: n
  complex(c_double_complex), intent(in) :: x(1:n)
  complex(c_double_complex), intent(in) :: y(1:n)
  complex(c_double_complex), intent(out) :: val

  integer(c_long) :: i
  complex(c_double_complex) :: tmp
  complex(c_double_complex) :: sum

! use blas
! call zaxpy(n, a, x, 1, y, 1)

! in case inc = 0
! do i = 1, n
!    y(i) = y(i) + a * x(1)
! end do

  sum = (0.d+0, 0.d+0)
  !$omp parallel default(shared) private(i, tmp) reduction(+:sum)
  tmp = (0.d+0, 0.d+0)
  !$omp do
  do i = 1, n
     tmp = tmp + conjg(x(i)) * y(i)
  end do
  !$omp end do
  sum = sum + tmp;
  !$omp end parallel

  val = sum

end subroutine zdotc_omp
!######################################################################
