!######################################################################
subroutine util_zcopy(n, x, incx, y, incy)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero

  implicit none
  integer(c_int), intent(in) :: n, incx, incy
  complex(c_double_complex), intent(in) :: x(1:*)
  complex(c_double_complex), intent(out) :: y(1:n)

  integer(c_int) :: i

  if (incx <  0) stop 'util_zcopy: bad incx.'
  if (incy <= 0) stop 'util_zcopy: bad incy.'

  if (incx == 0) then
     do i = 1, n
        y(i) = x(1)
     end do
#ifdef USE_BLAS
  else
     call zcopy(n, x, incx, y, incy)
  end if
#else
  else
     do i = 1, n
        y(i) = x(i)
     end do
  end if
#endif

  return

end subroutine util_zcopy
!######################################################################
!######################################################################
subroutine util_zscal(n, a, x, incx)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero

  implicit none
  integer(c_int), intent(in) :: n, incx
  complex(c_double_complex), intent(in) :: a
  complex(c_double_complex), intent(inout) :: x(1:*)

  integer(c_int) :: i

  if (incx <= 0) stop 'util_zscal: bad incx.'

#ifdef USE_BLAS
  call zscal(n, a, x, incx)
#else
  do i = 1, n
     x(i) = a * x(i)
  end do
#endif

  return

end subroutine util_zscal
!######################################################################
!######################################################################
subroutine util_zaxpy(n, a, x, incx, y, incy)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero

  implicit none
  integer(c_int), intent(in) :: n, incx, incy
  complex(c_double_complex), intent(in) :: a
  complex(c_double_complex), intent(in) :: x(1:*)
  complex(c_double_complex), intent(inout) :: y(1:n)

  integer(c_int) :: i

  if (incx <  0) stop 'util_zaxpy: bad incx.'
  if (incy <= 0) stop 'util_zaxpy: bad incy.'

  if (incx == 0) then
     do i = 1, n
        y(i) = y(i) + a * x(1)
     end do
#ifdef USE_BLAS
  else
     call zaxpy(n, a, x, incx, y, incy)
  end if
#else
  else
     do i = 1, n
        y(i) = y(i) + a * x(i)
     end do
  end if
#endif

  return

end subroutine util_zaxpy
!######################################################################
complex(c_double_complex) function util_zdotc(n, x, incx, y, incy)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero

  implicit none
  integer(c_int), intent(in) :: n, incx, incy
  complex(c_double_complex), intent(in) :: x(1:n)
  complex(c_double_complex), intent(in) :: y(1:n)

  integer(c_int) :: i
  complex(c_double_complex) :: tmp

  if (incx /= 1) stop 'util_zdotc: bad incx.'
  if (incy /= 1) stop 'util_zdotc: bad incy.'

  tmp = czero
  do i = 1, n
     tmp = tmp + conjg(x(i)) * y(i)
  end do

  util_zdotc = tmp
  return

end function util_zdotc
!######################################################################
!######################################################################
complex(c_double_complex) function util_zdotcp(n, a, b)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero

  implicit none
  integer(c_int), intent(in) :: n
  complex(c_double_complex), intent(in) :: a(1:n)
  complex(c_double_complex), intent(in) :: b(1:n)

  integer(c_int) :: i
  complex(c_double_complex) :: tmp

  tmp = czero
!$omp parallel default(shared) reduction(+:tmp)
!$omp do
  do i = 1, n
     tmp = tmp + conjg(a(i)) * b(i)
  end do
!$omp end do
!$omp end parallel

  util_zdotcp = tmp
  return

end function util_zdotcp
!######################################################################
!######################################################################
complex(c_double_complex) function util_zdotu(n, x, incx, y, incy)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero

  implicit none
  integer(c_int), intent(in) :: n, incx, incy
  complex(c_double_complex), intent(in) :: x(1:n)
  complex(c_double_complex), intent(in) :: y(1:n)

  integer(c_int) :: i
  complex(c_double_complex) :: tmp

  if (incx /= 1) stop 'util_zdotu: bad incx.'
  if (incy /= 1) stop 'util_zdotu: bad incy.'

  tmp = czero
  do i = 1, n
     tmp = tmp + x(i) * y(i)
  end do

  util_zdotu = tmp
  return

end function util_zdotu
!######################################################################
!######################################################################
complex(c_double_complex) function util_zdotup(n, a, b)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero

  implicit none
  integer(c_int), intent(in) :: n
  complex(c_double_complex), intent(in) :: a(1:n)
  complex(c_double_complex), intent(in) :: b(1:n)

  integer(c_int) :: i
  complex(c_double_complex) :: tmp

  tmp = czero
!$omp parallel default(shared) reduction(+:tmp)
!$omp do
  do i = 1, n
     tmp = tmp + a(i) * b(i)
  end do
!$omp end do
!$omp end parallel

  util_zdotup = tmp
  return

end function util_zdotup
!######################################################################
