!######################################################################
subroutine test_zaxpy_omp()

  use, intrinsic :: iso_c_binding

  implicit none

  integer(c_int), parameter :: n = 100
  complex(c_double_complex), parameter :: a = 1d0
  complex(c_double_complex) :: x(1:n)
  complex(c_double_complex) :: y(1:n)
  integer(c_int) :: i

  integer(c_int), external :: omp_get_num_threads
  integer(c_int), external :: omp_get_thread_num

  !$omp parallel default(shared) private(i)
  i = omp_get_num_threads()
  !$omp end parallel

!  write(6,"('hello!':,i5)") 
  write(6,"('hello!')")

!  x = 1d0
!  y = 0d0
!  call zaxpy_omp(n,a,x,y)
!  do i = 1, n
!     write(6,"('hello!')")
!     write(6,"(i5)") i
!     write(6,"(4f10.5)") x(i),y(i)
!  end do

end subroutine test_zaxpy_omp
!######################################################################
subroutine zaxpy_omp(n, a, x, y)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_int), intent(in) :: n
  complex(c_double_complex), intent(in) :: a
  complex(c_double_complex), intent(in) :: x(1:n)
  complex(c_double_complex), intent(inout) :: y(1:n)

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
     y(i) = y(i) + a * x(i)
  end do
  !$omp end do
  !$omp end parallel

end subroutine zaxpy_omp
!######################################################################
!subroutine zaxpy_omp(n, a, x, y)
!
!  implicit none
!  integer(c_int), intent(in) :: n
!  complex(kind(0d0)), intent(in) :: a
!  complex(kind(0d0)), intent(in) :: x(1:n)
!  complex(kind(0d0)), intent(inout) :: y(1:n)
!
!  integer(c_int) :: i
!
!  !$omp parallel default(shared) private(i)
!  !$omp do
!  do i = 1, n
!     y(i) = y(i) + a * x(i)
!  end do
!  !$omp end do
!  !$omp end parallel
!
!end subroutine zaxpy_omp
!######################################################################
