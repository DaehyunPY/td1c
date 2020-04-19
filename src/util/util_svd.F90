!################################################################################
subroutine util_svd_real(m, n, a, z, u, v)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_long), intent(in) :: m, n
  real(c_double), intent(in) :: a(1:m, 1:n)
  real(c_double), intent(out) :: z(1:*)
  real(c_double), intent(out) :: u(1:m, 1:m)
  real(c_double), intent(out) :: v(1:n, 1:n)

  integer(c_long) :: info, i, j, len
  integer(c_long), allocatable :: iwork(:)
  real(c_double), allocatable :: rwork(:)
  real(c_double), allocatable :: tmpa(:,:)
  real(c_double), allocatable :: vt(:,:)

  len = 8*min(m, n)
!debug
!  write(6, "('util_svd_real: len = ', i10)") len
!debug
  allocate(iwork(1:len))
  allocate(rwork(1:1))
  allocate(tmpa(1:m, 1:n))
  allocate(vt(1:n, 1:n))
  call dgesdd('A', m, n, tmpa, m, z, u, m, vt, n, rwork, -1, iwork, info)
  len = int(rwork(1))
!debug
!  write(6, "('util_svd_real: len = ', i10)") len
!debug
  deallocate(rwork)
  allocate(rwork(1:len))


  tmpa(1:m, 1:n) = a(1:m, 1:n)
  call dgesdd('A', m, n, tmpa, m, z, u, m, vt, n, rwork, len, iwork, info)

  if (info /= 0) then
     write(6, "('util_svd_real: info = ', i20)") info
     stop 'error in util_svd_real.'
  else
     do i = 1, n
        do j = 1, n
           v(j, i) = vt(i, j)
        end do
     end do
  end if

!debug  write(6, "('util_svd_real. singular values:')")
!debug  do i = 1, min(m, n)
!debug     write(6, "(i10, f20.10)") i, z(i)
!debug  end do

  deallocate(vt)
  deallocate(tmpa)
  deallocate(rwork)
  deallocate(iwork)

end subroutine util_svd_real
!################################################################################
