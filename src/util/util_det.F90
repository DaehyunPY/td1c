!################################################################################
complex(c_double_complex) function util_det(n, thresh, a)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_long), intent(in) :: n
  real(c_double), intent(in) :: thresh
  complex(c_double_complex), intent(in) :: a(1:n, 1:n)

  real(c_double), parameter :: one = 1.d+0
  complex(c_double_complex), parameter :: czero = (0.d+0, 0.d+0)
  complex(c_double_complex), parameter :: runit = (1.d+0, 0.d+0)

  logical :: zero_triv
  integer(c_long) :: npiv, i, j
  integer(c_long) :: nsys
  complex(c_double_complex) :: test
! complex(c_double_complex), external :: zdotu
  complex(c_double_complex), external :: util_zdotu
  complex(c_double_complex), allocatable :: tmpa(:,:), norm(:)

  nsys = n

  ! check the trivial zero determinant
  ! **********************************
  zero_triv = .false.

  allocate(norm(1:n))
  do i = 1, n
!    test = zdotu(nsys, a(1, i), 1, a(1, i), 1)
     test = util_zdotu(nsys, a(1, i), 1, a(1, i), 1)
     norm(i) = sqrt(test)
     if (abs(test) < thresh) then
        zero_triv = .true.
        exit
     end if
  end do

  do i = 1, n
     if (zero_triv) exit
     do j = 1, i - 1
!       test = zdotu(nsys, a(1, j), 1, a(1, i), 1) / (norm(j) * norm(i))
        test = util_zdotu(nsys, a(1, j), 1, a(1, i), 1) / (norm(j) * norm(i))
        if (abs(runit - test) < thresh .or. abs(-runit - test) < thresh) then
           zero_triv = .true.
           exit
        end if
     end do
  end do
  deallocate(norm)

  if (zero_triv) then
     util_det = czero
     return
  end if
  ! **********************************

  allocate(tmpa(1:n, 1:n))
  tmpa(1:n, 1:n) = a(1:n, 1:n)

!debug
!  write(6, "('before:')")
!  do i = 1, n
!     write(6, "(3f20.10)") real(tmpa(i, 1:n))
!  end do
!debug

  call util_trimat(n, thresh, npiv, tmpa)
  util_det = (-one) ** npiv
  do i = 1, n
     util_det = util_det * tmpa(i, i)
  end do

!debug
!  write(6, "('after:')")
!  do i = 1, n
!     write(6, "(3f20.10)") real(tmpa(i, 1:n))
!  end do
!debug

  deallocate(tmpa)

end function util_det
!################################################################################
