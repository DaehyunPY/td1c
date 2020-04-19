!################################################################################
subroutine util_trimat(n, thresh, npiv, a)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_long), intent(in) :: n
  real(c_double), intent(in) :: thresh
  integer(c_long), intent(out) :: npiv
  complex(c_double_complex), intent(inout) :: a(1:n, 1:n)

  integer(c_long) :: i, j
  complex(c_double_complex) :: diag, fac
  complex(c_double_complex), allocatable :: tcol(:)

  allocate(tcol(1:n))
  npiv = 0

  do i = 1, n
     diag = a(i, i)
     
     ! pivot
     if (abs(diag) < thresh) then
        npiv = npiv + 1
        do j = i + 1, n
           diag = a(i, j)
           if (abs(diag) > thresh) then
              tcol(1:n) = a(1:n, i)
              a(1:n, i) = a(1:n, j)
              a(1:n, j) = tcol(1:n)
              exit
           end if
        end do
     end if

!debug
!debug     write(6, "('debug: util_trimat ', i5, e12.5)") i, abs(diag)
!debug
     ! transform to upper trianguler matrix
     do j = i + 1, n
        fac = a(i, j) / diag
        a(1:n, j) = a(1:n, j) - a(1:n, i) * fac
     end do

  end do
  deallocate(tcol)

end subroutine util_trimat
!################################################################################
