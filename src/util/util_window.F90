!################################################################################
subroutine util_window1(n, i, j, a, wa)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_long), intent(in) :: n, i, j
  complex(c_double_complex), intent(in) :: a(1:n, 1:n)
  complex(c_double_complex), intent(out) :: wa(1:(n-1), 1:(n-1))

  integer(c_long) :: p, q

  do q = 1, j - 1
     do p = 1, i - 1
        wa(p, q) = a(p, q)
     end do
     do p = i + 1, n
        wa(p-1, q) = a(p, q)
     end do
  end do

  do q = j + 1, n
     do p = 1, i - 1
        wa(p, q-1) = a(p, q)
     end do
     do p = i + 1, n
        wa(p-1, q-1) = a(p, q)
     end do
  end do

!bug  do q = 1, j - 1
!bug     do p = 1, i - 1
!bug        wa(p, q) = a(p, q)
!bug     end do
!bug     do p = i, n - 1
!bug        wa(p, q) = a(p+1, q)
!bug     end do
!bug  end do
!bug
!bug  do q = j, n - 1
!bug     do p = 1, i - 1
!bug        wa(p, q) = a(p, q+1)
!bug     end do
!bug     do p = i, n - 1
!bug        wa(p, q) = a(p+1, q+1)
!bug     end do
!bug  end do

end subroutine util_window1
!################################################################################
!################################################################################
subroutine util_window2(n, i, j, k, l, a, wa)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_long), intent(in) :: n, i, j, k, l
  complex(c_double_complex), intent(in) :: a(1:n, 1:n)
  complex(c_double_complex), intent(out) :: wa(1:(n-2), 1:(n-2))

  integer(c_long) :: p, q

  do q = 1, k - 1
     do p = 1, i - 1
        wa(p, q) = a(p, q)
     end do
     do p = i + 1, j - 1
        wa(p-1, q) = a(p, q)
     end do
     do p = j + 1, n
        wa(p-2, q) = a(p, q)
     end do
  end do

  do q = k + 1, l - 1
     do p = 1, i - 1
        wa(p, q-1) = a(p, q)
     end do
     do p = i + 1, j - 1
        wa(p-1, q-1) = a(p, q)
     end do
     do p = j + 1, n
        wa(p-2, q-1) = a(p, q)
     end do
  end do

  do q = l + 1, n
     do p = 1, i - 1
        wa(p, q-2) = a(p, q)
     end do
     do p = i + 1, j - 1
        wa(p-1, q-2) = a(p, q)
     end do
     do p = j + 1, n
        wa(p-2, q-2) = a(p, q)
     end do
  end do

!bug  do q = 1, k - 1
!bug     do p = 1, i - 1
!bug        wa(p, q) = a(p, q)
!bug     end do
!bug     do p = i, j - 1
!bug        wa(p, q) = a(p+1, q)
!bug     end do
!bug     do p = j, n - 2
!bug        wa(p, q) = a(p+2, q)
!bug     end do
!bug  end do
!bug
!bug  do q = k, l - 1
!bug     do p = 1, i - 1
!bug        wa(p, q) = a(p, q+1)
!bug     end do
!bug     do p = i, j - 1
!bug        wa(p, q) = a(p+1, q+1)
!bug     end do
!bug     do p = j, n - 2
!bug        wa(p, q) = a(p+2, q+1)
!bug     end do
!bug  end do
!bug
!bug  do q = l, n - 2
!bug     do p = 1, i - 1
!bug        wa(p, q) = a(p, q+2)
!bug     end do
!bug     do p = i, j - 1
!bug        wa(p, q) = a(p+1, q+2)
!bug     end do
!bug     do p = j, n - 2
!bug        wa(p, q) = a(p+2, q+2)
!bug     end do
!bug  end do

end subroutine util_window2
!################################################################################
