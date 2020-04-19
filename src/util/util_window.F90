!################################################################################
subroutine util_window1(n, i, j, a, wa)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_int), intent(in) :: n, i, j
  complex(c_double_complex), intent(in) :: a(1:n, 1:n)
  complex(c_double_complex), intent(out) :: wa(1:(n-1), 1:(n-1))

  integer(c_int) :: p, q

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
  integer(c_int), intent(in) :: n, i, j, k, l
  complex(c_double_complex), intent(in) :: a(1:n, 1:n)
  complex(c_double_complex), intent(out) :: wa(1:(n-2), 1:(n-2))

  integer(c_int) :: p, q

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
subroutine util_window3(n, i1, i2, i3, j1, j2, j3, a, wa)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_int), intent(in) :: n, i1,i2,i3, j1,j2,j3
  complex(c_double_complex), intent(in) :: a(1:n, 1:n)
  complex(c_double_complex), intent(out) :: wa(1:(n-3), 1:(n-3))

  integer(c_int) :: p, q

  do q = 1, j1 - 1
     do p = 1, i1 - 1
        wa(p, q) = a(p, q)
     end do
     do p = i1 + 1, i2 - 1
        wa(p-1, q) = a(p, q)
     end do
     do p = i2 + 1, i3 - 1
        wa(p-2, q) = a(p, q)
     end do
     do p = i3 + 1, n
        wa(p-3, q) = a(p, q)
     end do
  end do

  do q = j1 + 1, j2 - 1
     do p = 1, i1 - 1
        wa(p, q-1) = a(p, q)
     end do
     do p = i1 + 1, i2 - 1
        wa(p-1, q-1) = a(p, q)
     end do
     do p = i2 + 1, i3 - 1
        wa(p-2, q-1) = a(p, q)
     end do
     do p = i3 + 1, n
        wa(p-3, q-1) = a(p, q)
     end do
  end do

  do q = j2 + 1, j3 - 1
     do p = 1, i1 - 1
        wa(p, q-2) = a(p, q)
     end do
     do p = i1 + 1, i2 - 1
        wa(p-1, q-2) = a(p, q)
     end do
     do p = i2 + 1, i3 - 1
        wa(p-2, q-2) = a(p, q)
     end do
     do p = i3 + 1, n
        wa(p-3, q-2) = a(p, q)
     end do
  end do

  do q = j3 + 1, n
     do p = 1, i1 - 1
        wa(p, q-3) = a(p, q)
     end do
     do p = i1 + 1, i2 - 1
        wa(p-1, q-3) = a(p, q)
     end do
     do p = i2 + 1, i3 - 1
        wa(p-2, q-3) = a(p, q)
     end do
     do p = i3 + 1, n
        wa(p-3, q-3) = a(p, q)
     end do
  end do

end subroutine util_window3
!################################################################################
subroutine util_window4(n, i1, i2, i3, i4, j1, j2, j3, j4, a, wa)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_int), intent(in) :: n, i1,i2,i3,i4, j1,j2,j3,j4
  complex(c_double_complex), intent(in) :: a(1:n, 1:n)
  complex(c_double_complex), intent(out) :: wa(1:(n-4), 1:(n-4))

  integer(c_int) :: p, q

  do q = 1, j1 - 1
     do p = 1, i1 - 1
        wa(p, q) = a(p, q)
     end do
     do p = i1 + 1, i2 - 1
        wa(p-1, q) = a(p, q)
     end do
     do p = i2 + 1, i3 - 1
        wa(p-2, q) = a(p, q)
     end do
     do p = i3 + 1, i4 - 1
        wa(p-3, q) = a(p, q)
     end do
     do p = i4 + 1, n
        wa(p-4, q) = a(p, q)
     end do
  end do

  do q = j1 + 1, j2 - 1
     do p = 1, i1 - 1
        wa(p, q-1) = a(p, q)
     end do
     do p = i1 + 1, i2 - 1
        wa(p-1, q-1) = a(p, q)
     end do
     do p = i2 + 1, i3 - 1
        wa(p-2, q-1) = a(p, q)
     end do
     do p = i3 + 1, i4 - 1
        wa(p-3, q-1) = a(p, q)
     end do
     do p = i4 + 1, n
        wa(p-4, q-1) = a(p, q)
     end do
  end do

  do q = j2 + 1, j3 - 1
     do p = 1, i1 - 1
        wa(p, q-2) = a(p, q)
     end do
     do p = i1 + 1, i2 - 1
        wa(p-1, q-2) = a(p, q)
     end do
     do p = i2 + 1, i3 - 1
        wa(p-2, q-2) = a(p, q)
     end do
     do p = i3 + 1, i4 - 1
        wa(p-3, q-2) = a(p, q)
     end do
     do p = i4 + 1, n
        wa(p-4, q-2) = a(p, q)
     end do
  end do

  do q = j3 + 1, j4 - 1
     do p = 1, i1 - 1
        wa(p, q-3) = a(p, q)
     end do
     do p = i1 + 1, i2 - 1
        wa(p-1, q-3) = a(p, q)
     end do
     do p = i2 + 1, i3 - 1
        wa(p-2, q-3) = a(p, q)
     end do
     do p = i3 + 1, i4 - 1
        wa(p-3, q-3) = a(p, q)
     end do
     do p = i4 + 1, n
        wa(p-4, q-3) = a(p, q)
     end do
  end do

  do q = j4 + 1, n
     do p = 1, i1 - 1
        wa(p, q-4) = a(p, q)
     end do
     do p = i1 + 1, i2 - 1
        wa(p-1, q-4) = a(p, q)
     end do
     do p = i2 + 1, i3 - 1
        wa(p-2, q-4) = a(p, q)
     end do
     do p = i3 + 1, i4 - 1
        wa(p-3, q-4) = a(p, q)
     end do
     do p = i4 + 1, n
        wa(p-4, q-4) = a(p, q)
     end do
  end do

end subroutine util_window4
!################################################################################
subroutine util_window5(n, i1, i2, i3, i4, i5, j1, j2, j3, j4, j5, a, wa)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_int), intent(in) :: n, i1,i2,i3,i4,i5, j1,j2,j3,j4,j5
  complex(c_double_complex), intent(in) :: a(1:n, 1:n)
  complex(c_double_complex), intent(out) :: wa(1:(n-5), 1:(n-5))

  integer(c_int) :: p, q

  do q = 1, j1 - 1
     do p = 1, i1 - 1
        wa(p, q) = a(p, q)
     end do
     do p = i1 + 1, i2 - 1
        wa(p-1, q) = a(p, q)
     end do
     do p = i2 + 1, i3 - 1
        wa(p-2, q) = a(p, q)
     end do
     do p = i3 + 1, i4 - 1
        wa(p-3, q) = a(p, q)
     end do
     do p = i4 + 1, i5 - 1
        wa(p-4, q) = a(p, q)
     end do
     do p = i5 + 1, n
        wa(p-5, q) = a(p, q)
     end do
  end do

  do q = j1 + 1, j2 - 1
     do p = 1, i1 - 1
        wa(p, q-1) = a(p, q)
     end do
     do p = i1 + 1, i2 - 1
        wa(p-1, q-1) = a(p, q)
     end do
     do p = i2 + 1, i3 - 1
        wa(p-2, q-1) = a(p, q)
     end do
     do p = i3 + 1, i4 - 1
        wa(p-3, q-1) = a(p, q)
     end do
     do p = i4 + 1, i5 - 1
        wa(p-4, q-1) = a(p, q)
     end do
     do p = i5 + 1, n
        wa(p-5, q-1) = a(p, q)
     end do
  end do

  do q = j2 + 1, j3 - 1
     do p = 1, i1 - 1
        wa(p, q-2) = a(p, q)
     end do
     do p = i1 + 1, i2 - 1
        wa(p-1, q-2) = a(p, q)
     end do
     do p = i2 + 1, i3 - 1
        wa(p-2, q-2) = a(p, q)
     end do
     do p = i3 + 1, i4 - 1
        wa(p-3, q-2) = a(p, q)
     end do
     do p = i4 + 1, i5 - 1
        wa(p-4, q-2) = a(p, q)
     end do
     do p = i5 + 1, n
        wa(p-5, q-2) = a(p, q)
     end do
  end do

  do q = j3 + 1, j4 - 1
     do p = 1, i1 - 1
        wa(p, q-3) = a(p, q)
     end do
     do p = i1 + 1, i2 - 1
        wa(p-1, q-3) = a(p, q)
     end do
     do p = i2 + 1, i3 - 1
        wa(p-2, q-3) = a(p, q)
     end do
     do p = i3 + 1, i4 - 1
        wa(p-3, q-3) = a(p, q)
     end do
     do p = i4 + 1, i5 - 1
        wa(p-4, q-3) = a(p, q)
     end do
     do p = i5 + 1, n
        wa(p-5, q-3) = a(p, q)
     end do
  end do

  do q = j4 + 1, j5 - 1
     do p = 1, i1 - 1
        wa(p, q-4) = a(p, q)
     end do
     do p = i1 + 1, i2 - 1
        wa(p-1, q-4) = a(p, q)
     end do
     do p = i2 + 1, i3 - 1
        wa(p-2, q-4) = a(p, q)
     end do
     do p = i3 + 1, i4 - 1
        wa(p-3, q-4) = a(p, q)
     end do
     do p = i4 + 1, i5 - 1
        wa(p-4, q-4) = a(p, q)
     end do
     do p = i5 + 1, n
        wa(p-5, q-4) = a(p, q)
     end do
  end do

  do q = j5 + 1, n
     do p = 1, i1 - 1
        wa(p, q-5) = a(p, q)
     end do
     do p = i1 + 1, i2 - 1
        wa(p-1, q-5) = a(p, q)
     end do
     do p = i2 + 1, i3 - 1
        wa(p-2, q-5) = a(p, q)
     end do
     do p = i3 + 1, i4 - 1
        wa(p-3, q-5) = a(p, q)
     end do
     do p = i4 + 1, i5 - 1
        wa(p-4, q-5) = a(p, q)
     end do
     do p = i5 + 1, n
        wa(p-5, q-5) = a(p, q)
     end do
  end do

end subroutine util_window5
!################################################################################
