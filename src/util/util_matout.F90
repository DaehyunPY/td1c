!################################################################################
subroutine util_matoutc(n, mat)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_long), intent(in) :: n
  complex(c_double_complex), intent(in) :: mat(1:n, 1:n)
  !--------------------------------------------------------------------
  integer(c_long) :: i, j


!old  do i = 1, n
!old     do j = 1, n
!old        write(6, "(e12.4)", advance='no') dble(mat(j, i))
!old     end do
!old     write(6, "(' | ')", advance='no')
!old     do j = 1, n
!old        write(6, "(e12.4)", advance='no') aimag(mat(j, i))
!old     end do
!old     write(6, *)
!old  end do

  write(6, "(' mat-real')")
  do i = 1, n
     do j = 1, n
        write(6, "(e12.4)", advance='no') dble(mat(j, i))
     end do
     write(6, *)
  end do


  write(6, "(' mat-imag')")
  do i = 1, n
     do j = 1, n
        write(6, "(e12.4)", advance='no') aimag(mat(j, i))
     end do
     write(6, *)
  end do

end subroutine util_matoutc
!################################################################################
!################################################################################
subroutine util_matoutc2(n, m, mat)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_long), intent(in) :: n, m
  complex(c_double_complex), intent(in) :: mat(1:n, 1:m)
  !--------------------------------------------------------------------
  integer(c_long) :: i, j


  write(6, "(' mat-real')")
  do i = 1, m
     do j = 1, n
        write(6, "(e12.4)", advance='no') dble(mat(j, i))
     end do
     write(6, *)
  end do


  write(6, "(' mat-imag')")
  do i = 1, m
     do j = 1, n
        write(6, "(e12.4)", advance='no') aimag(mat(j, i))
     end do
     write(6, *)
  end do

end subroutine util_matoutc2
!################################################################################
!################################################################################
subroutine util_linoutc(n, mat)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_long), intent(in) :: n
  complex(c_double_complex), intent(in) :: mat(1:*)
  !--------------------------------------------------------------------
  integer(c_long) :: i, j, ji


  do i = 1, n
     do j = 1, i
        ji = (i * (i - 1)) / 2 + j
        write(6, "(e12.4)", advance='no') dble(mat(ji))
     end do
     do j = i + 1, n
        write(6, "(12x)", advance='no')
     end do
     write(6, "(' | ')", advance='no')
     do j = 1, i
        ji = (i * (i - 1)) / 2 + j
        write(6, "(e12.4)", advance='no') aimag(mat(ji))
     end do
     write(6, *)
  end do

end subroutine util_linoutc
!################################################################################
