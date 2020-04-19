!################################################################################
subroutine util_print_vec(n, vec, fname)

  use, intrinsic :: iso_c_binding

  implicit none
  character(*), intent(in) :: fname
  integer(c_long), intent(in) :: n
  complex(c_double_complex), intent(in) :: vec(1:n)
  integer(c_long) :: i

  open(unit=1, file=trim(fname), status='unknown', form='formatted')
  do i = 1, n
     write(1, "(i10,2f20.10)") i, vec(i)
  end do
  close(unit=1)

end subroutine util_print_vec
!################################################################################
