!################################################################################
real(c_long_double) function util_fact(n)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_long), intent(in) :: n

  integer(c_long) i
  real(c_long_double) d4m, d4i

  d4m = 1
  do i = 1, n
     d4i = i
     d4m = d4m * d4i
  end do
  
  util_fact = d4m

end function util_fact
!################################################################################
