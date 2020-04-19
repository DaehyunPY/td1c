!################################################################################
integer(c_int) function util_ifact(n, k)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_int), intent(in) :: n, k
  integer(c_int) i, m

  m = 1
  do i = n, n - k + 1, -1
     m = m * i
  end do
  
  util_ifact = m

end function util_ifact
!################################################################################
