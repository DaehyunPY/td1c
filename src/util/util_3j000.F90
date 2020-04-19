!################################################################################
real(c_double) function util_3j000(l1, l2, l3)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_int), intent(in) :: l1, l2, l3

  integer(c_int) :: g, j, tmp
  real(c_long_double) :: w3j
  real(c_long_double), external :: util_fact

  util_3j000 = 0.d+0

  if (l1 < 0 .or. l2 < 0 .or. l3 < 0) return
  if (abs(l1 - l2) > l3 .or. l3 > l1 + l2) return

  j = l1 + l2 + l3
  if (mod(j, 2) .ne. 0) return

  g = j / 2
  tmp = j - 2 * l1
  w3j = util_fact(tmp)
  tmp = j - 2 * l2
  w3j = w3j * util_fact(tmp)
  tmp = j + 1
  w3j = w3j / util_fact(tmp)
  tmp = j - 2 * l3
  w3j = w3j * util_fact(tmp)
  
  w3j = sqrt(w3j)

  tmp = g - l1
  w3j = w3j / util_fact(tmp)
  tmp = g
  w3j = w3j * util_fact(tmp)
  tmp = g - l2
  w3j = w3j / util_fact(tmp)
  tmp = g - l3
  w3j = w3j / util_fact(tmp)
  
  w3j = w3j * (-1.d+0) ** g
  util_3j000 = dble(w3j)

end function util_3j000
!################################################################################
