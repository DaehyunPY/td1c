!################################################################################
integer(c_long) function util_bicoeff(n, k)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_long), intent(in) :: n, k
  integer(c_long), external :: util_ifact
  integer(c_long) :: kx, num, denom

!old  nk = n - k
!old  num = util_ifact(n)
!old  denom1 = util_ifact(nk)
!old  denom2 = util_ifact(k)
!old  util_bicoeff = num / (denom1 * denom2)
  kx = min(k, n - k)
  num = util_ifact(n, kx)
  denom = util_ifact(kx, kx)
  util_bicoeff = num / denom

!debugwrite(6, "('util_bicoeff: n = ', i5)"), n
!debugwrite(6, "('util_bicoeff: k = ', i5)"), kx
!debugwrite(6, "('util_bicoeff: num = ', i20)"), num
!debugwrite(6, "('util_bicoeff: de1 = ', i20)"), denom

end function util_bicoeff
!################################################################################
