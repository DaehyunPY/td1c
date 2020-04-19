!#######################################################################
integer(c_int) function ormas_sgn1x(nact, ifun, jfun, occvec)

  use, intrinsic :: iso_c_binding

  implicit none
  integer(c_int), intent(in) :: nact
  integer(c_int), intent(in) :: ifun, jfun
  integer(c_int), intent(in) :: occvec(1:nact)

  integer(c_int) :: kfun, eps

  eps = 0
  do kfun = min(ifun,jfun) + 1, max(ifun,jfun) - 1
     if (occvec(kfun) /= 0) eps = eps + 1
  end do

  ormas_sgn1x = (-1) ** eps

end function ormas_sgn1x
!#######################################################################
