!################################################################################
complex(c_double_complex) function ormas_dett_1(istr, jstr, ncore, nact, nfun, n1x, p1x, h1x, &
     & eq1x, sgn1x, ovlp)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : thrdet
  use mod_const, only : czero, runit

  implicit none
  !--------------------------------------------------------------------
  integer(c_int), intent(in) :: istr, jstr, ncore, nact, nfun, n1x(-3:1, 1:*)
  integer(c_int), intent(in) :: p1x(1:nact*nact, 1:*)
  integer(c_int), intent(in) :: h1x(1:nact*nact, 1:*)
  integer(c_int), intent(in) :: eq1x(1:nact*nact, 1:*)
  integer(c_int), intent(in) :: sgn1x(1:nact*nact, 1:*)
  complex(c_double_complex), intent(in) :: ovlp(1:nfun, 1:nfun)
  !--------------------------------------------------------------------
  integer(c_int) :: i1x, ifun, jfun

  ormas_dett_1 = czero

!2018/07/13  if (istr == jstr) then
  do ifun = 1, ncore
     ormas_dett_1 = ormas_dett_1 + ovlp(ifun, ifun)
  end do
!2018/07/13  end if
  if (nact == 0) return

  do i1x = 1, n1x(0, istr)
     if (eq1x(i1x, istr) == jstr) then
        ifun = ncore + h1x(i1x, istr)
        jfun = ncore + p1x(i1x, istr)
        ormas_dett_1 = ormas_dett_1 + ovlp(ifun, jfun) * sgn1x(i1x, istr)
     end if
  end do  

end function ormas_dett_1
!################################################################################
