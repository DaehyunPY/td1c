!################################################################################
complex(c_double_complex) function ormas_dets_single(ifun, jfun, nel, occi, occj, s0, work)

  use, intrinsic :: iso_c_binding
  use mod_const, only : one, czero, runit
  use mod_ormas, only : thrdet, ncore, nact

  implicit none
  !--------------------------------------------------------------------
  integer(c_int), intent(in) :: ifun, jfun, nel, occi(1:nact), occj(1:nact)
  complex(c_double_complex), intent(in) :: s0(1:nel, 1:nel)
  complex(c_double_complex), intent(out) :: work(1:nel-1, 1:nel-1)
  !--------------------------------------------------------------------
  integer(c_int) :: iel, jel, iorb, jorb, iact, jact
  complex(c_double_complex), external :: util_det

  ormas_dets_single = czero
  if (nel <= 1) return
  if (occi(ifun) == 0) return
  if (occj(jfun) == 0) return

  if (ifun <= ncore) then
     iorb = ifun
  else
     iorb = ncore
     do iact = 1, ifun - ncore
        iorb = iorb + occi(iact)
     end do
  end if
  if (jfun <= ncore) then
     jorb = jfun
  else
     jorb = ncore
     do jact = 1, jfun - ncore
        jorb = jorb + occj(jact)
     end do
  end if
  do iel = 1, iorb - 1
     do jel = 1, jorb - 1
        work(iel, jel) = s0(iel, jel)
     end do
  end do
  do iel = iorb + 1, nel
     do jel = 1, jorb - 1
        work(iel-1, jel) = s0(iel, jel)
     end do
  end do
  do iel = 1, iorb - 1
     do jel = jorb + 1, nel
        work(iel, jel-1) = s0(iel, jel)
     end do
  end do
  do iel = iorb + 1, nel
     do jel = jorb + 1, nel
        work(iel-1, jel-1) = s0(iel, jel)
     end do
  end do
  ormas_dets_single = util_det(nel-1, thrdet, work) * (-one) ** (iorb + jorb)
!debug
!  write(6, "('ormas_dets_single: passing through util_det.')")
!debug

  return

end function ormas_dets_single
!################################################################################
