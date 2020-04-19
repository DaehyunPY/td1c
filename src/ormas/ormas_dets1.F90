!################################################################################
complex(c_double_complex) function ormas_dets1(istr, jstr, ncore, nel, nela, orb, s0, work)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : thrdet
  use mod_const, only : czero, runit

  implicit none
  !--------------------------------------------------------------------
  integer(c_long), intent(in) :: istr, jstr, ncore, nel, nela, orb(0:nela, *)
  complex(c_double_complex), intent(in) :: s0(1:nel, 1:nel)
  complex(c_double_complex), intent(out) :: work(1:nel, 1:nel)
  !--------------------------------------------------------------------
  integer(c_long) :: iel, jel, iela, jela
  complex(c_double_complex), external :: util_det

  if (nel <= 0) then
     write(6, "('ormas_dets1: this is not fully tested...')")
     ormas_dets1 = czero
     return
  end if

  ormas_dets1 = czero
  do iel = 1, nel
     work(1:nel, 1:nel) = s0(1:nel, 1:nel)
     work(iel, 1:nel) = -s0(iel, 1:nel)
     if (iel <= ncore) then
        work(iel, iel) = work(iel, iel) + runit
     else
        iela = iel - ncore
        do jel = ncore + 1, nel
           jela = jel - ncore
!debug
!  write(6, "('ormas_dets1: ', 6i5)") istr, jstr, iela, jela, orb(iela, istr), orb(jela, jstr)
!debug
           if (orb(iela, istr) == orb(jela, jstr)) then
              work(iel, jel) = work(iel, jel) + runit
!debug
!  write(6, "('ormas_dets1: ', 6i5)") istr, jstr, iela, jela, orb(iela, istr), orb(jela, jstr)
!debug
              exit
           end if
        end do
     end if
     ormas_dets1 = ormas_dets1 + util_det(nel, thrdet, work)
!debug
!  write(6, "('ormas_dets1: passing through util_det.')")
!debug
  end do

  return

end function ormas_dets1
!################################################################################
