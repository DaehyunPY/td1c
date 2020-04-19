!################################################################################
subroutine ormas_sij(istr, jstr, ncore, nel, nela, nfun, orb, ovlp, sij)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : thrdet

  implicit none
  !--------------------------------------------------------------------
  integer(c_int), intent(in) :: istr, jstr, ncore, nel, nela, nfun, orb(0:nela, *)
  complex(c_double_complex), intent(in) :: ovlp(1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: sij(1:nel, 1:nel)
  !--------------------------------------------------------------------
  integer(c_int) :: iel, jel, iact, jact, ifun, jfun
  complex(c_double_complex), external :: util_det

  ! active - active
  do iel = 1, nela
     iact = orb(iel, istr)
     ifun = ncore + iact
     do jel = 1, nela
        jact = orb(jel, jstr)
        jfun = ncore + jact
        sij(ncore + iel, ncore + jel) = ovlp(ifun, jfun)
     end do
  end do

  ! active - core
  do iel = 1, nela
     iact = orb(iel, istr)
     ifun = ncore + iact
     do jel = 1, ncore
        sij(ncore + iel, jel) = ovlp(ifun, jel)
     end do
  end do

  ! core - active
  do iel = 1, ncore
     do jel = 1, nela
        jact = orb(jel, jstr)
        jfun = ncore + jact
        sij(iel, ncore + jel) = ovlp(iel, jfun)
     end do
  end do

  ! core - core
  do iel = 1, ncore
     do jel = 1, ncore
        sij(iel, jel) = ovlp(iel, jel)
     end do
  end do

end subroutine ormas_sij
!################################################################################
