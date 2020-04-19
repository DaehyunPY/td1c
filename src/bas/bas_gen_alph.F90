!######################################################################
subroutine bas_gen_alph()

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, xrad
  use mod_sph, only : lmax1, mmax1
  use mod_const, only : one, two, three
  use mod_bas, only : alph_lm

  implicit none
  integer(c_int) :: l, m, irad
  real(c_double) :: dl, dm, tmp

  do m = -mmax1, mmax1
     do l = abs(m), lmax1
        dl = dble(l)
        dm = dble(m)
        alph_lm(l, m) = sqrt(((dl + one)**2 - dm**2) / ((two * dl + one)*(two * dl + three)))
        !write(6, "('alph_lm: ', 2i5, f20.10)") l, m, alph_lm(l, m)
!        tmp = alph_lm(l, m) * (dl + one)
!        do irad = 1, nrad - 1
!           alph_rlm(irad, l, m) = tmp / xrad(irad)
!        end do
     end do
  end do

end subroutine bas_gen_alph
!######################################################################
