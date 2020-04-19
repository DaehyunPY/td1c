!######################################################################
subroutine bas_gen_zfac()

  use, intrinsic :: iso_c_binding
  use mod_sph, only : lmax1, mmax1
  use mod_rad, only : nrad, xrad, wrad
  use mod_bas, only : alph_lm, bas_zfac
  use mod_control, only : fedvr_normalized

  implicit none
  complex(c_double_complex) :: azfac, tmp
  integer(c_int) :: l, m, irad

  if (fedvr_normalized) then
     do m = -mmax1, mmax1
        do l = abs(m), lmax1
           azfac = alph_lm(l, m)
           do irad = 1, nrad - 1
              bas_zfac(irad, l, m) = xrad(irad) * azfac
           end do
        end do
     end do
  else
     do m = -mmax1, mmax1
        do l = abs(m), lmax1
           azfac = alph_lm(l, m)
           do irad = 1, nrad - 1
              bas_zfac(irad, l, m) = xrad(irad) * azfac * wrad(irad)
           end do
        end do
     end do
  end if

end subroutine bas_gen_zfac
!######################################################################
