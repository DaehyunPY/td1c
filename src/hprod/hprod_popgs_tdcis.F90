!######################################################################
subroutine hprod_popgs_tdcis(chi)
  !population of the GS within rVG
  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, xrad
  use mod_sph, only : lmax1
  use mod_ormas, only : nfun, thrdet
  use mod_hprod, only : orb, orb0, orb0rot
  use mod_const, only : runit
  use mod_control, only : tdcis_rvg

  implicit none
  complex(c_double_complex), intent(in) :: chi(1:(nrad - 1), 0:lmax1, 1:nfun)
  complex(c_double_complex) :: orb0chi(1:nfun, 1:nfun)
  complex(c_double_complex) :: orb0chi_i(1:nfun, 1:nfun)
  complex(c_double_complex) :: det1
  complex(c_double_complex) :: deti
  complex(c_double_complex), external :: util_det

  complex(c_double_complex) :: pop
  integer(c_int) :: ifun, jfun

  ifun = 1
  ifun = 1
  pop = 0d0
  call hprod_mkovlp(xrad(nrad), orb0, orb0rot, orb0chi)
  call hprod_mkovlp(xrad(nrad), orb0, chi, orb0chi)

  orb0chi_i = orb0chi
  det1 = util_det(nfun, thrdet, orb0chi)
  do ifun = 1, nfun
     do jfun = 1, nfun
        orb0chi_i(ifun, jfun) = orb0chi(ifun, jfun)
     end do
     deti = util_det(nfun, thrdet, orb0chi_i)
     orb0chi_i = orb0chi
  end do

  write(6, "('hprod_popgs_tdcis:   ', 4f20.10)") det1, deti

end subroutine hprod_popgs_tdcis
!######################################################################
