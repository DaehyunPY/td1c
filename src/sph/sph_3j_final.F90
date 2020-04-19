!#######################################################################
subroutine sph_3j_final()
  use, intrinsic :: iso_c_binding
  use mod_sph, only : sph_gaunt
  implicit none
  deallocate(sph_gaunt)
end subroutine sph_3j_final
!#######################################################################
