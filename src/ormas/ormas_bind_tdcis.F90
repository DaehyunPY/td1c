!################################################################################
subroutine ormas_bind_tdcis(nfcore_tdcis_)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nfcore_tdcis

  implicit none
  integer(c_int), target, intent(in) :: nfcore_tdcis_

  nfcore_tdcis =>   nfcore_tdcis_

end subroutine ormas_bind_tdcis
!################################################################################
