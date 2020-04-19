!#######################################################################
subroutine pes_bind(pes_numk_, pes_llr_, pes_ulr_, pes_k_min_, pes_k_max_, &
     pes_k_step_, pes_bess_, pes_psik_, pes_rhok_)

  use, intrinsic :: iso_c_binding
  use mod_pes
  use mod_rad, only : nrad
  use mod_sph, only : lmax1
  use mod_ormas, only : nfun

  implicit none
  integer(c_long), target, intent(in) :: pes_numk_, pes_llr_, pes_ulr_
  real(c_double), target, intent(in) :: pes_k_min_, pes_k_max_, pes_k_step_
  real(c_double), target, intent(in) :: pes_bess_(1:(nrad-1), 0:pes_numk_, 0:lmax1)
  complex(c_double_complex), target, intent(in) :: pes_psik_(0:pes_numk_, 0:lmax1, 1:nfun)
  complex(c_double_complex), target, intent(in) :: pes_rhok_(0:pes_numk_, 1:nfun, 1:nfun)

  pes_numk => pes_numk_
  pes_llr => pes_llr_
  pes_ulr => pes_ulr_
  pes_k_min => pes_k_min_
  pes_k_max => pes_k_max_
  pes_k_step => pes_k_step_
  pes_bess => pes_bess_
  pes_psik => pes_psik_
  pes_rhok => pes_rhok_

end subroutine pes_bind
!#######################################################################
