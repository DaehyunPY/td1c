!################################################################################
module mod_pes

  use, intrinsic :: iso_c_binding

  implicit none

  integer(c_int), pointer :: pes_numk, pes_llr, pes_ulr
  real(c_double), pointer :: pes_k_min, pes_k_max, pes_k_step
  real(c_double), pointer :: pes_bess(:,:,:)
  complex(c_double_complex), pointer :: pes_psik(:,:,:)
  complex(c_double_complex), pointer :: pes_rhok(:,:,:)

end module mod_pes
!################################################################################
