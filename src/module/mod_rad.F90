!2015/10/22 Yuki Orimo Changed
!################################################################################
module mod_rad

  use, intrinsic :: iso_c_binding

  implicit none

  integer(c_long), pointer :: nfe, ndvr, nrad, nradfc
  integer(c_long), pointer :: mapf(:), mapb(:)
  real(c_double), pointer :: xrad(:), wrad(:)
  real(c_double), pointer :: radk(:,:), radk0(:,:)
  real(c_double), pointer :: radp(:,:)
  real(c_double), pointer :: rmask, mask(:)

  integer(c_long), pointer :: ecs_flag, irad_ecs, irad_inf
  real(c_double), pointer :: theta, recs, bra_wrad(:) 
  complex(c_double_complex), pointer :: cwrad(:), radkI_ecs, cxrad(:), rdr(:), wdw(:)
  integer(c_long), pointer :: type_mkint1_sph, type_mkint2_sph, type_mkv2mf, type_mkxmat_aa, switchoff, irad_sw, pot_type
  logical(c_bool), pointer :: inf_range
  real(c_double), pointer :: exp_factor
  integer(c_long), pointer :: trunc_irad

end module mod_rad
!################################################################################
