!################################################################################
module mod_sph

  use, intrinsic :: iso_c_binding

  implicit none

  integer(c_long), pointer :: nsph, lmax1, lmax2, mmax1, mmax2
  integer(c_long), pointer :: nang, nlat, nphi
  real(c_double), pointer :: wlat(:), wang(:), wphi, cost(:), sint(:)
  real(c_double), pointer :: legf1(:,:,:), legb1(:,:,:)
  real(c_double), pointer :: legf2(:,:,:), legb2(:,:,:)
  real(c_double), allocatable :: sph_gaunt(:,:,:,:,:)

end module mod_sph
!################################################################################
