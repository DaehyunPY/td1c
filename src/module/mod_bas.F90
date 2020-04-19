!2015/10/22 Yuki Orimo Changed
!################################################################################
module mod_bas

  use, intrinsic :: iso_c_binding

  implicit none

  integer(c_int), pointer :: znuc, ltot, mtot, smul, nbas, nbas2, ngrid, nval(:), lval(:), mval(:), psp_label
  real(c_double), pointer :: grid(:,:), wgt(:)
  real(c_double), pointer :: alph_lm(:,:)
  complex(c_double_complex), pointer :: pmat(:,:)
  complex(c_double_complex), pointer :: kmat(:,:,:)
  complex(c_double_complex), pointer :: tmat(:,:,:)

  complex(c_double_complex), pointer :: bas_zfac(:,:,:)
  complex(c_double_complex), pointer :: bas_pzfac1(:,:)
  complex(c_double_complex), pointer :: bas_pzfac2(:,:,:)
  complex(c_double_complex), pointer :: bas_azfac(:,:,:)

  real(c_double), pointer :: d2ll(:,:,:)
  real(c_double), pointer :: bas_d2fac1(:)
  real(c_double), pointer :: bas_d2fac2(:,:)
  real(c_double), pointer :: bas_d2invr(:,:)
  real(c_double), pointer :: bas_d2rpl0(:,:)
  real(c_double), pointer :: bas_d2rpl1(:,:)

  real(c_double), save :: pp_znuc, pp_rloc, pp_cloc(1:4)
  real(c_double), save :: pp_rproj(0:3)
  real(c_double), save :: pp_hproj(1:3,1:3,0:3)
  integer(c_int), save :: pp_maxl
  integer(c_int), save :: pp_nump(0:3), pp_irmax(1:3,0:3)
  real(c_double), save, allocatable :: pp_vloc(:)
  real(c_double), save, allocatable :: pp_vlocHF(:,:)
  real(c_double), save, allocatable :: pp_pproj(:,:,:)
  real(c_double), save, allocatable :: pp_fproj(:,:,:)
  real(c_double), save, allocatable :: pp_gproj(:,:,:)
!nyi  real(c_double), save, allocatable :: pp_fprojx(:,:,:)
!nyi  real(c_double), save, allocatable :: pp_gprojx(:,:,:)

! Orimo_ECS
  complex(c_double_complex), pointer :: d2ll_ecs(:,:,:)
  integer(c_int), pointer :: ipiv_ecs(:,:)
  complex(c_double_complex), allocatable :: d2ll_herm(:,:,:)
  complex(c_double_complex), pointer :: d1mat(:,:)
  complex(c_double_complex), pointer :: d2mat(:,:)
  real(c_double), pointer :: confd2ll(:,:,:)
  complex(c_double_complex), pointer :: bas_d2crpl1(:,:)
! Orimo_ECS

end module mod_bas
!################################################################################
