!################################################################################
module mod_hprod

  use, intrinsic :: iso_c_binding

  implicit none

  real(c_double), pointer :: ene_fcore, ene_dcore, ene_core, ene_act, ene_tot
  real(c_double), pointer :: dip_exp, vel_exp, acc_exp
  logical(c_bool), pointer :: projhigh
  real(c_double), pointer :: projhigh_cutoff

  complex(c_double_complex), save, allocatable :: int1e(:,:)
  complex(c_double_complex), save, allocatable :: int2e(:,:,:,:)
  complex(c_double_complex), save, allocatable :: den1(:,:)
  complex(c_double_complex), save, allocatable :: rden(:,:)
  complex(c_double_complex), save, allocatable :: rrden(:,:)
  complex(c_double_complex), save, allocatable :: den2(:,:,:,:)
  complex(c_double_complex), save, allocatable :: den2r(:,:,:,:)

  complex(c_double_complex), save, allocatable :: orb(:,:,:)
  complex(c_double_complex), save, allocatable :: torb(:,:,:)
  complex(c_double_complex), save, allocatable :: h0orb(:,:,:)
  complex(c_double_complex), save, allocatable :: h1orb(:,:,:)
  complex(c_double_complex), save, allocatable :: v2orb(:,:,:)
  complex(c_double_complex), save, allocatable :: gorb(:,:,:)
  complex(c_double_complex), save, allocatable :: orbg(:,:,:)
  complex(c_double_complex), save, allocatable :: torbg(:,:,:)
  complex(c_double_complex), save, allocatable :: h0orbg(:,:,:)
  complex(c_double_complex), save, allocatable :: h1orbg(:,:,:)
  complex(c_double_complex), save, allocatable :: v2orbg(:,:,:)
  complex(c_double_complex), save, allocatable :: gorbg(:,:,:)
  complex(c_double_complex), save, allocatable :: orb0(:,:,:)
  complex(c_double_complex), save, allocatable :: rhofc(:)
  complex(c_double_complex), save, allocatable :: v2jfc(:)
  complex(c_double_complex), save, allocatable :: v2xfc(:)
  complex(c_double_complex), save, allocatable :: rho1(:,:)
  complex(c_double_complex), save, allocatable :: rho2(:,:,:,:)
  complex(c_double_complex), save, allocatable :: v2sph(:,:,:,:)
  complex(c_double_complex), save, allocatable :: v2ang(:,:,:,:)
  complex(c_double_complex), save, allocatable :: cic(:)
  complex(c_double_complex), save, allocatable :: dcic(:)
  complex(c_double_complex), save, allocatable :: cic0(:)
! ##### 3j selection rule #####
  complex(c_double_complex), save, allocatable :: orbe(:,:,:), gorbe(:,:,:), v2orbe(:,:,:), v2ange(:,:,:,:)
  complex(c_double_complex), save, allocatable :: orbo(:,:,:), gorbo(:,:,:), v2orbo(:,:,:), v2ango(:,:,:,:)
  complex(c_double_complex), save, allocatable :: rho23j(:,:,:,:), v2sph3j(:,:,:,:)
  complex(c_double_complex), save, allocatable :: torb3j(:,:,:,:)
! ##### 3j selection rule #####

  complex(c_double_complex), save, allocatable :: ovlp(:,:)
  complex(c_double_complex), save, allocatable :: fmat(:,:)
  complex(c_double_complex), save, allocatable :: bmat(:,:)
  complex(c_double_complex), save, allocatable :: xmat(:,:)

  integer(c_long), save :: projhigh_nfun
  integer(c_long), save, allocatable :: projhigh_ncut(:)
  complex(c_double_complex), save, allocatable :: projhigh_orbs(:,:,:)
  real(c_double), save, allocatable :: projhigh_eigs(:,:)

! Sato_ECS
  complex(c_double_complex), save, allocatable :: h0orb_out(:,:,:)
  complex(c_double_complex), save, allocatable :: h1orb_out(:,:,:)
! Sato_ECS

end module mod_hprod
!################################################################################
