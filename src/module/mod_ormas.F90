!################################################################################
module mod_ormas

  use, intrinsic :: iso_c_binding

  implicit none

  logical, parameter :: ormas_nod2x = .false.
  logical, parameter :: ormas_nocp = .false.
  logical, parameter :: ormas_regd = .true.
  logical, parameter :: ormas_sd1 = .false.
  logical, parameter :: ormas = .true.
  logical, parameter :: ormas0 = .false.
  logical, parameter :: rasscf = .false.
  logical, parameter :: qcas = .false.
  logical, parameter :: read_allowed = .false.
  logical, parameter :: s2zero = .false.
  logical, parameter :: ormas_phase = .true.

  integer(c_int), save :: iprint
  integer(c_int), parameter :: MAX_NBLOCK = 4
  integer(c_int), parameter :: MAX_NSUB = 10

!!!TEST
  integer(c_int), parameter :: hcic_type = 0
  integer(c_int), parameter :: den1_type = 0
  integer(c_int), parameter :: den2_type = 0
!  integer(c_int), parameter :: hcic_type = 2
!  integer(c_int), parameter :: den1_type = 2
!  integer(c_int), parameter :: den2_type = 4
!!!TEST
  integer(c_int), parameter :: ormas_pene = 1
  integer(c_int), parameter :: ormas_ncut = 0
  integer(c_int), parameter :: ormas_reg = 1

  real(c_double), parameter :: thrcic = 1.0D-15
  real(c_double), parameter :: thrdet = 1.0D-15
  real(c_double), parameter :: ormas_smin = 1.0D-10
  real(c_double), parameter :: ormas_damp = 1.0D-12
  real(c_double), parameter :: ormas_dsv =  1.0D-15
  real(c_double), parameter :: ormas_dyreg0 = 1.0D-10
  real(c_double), parameter :: ormas_dyreg1 = 1.0D-10
  real(c_double), parameter :: ormas_dyreg2 = 1.0D-8
  real(c_double), parameter :: ormas_maxang = 1.0D-10
  real(c_double), save, pointer :: thradfc
  logical(c_bool), save, pointer :: ormas_donly

  logical(c_bool), save, pointer :: dplus
  integer(c_int), save, pointer :: nact1,act1_ll,act1_ul
  integer(c_int) :: dplus_type ! 0: H2.H1|P1.P2, 1: H1|P1.P2, 2: H2.H1|P1

  integer(c_int), save, pointer :: nblock
  integer(c_int), save, pointer :: type_block(:)
  integer(c_int), save, pointer :: nfun_block(:)
  integer(c_int), save, pointer :: nfcore2, nfcore1, nfcore, ndcore, ncore, nact, nocc, nvir, nfun, mval(:), froz(:)
  integer(c_int), save, pointer :: nelcore(:), nelact(:), neltot(:)
  integer(c_int), save, pointer :: nsub, norb_sub(:), lorb_sub(:,:)
  integer(c_int), save, pointer :: min_sub(:), max_sub(:)
  integer(c_int), save, pointer :: nstr_alph, nstr_beta
  integer(c_int), save, pointer :: ndet
! tdcis-teramura
  integer(c_int), save, pointer :: nfcore_tdcis
! tdcis-teramura

  integer(c_int), save :: min_sub_alph(1:MAX_NSUB), min_sub_beta(1:MAX_NSUB)
  integer(c_int), save :: max_sub_alph(1:MAX_NSUB), max_sub_beta(1:MAX_NSUB)
  integer(c_int), save :: ndist, ndist_alph, ndist_beta
  integer(c_int), save, allocatable :: sub_orb(:)
  integer(c_int), save, allocatable :: dist(:,:), dist_alph(:,:), dist_beta(:,:)
  integer(c_int), save, allocatable :: det_allowed(:,:), ndet_dist(:), mapf_det(:,:), mapr_det(:,:)

!FOR M-ADAPTATION
  integer(c_int), save, pointer :: lcic
  integer(c_int), save, pointer :: ndetx
  integer(c_int), save :: mmin_alph, mmax_alph
  integer(c_int), save :: mmin_beta, mmax_beta
  integer(c_int), save, allocatable :: mval_alph(:), nstr_m_alph(:), llstr_m_alph(:)
  integer(c_int), save, allocatable :: mval_beta(:), nstr_m_beta(:), llstr_m_beta(:)
  integer(c_int), save, allocatable :: nstr_dist_m_alph(:,:), llstr_dist_m_alph(:,:)
  integer(c_int), save, allocatable :: nstr_dist_m_beta(:,:), llstr_dist_m_beta(:,:)

  integer(c_int), allocatable :: nstr_alph_beta(:), llstr_alph_beta(:), ntot_alph_beta(:)
  integer(c_int), allocatable :: nstr_beta_alph(:), llstr_beta_alph(:), ntot_beta_alph(:)
  integer(c_int), allocatable :: map2to1_alph(:), map1to2_alph(:)
  integer(c_int), allocatable :: map2to1_beta(:), map1to2_beta(:)
  integer(c_int), save, allocatable :: mapr_detx(:,:)
!FOR M-ADAPTATION

  integer(c_int), save, allocatable :: nstr_alph_dist(:), lstr_alph_dist(:,:), nstr_alph_dist_sub(:,:)
  integer(c_int), save, allocatable :: nstr_beta_dist(:), lstr_beta_dist(:,:), nstr_beta_dist_sub(:,:)
  integer(c_int), save, allocatable :: arc_alph(:,:)
  integer(c_int), save, allocatable :: arc_beta(:,:)
  integer(c_int), save, allocatable :: dist_str_alph(:,:), substr_alph(:,:)
  integer(c_int), save, allocatable :: dist_str_beta(:,:), substr_beta(:,:)
  integer(c_int), save, allocatable :: onv_alph(:,:), orb_alph(:,:)
  integer(c_int), save, allocatable :: onv_beta(:,:), orb_beta(:,:)
  integer(c_int), save, allocatable :: n1x_alph(:,:), p1x_alph(:,:), h1x_alph(:,:), eq1x_alph(:,:), sgn1x_alph(:,:)
!OLD  integer(c_int), save, allocatable :: n1xr_alph(:,:), r1xr_alph(:,:,:), l1xr_alph(:,:,:), sgn1xr_alph(:,:,:)
  integer(c_int), save, allocatable :: n1x_beta(:,:), p1x_beta(:,:), h1x_beta(:,:), eq1x_beta(:,:), sgn1x_beta(:,:)
!OLD  integer(c_int), save, allocatable :: n1xr_beta(:,:), r1xr_beta(:,:,:), l1xr_beta(:,:,:), sgn1xr_beta(:,:,:)
  integer(c_int), save :: nrotoo, nrotca, nrotaa
  integer(c_int), save, allocatable :: rotoo_mapf(:,:), rotoo_mapb(:,:)
  integer(c_int), save, allocatable :: rotca_mapf(:,:), rotca_mapb(:,:)
  integer(c_int), save, allocatable :: rotaa_mapf(:,:), rotaa_mapb(:,:)
  integer(c_int), save, allocatable :: nfunpp(:), pp2fun(:), fun2pp(:)

  ! M adaptation
  integer(c_int), save :: max1x_alph, max1x_m_alph
  integer(c_int), save :: max1x_beta, max1x_m_beta
  integer(c_int), save, allocatable :: n1x_m_alph(:,:), map1x_m_alph(:,:,:)
  integer(c_int), save, allocatable :: n1x_m_beta(:,:), map1x_m_beta(:,:,:)

  !  integer(c_int), save :: nrotov
  !  integer(c_int), save, allocatable :: rotov_mapf(:,:), rotov_mapb(:,:)

  real(c_double), save :: max_grad(1:3)
  real(c_double), save :: rms_grad(1:3)

  ! tdcc flag
  logical(c_bool), save, pointer :: tdcc
  logical(c_bool), save, pointer :: den2_abonly

  ! arrays for mkden2 version 4 by Fabian  
  logical(c_bool), save, pointer :: fab_den2
  integer(c_int), allocatable :: fab_eq1x(:,:),fab_sgn1x(:,:),fab_h1x(:,:),fab_p1x(:,:)
  integer(c_int), allocatable :: fab_nr1x(:)
  integer(c_int) :: fab_n1x

  logical(c_bool), save :: cic_old

end module mod_ormas
!################################################################################
