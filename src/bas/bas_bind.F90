!#######################################################################
subroutine bas_rad_bind(nfe_, ndvr_, nrad_, nradfc_, mapf_, mapb_, xrad_, wrad_, &
     radp_, radk_, radk0_, rmask_, mask_, ecs_flag_, theta_, recs_, cwrad_, &
     radkI_ecs_, irad_ecs_, cxrad_, bra_wrad_, &
     rdr_, wdw_, type_mkint1_sph_, type_mkint2_sph_, type_mkv2mf_, type_mkxmat_aa_, &
     inf_range_, irad_inf_, exp_factor_, switchoff_, irad_sw_, pot_type_, trunc_irad_, &
     nradgs_)

  use, intrinsic :: iso_c_binding
  use mod_rad

  implicit none
  integer(c_int), target, intent(in) :: nfe_, ndvr_, nrad_, nradfc_, mapf_(0:(nfe_-1)), mapb_(0:nrad_)
  real(c_double), target, intent(in) :: xrad_(0:nrad_), wrad_(0:nrad_)
  real(c_double), target, intent(in) :: radp_(1:(2*ndvr_+1), 0:nrad_)
  real(c_double), target, intent(in) :: radk_(1:(2*ndvr_+1), 0:nrad_)
  real(c_double), target, intent(in) :: radk0_(1:(2*ndvr_+1), 0:nrad_)
  real(c_double), target, intent(in) :: rmask_, mask_(0:nrad_)
! Orimo_ECS
  integer(c_int), target, intent(in) :: ecs_flag_, irad_ecs_, irad_inf_
  real(c_double), target, intent(in) :: exp_factor_
  real(c_double), target, intent(in) :: theta_, recs_, bra_wrad_(0:nrad_)
  complex(c_double_complex), target, intent(in) :: cwrad_(0:nrad_), radkI_ecs_
  complex(c_double_complex), target, intent(in) :: cxrad_(0:nrad_), rdr_(0:nrad_), wdw_(0:nrad_)
  integer(c_int), target, intent(in) :: type_mkint1_sph_, type_mkint2_sph_, type_mkv2mf_ 
  integer(c_int), target, intent(in) :: type_mkxmat_aa_, switchoff_, irad_sw_, pot_type_
  logical(c_bool), target, intent(in) :: inf_range_
  integer(c_int), target, intent(in) :: trunc_irad_
! Orimo_ECS
! tdcis-teramura
  integer(c_int), target, intent(in) :: nradgs_
! tdcis-teramura

  nfe => nfe_
  ndvr => ndvr_
  nrad => nrad_
  nradfc => nradfc_
  mapf => mapf_
  mapb => mapb_
  xrad => xrad_
  wrad => wrad_
  radk => radk_
  radk0 => radk0_
  radp => radp_
  rmask => rmask_
  mask => mask_
! Orimo_ECS
  ecs_flag => ecs_flag_
  theta => theta_
  recs => recs_
  cwrad => cwrad_
  radkI_ecs => radkI_ecs_
  irad_ecs => irad_ecs_
  cxrad => cxrad_
  bra_wrad => bra_wrad_
  rdr => rdr_
  wdw => wdw_
  type_mkint1_sph => type_mkint1_sph_
  type_mkint2_sph => type_mkint2_sph_
  type_mkv2mf => type_mkv2mf_
  type_mkxmat_aa => type_mkxmat_aa_
  inf_range => inf_range_
  irad_inf => irad_inf_
  exp_factor => exp_factor_
  switchoff => switchoff_
  irad_sw => irad_sw_
  pot_type => pot_type_
  trunc_irad => trunc_irad_
! Orimo_ECS
! tdcis-teramura
  nradgs => nradgs_
! tdcis-teramura

end subroutine bas_rad_bind
!#######################################################################
subroutine bas_sph_bind(lmax1_, lmax2_, mmax1_, mmax2_, nsph_, nlat_, nphi_, nang_, &
     wlat_, wphi_, wang_, cost_, sint_, legf1_, legb1_, legf2_, legb2_, &
     lmax1gs_, lmax2gs_ )

  use, intrinsic :: iso_c_binding
  use mod_sph

  implicit none
  integer(c_int), target, intent(in) :: lmax1_, lmax2_, mmax1_, mmax2_, nsph_, nlat_, nphi_, nang_
  real(c_double), target, intent(in) :: wlat_(1:nlat_), wphi_, wang_(1:nang_), cost_(1:nlat_), sint_(1:nlat_)
  real(c_double), target, intent(in) :: legf1_(0:lmax1_, 1:nlat_, -mmax1_:mmax1_)
  real(c_double), target, intent(in) :: legb1_(1:nlat_, 0:lmax1_, -mmax1_:mmax1_)
  real(c_double), target, intent(in) :: legf2_(0:lmax2_, 1:nlat_, -mmax2_:mmax2_)
  real(c_double), target, intent(in) :: legb2_(1:nlat_, 0:lmax2_, -mmax2_:mmax2_)
! tdcis-teramura
  integer(c_int), target, intent(in) :: lmax1gs_, lmax2gs_
! tdcis-teramura

  lmax1 => lmax1_
  lmax2 => lmax2_
  mmax1 => mmax1_
  mmax2 => mmax2_
  nsph => nsph_
  nlat => nlat_
  nphi => nphi_
  nang => nang_
  wlat => wlat_
  wphi => wphi_
  wang => wang_
  cost => cost_
  sint => sint_
  legf1 => legf1_
  legb1 => legb1_
  legf2 => legf2_
  legb2 => legb2_
! tdcis-teramura
  lmax1gs => lmax1gs_
  lmax2gs => lmax2gs_
! tdcis-teramura

end subroutine bas_sph_bind
!#######################################################################
subroutine bas_bas_bind(znuc_, smul_, ltot_, mtot_, nbas_, nbas2_, ngrid_, &
     nfun_, nval_, lval_, mval_, grid_, wgt_, alph_lm_, pmat_, kmat_, tmat_, d2ll_, &
     bas_zfac_, bas_pzfac1_, bas_pzfac2_, bas_azfac_, bas_d2fac1_, bas_d2fac2_, bas_d2invr_, &
     bas_d2rpl0_, bas_d2rpl1_, d2ll_ecs_, ipiv_ecs_, d1mat_, d2mat_, confd2ll_, bas_d2crpl1_, &     
     psp_label_)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : lmax1, lmax2, mmax1
  use mod_rad, only : ndvr, nrad, irad_ecs
  use mod_bas

  implicit none
  integer(c_int), target, intent(in) :: znuc_, ltot_, mtot_, smul_, nbas_, &
       nbas2_, ngrid_, nfun_, nval_(1:nfun_), lval_(1:nfun_), mval_(1:nfun_)
  real(c_double), target, intent(in) :: grid_(1:ngrid_, 1:4), wgt_(1:ngrid_)
  real(c_double), target, intent(in) :: alph_lm_(0:lmax1, -mmax1:mmax1)
  complex(c_double_complex), target, intent(in) :: pmat_(1:(2*ndvr+1), 1:(nrad-1))
  complex(c_double_complex), target, intent(in) :: kmat_(1:(2*ndvr+1), 1:(nrad-1), 0:lmax1)
  complex(c_double_complex), target, intent(in) :: tmat_(1:(2*ndvr+1), 1:(nrad-1), 0:lmax1)
  complex(c_double_complex), target, intent(in) :: bas_zfac_  (1:(nrad-1), 0:lmax1, -mmax1:mmax1)
  complex(c_double_complex), target, intent(in) :: bas_pzfac1_(            0:lmax1, -mmax1:mmax1)
  complex(c_double_complex), target, intent(in) :: bas_pzfac2_(1:(nrad-1), 0:lmax1, -mmax1:mmax1)
  complex(c_double_complex), target, intent(in) :: bas_azfac_ (1:(nrad-1), 0:lmax1, -mmax1:mmax1)
  real(c_double), target, intent(in) :: d2ll_(1:(ndvr+1), 1:(nrad-1), 0:lmax2)
  real(c_double), target, intent(in) :: bas_d2fac1_(            0:lmax2) ! 1 / R^{2l+1}
  real(c_double), target, intent(in) :: bas_d2fac2_(1:(nrad-1), 0:lmax2) ! 4\pi / (2l+1) / r
  real(c_double), target, intent(in) :: bas_d2invr_(1:(nrad-1), 0:lmax2) ! (2l+1) / r
  real(c_double), target, intent(in) :: bas_d2rpl0_(1:(nrad-1), 0:lmax2) ! r^l
  real(c_double), target, intent(in) :: bas_d2rpl1_(1:(nrad-1), 0:lmax2) ! r^{l+1}
! Orimo_ECS
  complex(c_double_complex), target, intent(in) :: d2ll_ecs_(1:(3*ndvr+1), 1:(nrad-1), 0:lmax2)
  integer(c_int), target, intent(in) :: ipiv_ecs_(1:(nrad-1), 0:lmax2)
  complex(c_double_complex), target, intent(in) :: d1mat_(1:(2*ndvr+1), 1:(nrad-1))
  complex(c_double_complex), target, intent(in) :: d2mat_(1:(2*ndvr+1), 1:(nrad-1))
  real(c_double), target, intent(in) :: confd2ll_(1:(  ndvr+1), 1:(irad_ecs-1), 0:lmax2)
  complex(c_double_complex), target, intent(in) :: bas_d2crpl1_(1:(nrad-1), 0:lmax2) ! R(r)^{l+1} * 4 \pi / (2*l - 1)
! Orimo_ECS
  integer(c_int), target, intent(in) :: psp_label_

  znuc => znuc_
  smul => smul_
  ltot => ltot_
  mtot => mtot_
  nbas => nbas_
  nbas2 => nbas2_
  ngrid => ngrid_
  nval => nval_
  lval => lval_
  mval => mval_
  grid => grid_
  wgt => wgt_
  alph_lm => alph_lm_
  pmat => pmat_
  kmat => kmat_
  tmat => tmat_
  d2ll => d2ll_
  bas_zfac   => bas_zfac_
  bas_pzfac1 => bas_pzfac1_
  bas_pzfac2 => bas_pzfac2_
  bas_azfac  => bas_azfac_
  bas_d2fac1 => bas_d2fac1_
  bas_d2fac2 => bas_d2fac2_
  bas_d2invr => bas_d2invr_
  bas_d2rpl0 => bas_d2rpl0_
  bas_d2rpl1 => bas_d2rpl1_
! Orimo_ECS
  d2ll_ecs => d2ll_ecs_
  ipiv_ecs => ipiv_ecs_
  d1mat => d1mat_
  d2mat => d2mat_
  confd2ll => confd2ll_
  bas_d2crpl1 => bas_d2crpl1_
! Orimo_ECS
  psp_label => psp_label_

end subroutine bas_bas_bind
!#######################################################################
