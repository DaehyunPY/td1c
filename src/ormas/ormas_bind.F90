!################################################################################
subroutine ormas_bind(nblock_, type_block_, nfun_block_, nfcore2_, nfcore1_, &
     nfcore_, ndcore_, ncore_, nact_, nocc_, nvir_, nfun_, mval_, froz_, nelcore_, nelact_, &
     neltot_, nsub_, norb_sub_, lorb_sub_, min_sub_, max_sub_, nstr_alph_, nstr_beta_, &
     ndet_, lcic_, ndetx_, thradfc_, fab_den2_, donly_, den2_abonly_, dplus_, &
     nact1_, act1_ll_, act1_ul_, tdcc_, cc_type_, cc_code_)
     

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nblock, type_block, nfun_block
  use mod_ormas, only : nfcore2, nfcore1, nfcore, ndcore, ncore, nact, &
       nocc, nvir, nfun, mval, froz
  use mod_ormas, only : nelcore, nelact, neltot
  use mod_ormas, only : nsub, norb_sub, lorb_sub, min_sub, max_sub
  use mod_ormas, only : nstr_alph, nstr_beta, ndet, lcic, ndetx, thradfc
  use mod_ormas, only : fab_den2, ormas_donly, den2_abonly, dplus, nact1, act1_ll, act1_ul, tdcc
  use mod_cc, only : cc_type, cc_code

  implicit none
  integer(c_int), target, intent(in) :: nblock_
  integer(c_int), target, intent(in) :: type_block_(1:nblock_)
  integer(c_int), target, intent(in) :: nfun_block_(1:nblock_)
  integer(c_int), target, intent(in) :: nfcore2_, nfcore1_, nfcore_, ndcore_, ncore_, nact_, &
       nocc_, nvir_, nfun_, mval_(1:nfun_), froz_(1:nfun_)
  integer(c_int), target, intent(in) :: nelcore_(1:3), nelact_(1:3), neltot_(1:3)
  integer(c_int), target, intent(in) :: nsub_, norb_sub_(1:nsub_), lorb_sub_(1:2,1:nsub_)
  integer(c_int), target, intent(in) :: min_sub_(1:nsub_), max_sub_(1:nsub_)
  integer(c_int), target, intent(in) :: nstr_alph_, nstr_beta_, ndet_, lcic_, ndetx_
  real(c_double), target, intent(in) :: thradfc_
  logical(c_bool), target, intent(in) :: fab_den2_, donly_, den2_abonly_, dplus_
  integer(c_int), target, intent(in) :: nact1_, act1_ll_, act1_ul_
  logical(c_bool), target, intent(in) :: tdcc_
  integer(c_int), target, intent(in) :: cc_type_
  character(*), intent(in) :: cc_code_

  nblock => nblock_
  type_block => type_block_
  nfun_block => nfun_block_
  nfcore2 => nfcore2_
  nfcore1 => nfcore1_
  nfcore => nfcore_
  ndcore => ndcore_
  ncore => ncore_
  nact => nact_
  nocc => nocc_
  nvir => nvir_
  nfun => nfun_
  mval => mval_
  froz => froz_
  nelcore => nelcore_
  nelact => nelact_
  neltot => neltot_
  nsub => nsub_
  norb_sub => norb_sub_
  lorb_sub => lorb_sub_
  min_sub => min_sub_
  max_sub => max_sub_
  nstr_alph => nstr_alph_
  nstr_beta => nstr_beta_
  ndet => ndet_
  lcic => lcic_
  ndetx => ndetx_
  thradfc => thradfc_
  fab_den2 => fab_den2_
  ormas_donly => donly_
  den2_abonly => den2_abonly_
  dplus => dplus_
  nact1 => nact1_
  act1_ll => act1_ll_
  act1_ul => act1_ul_
  tdcc => tdcc_
  cc_type => cc_type_
  cc_code = cc_code_

end subroutine ormas_bind
!################################################################################
