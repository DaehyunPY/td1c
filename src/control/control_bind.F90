!#######################################################################
subroutine control_bind(icomp_, igauge_, ioorot_, isplit_, iprojfc_, &
     type_dcx_, jfc_implicit_, xfc_implicit_, h1rat_maxcyc_, h1rat_thresh_, &
     docs1_, docs2_, sae_, psp_, psp_type_, dft_type_, reg_type_, throcc1_, throcc2_, throcc3_, ncut_occ3_, exact3j_)

  use, intrinsic :: iso_c_binding
  use mod_control

  implicit none
  integer(c_long), target, intent(in) :: icomp_, igauge_, ioorot_, isplit_, iprojfc_, type_dcx_
  logical(c_bool), target, intent(in) :: jfc_implicit_, xfc_implicit_
  integer(c_long), target, intent(in) :: h1rat_maxcyc_
  real(c_double), target, intent(in) :: h1rat_thresh_
  logical(c_bool), target, intent(in) :: docs1_, docs2_, sae_, psp_
  integer(c_long), target, intent(in) :: psp_type_, dft_type_
  integer(c_long), target, intent(in) :: reg_type_
  real(c_double), target, intent(in) :: throcc1_, throcc2_, throcc3_
  integer(c_long), target, intent(in) :: ncut_occ3_
  logical(c_bool), target, intent(in) :: exact3j_

  icomp  => icomp_
  igauge => igauge_
  ioorot => ioorot_
  isplit => isplit_
  iprojfc  => iprojfc_
  type_dcx => type_dcx_
  jfc_implicit => jfc_implicit_
  xfc_implicit => xfc_implicit_
  h1rat_maxcyc => h1rat_maxcyc_
  h1rat_thresh => h1rat_thresh_
  DoCS1 => docs1_
  DoCS2 => docs2_
  SAE => sae_
  PSP => psp_
  PSP_Type => psp_type_
  DFT_Type => dft_type_
  reg_type => reg_type_
  throcc1 => throcc1_
  throcc2 => throcc2_
  throcc3 => throcc3_
  ncut_occ3 => ncut_occ3_
  exact3j => exact3j_

end subroutine control_bind
!#######################################################################
