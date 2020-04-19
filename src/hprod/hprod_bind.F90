!#######################################################################
subroutine hprod_bind(ene_fcore_, ene_dcore_, ene_core_, &
     ene_act_, ene_tot_, dip_exp_, vel_exp_, acc_exp_, projhigh_, projhigh_cutoff_)

  use, intrinsic :: iso_c_binding
  use mod_hprod

  implicit none
  real(c_double), target, intent(in) :: ene_fcore_, ene_dcore_, ene_core_, ene_act_, ene_tot_
  real(c_double), target, intent(in) :: dip_exp_, vel_exp_, acc_exp_
  logical(c_bool), target, intent(in) :: projhigh_
  real(c_double), target, intent(in) :: projhigh_cutoff_

  ene_fcore    => ene_fcore_
  ene_dcore    => ene_dcore_
  ene_core     => ene_core_
  ene_act      => ene_act_
  ene_tot      => ene_tot_
  dip_exp      => dip_exp_
  vel_exp      => vel_exp_
  acc_exp      => acc_exp_
  projhigh     => projhigh_
  projhigh_cutoff     => projhigh_cutoff_

end subroutine hprod_bind
!#######################################################################
