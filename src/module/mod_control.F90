!################################################################################
module mod_control

  use, intrinsic :: iso_c_binding
  use xc_f03_lib_m

  implicit none

  character(len = 256) :: name
  logical(c_bool), parameter :: fedvr_normalized = .true.
  integer(c_long), pointer :: icomp
  integer(c_long), pointer :: igauge
  integer(c_long), pointer :: ioorot
  integer(c_long), pointer :: isplit
  integer(c_long), pointer :: iprojfc
  integer(c_long), pointer :: type_dcx
  logical(c_bool), pointer :: jfc_implicit
  logical(c_bool), pointer :: xfc_implicit
  integer(c_long), pointer :: h1rat_maxcyc
  real(c_double), pointer :: h1rat_thresh

  logical(c_bool), pointer :: DoCS1
  logical(c_bool), pointer :: DoCS2
  logical(c_bool), pointer :: SAE, PSP
  integer(c_long), pointer :: PSP_Type
  integer(c_long), pointer :: DFT_Type

  integer(c_long), pointer :: reg_type
  real(c_double), pointer :: throcc1 ! for 1/D
  real(c_double), pointer :: throcc2 ! for 1/(2-D)
  real(c_double), pointer :: throcc3 ! for 1/A
  integer(c_long), pointer :: ncut_occ3
  logical(c_bool), pointer :: exact3j

  !libxc
  TYPE(xc_f03_func_t), save :: xc_funcx, xc_funcc
  TYPE(xc_f03_func_info_t), save :: xc_infox, xc_infoc
  integer(c_int), save :: func_idx
  integer(c_int), save :: func_idc
  !libxc

end module mod_control
!################################################################################
