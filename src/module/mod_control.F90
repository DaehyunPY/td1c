!################################################################################
module mod_control

  use, intrinsic :: iso_c_binding
  use xc_f03_lib_m

  implicit none

  character(len = 256) :: name
  logical(c_bool), parameter :: fedvr_normalized = .true.
  integer(c_int), pointer :: icomp
  integer(c_int), pointer :: igauge
  integer(c_int), pointer :: ioorot
  integer(c_int), pointer :: isplit
  integer(c_int), pointer :: iprojfc
  integer(c_int), pointer :: type_dcx
  logical(c_bool), pointer :: jfc_implicit
  logical(c_bool), pointer :: xfc_implicit
  integer(c_int), pointer :: h1rat_maxcyc
  real(c_double), pointer :: h1rat_thresh

  logical(c_bool), pointer :: DoCS1
  logical(c_bool), pointer :: DoCS2
  logical(c_bool), pointer :: SAE, PSP
  integer(c_int), pointer :: PSP_Type
  integer(c_int), pointer :: DFT_Type

  integer(c_int), pointer :: reg_type
  real(c_double), pointer :: throcc1 ! for 1/D
  real(c_double), pointer :: throcc2 ! for 1/(2-D)
  real(c_double), pointer :: throcc3 ! for 1/A
  integer(c_int), pointer :: xact2_type ! 0:direct inverse of A, >=1:Jacob iteration of A.
  integer(c_int), pointer :: xact2_maxitr ! maximum iteration number for iterative A solver
  real(c_double), pointer :: xact2_thresh ! error threshold for iterative A solver
  integer(c_int), pointer :: ncut_occ3
  logical(c_bool), pointer :: exact3j
  logical(c_bool), pointer :: cionly
! tdcis-teramura
  logical(c_bool), pointer :: istdcis
  logical(c_bool), pointer :: tdcis_rvg
! tdcis-teramura

  !libxc
  TYPE(xc_f03_func_t), save :: xc_funcx, xc_funcc
  TYPE(xc_f03_func_info_t), save :: xc_infox, xc_infoc
  integer(c_int), save :: func_idx
  integer(c_int), save :: func_idc
  !libxc

end module mod_control
!################################################################################
