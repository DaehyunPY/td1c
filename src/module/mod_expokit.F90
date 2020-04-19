!################################################################################
module mod_expokit

  use, intrinsic :: iso_c_binding

  implicit none

  integer(c_int) :: icomp_exp
  integer(c_int) :: oorot_exp
  integer(c_int) :: isplit_exp
  integer(c_int) :: igauge_exp
  real(c_double) :: lfield_exp(1:9)
  real(c_double_complex), allocatable :: work_exp(:,:)

end module mod_expokit
!################################################################################
